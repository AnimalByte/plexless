//! Header handling and deterministic bounded duplicate evaluation.
use crate::fastp_options::{FastpOptions, UmiLocation};
use std::collections::HashSet;

pub(crate) fn umi_bases(o: &FastpOptions, mate: usize) -> usize {
    if !o.umi {
        return 0;
    }
    if matches!(
        (o.umi_loc, mate),
        (Some(UmiLocation::Read1), 0)
            | (Some(UmiLocation::Read2), 1)
            | (Some(UmiLocation::PerRead), _)
    ) {
        o.umi_len + o.umi_skip
    } else {
        0
    }
}
pub(crate) fn index(id: &[u8], last: bool) -> &[u8] {
    let token = id.rsplit(|b| b.is_ascii_whitespace()).next().unwrap_or(id);
    let token = token.rsplit(|&b| b == b':').next().unwrap_or(token);
    if last {
        token.rsplit(|&b| b == b'+').next().unwrap_or(token)
    } else {
        token.split(|&b| b == b'+').next().unwrap_or(token)
    }
}
pub(crate) fn validated_index(id: &[u8], last: bool) -> Result<&[u8], String> {
    let token = id.rsplit(|b| b.is_ascii_whitespace()).next().unwrap_or(id);
    if !token.contains(&b':') || (last && !token.contains(&b'+')) {
        return Err("Index/UMI processing requires a CASAVA header with DNA index1 or index1+index2 after the final colon".into());
    }
    let bases = index(id, last);
    if bases.is_empty() || !bases.iter().all(|b| b"ACGTNacgtn".contains(b)) {
        return Err("FASTQ header index must be a nonempty A/C/G/T/N sequence".into());
    }
    Ok(bases)
}
pub(crate) fn header(
    o: &FastpOptions,
    id: &[u8],
    ids: [&[u8]; 2],
    sequences: [&[u8]; 2],
    starts: [usize; 2],
    paired: bool,
) -> Result<Vec<u8>, String> {
    let mut umi = Vec::new();
    if o.umi {
        match o.umi_loc {
            Some(UmiLocation::Index1) => umi.extend_from_slice(validated_index(ids[0], false)?),
            Some(UmiLocation::Index2) => umi.extend_from_slice(validated_index(ids[1], true)?),
            Some(UmiLocation::PerIndex) => {
                umi.extend_from_slice(validated_index(ids[0], false)?);
                if paired {
                    umi.push(b'_');
                    umi.extend_from_slice(validated_index(ids[1], true)?);
                }
            }
            Some(UmiLocation::Read1 | UmiLocation::Read2 | UmiLocation::PerRead) => {
                let mates: &[usize] = match o.umi_loc {
                    Some(UmiLocation::Read1) => &[0],
                    Some(UmiLocation::Read2) => &[1],
                    _ if paired => &[0, 1],
                    _ => &[0],
                };
                for (i, &mate) in mates.iter().enumerate() {
                    if i > 0 {
                        umi.push(b'_');
                    }
                    let start = starts[mate].min(sequences[mate].len());
                    if sequences[mate].len() - start < o.umi_len + o.umi_skip {
                        return Err(format!(
                            "R{} is too short for its proximal structure plus UMI length/skip",
                            mate + 1
                        ));
                    }
                    let end = (start + o.umi_len.min(sequences[mate].len() - start))
                        .min(sequences[mate].len());
                    umi.extend_from_slice(&sequences[mate][start..end]);
                }
            }
            None => {}
        }
    }
    let split = id
        .iter()
        .position(u8::is_ascii_whitespace)
        .unwrap_or(id.len());
    let mut main = id[..split].to_vec();
    if o.fix_mgi_id {
        for b in &mut main {
            if *b == b'/' {
                *b = b':';
            }
        }
    }
    if !umi.is_empty() {
        let suffix = if !o.fix_mgi_id && (main.ends_with(b"/1") || main.ends_with(b"/2")) {
            main.split_off(main.len() - 2)
        } else {
            Vec::new()
        };
        main.extend_from_slice(o.umi_delim.as_bytes());
        if !o.umi_prefix.is_empty() {
            main.extend_from_slice(o.umi_prefix.as_bytes());
            main.push(b'_');
        }
        main.extend_from_slice(&umi);
        main.extend_from_slice(&suffix);
    }
    main.extend_from_slice(&id[split..]);
    Ok(main)
}

#[derive(Debug)]
pub(crate) struct DuplicateCatalog {
    seen: HashSet<Box<[u8]>>,
    scratch: Vec<u8>,
    bytes: usize,
    limit: usize,
    disk: Option<rusqlite::Connection>,
    bloom: Option<Vec<u64>>,
    pub processed: u64,
    pub duplicates: u64,
}
impl DuplicateCatalog {
    pub(crate) fn configured(o: &FastpOptions) -> Result<Option<Self>, String> {
        if o.dedup {
            Self::new(o.duplication_memory_mib).map(Some)
        } else if o.evaluate_duplication || (o.fastp_defaults && !o.dont_eval_duplication) {
            Self::reporting(o.duplication_memory_mib).map(Some)
        } else {
            Ok(None)
        }
    }
    pub(crate) fn new(mib: usize) -> Result<Self, String> {
        Ok(Self {
            seen: HashSet::new(),
            scratch: Vec::new(),
            bytes: 0,
            limit: mib
                .checked_mul(1024 * 1024)
                .ok_or("Duplicate memory budget overflow")?,
            disk: None,
            bloom: None,
            processed: 0,
            duplicates: 0,
        })
    }
    pub(crate) fn reporting(mib: usize) -> Result<Self, String> {
        let mut catalog = Self::new(mib)?;
        catalog.bloom = Some(vec![0; catalog.limit.div_ceil(8)]);
        Ok(catalog)
    }
    fn spill(&mut self) -> Result<(), String> {
        let connection = rusqlite::Connection::open("")
            .map_err(|e| format!("Could not create temporary duplicate catalog: {e}"))?;
        connection
            .pragma_update(None, "journal_mode", "OFF")
            .map_err(|e| e.to_string())?;
        connection
            .pragma_update(None, "synchronous", "OFF")
            .map_err(|e| e.to_string())?;
        connection
            .pragma_update(None, "cache_size", -((self.limit / 1024) as i64))
            .map_err(|e| e.to_string())?;
        connection
            .execute_batch("CREATE TABLE seen (key BLOB PRIMARY KEY) WITHOUT ROWID; BEGIN;")
            .map_err(|e| e.to_string())?;
        {
            let mut insert = connection
                .prepare_cached("INSERT INTO seen (key) VALUES (?1)")
                .map_err(|e| e.to_string())?;
            for key in self.seen.drain() {
                insert.execute([key.as_ref()]).map_err(|e| e.to_string())?;
            }
        }
        self.seen = HashSet::new();
        self.bytes = 0;
        self.disk = Some(connection);
        Ok(())
    }
    pub(crate) fn observe(&mut self, a: &[u8], b: Option<&[u8]>) -> Result<bool, String> {
        self.scratch.clear();
        self.scratch
            .extend_from_slice(&(a.len() as u64).to_le_bytes());
        self.scratch.extend(a.iter().map(u8::to_ascii_uppercase));
        self.scratch.push(u8::from(b.is_some()));
        if let Some(b) = b {
            self.scratch.extend(b.iter().map(u8::to_ascii_uppercase));
        }
        self.processed += 1;
        if let Some(bits) = self.bloom.as_mut() {
            let mut first = 14695981039346656037u64;
            let mut second = 0x9e3779b97f4a7c15u64;
            for &byte in &self.scratch {
                first = (first ^ u64::from(byte)).wrapping_mul(1099511628211);
                second = second.rotate_left(7) ^ u64::from(byte);
                second = second.wrapping_mul(0x517cc1b727220a95);
            }
            let mut duplicate = true;
            for i in 0..4 {
                let bit = first.wrapping_add((i as u64).wrapping_mul(second | 1))
                    % (bits.len() as u64 * 64);
                let word = (bit / 64) as usize;
                let mask = 1 << (bit % 64);
                duplicate &= bits[word] & mask != 0;
                bits[word] |= mask;
            }
            self.duplicates += u64::from(duplicate);
            return Ok(duplicate);
        }
        if let Some(connection) = &self.disk {
            let inserted = connection
                .prepare_cached("INSERT OR IGNORE INTO seen (key) VALUES (?1)")
                .and_then(|mut statement| statement.execute([self.scratch.as_slice()]))
                .map_err(|e| format!("Duplicate catalog write failed: {e}"))?;
            let duplicate = inserted == 0;
            self.duplicates += u64::from(duplicate);
            return Ok(duplicate);
        }
        if self.seen.contains(self.scratch.as_slice()) {
            self.duplicates += 1;
            return Ok(true);
        }
        let next = self
            .bytes
            .checked_add(self.scratch.len() + 64)
            .ok_or("Duplicate memory usage overflow")?;
        if next > self.limit {
            self.spill()?;
            self.disk
                .as_ref()
                .expect("spill initialized catalog")
                .prepare_cached("INSERT INTO seen (key) VALUES (?1)")
                .and_then(|mut statement| statement.execute([self.scratch.as_slice()]))
                .map_err(|e| format!("Duplicate catalog write failed: {e}"))?;
            return Ok(false);
        }

        self.seen
            .try_reserve(1)
            .map_err(|e| format!("Could not grow duplicate catalog: {e}"))?;
        self.seen.insert(self.scratch.clone().into_boxed_slice());
        self.bytes = next;
        Ok(false)
    }
}

pub(crate) fn normalize_quality<'a>(
    raw: &'a [u8],
    phred64: bool,
    scratch: &'a mut Vec<u8>,
) -> Result<&'a [u8], String> {
    if !phred64 {
        return Ok(raw);
    }
    scratch.clear();
    for &q in raw {
        if !(64..=126).contains(&q) {
            return Err("Invalid Phred+64 quality byte".into());
        }
        scratch.push(q - 31);
    }
    Ok(scratch)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_dedup_preserves_pair_boundaries_and_ignores_case() {
        let mut d = DuplicateCatalog::new(1).unwrap();
        assert!(!d.observe(b"AC", Some(b"GT")).unwrap());
        assert!(d.observe(b"ac", Some(b"gt")).unwrap());
        assert!(!d.observe(b"ACG", Some(b"T")).unwrap());
        assert!(!d.observe(b"ACGT", None).unwrap());
    }
    #[test]
    fn exact_catalog_spills_without_losing_keys() {
        let mut d = DuplicateCatalog::new(1).unwrap();
        d.limit = 100;
        assert!(!d.observe(b"ACGTCAGG", None).unwrap());
        assert!(!d.observe(b"TGACTACG", None).unwrap());
        assert!(d.disk.is_some());
        assert!(d.observe(b"ACGTCAGG", None).unwrap());
        assert!(d.observe(b"TGACTACG", None).unwrap());
        assert!(!d.observe(b"GGTACCAT", None).unwrap());
        let mut reporting = DuplicateCatalog::reporting(1).unwrap();
        assert!(!reporting.observe(b"ACGTCAGG", None).unwrap());
        assert!(reporting.observe(b"ACGTCAGG", None).unwrap());
    }
    #[test]
    fn umi_from_mates_keeps_headers_and_mate_suffixes() {
        let o = FastpOptions {
            umi: true,
            umi_loc: Some(UmiLocation::PerRead),
            umi_len: 2,
            ..Default::default()
        };
        assert_eq!(
            header(
                &o,
                b"read/1 extra",
                [b"read/1 extra", b"read/2 extra"],
                [b"ACGATT", b"GTCAAA"],
                [0, 0],
                true
            )
            .unwrap(),
            b"read:AC_GT/1 extra"
        );
    }
}
