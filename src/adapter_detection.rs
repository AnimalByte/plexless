//! Bounded startup adapter inference inspired by fastp's frequent-seed extension.
//! No manufacturer sequences are compiled into Plexless.
use crate::fastq::InputReader;
use std::collections::{HashMap, HashSet};
use std::path::Path;
const SAMPLE_READS: usize = 5000;
const TAIL_BASES: usize = 200;
const SEED: usize = 12;
fn base(b: u8) -> Option<u8> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}
fn informative(s: &[u8]) -> bool {
    let mut c = [0; 4];
    let mut pairs = 0u16;
    for (i, &b) in s.iter().enumerate() {
        let Some(b) = base(b) else {
            return false;
        };
        c[b as usize] += 1;
        if i > 0 {
            pairs |= 1u16 << (base(s[i - 1]).unwrap() * 4 + b);
        }
    }
    c.iter().filter(|&&n| n > 0).count() >= 3
        && c.iter().all(|&n| n * 4 <= s.len() * 3)
        && pairs.count_ones() >= 6
}
pub(crate) fn infer(path: &Path, prefix: usize) -> Result<Option<Vec<u8>>, String> {
    infer_mate(path, prefix, None)
}
pub(crate) fn infer_mate(
    path: &Path,
    prefix: usize,
    mate: Option<usize>,
) -> Result<Option<Vec<u8>>, String> {
    if path == Path::new("-") {
        return Ok(None);
    }
    let mut reader = InputReader::try_open(path)?;
    let mut tails = Vec::new();
    let mut records = 0usize;
    while tails.len() < SAMPLE_READS {
        let Some(r) = reader.next_record() else {
            break;
        };
        let r = r.map_err(|e| format!("Adapter sampling parse error: {e}"))?;
        let selected = mate.is_none_or(|mate| records % 2 == mate);
        records += 1;
        if !selected {
            continue;
        }
        let s = r.seq();
        let start = prefix.min(s.len()).max(s.len().saturating_sub(TAIL_BASES));
        tails.push(s[start..].to_ascii_uppercase());
    }
    infer_sample(&tails)
}
fn infer_sample(tails: &[Vec<u8>]) -> Result<Option<Vec<u8>>, String> {
    let mut counts = HashMap::<u64, u32>::new();
    for tail in tails {
        let mut seen = HashSet::new();
        for window in tail.windows(SEED) {
            if !informative(window) {
                continue;
            }
            let key = window
                .iter()
                .fold(0u64, |key, &b| (key << 2) | u64::from(base(b).unwrap()));
            if seen.insert(key) {
                *counts.entry(key).or_default() += 1;
            }
        }
    }
    let required = 20.max(tails.len().div_ceil(20));
    let best = counts
        .into_iter()
        .filter(|(_, n)| *n as usize >= required)
        .max_by_key(|&(key, n)| (n, std::cmp::Reverse(key)));
    let Some((key, _)) = best else {
        return Ok(None);
    };
    let seed: Vec<u8> = (0..SEED)
        .map(|i| b"ACGT"[((key >> (2 * (SEED - 1 - i))) & 3) as usize])
        .collect();
    let mut positions: Vec<(usize, usize)> = tails
        .iter()
        .enumerate()
        .filter_map(|(i, s)| s.windows(SEED).position(|w| w == seed).map(|p| (i, p)))
        .collect();
    let mut adapter = seed;
    let mut offsets = HashSet::new();
    for &(_, p) in &positions {
        offsets.insert(p);
    }
    if offsets.len() < 3 {
        return Ok(None);
    }
    // Extend toward the physical adapter start until the flanking biological
    // bases cease to agree. All operations are bounded by the 200-base sample.
    for left in [true, false] {
        loop {
            if adapter.len() >= TAIL_BASES {
                break;
            }
            let mut c = [0usize; 4];
            for &(i, p) in &positions {
                let index = if left {
                    p.checked_sub(1)
                } else {
                    p.checked_add(adapter.len())
                };
                if let Some(b) = index.and_then(|j| tails[i].get(j)).and_then(|&b| base(b)) {
                    c[b as usize] += 1;
                }
            }
            let (b, n) = c
                .iter()
                .copied()
                .enumerate()
                .max_by_key(|&(b, n)| (n, std::cmp::Reverse(b)))
                .unwrap();
            if n < required || n * 10 < positions.len() * 9 {
                break;
            }
            positions.retain(|&(i, p)| {
                let index = if left {
                    p.checked_sub(1)
                } else {
                    p.checked_add(adapter.len())
                };
                index
                    .and_then(|j| tails[i].get(j))
                    .is_some_and(|&v| base(v) == Some(b as u8))
            });
            if left {
                for (_, p) in &mut positions {
                    *p -= 1;
                }
                adapter.insert(0, b"ACGT"[b]);
            } else {
                adapter.push(b"ACGT"[b]);
            }
        }
    }
    Ok((adapter.len() >= 20 && informative(&adapter)).then_some(adapter))
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn inferred_adapter_uses_variable_boundaries_and_no_known_catalog() {
        let adapter = b"ACGTTGCACTGATCGAGGTACCTATGCCAGTA";
        let mut tails = Vec::new();
        let mut x = 237u64;
        for i in 0..100 {
            let mut s: Vec<u8> = (0..50 + i % 30)
                .map(|_| {
                    x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGT"[(x >> 32) as usize % 4]
                })
                .collect();
            s.extend_from_slice(adapter);
            tails.push(s);
        }
        assert_eq!(infer_sample(&tails).unwrap().unwrap(), adapter);
    }
}
