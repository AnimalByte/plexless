//! Assignment-independent distal boundary detection and startup compilation.
use crate::barcodes::BarcodeCatalog;
use crate::cli::DemuxArgs;
use crate::decoder::hamming_distance;
use crate::encoder::{EncodedBarcode, encode};
use crate::samples::{Sample, SampleSheet};
use crate::structure::{BarcodeExtractionPlans, Orientation, ReadLayout, ReadMate, ReadStructure};
use crate::technical::{TechnicalSequences, dna, structures};

/// Original-read coordinates carried through every output adapter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct BiologicalBounds {
    pub start: usize,
    pub end: usize,
}

impl BiologicalBounds {
    pub(crate) fn full(len: usize) -> Self {
        Self { start: 0, end: len }
    }
    pub(crate) fn range(self) -> std::ops::Range<usize> {
        self.start..self.end
    }
    pub(crate) fn changed(self, len: usize) -> bool {
        self.start != 0 || self.end != len
    }
    pub(crate) fn len(self) -> usize {
        self.end - self.start
    }
}

#[derive(Debug)]
struct VerificationChunk {
    expected: EncodedBarcode,
    /// Two bits per known reference position. Unknown technical bases alone are masked.
    known_mask: u64,
}

#[derive(Debug)]
pub struct PassThroughSignature {
    anchor: u64,
    anchor_mask: u64,
    anchor_len: usize,
    anchor_offset: usize,
    verification: Box<[VerificationChunk]>,
    verification_len: usize,
    mismatch_limit: usize,
    partial: PartialPlan,
}

const fn base_table() -> [u8; 256] {
    let mut table = [255; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table
}
const BASE: [u8; 256] = base_table();

impl PassThroughSignature {
    // N here denotes an EXPECTED unknown technical base, never an observed wildcard.
    #[cfg(test)]
    fn compile(expected: &[u8], mismatch_limit: usize) -> Option<Self> {
        Self::compile_with_partial(expected, mismatch_limit, None)
    }

    #[cfg(test)]
    fn compile_with_partial(
        expected: &[u8],
        mismatch_limit: usize,
        minimum: Option<usize>,
    ) -> Option<Self> {
        Self::compile_with_policy(expected, mismatch_limit, minimum, 20)
    }

    fn compile_with_policy(
        expected: &[u8],
        mismatch_limit: usize,
        minimum: Option<usize>,
        min_full_known: usize,
    ) -> Option<Self> {
        let mut best = None;
        for length in (12..=16).rev() {
            if expected.len() < length {
                continue;
            }
            let mut best_diversity = 0;
            for (offset, window) in expected.windows(length).enumerate() {
                let mut counts = [0usize; 4];
                let mut pairs = 0u16;
                let mut previous = None;
                let mut valid = true;
                for &base in window {
                    let code = BASE[base as usize];
                    if code > 3 {
                        valid = false;
                        break;
                    }
                    counts[code as usize] += 1;
                    if let Some(prev) = previous {
                        pairs |= 1u16 << (prev * 4 + code);
                    }
                    previous = Some(code);
                }
                // Reject low-complexity anchors, even when they are long enough.
                let diversity = pairs.count_ones();
                if valid
                    && counts.iter().filter(|&&n| n > 0).count() >= 3
                    && counts.iter().all(|&n| n * 4 <= length * 3)
                    && diversity >= 6
                    && diversity > best_diversity
                {
                    best = Some((offset, length));
                    best_diversity = diversity;
                }
            }
            if best.is_some() {
                break;
            }
        }
        if expected.iter().filter(|&&b| BASE[b as usize] < 4).count() < min_full_known {
            best = None;
        }
        let partial = PartialPlan::compile(expected, minimum);
        if best.is_none() && partial.eligible.iter().all(|&e| !e) {
            return None;
        }
        let (anchor_offset, anchor_len) = best.unwrap_or((0, 0));
        let anchor = encode(&expected[anchor_offset..anchor_offset + anchor_len])
            .map_or(0, |encoded| encoded.value);
        let verification = expected
            .chunks(32)
            .map(|chunk| {
                let expected =
                    encode(chunk).expect("compiled structure contains DNA or technical N");
                let known_mask = chunk.iter().fold(0u64, |mask, &base| {
                    (mask << 2) | if base == b'N' { 0 } else { 3 }
                });
                VerificationChunk {
                    expected,
                    known_mask,
                }
            })
            .collect();
        Some(Self {
            anchor,
            anchor_mask: (1u64 << (2 * anchor_len)) - 1,
            anchor_len,
            anchor_offset,
            verification,
            verification_len: expected.len(),
            mismatch_limit,
            partial,
        })
    }

    pub(crate) fn full_at(&self, sequence: &[u8], start: usize) -> bool {
        if self.anchor_len == 0 {
            return false;
        }
        let Some(end) = start.checked_add(self.verification_len) else {
            return false;
        };
        let Some(candidate) = sequence.get(start..end) else {
            return false;
        };
        let anchor = &candidate[self.anchor_offset..self.anchor_offset + self.anchor_len];
        encode(anchor).is_ok_and(|observed| observed.n_mask == 0 && observed.value == self.anchor)
            && self.verify(candidate)
    }

    // Keep verification setup and register pressure off the common negative scan.
    #[inline(never)]
    fn verify(&self, candidate: &[u8]) -> bool {
        let mut mismatches = 0usize;
        for (observed, chunk) in candidate.chunks(32).zip(self.verification.iter()) {
            let mut value = 0u64;
            for (index, &base) in observed.iter().enumerate() {
                let code = BASE[base as usize];
                let shift = 2 * (observed.len() - index - 1);
                if code > 3 && (chunk.known_mask >> shift) & 3 != 0 {
                    return false; // Conservative: any unknown at a known position rejects.
                }
                value = (value << 2) | u64::from(code & 3);
            }
            mismatches += usize::from(hamming_distance(
                value & chunk.known_mask,
                chunk.expected.value & chunk.known_mask,
            ));
            if mismatches > self.mismatch_limit {
                return false;
            }
        }
        true
    }
}

/// Finds the complete distal structure's start, relative to `seq`.
/// No assignment, orientation, parsing, allocation, or output decisions occur here.
pub fn detect_pass_through(seq: &[u8], signature: &PassThroughSignature) -> Option<usize> {
    if signature.anchor_len == 0 {
        return None;
    }
    let last_start = seq.len().checked_sub(signature.verification_len)?;
    let scan_start = signature.anchor_offset;
    let scan_end = last_start + signature.anchor_offset + signature.anchor_len;
    let mut window = 0u64;
    let mut remaining = seq[scan_start..scan_end].iter();
    // Track validity by the iterator's remaining length. This avoids advancing
    // separate position and clean-base counters on every negative-path base.
    let mut valid_at_remaining = remaining.len() - signature.anchor_len;
    while let Some(&base) = remaining.next() {
        let code = BASE[base as usize];
        if code > 3 {
            window = 0;
            // No later complete clean anchor can fit if this subtraction fails.
            valid_at_remaining = remaining.len().checked_sub(signature.anchor_len)?;
            continue;
        }
        window = ((window << 2) | u64::from(code)) & signature.anchor_mask;
        if window == signature.anchor && remaining.len() <= valid_at_remaining {
            let anchor_start = scan_end - remaining.len() - signature.anchor_len;
            let distal_start = anchor_start - signature.anchor_offset;
            if signature.verify(&seq[distal_start..distal_start + signature.verification_len]) {
                return Some(distal_start);
            }
        }
    }
    None
}

#[derive(Debug, Default)]
pub(crate) struct PassThroughPlans {
    pub(crate) enabled: bool,
    samples: Vec<[Option<PassThroughSignature>; 2]>,
    shared: [Option<PassThroughSignature>; 2],
    additional: Vec<PassThroughSignature>,
}

impl PassThroughPlans {
    pub(crate) fn compile(
        layout: &ReadLayout,
        catalog: &BarcodeCatalog,
        samples: &SampleSheet,
        args: &DemuxArgs,
    ) -> Result<Self, String> {
        let technical =
            TechnicalSequences::load(args.preprocessing.technical_sequences.as_deref(), layout)?;
        validate_distal(args, matches!(layout, ReadLayout::Paired { .. }))?;
        if args.no_pass_through {
            return Ok(Self::default());
        }
        let (r1, r2, paired) = match layout {
            ReadLayout::Single { r1 } => (r1.as_ref(), r1.as_ref(), false),
            ReadLayout::Paired { r1, r2 } => (r1.as_ref(), r2.as_ref(), true),
        };
        if !paired && args.r2_distal_structure.is_some() {
            return Err("R2 distal override requires paired input".into());
        }
        let overrides = [
            args.r1_distal_structure.as_deref(),
            args.r2_distal_structure.as_deref(),
        ]
        .map(|value| value.map(ReadStructure::parse).transpose());
        let [override1, override2] = overrides;
        let overrides = [override1?, override2?];
        for (index, structure) in overrides.iter().enumerate() {
            if let Some(structure) = structure {
                let mate = if index == 0 {
                    ReadMate::R1
                } else {
                    ReadMate::R2
                };
                if structure.mate != mate {
                    return Err("Distal override mate must match its option".into());
                }
            }
        }
        // Explicit adapters are independent of assignment: compile once per mate,
        // even for large sample catalogs.
        let explicit = [
            &args.preprocessing.r1_distal_sequence,
            &args.preprocessing.r2_distal_sequence,
        ];
        let mut shared = [None, None];
        for mate in 0..if paired { 2 } else { 1 } {
            if let Some(sequence) = explicit[mate] {
                shared[mate] = PassThroughSignature::compile_with_policy(
                    &dna(sequence)?,
                    args.pass_through_mismatches,
                    args.preprocessing
                        .partial_pass_through
                        .then_some(args.preprocessing.min_partial_adapter_overlap),
                    args.preprocessing.min_full_adapter_known_bases,
                );
            }
        }
        let plans = layout.compile_extraction_plans()?;
        let mut compiled = Vec::with_capacity(samples.samples.len());
        let mut unavailable = [0usize; 2];
        for sample in &samples.samples {
            let mut signatures = [None, None];
            for mate in 0..if paired { 2 } else { 1 } {
                if explicit[mate].is_some() {
                    if shared[mate].is_none() {
                        unavailable[mate] += 1;
                    }
                    continue;
                }
                let expected = if let Some(structure) = &overrides[mate] {
                    compile_override(structure, &plans, catalog, sample, &technical)?
                } else if let Some(source) = if mate == 0 { r2 } else { r1 } {
                    compile_opposite(
                        source,
                        &plans,
                        catalog,
                        sample,
                        &technical,
                        crate::fastp_runtime::umi_bases(
                            &args.preprocessing.fastp,
                            usize::from(source.mate == ReadMate::R2),
                        ),
                    )?
                } else {
                    Vec::new()
                };
                signatures[mate] = PassThroughSignature::compile_with_policy(
                    &expected,
                    args.pass_through_mismatches,
                    args.preprocessing
                        .partial_pass_through
                        .then_some(args.preprocessing.min_partial_adapter_overlap),
                    args.preprocessing.min_full_adapter_known_bases,
                );
                if signatures[mate].is_none() {
                    unavailable[mate] += 1;
                }
            }
            compiled.push(signatures);
        }
        for (mate, count) in unavailable
            .iter()
            .enumerate()
            .take(if paired { 2 } else { 1 })
        {
            if *count > 0 {
                eprintln!(
                    "Pass-through: R{} unavailable for {count}/{} samples: insufficient known-base evidence or no sufficiently varied contiguous known anchor of at least 12 bases",
                    mate + 1,
                    samples.samples.len()
                );
            }
        }
        let mut result = Self {
            enabled: true,
            samples: compiled,
            shared,
            additional: Vec::new(),
        };
        result.load_fasta(args)?;
        Ok(result)
    }

    pub(crate) fn compile_trim(layout: &ReadLayout, args: &DemuxArgs) -> Result<Self, String> {
        let technical =
            TechnicalSequences::load(args.preprocessing.technical_sequences.as_deref(), layout)?;
        let paired = matches!(layout, ReadLayout::Paired { .. });
        validate_distal(args, paired)?;
        if args.no_pass_through {
            return Ok(Self::default());
        }
        let proximal = structures(layout);
        let mut signatures = [None, None];
        for mate in 0..if paired { 2 } else { 1 } {
            let explicit = if mate == 0 {
                &args.preprocessing.r1_distal_sequence
            } else {
                &args.preprocessing.r2_distal_sequence
            };
            let distal = if mate == 0 {
                &args.r1_distal_structure
            } else {
                &args.r2_distal_structure
            };
            let expected = if let Some(sequence) = explicit {
                dna(sequence)?
            } else if let Some(text) = distal {
                let structure = ReadStructure::parse(text)?;
                if usize::from(structure.mate == ReadMate::R2) != mate {
                    return Err("Distal override mate must match its option".into());
                }
                if !structure.barcode_symbols().is_empty() {
                    return Err("Trim distal structures cannot reference sample barcodes; use --r1-distal-sequence/--r2-distal-sequence".into());
                }
                let bases = technical.0[mate].as_deref().ok_or("Trim distal T structures require known --technical-sequences; use an explicit distal sequence without proximal configuration")?;
                if bases.len() != structure.prefix_len {
                    return Err("Distal T length differs from known technical sequence".into());
                }
                bases.to_vec()
            } else if let Some(source) = if paired {
                proximal[1 - mate]
            } else {
                proximal[0]
            } {
                // Only a completely known technical proximal structure can infer
                // an assignment-independent opposite adapter in trim mode.
                if !source.barcode_symbols().is_empty() {
                    Vec::new()
                } else if let Some(bases) = technical.mate(source.mate) {
                    let mut out = distal_buffer(
                        bases.len(),
                        crate::fastp_runtime::umi_bases(
                            &args.preprocessing.fastp,
                            usize::from(source.mate == ReadMate::R2),
                        ),
                    )?;
                    append_oriented(&mut out, bases, true);
                    out
                } else {
                    Vec::new()
                }
            } else {
                Vec::new()
            };
            signatures[mate] = PassThroughSignature::compile_with_policy(
                &expected,
                args.pass_through_mismatches,
                args.preprocessing
                    .partial_pass_through
                    .then_some(args.preprocessing.min_partial_adapter_overlap),
                args.preprocessing.min_full_adapter_known_bases,
            );
        }
        let mut result = Self {
            enabled: true,
            samples: Vec::new(),
            shared: signatures,
            additional: Vec::new(),
        };
        result.load_fasta(args)?;
        Ok(result)
    }

    fn load_fasta(&mut self, args: &DemuxArgs) -> Result<(), String> {
        let Some(path) = &args.preprocessing.fastp.adapter_fasta else {
            return Ok(());
        };
        use std::io::BufRead;
        let file =
            std::fs::File::open(path).map_err(|e| format!("Could not open adapter FASTA: {e}"))?;
        let mut sequences = Vec::new();
        let mut sequence = String::new();
        let mut header = false;
        for line in std::io::BufReader::new(file).lines() {
            let line = line.map_err(|e| e.to_string())?;
            if line.starts_with('>') {
                if header {
                    if sequence.is_empty() {
                        return Err("Empty adapter FASTA sequence".into());
                    }
                    sequences.push(std::mem::take(&mut sequence));
                }
                header = true;
            } else {
                if !header || line.is_empty() {
                    return Err("Malformed adapter FASTA".into());
                }
                sequence.push_str(&line);
            }
        }
        if !header || sequence.is_empty() {
            return Err("Adapter FASTA must contain nonempty sequences".into());
        }
        sequences.push(sequence);
        for sequence in sequences {
            let sequence = dna(&sequence)?;
            let signature = PassThroughSignature::compile_with_policy(
                &sequence,
                args.pass_through_mismatches,
                args.preprocessing
                    .partial_pass_through
                    .then_some(args.preprocessing.min_partial_adapter_overlap),
                args.preprocessing.min_full_adapter_known_bases,
            )
            .ok_or("Adapter FASTA sequence lacks the required known-base evidence or complexity")?;
            self.additional.push(signature);
        }
        Ok(())
    }
    pub(crate) fn signatures(
        &self,
        sample: u32,
        mate: usize,
    ) -> impl Iterator<Item = &PassThroughSignature> {
        self.signature(sample, mate)
            .into_iter()
            .chain(self.additional.iter())
    }
    pub(crate) fn signature(&self, sample: u32, mate: usize) -> Option<&PassThroughSignature> {
        self.shared
            .get(mate)?
            .as_ref()
            .or_else(|| self.samples.get(sample as usize)?.get(mate)?.as_ref())
    }
}

fn barcode<'a>(
    symbol: u8,
    plans: &BarcodeExtractionPlans,
    catalog: &'a BarcodeCatalog,
    sample: &Sample,
) -> Result<&'a [u8], String> {
    let index = plans.plan_index(symbol).ok_or_else(|| {
        format!(
            "Distal structure references undeclared barcode {}",
            symbol as char
        )
    })?;
    Ok(catalog
        .set(symbol)
        .and_then(|set| set.barcode(sample.barcode_ids[index]))
        .ok_or("Distal structure references missing assigned barcode")?
        .sequence
        .as_bytes())
}

fn append_oriented(output: &mut Vec<u8>, bases: &[u8], reverse: bool) {
    if reverse {
        output.extend(
            bases
                .iter()
                .rev()
                .map(|&base| b"TGCA"[BASE[base as usize] as usize]),
        );
    } else {
        output.extend_from_slice(bases);
    }
}

// Read UMIs/skipped cycles follow the proximal prefix. Their reverse
// complement precedes the opposite observed distal prefix. Mask their variable
// sequence without counting it as known adapter evidence.
fn distal_buffer(prefix: usize, umi: usize) -> Result<Vec<u8>, String> {
    let length = prefix
        .checked_add(umi)
        .ok_or("Distal structure length overflow")?;
    let mut output = Vec::new();
    output
        .try_reserve_exact(length)
        .map_err(|_| "Distal structure is too large")?;
    output.resize(umi, b'N');
    Ok(output)
}
fn compile_opposite(
    source: &ReadStructure,
    plans: &BarcodeExtractionPlans,
    catalog: &BarcodeCatalog,
    sample: &Sample,
    technical: &TechnicalSequences,
    umi: usize,
) -> Result<Vec<u8>, String> {
    let mut output = distal_buffer(source.prefix_len, umi)?;
    let known = technical.mate(source.mate);
    let mut technical_end = known.map_or(0, <[u8]>::len);
    for segment in source.segments.iter().rev() {
        let length = segment.end - segment.start;
        if segment.symbol == b'T' {
            if let Some(known) = known {
                append_oriented(
                    &mut output,
                    &known[technical_end - length..technical_end],
                    true,
                );
                technical_end -= length;
            } else {
                output.resize(output.len() + length, b'N');
            }
            continue;
        }
        let plan = plans
            .plan(plans.plan_index(segment.symbol).expect("compiled symbol"))
            .expect("compiled plan");
        let mut offset = 0;
        for piece in &plan.pieces {
            if piece.mate == source.mate && piece.start == segment.start {
                break;
            }
            offset += piece.end - piece.start;
        }
        let bases = barcode(segment.symbol, plans, catalog, sample)?;
        // Invert physical orientation while reversing segment order. An rc piece
        // was already opposite to canonical, so the two transforms cancel.
        append_oriented(
            &mut output,
            &bases[offset..offset + length],
            segment.orientation == Orientation::Forward,
        );
    }
    Ok(output)
}

fn compile_override(
    structure: &ReadStructure,
    plans: &BarcodeExtractionPlans,
    catalog: &BarcodeCatalog,
    sample: &Sample,
    technical: &TechnicalSequences,
) -> Result<Vec<u8>, String> {
    let mut offsets = [0usize; 26];
    let mut output = Vec::new();
    output
        .try_reserve_exact(structure.prefix_len)
        .map_err(|_| "Distal override is too large")?;
    let known = technical.mate(structure.mate);
    let technical_length: usize = structure
        .segments
        .iter()
        .filter(|s| s.symbol == b'T')
        .map(|s| s.end - s.start)
        .sum();
    if known.is_some_and(|k| k.len() != technical_length) {
        return Err("Known technical sequence length differs from distal override T length; use an explicit distal sequence".into());
    }
    let mut technical_offset = 0;
    for segment in &structure.segments {
        let length = segment.end - segment.start;
        if segment.symbol == b'T' {
            if let Some(known) = known {
                output.extend_from_slice(&known[technical_offset..technical_offset + length]);
                technical_offset += length;
            } else {
                output.resize(output.len() + length, b'N');
            }
            continue;
        }
        let offset = &mut offsets[(segment.symbol - b'A') as usize];
        let bases = barcode(segment.symbol, plans, catalog, sample)?;
        let piece = bases
            .get(*offset..*offset + length)
            .ok_or("Distal override exceeds assigned logical barcode length")?;
        append_oriented(
            &mut output,
            piece,
            segment.orientation == Orientation::ReverseComplement,
        );
        *offset += length;
    }
    for symbol in structure.barcode_symbols() {
        if offsets[(symbol - b'A') as usize] != barcode(symbol, plans, catalog, sample)?.len() {
            return Err(
                "Distal override must consume each referenced logical barcode completely".into(),
            );
        }
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EXPECTED: &[u8] = b"NNNACGTTGCACTGATCGAGGTACCTTACGA";

    fn signature() -> PassThroughSignature {
        PassThroughSignature::compile(EXPECTED, 1).unwrap()
    }

    #[test]
    fn boundaries_and_sequence_quality_coordinates() {
        let signature = signature();
        for insert in [b"".as_slice(), b"G", b"GATTACA"] {
            for tail in [b"".as_slice(), b"AAAAAAA"] {
                let seq = [insert, EXPECTED, tail].concat();
                let start = detect_pass_through(&seq, &signature).unwrap();
                let anchor_start = insert.len() + signature.anchor_offset;
                assert_eq!(start, anchor_start - signature.anchor_offset);
                assert_eq!(&seq[..start], insert);
                let qual: Vec<_> = (0..seq.len()).map(|i| 33 + (i % 60) as u8).collect();
                assert_eq!(qual[..start].len(), insert.len());
            }
        }
        assert_eq!(detect_pass_through(b"", &signature), None);
        assert_eq!(detect_pass_through(b"GATTACA", &signature), None);
        assert_eq!(detect_pass_through(&EXPECTED[3..], &signature), None);
        assert_eq!(
            detect_pass_through(&EXPECTED[..EXPECTED.len() - 1], &signature),
            None
        );
        assert_eq!(detect_pass_through(&[b'A'; 200], &signature), None);
    }

    #[test]
    fn failed_anchor_verification_continues_to_real_distal_start() {
        let signature = signature();
        let mut fake = EXPECTED.to_vec();
        for (index, base) in fake.iter_mut().enumerate().skip(3) {
            if !(signature.anchor_offset..signature.anchor_offset + signature.anchor_len)
                .contains(&index)
            {
                *base = if *base == b'A' { b'C' } else { b'A' };
            }
        }
        let bio = [b"GATTACA".as_slice(), &fake, b"GGGGGG"].concat();
        assert_eq!(detect_pass_through(&bio, &signature), None);
        let seq = [bio.as_slice(), EXPECTED, b"AAAA"].concat();
        let start = detect_pass_through(&seq, &signature).unwrap();
        assert_eq!(start, bio.len());
        assert_eq!(&seq[..start], bio);
    }

    #[test]
    fn mismatch_budget_and_conservative_observed_unknowns() {
        let signature = signature();
        let positions: Vec<_> = (3..EXPECTED.len())
            .filter(|index| {
                !(signature.anchor_offset..signature.anchor_offset + signature.anchor_len)
                    .contains(index)
            })
            .collect();
        let mut candidate = EXPECTED.to_vec();
        candidate[positions[0]] = if candidate[positions[0]] == b'A' {
            b'C'
        } else {
            b'A'
        };
        assert_eq!(detect_pass_through(&candidate, &signature), Some(0));
        candidate[positions[1]] = if candidate[positions[1]] == b'A' {
            b'C'
        } else {
            b'A'
        };
        assert_eq!(detect_pass_through(&candidate, &signature), None);
        for index in 3..EXPECTED.len() {
            for invalid in *b"Nn?-" {
                let mut candidate = EXPECTED.to_vec();
                candidate[index] = invalid;
                assert_eq!(detect_pass_through(&candidate, &signature), None);
                candidate.extend_from_slice(EXPECTED);
                assert_eq!(
                    detect_pass_through(&candidate, &signature),
                    Some(EXPECTED.len())
                );
            }
        }
        let mut candidate = EXPECTED.to_ascii_lowercase();
        candidate[..3].copy_from_slice(b"?-N");
        assert_eq!(detect_pass_through(&candidate, &signature), Some(0));
    }

    #[test]
    fn long_chunked_verification_checks_all_known_bases() {
        let expected = [EXPECTED, b"NN", EXPECTED, EXPECTED].concat();
        let signature = PassThroughSignature::compile(&expected, 0).unwrap();
        assert!(signature.verification.len() > 2);
        assert_eq!(detect_pass_through(&expected, &signature), Some(0));
        for index in 0..expected.len() {
            if expected[index] == b'N' {
                continue;
            }
            let mut candidate = expected.clone();
            candidate[index] = if candidate[index] == b'A' { b'C' } else { b'A' };
            assert_eq!(
                detect_pass_through(&candidate, &signature),
                None,
                "position {index}"
            );
        }
    }

    #[test]
    fn weak_or_short_anchors_are_unavailable() {
        for expected in [
            b"ACGTTGCACTG".as_slice(),
            b"AAAAAAAAAAAAAAAAAAAA",
            b"ACACACACACACACAC",
            b"ACGTTGCANACTGATCGA",
        ] {
            assert!(PassThroughSignature::compile(expected, 1).is_none());
        }
    }

    #[test]
    fn rolling_scan_agrees_with_scalar_reference() {
        let signature = signature();
        let mut seed = 12345u64;
        for length in 0..240 {
            let mut seq: Vec<u8> = (0..length)
                .map(|_| {
                    seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGTN"[((seed >> 32) % 5) as usize]
                })
                .collect();
            if length >= EXPECTED.len() && length % 2 == 0 {
                let start = (seed as usize) % (length - EXPECTED.len() + 1);
                seq[start..start + EXPECTED.len()].copy_from_slice(EXPECTED);
            }
            let reference = (0..=length.saturating_sub(EXPECTED.len())).find(|&start| {
                let Some(candidate) = seq.get(start..start + EXPECTED.len()) else {
                    return false;
                };
                let anchor =
                    signature.anchor_offset..signature.anchor_offset + signature.anchor_len;
                if candidate[anchor.clone()] != EXPECTED[anchor] {
                    return false;
                }
                let mut mismatches = 0;
                for (&observed, &expected) in candidate.iter().zip(EXPECTED) {
                    if expected == b'N' {
                        continue;
                    }
                    if !b"ACGT".contains(&observed) {
                        return false;
                    }
                    mismatches += usize::from(observed != expected);
                }
                mismatches <= 1
            });
            assert_eq!(detect_pass_through(&seq, &signature), reference);
        }
    }
}

fn validate_distal(args: &DemuxArgs, paired: bool) -> Result<(), String> {
    let p = &args.preprocessing;
    if !paired && (p.r2_distal_sequence.is_some() || args.r2_distal_structure.is_some()) {
        return Err("R2 distal configuration requires paired input".into());
    }
    for (sequence, structure) in [
        (&p.r1_distal_sequence, &args.r1_distal_structure),
        (&p.r2_distal_sequence, &args.r2_distal_structure),
    ] {
        if sequence.is_some() && structure.is_some() {
            return Err("Conflicting distal sequence and structure for the same mate".into());
        }
        if let Some(sequence) = sequence {
            dna(sequence)?;
        }
        if args.no_pass_through
            && (sequence.is_some() || structure.is_some() || p.partial_pass_through)
        {
            return Err("Distal configuration conflicts with --no-pass-through".into());
        }
    }
    if p.partial_pass_through && p.min_partial_adapter_overlap == 0 {
        return Err("Minimum partial adapter overlap must be positive".into());
    }
    Ok(())
}

/// Shift-And automaton for all adapter prefixes. Unknown expected positions
/// accept any observed byte, but never add evidence. No full-match anchor is used.
#[derive(Debug)]
struct PartialPlan {
    masks: Vec<[u64; 5]>,
    known: Vec<usize>,
    eligible: Vec<bool>,
    minimum_overlap: usize,
}
impl PartialPlan {
    fn compile(expected: &[u8], minimum: Option<usize>) -> Self {
        let Some(minimum) = minimum else {
            return Self {
                masks: Vec::new(),
                known: Vec::new(),
                eligible: Vec::new(),
                minimum_overlap: usize::MAX,
            };
        };
        let mut masks = vec![[0u64; 5]; expected.len().div_ceil(64)];
        let mut known = Vec::with_capacity(expected.len());
        let mut eligible = Vec::with_capacity(expected.len());
        let mut counts = [0usize; 4];
        let mut pairs = 0u16;
        let mut previous = None;
        let mut evidence = 0;
        for (i, &base) in expected.iter().enumerate() {
            let code = BASE[base as usize];
            if code < 4 {
                masks[i / 64][code as usize] |= 1u64 << (i % 64);
                counts[code as usize] += 1;
                evidence += 1;
                if let Some(prev) = previous {
                    pairs |= 1u16 << (prev * 4 + code);
                }
                previous = Some(code);
            } else {
                for mask in &mut masks[i / 64] {
                    *mask |= 1u64 << (i % 64);
                }
                previous = None;
            }
            known.push(evidence);
            eligible.push(
                evidence >= minimum
                    && counts.iter().filter(|&&n| n > 0).count() >= 3
                    && counts.iter().all(|&n| n * 4 <= evidence * 3)
                    && pairs.count_ones() >= 6,
            );
        }
        let minimum_overlap = eligible
            .iter()
            .position(|&value| value)
            .map_or(usize::MAX, |i| i + 1);
        Self {
            masks,
            known,
            eligible,
            minimum_overlap,
        }
    }
}

impl PassThroughSignature {
    pub(crate) fn partial_evidence(&self, sequence_len: usize, start: usize) -> usize {
        sequence_len
            .checked_sub(start)
            .and_then(|n| n.checked_sub(1))
            .and_then(|i| self.partial.known.get(i))
            .copied()
            .unwrap_or(0)
    }
    pub(crate) fn partial_at(&self, sequence: &[u8], start: usize) -> bool {
        let Some(overlap) = sequence.len().checked_sub(start) else {
            return false;
        };
        if overlap == 0
            || overlap >= self.verification_len
            || self.partial.masks.is_empty()
            || !self.partial.eligible[overlap - 1]
        {
            return false;
        }
        sequence[start..].iter().enumerate().all(|(i, &base)| {
            let code = usize::from(BASE[base as usize].min(4));
            self.partial.masks[i / 64][code] & (1 << (i % 64)) != 0
        })
    }

    /// Return the strongest suffix/prefix boundary only after full detection fails.
    /// Maximize known evidence; ties use the shorter overlap to avoid extra clipping.
    pub(crate) fn detect_partial(&self, seq: &[u8], scratch: &mut Vec<u64>) -> Option<usize> {
        let plan = &self.partial;
        if plan.masks.is_empty() {
            return None;
        }
        let limit = seq.len().min(self.verification_len.saturating_sub(1));
        if limit < plan.minimum_overlap {
            return None;
        }
        let words = limit.div_ceil(64);
        let mut inline = [0u64; 4];
        let states = if words <= inline.len() {
            &mut inline[..words]
        } else {
            scratch.resize(words, 0);
            scratch.fill(0);
            scratch.as_mut_slice()
        };
        for &base in &seq[seq.len() - limit..] {
            let code = usize::from(BASE[base as usize].min(4));
            let mut carry = 1u64;
            for (state, masks) in states.iter_mut().zip(&plan.masks) {
                let next = *state >> 63;
                *state = ((*state << 1) | carry) & masks[code];
                carry = next;
            }
        }
        let mut best = None;
        let mut evidence = 0;
        for (word_index, &word) in states.iter().enumerate() {
            let mut bits = word;
            while bits != 0 {
                let i = word_index * 64 + bits.trailing_zeros() as usize;
                bits &= bits - 1;
                let overlap = i + 1;
                if overlap <= limit && plan.eligible[i] && plan.known[i] > evidence {
                    best = Some(seq.len() - overlap);
                    evidence = plan.known[i];
                }
            }
        }
        best
    }
}

#[cfg(test)]
mod partial_tests {
    use super::*;
    const ADAPTER: &[u8] = b"ACGTTGCACTGATCGAGGTACCTATGCCAGTATCGGATCAACGTGTCACAGTCTAGGACTTCAGATGCCTAGCATCGTACGATTCGACCTGATGACCTGACGTTGCACTGATCGAGGTACCTATGCCAGTA";
    fn partial(expected: &[u8], seq: &[u8], minimum: usize) -> Option<usize> {
        PassThroughSignature::compile_with_partial(expected, 8, Some(minimum))?
            .detect_partial(seq, &mut Vec::new())
    }
    #[test]
    fn bounded_suffix_prefix_matching_lengths_masks_and_invalid_bases() {
        for length in [42, 64, 100, 128] {
            let adapter = &ADAPTER[..length];
            let disabled = PassThroughSignature::compile(adapter, 8).unwrap();
            for overlap in [8, 20, 30] {
                let seq = [b"GATTACA".as_slice(), &adapter[..overlap]].concat();
                assert_eq!(detect_pass_through(&seq, &disabled), None);
                assert_eq!(disabled.detect_partial(&seq, &mut Vec::new()), None);
                assert_eq!(partial(adapter, &seq, 20), (overlap >= 20).then_some(7));
            }
            let internal = [b"GATTACA".as_slice(), &adapter[..30], b"NNNN"].concat();
            assert_eq!(partial(adapter, &internal, 20), None);
            let fragment = [b"GATTACA".as_slice(), &adapter[5..35]].concat();
            assert_eq!(partial(adapter, &fragment, 20), None);
        }
        let mut adapter = ADAPTER[..100].to_vec();
        adapter[..10].fill(b'N');
        assert_eq!(partial(&adapter, &ADAPTER[..29], 20), None);
        assert_eq!(partial(&adapter, &ADAPTER[..30], 20), Some(0));
        for invalid in *b"N?-" {
            let mut read = ADAPTER[..30].to_vec();
            read[15] = invalid;
            assert_eq!(partial(ADAPTER, &read, 20), None);
            read[15] = ADAPTER[15];
            read[0] = invalid;
            assert_eq!(partial(&adapter, &read, 20), Some(0));
        }
        let mut mismatch = ADAPTER[..30].to_vec();
        mismatch[17] = b'A';
        assert_ne!(mismatch[17], ADAPTER[17]);
        assert_eq!(partial(ADAPTER, &mismatch, 20), None); // full mismatch allowance is irrelevant
        assert_eq!(partial(ADAPTER, &ADAPTER[..8], 8), Some(0));
        assert_eq!(partial(ADAPTER, &ADAPTER[..30], 31), None);
        assert_eq!(partial(ADAPTER, &ADAPTER[..30], 30), Some(0));
    }
    #[test]
    fn partial_does_not_depend_on_full_anchor_and_rejects_low_complexity() {
        let expected = [b"ACGTTGCACTNGATCGAGGTANACCTATGCCAN".as_slice(), ADAPTER].concat();
        let signature = PassThroughSignature::compile_with_partial(&expected, 1, Some(20)).unwrap();
        assert!(signature.anchor_offset >= 33);
        assert_eq!(
            signature.detect_partial(&expected[..33], &mut Vec::new()),
            Some(0)
        );
        for prefix in [
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".as_slice(),
            b"ACACACACACACACACACACACACACACAC",
            b"ACGTACGTACGTACGTACGTACGTACGTAC",
        ] {
            let expected = [prefix, ADAPTER].concat();
            assert_eq!(partial(&expected, prefix, 20), None);
        }
    }
    #[test]
    fn packed_partial_agrees_with_scalar_on_negatives_and_planted_matches() {
        let mut seed = 4849201u64;
        for i in 0..4000 {
            let length = 30 + i % 180;
            let seq: Vec<u8> = (0..length)
                .map(|_| {
                    seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGT"[(seed >> 32) as usize % 4]
                })
                .collect();
            assert_eq!(partial(ADAPTER, &seq, 20), None, "negative {i}");
        }
        let mut expected = ADAPTER.to_vec();
        for position in [1, 15, 31, 64, 99] {
            expected[position] = b'N';
        }
        let signature = PassThroughSignature::compile_with_partial(&expected, 0, Some(20)).unwrap();
        let mut scratch = Vec::new();
        for length in 1..expected.len() {
            let seq = [b"GGTTAG".as_slice(), &ADAPTER[..length]].concat();
            let reference = (1..=seq.len().min(expected.len() - 1))
                .filter(|&n| {
                    signature.partial.eligible[n - 1]
                        && seq[seq.len() - n..]
                            .iter()
                            .zip(&expected[..n])
                            .all(|(&a, &b)| b == b'N' || a == b)
                })
                .max_by_key(|&n| (signature.partial.known[n - 1], std::cmp::Reverse(n)))
                .map(|n| seq.len() - n);
            assert_eq!(
                signature.detect_partial(&seq, &mut scratch),
                reference,
                "length {length}"
            );
        }
    }
    #[test]
    fn long_full_signature_mismatches_have_no_two_base_ceiling() {
        let signature = PassThroughSignature::compile(ADAPTER, 5).unwrap();
        let mut observed = ADAPTER.to_vec();
        let positions: Vec<_> = (0..ADAPTER.len())
            .filter(|i| {
                !(signature.anchor_offset..signature.anchor_offset + signature.anchor_len)
                    .contains(i)
            })
            .collect();
        for (i, &position) in positions.iter().take(6).enumerate() {
            observed[position] = if observed[position] == b'A' {
                b'C'
            } else {
                b'A'
            };
            assert_eq!(
                detect_pass_through(&observed, &signature),
                if i < 5 { Some(0) } else { None }
            );
        }
        let signature = PassThroughSignature::compile(ADAPTER, usize::MAX).unwrap();
        assert_eq!(detect_pass_through(ADAPTER, &signature), Some(0));
    }
}

#[cfg(test)]
mod wide_budget_tests {
    use super::*;
    #[test]
    fn full_mismatch_counter_and_chunks_support_more_than_255_bases() {
        let expected =
            b"ACGTTGCACTGATCGAGGTACCTATGCCAGTATCGGATCAACGTGTCACAGTCTAGGACTTCAGAT".repeat(8);
        let signature = PassThroughSignature::compile(&expected, 256).unwrap();
        let mut observed = expected.clone();
        for i in (0..expected.len())
            .filter(|i| {
                !(signature.anchor_offset..signature.anchor_offset + signature.anchor_len)
                    .contains(i)
            })
            .take(256)
        {
            observed[i] = if observed[i] == b'A' { b'C' } else { b'A' };
        }
        assert_eq!(detect_pass_through(&observed, &signature), Some(0));
        let strict = PassThroughSignature::compile(&expected, 255).unwrap();
        assert_eq!(detect_pass_through(&observed, &strict), None);
    }
}

#[cfg(test)]
mod genomic_negative_tests {
    use super::*;
    #[test]
    fn partial_and_full_preserve_human_and_bacterial_genomic_windows_both_strands() {
        let adapter=b"ACGTTGCACTGATCGAGGTACCTATGCCAGTATCGGATCAACGTGTCACAGTCTAGGACTTCAGATGCCTAGCATCGTACGATTCGACCTGATGACCTGACGTTGCACTGATCGAGGTACCTATGCCAGTA";
        for fasta in [
            include_str!("../tests/fixtures/ecoli_U00096_3_1_10000.fa"),
            include_str!("../tests/fixtures/human_NC_000022_11_20000001_20010000.fa"),
        ] {
            let sequence: Vec<u8> = fasta.lines().skip(1).flat_map(str::bytes).collect();
            assert_eq!(sequence.len(), 10000);
            let mut opposite = Vec::new();
            append_oriented(&mut opposite, &sequence, true);
            for length in [42, 64, 100, 128] {
                let signature =
                    PassThroughSignature::compile_with_partial(&adapter[..length], 1, Some(20))
                        .unwrap();
                let mut scratch = Vec::new();
                for strand in [&sequence, &opposite] {
                    for start in (0..=strand.len() - 150).step_by(13) {
                        let read = &strand[start..start + 150];
                        assert_eq!(
                            detect_pass_through(read, &signature),
                            None,
                            "full {length}, {start}"
                        );
                        assert_eq!(
                            signature.detect_partial(read, &mut scratch),
                            None,
                            "partial {length}, {start}"
                        );
                    }
                }
            }
        }
    }
}
