//! Assignment-independent distal boundary detection and startup compilation.
use crate::barcodes::BarcodeCatalog;
use crate::cli::DemuxArgs;
use crate::decoder::hamming_distance;
use crate::encoder::{EncodedBarcode, encode};
use crate::samples::{Sample, SampleSheet};
use crate::structure::{BarcodeExtractionPlans, Orientation, ReadLayout, ReadMate, ReadStructure};

/// Original-read coordinates carried through every output adapter.
#[derive(Debug, Clone, Copy)]
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
    mismatch_limit: u8,
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
    fn compile(expected: &[u8], mismatch_limit: u8) -> Option<Self> {
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
        let (anchor_offset, anchor_len) = best?;
        let anchor = encode(&expected[anchor_offset..anchor_offset + anchor_len])
            .ok()?
            .value;
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
        })
    }

    fn verify(&self, candidate: &[u8]) -> bool {
        let mut mismatches = 0u32;
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
            mismatches += u32::from(hamming_distance(
                value & chunk.known_mask,
                chunk.expected.value & chunk.known_mask,
            ));
            if mismatches > u32::from(self.mismatch_limit) {
                return false;
            }
        }
        true
    }
}

/// Finds the complete distal structure's start, relative to `seq`.
/// No assignment, orientation, parsing, allocation, or output decisions occur here.
pub fn detect_pass_through(seq: &[u8], signature: &PassThroughSignature) -> Option<usize> {
    let last_start = seq.len().checked_sub(signature.verification_len)?;
    let scan_start = signature.anchor_offset;
    let scan_end = last_start + signature.anchor_offset + signature.anchor_len;
    let mut window = 0u64;
    let mut valid = 0usize;
    for (index, &base) in seq[scan_start..scan_end].iter().enumerate() {
        let code = BASE[base as usize];
        if code > 3 {
            window = 0;
            valid = 0;
            continue;
        }
        window = ((window << 2) | u64::from(code)) & signature.anchor_mask;
        valid += 1;
        if valid >= signature.anchor_len && window == signature.anchor {
            let anchor_start = scan_start + index + 1 - signature.anchor_len;
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
    samples: Vec<[Option<PassThroughSignature>; 2]>,
}

impl PassThroughPlans {
    pub(crate) fn compile(
        layout: &ReadLayout,
        catalog: &BarcodeCatalog,
        samples: &SampleSheet,
        args: &DemuxArgs,
    ) -> Result<Self, String> {
        if args.no_pass_through {
            return Ok(Self::default());
        }
        if args.pass_through_mismatches > 2 {
            return Err("Pass-through mismatches must be between 0 and 2".into());
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
        let plans = layout.compile_extraction_plans()?;
        let mut compiled = Vec::with_capacity(samples.samples.len());
        let mut unavailable = [0usize; 2];
        for sample in &samples.samples {
            let mut signatures = [None, None];
            for mate in 0..if paired { 2 } else { 1 } {
                let expected = if let Some(structure) = &overrides[mate] {
                    compile_override(structure, &plans, catalog, sample)?
                } else if let Some(source) = if mate == 0 { r2 } else { r1 } {
                    compile_opposite(source, &plans, catalog, sample)?
                } else {
                    Vec::new()
                };
                signatures[mate] =
                    PassThroughSignature::compile(&expected, args.pass_through_mismatches);
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
                    "Pass-through: R{} unavailable for {count}/{} samples: no sufficiently varied contiguous known anchor of at least 12 bases",
                    mate + 1,
                    samples.samples.len()
                );
            }
        }
        Ok(Self { samples: compiled })
    }

    pub(crate) fn signature(&self, sample: u32, mate: usize) -> Option<&PassThroughSignature> {
        self.samples.get(sample as usize)?.get(mate)?.as_ref()
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

fn compile_opposite(
    source: &ReadStructure,
    plans: &BarcodeExtractionPlans,
    catalog: &BarcodeCatalog,
    sample: &Sample,
) -> Result<Vec<u8>, String> {
    let mut output = Vec::new();
    output
        .try_reserve_exact(source.prefix_len)
        .map_err(|_| "Distal structure is too large")?;
    for segment in source.segments.iter().rev() {
        let length = segment.end - segment.start;
        if segment.symbol == b'T' {
            output.resize(output.len() + length, b'N');
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
) -> Result<Vec<u8>, String> {
    let mut offsets = [0usize; 26];
    let mut output = Vec::new();
    output
        .try_reserve_exact(structure.prefix_len)
        .map_err(|_| "Distal override is too large")?;
    for segment in &structure.segments {
        let length = segment.end - segment.start;
        if segment.symbol == b'T' {
            output.resize(output.len() + length, b'N');
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
