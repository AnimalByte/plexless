//! Ungapped paired-read overlap analysis adapted from OpenGene/fastp.
//! See THIRD_PARTY.md and licenses/fastp-MIT.txt. Unlike the reference protected
//! prefix shortcut, every position must satisfy the global mismatch policy.
use crate::decoder::hamming_distance;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Overlap {
    pub insert_length: usize,
    pub length: usize,
    pub mismatches: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum OverlapResult {
    Unique(Overlap),
    Ambiguous,
    None,
}

#[derive(Debug, Default)]
struct PackedRead {
    words: Vec<u64>,
    invalid: Vec<u32>,
    len: usize,
}

fn code(base: u8) -> u8 {
    match base {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 255,
    }
}
impl PackedRead {
    fn load(&mut self, sequence: &[u8], reverse: bool) {
        self.len = sequence.len();
        let words = sequence.len().div_ceil(32);
        self.words.resize(words, 0);
        self.words.fill(0);
        self.invalid.resize(words, 0);
        self.invalid.fill(0);
        for i in 0..sequence.len() {
            let mut base = code(sequence[if reverse { sequence.len() - 1 - i } else { i }]);
            if base > 3 {
                self.invalid[i / 32] |= 1 << (i % 32);
            } else if reverse {
                base ^= 3;
            }
            self.words[i / 32] |= u64::from(base & 3) << (2 * (i % 32));
        }
    }
    fn window(&self, start: usize, length: usize) -> (u64, u32) {
        let index = start / 32;
        let shift = start % 32;
        let mut bases = self.words[index] >> (shift * 2);
        let mut invalid = self.invalid[index] >> shift;
        if shift != 0 && shift + length > 32 {
            bases |= self.words[index + 1] << (64 - shift * 2);
            invalid |= self.invalid[index + 1] << (32 - shift);
        }
        let mask = if length == 32 {
            u64::MAX
        } else {
            (1u64 << (2 * length)) - 1
        };
        let nmask = if length == 32 {
            u32::MAX
        } else {
            (1u32 << length) - 1
        };
        (bases & mask, invalid & nmask)
    }
}

#[derive(Debug, Default)]
pub(crate) struct OverlapScratch {
    r1: PackedRead,
    r2: PackedRead,
    reverse: Vec<u8>,
    prefix: Vec<usize>,
}

/// Search legal relative offsets of R1 and reverse-complement R2. A unique
/// alignment is required; multiple qualifying lengths preserve the pair.
/// Scratch is reused across records in a batch; memory is O(read length).
pub(crate) fn analyze(
    r1: &[u8],
    r2: &[u8],
    minimum: usize,
    limit: usize,
    percent: f64,
    scratch: &mut OverlapScratch,
) -> OverlapResult {
    if minimum == 0 || r1.len() < minimum || r2.len() < minimum {
        return OverlapResult::None;
    }
    let Some(total) = r1.len().checked_add(r2.len()) else {
        return OverlapResult::None;
    };
    scratch.r1.load(r1, false);
    scratch.r2.load(r2, true);
    let mut result = None;
    // Equivalent to scanning positive and negative offsets, expressed as an
    // unsigned insert length so unequal mates and overflow are handled safely.
    for insert in minimum..=total - minimum {
        let start1 = insert.saturating_sub(r2.len());
        let start2 = r2.len().saturating_sub(insert);
        let length = (r1.len() - start1).min(r2.len() - start2);
        if length < minimum {
            continue;
        }
        let budget = limit.min((length as f64 * percent / 100.0).floor() as usize);
        let mut differences = 0usize;
        let mut accepted = true;
        for offset in (0..length).step_by(32) {
            let width = (length - offset).min(32);
            let (a, na) = scratch.r1.window(start1 + offset, width);
            let (b, nb) = scratch.r2.window(start2 + offset, width);
            if na | nb != 0 {
                accepted = false;
                break;
            }
            differences += usize::from(hamming_distance(a, b));
            if differences > budget {
                accepted = false;
                break;
            }
        }
        if !accepted || !informative(&r1[start1..start1 + length]) {
            continue;
        }
        if result.is_some() {
            return OverlapResult::Ambiguous;
        }
        result = Some(Overlap {
            insert_length: insert,
            length,
            mismatches: differences,
        });
    }
    result.map_or(OverlapResult::None, OverlapResult::Unique)
}

fn informative(sequence: &[u8]) -> bool {
    let mut counts = [0usize; 4];
    let mut pairs = 0u16;
    let mut previous = None;
    for &base in sequence {
        let base = code(base);
        if base > 3 {
            return false;
        }
        counts[base as usize] += 1;
        if let Some(prior) = previous {
            pairs |= 1u16 << (prior * 4 + base);
        }
        previous = Some(base);
    }
    counts.iter().filter(|&&n| n > 0).count() >= 3
        && counts.iter().all(|&n| n * 4 <= sequence.len() * 3)
        && pairs.count_ones() >= 6
}

/// One-gap, short-insert fallback. Only adapter boundaries are returned: a gap
/// does not authorize base correction or merging. Compare the complete inferred
/// shared insert; gaps near either edge and non-unique boundaries are rejected.
pub(crate) fn analyze_gap(
    reads: [&[u8]; 2],
    fronts: [usize; 2],
    minimum: usize,
    limit: usize,
    percent: f64,
    scratch: &mut OverlapScratch,
) -> Result<Option<[usize; 2]>, ()> {
    let maximum = reads[0]
        .len()
        .saturating_sub(fronts[1])
        .min(reads[1].len().saturating_sub(fronts[0]));
    if minimum < 16 || maximum < minimum {
        return Ok(None);
    }
    let mut result = None;
    for length in minimum..=maximum {
        for longer in 0..2 {
            let ends = if longer == 0 {
                [length + 1 + fronts[1], length + fronts[0]]
            } else {
                [length + fronts[1], length + 1 + fronts[0]]
            };
            if ends[0] > reads[0].len() || ends[1] > reads[1].len() {
                continue;
            }
            let a = &reads[0][..ends[0] - fronts[1]];
            let b = &reads[1][..ends[1] - fronts[0]];
            if !informative(a) || !informative(b) {
                continue;
            }
            scratch.reverse.clear();
            scratch
                .reverse
                .extend(b.iter().rev().map(|&b| b"TGCA"[usize::from(code(b))]));
            let (short, long) = if longer == 0 {
                (scratch.reverse.as_slice(), a)
            } else {
                (a, scratch.reverse.as_slice())
            };
            scratch.prefix.clear();
            scratch.prefix.push(0);
            for i in 0..length {
                scratch.prefix.push(
                    scratch.prefix[i] + usize::from(!short[i].eq_ignore_ascii_case(&long[i])),
                );
            }
            let budget = limit.min((length as f64 * percent / 100.0).floor() as usize);
            let mut suffix = 0;
            let mut valid = false;
            for gap in (0..=length).rev() {
                if gap >= 8 && length - gap >= 8 && scratch.prefix[gap] + suffix <= budget {
                    valid = true;
                    break;
                }
                if gap > 0 {
                    suffix += usize::from(!short[gap - 1].eq_ignore_ascii_case(&long[gap]));
                }
            }
            if valid {
                if result.is_some_and(|prior| prior != ends) {
                    return Err(());
                }
                result = Some(ends);
            }
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn rc(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => b'N',
            })
            .collect()
    }
    fn dna(length: usize) -> Vec<u8> {
        let mut x = 87312u64;
        (0..length)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(x >> 32) as usize % 4]
            })
            .collect()
    }
    fn analyze_default(a: &[u8], b: &[u8]) -> OverlapResult {
        analyze(a, b, 30, 5, 20.0, &mut OverlapScratch::default())
    }
    #[test]
    fn fastp_reference_vector() {
        let a=b"CAGCGCCTACGGGCCCCTTTTTCTGCGCGACCGCGTGGCTGTGGGCGCGGATGCCTTTGAGCGCGGTGACTTCTCACTGCGTATCGAGC";
        let b=b"ACCTCCAGCGGCTCGATACGCAGTGAGAAGTCACCGCGCTCAAAGGCATCCGCGCCCACAGCCACGCGGTCGCGCAGAAAAAGGGGTCC";
        assert_eq!(
            analyze(a, b, 30, 2, 20.0, &mut OverlapScratch::default()),
            OverlapResult::Unique(Overlap {
                insert_length: b.len() + 10,
                length: 79,
                mismatches: 1
            })
        );
    }
    #[test]
    fn long_insert_overlap_does_not_imply_adapter_and_short_insert_is_inferred() {
        let insert = dna(220);
        assert_eq!(
            analyze_default(&insert[..150], &rc(&insert[70..])),
            OverlapResult::Unique(Overlap {
                insert_length: 220,
                length: 80,
                mismatches: 0
            })
        );
        let insert = dna(120);
        let a = [insert.as_slice(), b"ACCTGAGTTACGAGCTGCTAGGCTTCAACT"].concat();
        let b = [rc(&insert).as_slice(), b"TGCTACCTGAGCTTGATCGACTGACGTAGG"].concat();
        assert_eq!(
            analyze_default(&a, &b),
            OverlapResult::Unique(Overlap {
                insert_length: 120,
                length: 120,
                mismatches: 0
            })
        );
        assert_eq!(
            analyze_default(&a, &rc(&insert[20..])),
            OverlapResult::Unique(Overlap {
                insert_length: 120,
                length: 100,
                mismatches: 0
            })
        );
    }
    #[test]
    fn exact_minimum_ambiguous_low_complexity_ns_and_late_mismatches() {
        let a = dna(80);
        let mut b = rc(&a);
        assert_eq!(
            analyze_default(&a[..30], &rc(&a[..30])),
            OverlapResult::Unique(Overlap {
                insert_length: 30,
                length: 30,
                mismatches: 0
            })
        );
        for base in b.iter_mut().take(20) {
            *base = if *base == b'A' { b'C' } else { b'A' };
        }
        assert_eq!(analyze_default(&a, &b), OverlapResult::None);
        assert_eq!(
            analyze_default(&[b'A'; 150], &[b'T'; 150]),
            OverlapResult::None
        );
        let repeat = b"ACGTTGCACTGATCGAGGTACCTA".repeat(8);
        assert_eq!(
            analyze_default(&repeat, &rc(&repeat)),
            OverlapResult::Ambiguous
        );
        let mut n = a.clone();
        n[40] = b'N';
        assert_eq!(analyze_default(&n, &rc(&a)), OverlapResult::None);
    }
    #[test]
    fn packed_windows_equal_scalar_across_chunk_boundaries() {
        let sequence = dna(200);
        let mut packed = PackedRead::default();
        packed.load(&sequence, false);
        for start in 0..sequence.len() {
            for len in 1..=32.min(sequence.len() - start) {
                let (word, invalid) = packed.window(start, len);
                assert_eq!(invalid, 0);
                for i in 0..len {
                    assert_eq!((word >> (2 * i)) & 3, u64::from(code(sequence[start + i])));
                }
            }
        }
    }
}
