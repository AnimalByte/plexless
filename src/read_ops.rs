//! Rust implementations of fastp-style overlap correction and merging.
use crate::overlap::Overlap;
use crate::parallel::OwnedFastqRecord;
use crate::pass_through::BiologicalBounds;
fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => b'N',
    }
}
pub(crate) fn correct(
    a: &mut OwnedFastqRecord,
    b: &mut OwnedFastqRecord,
    bounds: [BiologicalBounds; 2],
    overlap: Overlap,
) -> u64 {
    let start = overlap.insert_length.saturating_sub(bounds[1].len());
    let mut corrected = 0;
    for i in start..bounds[0].len().min(overlap.insert_length) {
        let p1 = bounds[0].start + i;
        let p2 = bounds[1].start + overlap.insert_length - 1 - i;
        if a.seq[p1].to_ascii_uppercase() == complement(b.seq[p2]) {
            continue;
        }
        if a.qual[p1] >= b'!' + 30 && b.qual[p2] <= b'!' + 14 {
            b.seq[p2] = complement(a.seq[p1]);
            b.qual[p2] = a.qual[p1];
            corrected += 1;
        } else if b.qual[p2] >= b'!' + 30 && a.qual[p1] <= b'!' + 14 {
            a.seq[p1] = complement(b.seq[p2]);
            a.qual[p1] = b.qual[p2];
            corrected += 1;
        }
    }
    corrected
}
pub(crate) fn merge(
    a: &OwnedFastqRecord,
    b: &OwnedFastqRecord,
    bounds: [BiologicalBounds; 2],
    overlap: Overlap,
    fronts: [usize; 2],
) -> OwnedFastqRecord {
    // The overlap defines relative coordinates, not an instruction to discard
    // biological overhangs exposed by front clipping. Preserve their union when
    // the other mate covers them, bounded by the inferred physical insert ends.
    let left = fronts[0].min(bounds[1].len().saturating_sub(overlap.insert_length));
    let right = overlap.insert_length.max(
        bounds[0]
            .len()
            .min(overlap.insert_length.saturating_add(fronts[1])),
    );
    let length = left + right;
    let a_length = bounds[0]
        .len()
        .min(overlap.insert_length.saturating_add(fronts[1]));
    let b_length = bounds[1]
        .len()
        .min(overlap.insert_length.saturating_add(fronts[0]));
    let mut seq = Vec::with_capacity(length);
    let mut qual = Vec::with_capacity(length);
    for coordinate in 0..length {
        let a_index = coordinate
            .checked_sub(left)
            .filter(|&i| i < a_length)
            .map(|i| bounds[0].start + i);
        let b_relative = if coordinate < left {
            Some(overlap.insert_length + left - coordinate - 1)
        } else {
            overlap.insert_length.checked_sub(coordinate - left + 1)
        };
        let b_index = b_relative
            .filter(|&j| j < b_length)
            .map(|j| bounds[1].start + j);
        let (base, q) = match (a_index, b_index) {
            (Some(p1), Some(p2)) if b.qual[p2] > a.qual[p1] => (complement(b.seq[p2]), b.qual[p2]),
            (Some(p1), _) => (a.seq[p1], a.qual[p1]),
            (None, Some(p2)) => (complement(b.seq[p2]), b.qual[p2]),
            (None, None) => unreachable!("accepted overlap has no uncovered union positions"),
        };
        seq.push(base);
        qual.push(q);
    }
    let split =
        a.id.iter()
            .position(u8::is_ascii_whitespace)
            .unwrap_or(a.id.len());
    let mut id = a.id.clone();
    if id[..split].ends_with(b"/1") {
        id.drain(split - 2..split);
    }
    id.extend_from_slice(b" merged");
    OwnedFastqRecord { id, seq, qual }
}
