//! Quality-window/global and polyG/polyX trimming adapted from fastp.
//! See THIRD_PARTY.md. All coordinates refer to the original input read.
use crate::fastp_options::FastpOptions;
use crate::pass_through::BiologicalBounds;

pub(crate) fn global_and_quality(
    seq: &[u8],
    qual: &[u8],
    mut bounds: BiologicalBounds,
    mate: usize,
    o: &FastpOptions,
) -> BiologicalBounds {
    let front = if mate == 0 {
        o.trim_front1
    } else {
        o.trim_front2.unwrap_or(o.trim_front1)
    };
    let tail = if mate == 0 {
        o.trim_tail1
    } else {
        o.trim_tail2.unwrap_or(o.trim_tail1)
    };
    bounds.start = bounds.start.saturating_add(front).min(bounds.end);
    bounds.end = bounds.end.saturating_sub(tail).max(bounds.start);
    if o.cut_front {
        let w = usize::from(o.cut_front_window_size.unwrap_or(o.cut_window_size));
        let q = o.cut_front_mean_quality.unwrap_or(o.cut_mean_quality);
        let mut start = bounds.start;
        let mut sum = window_sum(qual, start, w, bounds.end);
        while start + w <= bounds.end {
            if sum >= w * u8::saturating_add(q, 33) as usize {
                break;
            }
            if start + w == bounds.end {
                start = bounds.end;
                break;
            }
            sum = sum - usize::from(qual[start]) + usize::from(qual[start + w]);
            start += 1;
        }
        if start > bounds.start && start < bounds.end {
            start = (start + w - 1).min(bounds.end);
        }
        while start < bounds.end && seq[start].eq_ignore_ascii_case(&b'N') {
            start += 1;
        }
        bounds.start = start.min(bounds.end);
    }
    if o.cut_right {
        let w = usize::from(o.cut_right_window_size.unwrap_or(o.cut_window_size));
        let q = o.cut_right_mean_quality.unwrap_or(o.cut_mean_quality);
        let mut start = bounds.start;
        let mut sum = window_sum(qual, start, w, bounds.end);
        while start + w <= bounds.end {
            if sum < w * usize::from(q + 33) {
                while start < bounds.end && qual[start] >= q + 33 {
                    start += 1;
                }
                bounds.end = start;
                break;
            }
            if start + w == bounds.end {
                break;
            }
            sum = sum - usize::from(qual[start]) + usize::from(qual[start + w]);
            start += 1;
        }
    } else if o.cut_tail {
        let w = usize::from(o.cut_tail_window_size.unwrap_or(o.cut_window_size));
        let q = o.cut_tail_mean_quality.unwrap_or(o.cut_mean_quality);
        let mut end = bounds.end;
        let mut sum = if end - bounds.start >= w {
            qual[end - w..end].iter().map(|&q| usize::from(q)).sum()
        } else {
            0
        };
        while end - bounds.start >= w {
            if sum >= w * usize::from(q + 33) {
                break;
            }
            if end - bounds.start == w {
                end = bounds.start;
                break;
            }
            sum = sum - usize::from(qual[end - 1]) + usize::from(qual[end - w - 1]);
            end -= 1;
        }
        if end < bounds.end && end > bounds.start {
            end = end.saturating_sub(w - 1).max(bounds.start);
        }
        while end > bounds.start && seq[end - 1].eq_ignore_ascii_case(&b'N') {
            end -= 1;
        }
        bounds.end = end;
    }
    bounds
}
fn window_sum(qual: &[u8], start: usize, w: usize, end: usize) -> usize {
    if end - start < w {
        0
    } else {
        qual[start..start + w].iter().map(|&q| usize::from(q)).sum()
    }
}

pub(crate) fn poly_g(seq: &[u8], mut bounds: BiologicalBounds, minimum: usize) -> BiologicalBounds {
    let mut mismatches = 0;
    let mut seen = 0;
    let mut first_g = bounds.end;
    for i in 0..bounds.len() {
        let position = bounds.end - i - 1;
        if seq[position].eq_ignore_ascii_case(&b'G') {
            first_g = position;
        } else {
            mismatches += 1;
        }
        let allowed = (i + 1) / 8;
        if mismatches > 5 || (mismatches > allowed && i + 1 >= minimum) {
            break;
        }
        seen = i + 1;
    }
    if seen >= minimum {
        bounds.end = first_g;
    }
    bounds
}
pub(crate) fn poly_x(seq: &[u8], mut bounds: BiologicalBounds, minimum: usize) -> BiologicalBounds {
    let mut counts = [0usize; 4];
    let mut best = None;
    for width in 1..=bounds.len() {
        let base = seq[bounds.end - width].to_ascii_uppercase();
        for (index, expected) in b"ATCG".iter().enumerate() {
            if base == *expected || base == b'N' {
                counts[index] += 1;
            }
        }
        let allowed = 5.min(width / 8);
        let maximum = *counts.iter().max().unwrap();
        if width - maximum > allowed && width >= 8.min(minimum.saturating_sub(1)) {
            break;
        }
        if width >= minimum {
            best = Some((width, counts.iter().position(|&n| n == maximum).unwrap()));
        }
    }
    if let Some((mut width, index)) = best {
        while width > 0 && !seq[bounds.end - width].eq_ignore_ascii_case(&b"ATCG"[index]) {
            width -= 1;
        }
        bounds.end -= width;
    }
    bounds
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn tail_right_front_and_global_keep_sequence_quality_coordinates() {
        let s = b"ACGTTGCACTGATCGAGGTACCTA";
        let mut q = vec![b'I'; s.len()];
        q[..6].fill(b'!');
        q[s.len() - 6..].fill(b'!');
        let o = FastpOptions {
            cut_front: true,
            cut_tail: true,
            ..Default::default()
        };
        let b = global_and_quality(s, &q, BiologicalBounds::full(s.len()), 0, &o);
        assert!(b.start >= 6 && b.end <= s.len() - 6);
        assert!(b.start <= b.end);
        let b = global_and_quality(
            s,
            &q,
            BiologicalBounds::full(s.len()),
            0,
            &FastpOptions {
                trim_front1: 100,
                trim_tail1: 100,
                ..Default::default()
            },
        );
        assert_eq!(b.len(), 0);
    }
    #[test]
    fn polyx_reference_and_empty_edges() {
        let s = b"ATTTTAAAAAAAAAATAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAT";
        assert_eq!(
            &s[poly_x(s, BiologicalBounds::full(s.len()), 10).range()],
            b"ATTTT"
        );
        assert_eq!(poly_x(b"", BiologicalBounds::full(0), 10).len(), 0);
        assert_eq!(
            poly_g(b"ACCTGGGGGGGGGG", BiologicalBounds::full(14), 10).len(),
            4
        );
    }
}
