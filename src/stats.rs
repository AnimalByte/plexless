use crate::structure::{Orientation, ReadLayout, ReadStructure, is_barcode_symbol};

#[derive(Debug, Default)]
pub struct FastqStats {
    pub reads: u64,
    pub bases: u64,

    pub min_length: Option<usize>,
    pub max_length: usize,

    pub gc_bases: u64,
    pub n_bases: u64,

    pub quality_sum: u64,
    pub q20_bases: u64,
    pub q30_bases: u64,
}

#[derive(Debug)]
pub(crate) struct BarcodeSegmentStats {
    pub(crate) symbol: u8,
    pub(crate) piece: usize,
    pub(crate) start: usize,
    pub(crate) end: usize,
    pub(crate) orientation: Orientation,
    pub(crate) stats: FastqStats,
}

/// Statistics projected onto the biological suffix and physical barcode
/// segments of one mate. Technical (`T`) segments are intentionally omitted.
#[derive(Debug)]
pub(crate) struct MateQcStats {
    prefix_len: usize,
    biological: FastqStats,
    barcode_segments: Vec<BarcodeSegmentStats>,
}

impl MateQcStats {
    fn new(structure: Option<&ReadStructure>) -> Self {
        let prefix_len = structure.map_or(0, |structure| structure.prefix_len);
        let mut pieces_per_symbol = [0usize; 26];
        let barcode_segments = structure
            .into_iter()
            .flat_map(|structure| &structure.segments)
            .filter(|segment| is_barcode_symbol(segment.symbol))
            .map(|segment| {
                let symbol_index = usize::from(segment.symbol - b'A');
                pieces_per_symbol[symbol_index] += 1;
                BarcodeSegmentStats {
                    symbol: segment.symbol,
                    piece: pieces_per_symbol[symbol_index],
                    start: segment.start,
                    end: segment.end,
                    orientation: segment.orientation,
                    stats: FastqStats::default(),
                }
            })
            .collect();

        Self {
            prefix_len,
            biological: FastqStats::default(),
            barcode_segments,
        }
    }

    pub(crate) fn for_layout(layout: &ReadLayout) -> (Self, Option<Self>) {
        match layout {
            ReadLayout::Single { r1 } => (Self::new(r1.as_ref()), None),
            ReadLayout::Paired { r1, r2 } => (Self::new(r1.as_ref()), Some(Self::new(r2.as_ref()))),
        }
    }

    pub(crate) fn update(&mut self, seq: &[u8], qual: &[u8]) -> Result<(), String> {
        if seq.len() != qual.len() {
            return Err("Sequence and quality lengths do not match".into());
        }
        validate_quality(qual)?;

        // A short structured read has no biological suffix. Clipping barcode
        // coordinates still records the observed bases and makes truncation
        // visible through the segment length statistics.
        let biological_start = self.prefix_len.min(seq.len());
        self.biological
            .update_validated(&seq[biological_start..], &qual[biological_start..]);

        for segment in &mut self.barcode_segments {
            let start = segment.start.min(seq.len());
            let end = segment.end.min(seq.len());
            segment
                .stats
                .update_validated(&seq[start..end], &qual[start..end]);
        }

        Ok(())
    }

    pub(crate) fn biological(&self) -> &FastqStats {
        &self.biological
    }

    pub(crate) fn barcode_segments(&self) -> &[BarcodeSegmentStats] {
        &self.barcode_segments
    }
}

impl FastqStats {
    pub fn update(&mut self, seq: &[u8], qual: &[u8]) -> Result<(), String> {
        if seq.len() != qual.len() {
            return Err("Sequence and quality lengths do not match".into());
        }
        validate_quality(qual)?;
        self.update_validated(seq, qual);
        Ok(())
    }

    fn update_validated(&mut self, seq: &[u8], qual: &[u8]) {
        let length = seq.len();

        // Cheap values: no sequence scan required.
        self.reads += 1;
        self.bases += length as u64;

        self.min_length = Some(self.min_length.map_or(length, |min| min.min(length)));

        self.max_length = self.max_length.max(length);

        // Only one full scan.
        for (&base, &quality) in seq.iter().zip(qual) {
            match base {
                b'G' | b'g' | b'C' | b'c' => self.gc_bases += 1,
                b'N' | b'n' => self.n_bases += 1,
                _ => {}
            }

            let phred = quality - b'!';

            self.quality_sum += u64::from(phred);

            if phred >= 20 {
                self.q20_bases += 1;
            }

            if phred >= 30 {
                self.q30_bases += 1;
            }
        }
    }

    pub fn mean_length(&self) -> f64 {
        if self.reads == 0 {
            0.0
        } else {
            self.bases as f64 / self.reads as f64
        }
    }

    pub fn gc_percent(&self) -> f64 {
        percent(self.gc_bases, self.bases)
    }

    pub fn n_percent(&self) -> f64 {
        percent(self.n_bases, self.bases)
    }

    pub fn mean_quality(&self) -> f64 {
        if self.bases == 0 {
            0.0
        } else {
            self.quality_sum as f64 / self.bases as f64
        }
    }

    pub fn q20_percent(&self) -> f64 {
        percent(self.q20_bases, self.bases)
    }

    pub fn q30_percent(&self) -> f64 {
        percent(self.q30_bases, self.bases)
    }
}

fn validate_quality(qual: &[u8]) -> Result<(), String> {
    if qual.iter().any(|quality| *quality < b'!') {
        Err("Invalid FASTQ quality score".into())
    } else {
        Ok(())
    }
}

fn percent(count: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        count as f64 / total as f64 * 100.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn projects_biological_and_physical_barcode_regions() {
        let layout = ReadLayout::paired(Some("R1_2A1T2A(rc)"), None).unwrap();
        let (mut r1, r2) = MateQcStats::for_layout(&layout);
        let mut r2 = r2.unwrap();

        r1.update(b"ACGTTGAT", b"II!55III").unwrap();
        r2.update(b"GATTACA", b"IIIIIII").unwrap();

        assert_eq!(r1.biological().reads, 1);
        assert_eq!(r1.biological().bases, 3);
        assert_eq!(r1.barcode_segments().len(), 2);
        assert_eq!(r1.barcode_segments()[0].piece, 1);
        assert_eq!(r1.barcode_segments()[0].stats.bases, 2);
        assert_eq!(r1.barcode_segments()[1].piece, 2);
        assert_eq!(
            r1.barcode_segments()[1].orientation,
            Orientation::ReverseComplement
        );
        assert_eq!(r1.barcode_segments()[1].stats.bases, 2);
        assert_eq!(r2.biological().bases, 7);
        assert!(r2.barcode_segments().is_empty());
    }

    #[test]
    fn clips_short_structured_reads_without_panicking() {
        let layout = ReadLayout::single(Some("R1_2A2B2T")).unwrap();
        let (mut r1, r2) = MateQcStats::for_layout(&layout);

        r1.update(b"AC", b"II").unwrap();

        assert!(r2.is_none());
        assert_eq!(r1.biological().reads, 1);
        assert_eq!(r1.biological().bases, 0);
        assert_eq!(r1.biological().min_length, Some(0));
        assert_eq!(r1.barcode_segments()[0].stats.reads, 1);
        assert_eq!(r1.barcode_segments()[0].stats.bases, 2);
        assert_eq!(r1.barcode_segments()[0].stats.min_length, Some(2));
        assert_eq!(r1.barcode_segments()[1].stats.reads, 1);
        assert_eq!(r1.barcode_segments()[1].stats.bases, 0);
        assert_eq!(r1.barcode_segments()[1].stats.min_length, Some(0));
    }
}
