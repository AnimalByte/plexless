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

impl FastqStats {
    pub fn update(&mut self, seq: &[u8], qual: &[u8]) -> Result<(), String> {
        let length = seq.len();

        if length != qual.len() {
            return Err("Sequence and quality lengths do not match".into());
        }

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

            let phred = quality
                .checked_sub(b'!')
                .ok_or("Invalid FASTQ quality score")?;

            self.quality_sum += u64::from(phred);

            if phred >= 20 {
                self.q20_bases += 1;
            }

            if phred >= 30 {
                self.q30_bases += 1;
            }
        }

        Ok(())
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

fn percent(count: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        count as f64 / total as f64 * 100.0
    }
}
