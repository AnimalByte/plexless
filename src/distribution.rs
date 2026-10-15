//! Detailed per-cycle QC used only when JSON/HTML reporting is requested.
use serde_json::{Value, json};
use std::collections::BTreeMap;
#[derive(Debug, Default, Clone)]
struct Cycle {
    bases: [u64; 6],
    qualities: [u64; 6],
}
#[derive(Debug, Default)]
pub(crate) struct Distribution {
    cycles: Vec<Cycle>,
    lengths: BTreeMap<usize, u64>,
    quality_histogram: Vec<u64>,
}
fn index(base: u8) -> usize {
    match base {
        b'A' | b'a' => 0,
        b'T' | b't' => 1,
        b'C' | b'c' => 2,
        b'G' | b'g' => 3,
        b'N' | b'n' => 4,
        _ => 5,
    }
}
impl Distribution {
    pub(crate) fn update(&mut self, seq: &[u8], qual: &[u8]) {
        self.cycles
            .resize_with(self.cycles.len().max(seq.len()), Cycle::default);
        if self.quality_histogram.is_empty() {
            self.quality_histogram.resize(94, 0);
        }
        *self.lengths.entry(seq.len()).or_default() += 1;
        for (i, (&base, &quality)) in seq.iter().zip(qual).enumerate() {
            let b = index(base);
            let q = quality.saturating_sub(b'!').min(93);
            self.cycles[i].bases[b] += 1;
            self.cycles[i].qualities[b] += u64::from(q);
            self.quality_histogram[q as usize] += 1;
        }
    }
    pub(crate) fn merge(&mut self, other: &Self) {
        self.cycles
            .resize_with(self.cycles.len().max(other.cycles.len()), Cycle::default);
        for (a, b) in self.cycles.iter_mut().zip(&other.cycles) {
            for i in 0..6 {
                a.bases[i] += b.bases[i];
                a.qualities[i] += b.qualities[i];
            }
        }
        for (&length, &count) in &other.lengths {
            *self.lengths.entry(length).or_default() += count;
        }
        if !other.quality_histogram.is_empty() {
            self.quality_histogram.resize(94, 0);
            for (a, &b) in self
                .quality_histogram
                .iter_mut()
                .zip(&other.quality_histogram)
            {
                *a += b;
            }
        }
    }
    pub(crate) fn json(&self) -> Value {
        let mut quality = serde_json::Map::new();
        let mut content = serde_json::Map::new();
        for (b, name) in ["A", "T", "C", "G", "N", "other"].into_iter().enumerate() {
            quality.insert(
                name.into(),
                json!(
                    self.cycles
                        .iter()
                        .map(|c| ratio(c.qualities[b], c.bases[b]))
                        .collect::<Vec<_>>()
                ),
            );
            content.insert(
                name.into(),
                json!(
                    self.cycles
                        .iter()
                        .map(|c| ratio(c.bases[b], c.bases.iter().sum()))
                        .collect::<Vec<_>>()
                ),
            );
        }
        quality.insert(
            "mean".into(),
            json!(
                self.cycles
                    .iter()
                    .map(|c| ratio(c.qualities.iter().sum(), c.bases.iter().sum()))
                    .collect::<Vec<_>>()
            ),
        );
        json!({"cycles":self.cycles.len(),"quality_curves":quality,"content_curves":content,"quality_histogram":self.quality_histogram,"length_distribution":self.lengths.iter().map(|(&length,&reads)|json!({"length":length,"reads":reads})).collect::<Vec<_>>()})
    }
}
fn ratio(a: u64, b: u64) -> f64 {
    if b == 0 { 0.0 } else { a as f64 / b as f64 }
}
