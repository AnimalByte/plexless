//! Bounded, sampled Space-Saving heavy hitters. Counts include substring
//! occurrences (not distinct reads); replacement error bounds are reported.
use serde_json::{Value, json};
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};
const CAPACITY: usize = 4096;
#[derive(Debug)]
pub(crate) struct Overrepresentation {
    sampling: u32,
    reads: u64,
    sampled: u64,
    windows: u64,
    candidates: HashMap<Vec<u8>, (u64, u64)>,
    minimum: BinaryHeap<Reverse<(u64, Vec<u8>)>>,
}
impl Overrepresentation {
    pub(crate) fn new(sampling: u32) -> Self {
        Self {
            sampling,
            reads: 0,
            sampled: 0,
            windows: 0,
            candidates: HashMap::new(),
            minimum: BinaryHeap::new(),
        }
    }
    fn add(&mut self, sequence: &[u8], weight: u64, error: u64) {
        if let Some(value) = self.candidates.get_mut(sequence) {
            value.0 += weight;
            value.1 += error;
            self.minimum.push(Reverse((value.0, sequence.to_vec())));
        } else {
            let mut offset = 0;
            if self.candidates.len() == CAPACITY {
                while let Some(Reverse((count, key))) = self.minimum.pop() {
                    if self.candidates.get(&key).is_some_and(|v| v.0 == count) {
                        self.candidates.remove(&key);
                        offset = count;
                        break;
                    }
                }
            }
            self.candidates
                .insert(sequence.to_vec(), (offset + weight, offset + error));
            self.minimum
                .push(Reverse((offset + weight, sequence.to_vec())));
        }
        if self.minimum.len() > CAPACITY * 4 {
            self.minimum = self
                .candidates
                .iter()
                .map(|(key, value)| Reverse((value.0, key.clone())))
                .collect();
        }
    }
    pub(crate) fn update(&mut self, sequence: &[u8]) {
        self.reads += 1;
        if !(self.reads - 1).is_multiple_of(u64::from(self.sampling)) {
            return;
        }
        self.sampled += 1;
        let upper = sequence
            .iter()
            .any(u8::is_ascii_lowercase)
            .then(|| sequence.to_ascii_uppercase());
        let sequence = upper.as_deref().unwrap_or(sequence);
        // Bound full-read candidate storage for long-read input. All substring
        // windows are counted, preserving position-independent discovery.
        for length in [10, 20, 40, 100, 150] {
            for window in sequence.windows(length) {
                if window.iter().all(u8::is_ascii_uppercase)
                    && window.iter().all(|b| b"ACGT".contains(b))
                {
                    self.windows += 1;
                    self.add(window, 1, 0);
                }
            }
        }
        if sequence.len() <= 512
            && ![10, 20, 40, 100, 150].contains(&sequence.len())
            && sequence.iter().all(|b| b"ACGT".contains(b))
        {
            self.windows += 1;
            self.add(sequence, 1, 0);
        }
    }
    pub(crate) fn merge(&mut self, other: &Self) {
        self.sampling = other.sampling;
        self.reads += other.reads;
        self.sampled += other.sampled;
        self.windows += other.windows;
        let own_min = if self.candidates.len() == CAPACITY {
            self.candidates.values().map(|v| v.0).min().unwrap_or(0)
        } else {
            0
        };
        let other_min = if other.candidates.len() == CAPACITY {
            other.candidates.values().map(|v| v.0).min().unwrap_or(0)
        } else {
            0
        };
        for (key, value) in self.candidates.iter_mut() {
            let (count, error) = other
                .candidates
                .get(key)
                .copied()
                .unwrap_or((other_min, other_min));
            value.0 += count;
            value.1 += error;
        }
        for (key, &(count, error)) in &other.candidates {
            if !self.candidates.contains_key(key) {
                self.candidates
                    .insert(key.clone(), (count + own_min, error + own_min));
            }
        }
        if self.candidates.len() > CAPACITY {
            let mut entries: Vec<_> = self.candidates.drain().collect();
            entries.sort_unstable_by(|a, b| b.1.0.cmp(&a.1.0).then_with(|| a.0.cmp(&b.0)));
            entries.truncate(CAPACITY);
            self.candidates.extend(entries);
        }
        self.minimum = self
            .candidates
            .iter()
            .map(|(key, value)| Reverse((value.0, key.clone())))
            .collect();
    }
    pub(crate) fn json(&self) -> Value {
        let mut rows: Vec<_> = self
            .candidates
            .iter()
            .filter(|(_, v)| v.0.saturating_sub(v.1) >= 2)
            .collect();
        rows.sort_unstable_by(|a, b| {
            b.1.0
                .saturating_sub(b.1.1)
                .cmp(&a.1.0.saturating_sub(a.1.1))
                .then_with(|| a.0.cmp(b.0))
        });
        json!({"sampling":self.sampling,"reads":self.reads,"sampled_reads":self.sampled,"sampled_windows":self.windows,"candidate_capacity":CAPACITY,"method":"bounded Space-Saving substring occurrence counts; lower bounds exclude replacement error","sequences":rows.into_iter().take(100).map(|(key,&(count,error))|json!({"sequence":String::from_utf8_lossy(key),"count_upper":count,"count_lower":count-error,"replacement_error":error})).collect::<Vec<_>>()})
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn merge_preserves_occurrence_bounds_after_candidate_replacement() {
        let mut left = Overrepresentation::new(1);
        let mut right = Overrepresentation::new(1);
        let mut truth = HashMap::<Vec<u8>, u64>::new();
        for round in 0..4 {
            for key in 0usize..CAPACITY * 2 {
                let sequence = key.to_le_bytes().to_vec();
                let weight = if key % 17 == 0 { 5 } else { 1 };
                *truth.entry(sequence.clone()).or_default() += weight;
                if (key + round) % 2 == 0 {
                    left.add(&sequence, weight, 0);
                } else {
                    right.add(&sequence, weight, 0);
                }
            }
        }
        left.merge(&right);
        assert_eq!(left.candidates.len(), CAPACITY);
        for (key, &(count, error)) in &left.candidates {
            let exact = truth[key];
            assert!(
                count - error <= exact && exact <= count,
                "{count}-{error} <= {exact} <= {count}"
            );
        }
    }
    #[test]
    fn bounds_survive_capacity_and_report_repeated_sequences() {
        let mut c = Overrepresentation::new(1);
        for _ in 0..5 {
            c.update(b"ACGTTGCAGTCAGGCTACGATCGTACTGAC");
        }
        let value = c.json();
        assert_eq!(value["sampled_reads"], 5);
        assert!(
            value["sequences"]
                .as_array()
                .unwrap()
                .iter()
                .any(|r| r["count_lower"] == 5)
        );
        for i in 0..CAPACITY * 2 {
            c.add(&i.to_le_bytes(), 1, 0);
        }
        assert!(c.candidates.len() <= CAPACITY);
        assert!(c.minimum.len() <= CAPACITY * 4);
    }
}
