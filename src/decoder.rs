use std::collections::HashMap;

use crate::barcodes::BarcodeSet;
use crate::encoder::{EncodedBarcode, encode};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DecodeResult {
    /// `barcode_id` is the decoder-local candidate index. Routing nodes use it
    /// directly to index their contiguous target arrays.
    Exact {
        barcode_id: u32,
    },
    Corrected {
        barcode_id: u32,
        mismatches: u8,
    },
    Ambiguous,
    Unmatched,
}

#[derive(Debug, Clone, Copy)]
struct CorrectionTarget {
    barcode_id: u32,
    mismatches: u8,
}

#[derive(Debug)]
pub struct Decoder {
    barcodes: Vec<EncodedBarcode>,
    exact: HashMap<u64, u32>,
    corrected: HashMap<u64, CorrectionTarget>,
    length: u8,
    max_mismatches: u8,
    indexed: bool,
}

struct NeighborBuilder<'a> {
    original: u64,
    length: u8,
    distance: u8,
    barcode_id: u32,
    exact: &'a HashMap<u64, u32>,
    corrected: &'a mut HashMap<u64, CorrectionTarget>,
}

impl Decoder {
    /// Builds a decoder over an entire set. This remains useful for callers
    /// that do not need route-local namespaces; hierarchical routing uses
    /// `new_subset` instead.
    pub fn new(set: &BarcodeSet, max_mismatches: u8) -> Result<Self, String> {
        let candidate_ids: Result<Vec<u32>, _> =
            (0..set.barcodes.len()).map(u32::try_from).collect();
        Self::new_subset(
            set,
            &candidate_ids.map_err(|_| "Too many barcodes")?,
            max_mismatches,
            "routing root",
        )
    }

    pub(crate) fn new_subset(
        set: &BarcodeSet,
        candidate_ids: &[u32],
        max_mismatches: u8,
        parent_path: &str,
    ) -> Result<Self, String> {
        if candidate_ids.is_empty() {
            return Err(format!(
                "Barcode candidate set {} is empty under {parent_path}",
                set.symbol as char
            ));
        }

        let mut barcodes = Vec::with_capacity(candidate_ids.len());
        let mut names = Vec::with_capacity(candidate_ids.len());
        let mut exact = HashMap::with_capacity(candidate_ids.len());
        let first_barcode = set
            .barcode(candidate_ids[0])
            .ok_or("Internal error: invalid barcode candidate ID")?;
        let length = encode(first_barcode.sequence.as_bytes())?.length;

        if max_mismatches > length {
            return Err(format!(
                "Maximum mismatches cannot exceed barcode length for set {} under {parent_path}",
                set.symbol as char
            ));
        }

        for &global_id in candidate_ids {
            let barcode = set
                .barcode(global_id)
                .ok_or("Internal error: invalid barcode candidate ID")?;
            let encoded = encode(barcode.sequence.as_bytes())?;

            if encoded.length != length {
                return Err(format!(
                    "All barcodes in set {} must have the same length",
                    set.symbol as char
                ));
            }
            if encoded.n_mask != 0 {
                return Err(format!("Whitelist barcode '{}' contains N", barcode.id));
            }

            let local_id = u32::try_from(barcodes.len()).map_err(|_| "Too many barcodes")?;
            if let Some(existing_id) = exact.insert(encoded.value, local_id) {
                let existing_name = names
                    .get(usize::try_from(existing_id).expect("local barcode ID fits usize"))
                    .expect("exact index refers to an existing barcode name");
                return Err(format!(
                    "Duplicate barcode sequence in routing node under {parent_path}: \
                     {} barcodes '{}' and '{}' have the same sequence",
                    set.symbol as char, existing_name, barcode.id
                ));
            }

            names.push(barcode.id.as_str());
            barcodes.push(encoded);
        }

        validate_barcode_distance(&barcodes, &names, set.symbol, max_mismatches, parent_path)?;

        let indexed = correction_index_size(length, max_mismatches).saturating_mul(barcodes.len())
            <= 1_000_000;
        let corrected = build_correction_index(
            &barcodes,
            &exact,
            length,
            if indexed { max_mismatches } else { 0 },
        )
        .map_err(|error| format!("{error} under {parent_path}"))?;

        Ok(Self {
            barcodes,
            exact,
            corrected,
            length,
            max_mismatches,
            indexed,
        })
    }

    pub fn decode(&self, observed: EncodedBarcode) -> DecodeResult {
        if observed.length != self.length {
            return DecodeResult::Unmatched;
        }
        if observed.n_mask != 0 {
            return self.decode_with_n(observed);
        }
        if let Some(&barcode_id) = self.exact.get(&observed.value) {
            return DecodeResult::Exact { barcode_id };
        }
        if let Some(target) = self.corrected.get(&observed.value) {
            return DecodeResult::Corrected {
                barcode_id: target.barcode_id,
                mismatches: target.mismatches,
            };
        }
        if self.indexed {
            DecodeResult::Unmatched
        } else {
            self.decode_with_n(observed)
        }
    }

    fn decode_with_n(&self, observed: EncodedBarcode) -> DecodeResult {
        let n_count = u8::try_from(observed.n_mask.count_ones())
            .expect("N mask cannot contain more than 32 set bits");

        if n_count > self.max_mismatches {
            return DecodeResult::Unmatched;
        }

        let ignored_positions = expand_n_mask(observed.n_mask);
        let mut best_id = 0u32;
        let mut best_distance = u8::MAX;
        let mut tied = false;

        // This scans only the candidates in this routing node. Barcodes in
        // unrelated parent namespaces are never considered.
        for (index, expected) in self.barcodes.iter().enumerate() {
            let diff = observed.value ^ expected.value;
            let mismatch_bits = (diff | (diff >> 1)) & 0x5555_5555_5555_5555;
            let known_mismatches = u8::try_from((mismatch_bits & !ignored_positions).count_ones())
                .expect("Barcode cannot contain more than 32 mismatches");
            let distance = known_mismatches + n_count;

            if distance > self.max_mismatches {
                continue;
            }
            if distance < best_distance {
                best_distance = distance;
                best_id = u32::try_from(index).expect("Barcode index must fit in u32");
                tied = false;
            } else if distance == best_distance {
                tied = true;
            }
        }

        if best_distance == u8::MAX {
            DecodeResult::Unmatched
        } else if tied {
            DecodeResult::Ambiguous
        } else {
            DecodeResult::Corrected {
                barcode_id: best_id,
                mismatches: best_distance,
            }
        }
    }
}

impl NeighborBuilder<'_> {
    fn add(&mut self, current: u64, start_position: u8, remaining: u8) -> Result<(), String> {
        if remaining == 0 {
            if self.exact.contains_key(&current) {
                return Err("Unsafe barcode correction collision detected".into());
            }
            if let Some(existing) = self.corrected.get(&current) {
                if existing.barcode_id != self.barcode_id {
                    return Err("Unsafe barcode correction collision detected".into());
                }
                return Ok(());
            }

            self.corrected.insert(
                current,
                CorrectionTarget {
                    barcode_id: self.barcode_id,
                    mismatches: self.distance,
                },
            );
            return Ok(());
        }

        for position in start_position..self.length {
            let shift = u32::from(position) * 2;
            let mask = 0b11u64 << shift;
            let original_base = (self.original >> shift) & 0b11;

            for replacement in 0u64..4 {
                if replacement == original_base {
                    continue;
                }
                let mutated = (current & !mask) | (replacement << shift);
                self.add(mutated, position + 1, remaining - 1)?;
            }
        }
        Ok(())
    }
}

fn validate_barcode_distance(
    barcodes: &[EncodedBarcode],
    names: &[&str],
    symbol: u8,
    max_mismatches: u8,
    parent_path: &str,
) -> Result<(), String> {
    if max_mismatches == 0 {
        return Ok(());
    }

    let required_distance = max_mismatches
        .checked_mul(2)
        .and_then(|value| value.checked_add(1))
        .ok_or("Mismatch distance overflow")?;

    for i in 0..barcodes.len() {
        for j in (i + 1)..barcodes.len() {
            let distance = hamming_distance(barcodes[i].value, barcodes[j].value);
            if distance < required_distance {
                return Err(format!(
                    "Unsafe barcode correction collision under {parent_path}:\n\
                     {} barcodes '{}' and '{}' require minimum Hamming distance {} \
                     for max_mismatches={}, but observed distance is {}.",
                    symbol as char, names[i], names[j], required_distance, max_mismatches, distance
                ));
            }
        }
    }
    Ok(())
}

// Bound neighborhood construction before allocating. Above this cap, decode by
// packed distance over route-local candidates, with O(catalog) memory.
fn correction_index_size(length: u8, budget: u8) -> usize {
    let mut term = 1usize;
    let mut total = 0usize;
    for distance in 1..=usize::from(budget) {
        term = term
            .saturating_mul(usize::from(length) + 1 - distance)
            .saturating_mul(3)
            / distance;
        total = total.saturating_add(term);
        if total > 1_000_000 {
            return total;
        }
    }
    total
}

fn build_correction_index(
    barcodes: &[EncodedBarcode],
    exact: &HashMap<u64, u32>,
    length: u8,
    max_mismatches: u8,
) -> Result<HashMap<u64, CorrectionTarget>, String> {
    let mut corrected = HashMap::new();
    if max_mismatches == 0 {
        return Ok(corrected);
    }

    for (index, barcode) in barcodes.iter().enumerate() {
        let barcode_id = u32::try_from(index).map_err(|_| "Too many barcodes")?;
        for distance in 1..=max_mismatches {
            let mut builder = NeighborBuilder {
                original: barcode.value,
                length,
                distance,
                barcode_id,
                exact,
                corrected: &mut corrected,
            };
            builder.add(barcode.value, 0, distance)?;
        }
    }
    Ok(corrected)
}

pub(crate) fn hamming_distance(a: u64, b: u64) -> u8 {
    let diff = a ^ b;
    let mismatch_bits = (diff | (diff >> 1)) & 0x5555_5555_5555_5555;
    u8::try_from(mismatch_bits.count_ones())
        .expect("Barcode cannot contain more than 32 mismatches")
}

fn expand_n_mask(mask: u32) -> u64 {
    let mut input = mask;
    let mut output = 0u64;
    let mut position = 0u32;

    while input != 0 {
        if input & 1 != 0 {
            output |= 1u64 << (position * 2);
        }
        input >>= 1;
        position += 1;
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Production constructors reject this deliberately unsafe two-barcode
    /// geometry. Building it directly here exercises the conservative runtime
    /// ambiguity branch as a defense-in-depth behavior.
    pub(super) fn unsafe_n_test_decoder() -> Decoder {
        let barcodes = vec![encode(b"AAAA").unwrap(), encode(b"AAAT").unwrap()];
        let exact = HashMap::from([(barcodes[0].value, 0), (barcodes[1].value, 1)]);
        Decoder {
            barcodes,
            exact,
            corrected: HashMap::new(),
            length: 4,
            max_mismatches: 1,
            indexed: true,
        }
    }

    #[test]
    fn n_aware_decode_covers_exact_corrected_ambiguous_and_unmatched() {
        let decoder = unsafe_n_test_decoder();

        assert_eq!(
            decoder.decode(encode(b"AAAA").unwrap()),
            DecodeResult::Exact { barcode_id: 0 }
        );
        assert_eq!(
            decoder.decode(encode(b"ANAA").unwrap()),
            DecodeResult::Corrected {
                barcode_id: 0,
                mismatches: 1
            }
        );
        assert_eq!(
            decoder.decode(encode(b"AAAN").unwrap()),
            DecodeResult::Ambiguous
        );
        assert_eq!(
            decoder.decode(encode(b"AANN").unwrap()),
            DecodeResult::Unmatched
        );
    }
}

#[cfg(test)]
mod scaling_tests {
    use super::*;
    use crate::barcodes::BarcodeCatalog;
    use crate::structure::ReadLayout;

    fn catalog() -> BarcodeCatalog {
        let path = std::env::temp_dir().join(format!(
            "plexless-decoder-scaling-{}.tsv",
            std::process::id()
        ));
        std::fs::write(
            &path,
            "Set\tID\tSequence\nA\ta\tACGTTGCACTGATCGAGGTACCTATGCCAGTA\n",
        )
        .unwrap();
        let result =
            BarcodeCatalog::load(&path, &ReadLayout::single(Some("R1_32A")).unwrap()).unwrap();
        std::fs::remove_file(path).unwrap();
        result
    }
    #[test]
    fn large_budgets_use_bounded_fallback_and_n_aware_distance() {
        let catalog = catalog();
        let set = catalog.set(b'A').unwrap();
        for budget in [0, 1, 2, 3, 4, 8, 16, 32] {
            let decoder = Decoder::new(set, budget).unwrap();
            assert!(decoder.corrected.len() <= 1_000_000);
            if budget >= 4 {
                assert!(!decoder.indexed);
                assert!(decoder.corrected.is_empty());
            }
            let mut observed = set.barcodes[0].sequence.as_bytes().to_vec();
            assert_eq!(
                decoder.decode(encode(&observed).unwrap()),
                DecodeResult::Exact { barcode_id: 0 }
            );
            for (i, base) in observed.iter_mut().enumerate().take(usize::from(budget)) {
                *base = if i % 2 == 0 {
                    b'N'
                } else if *base == b'A' {
                    b'C'
                } else {
                    b'A'
                };
            }
            if budget > 0 {
                assert_eq!(
                    decoder.decode(encode(&observed).unwrap()),
                    DecodeResult::Corrected {
                        barcode_id: 0,
                        mismatches: budget
                    }
                );
            }
            if budget < 32 {
                observed[usize::from(budget)] = b'N';
                assert_eq!(
                    decoder.decode(encode(&observed).unwrap()),
                    DecodeResult::Unmatched
                );
            }
        }
        assert!(
            Decoder::new(set, 33)
                .unwrap_err()
                .contains("barcode length")
        );
    }
    #[test]
    fn fallback_ambiguity_and_exact_precedence() {
        let mut decoder = super::tests::unsafe_n_test_decoder();
        decoder.indexed = false;
        assert_eq!(
            decoder.decode(encode(b"AAAC").unwrap()),
            DecodeResult::Ambiguous
        );
        assert_eq!(
            decoder.decode(encode(b"AAAA").unwrap()),
            DecodeResult::Exact { barcode_id: 0 }
        );
    }
}
