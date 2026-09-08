use std::collections::HashMap;

use crate::barcodes::BarcodeSet;
use crate::encoder::{EncodedBarcode, encode};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DecodeResult {
    Exact { barcode_id: u32 },
    Corrected { barcode_id: u32, mismatches: u8 },
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
    pub fn new(set: &BarcodeSet, max_mismatches: u8) -> Result<Self, String> {
        if set.barcodes.is_empty() {
            return Err(format!("Barcode set {} is empty", set.symbol as char));
        }

        if max_mismatches > 2 {
            return Err("Maximum mismatches greater than 2 are not supported".into());
        }

        let mut barcodes = Vec::with_capacity(set.barcodes.len());
        let mut exact = HashMap::with_capacity(set.barcodes.len());

        let first = encode(set.barcodes[0].sequence.as_bytes())?;
        let length = first.length;

        if max_mismatches > length {
            return Err("Maximum mismatches cannot exceed barcode length".into());
        }

        for (index, barcode) in set.barcodes.iter().enumerate() {
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

            let barcode_id = u32::try_from(index).map_err(|_| "Too many barcodes")?;

            if exact.insert(encoded.value, barcode_id).is_some() {
                return Err(format!(
                    "Duplicate barcode sequence '{}' in set {}",
                    barcode.sequence, set.symbol as char
                ));
            }

            barcodes.push(encoded);
        }

        validate_barcode_distance(&barcodes, set, max_mismatches)?;

        let corrected = build_correction_index(&barcodes, &exact, length, max_mismatches)?;

        Ok(Self {
            barcodes,
            exact,
            corrected,
            length,
            max_mismatches,
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

        DecodeResult::Unmatched
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
    set: &BarcodeSet,
    max_mismatches: u8,
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
                    "Barcodes '{}' and '{}' in set {} have Hamming \
                     distance {}, but at least {} is required for \
                     {}-mismatch correction",
                    set.barcodes[i].id,
                    set.barcodes[j].id,
                    set.symbol as char,
                    distance,
                    required_distance,
                    max_mismatches
                ));
            }
        }
    }

    Ok(())
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

fn hamming_distance(a: u64, b: u64) -> u8 {
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
