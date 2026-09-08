use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::encoder::encode;
use crate::structure::{ReadLayout, is_barcode_symbol};

#[derive(Debug)]
pub struct Barcode {
    pub id: String,
    pub sequence: String,
}

#[derive(Debug)]
pub struct BarcodeSet {
    pub symbol: u8,
    pub barcodes: Vec<Barcode>,
    id_lookup: HashMap<String, u32>,
}

#[derive(Debug)]
pub struct BarcodeCatalog {
    sets: HashMap<u8, BarcodeSet>,
}

impl BarcodeCatalog {
    pub fn load(path: &Path, layout: &ReadLayout) -> Result<Self, String> {
        let expected_lengths = barcode_lengths(layout)?;
        let root_symbol = layout.barcode_symbols().first().copied();

        let file = File::open(path).map_err(|e| format!("Could not open barcode file: {e}"))?;

        let mut lines = BufReader::new(file).lines();

        let header = lines
            .next()
            .ok_or("Barcode file is empty")?
            .map_err(|e| format!("Could not read barcode file: {e}"))?;

        let header = header.trim_end_matches('\r');

        if header != "Set\tID\tSequence" {
            return Err("Barcode file header must be: Set<TAB>ID<TAB>Sequence".into());
        }

        let mut sets: HashMap<u8, BarcodeSet> = HashMap::new();
        let mut seen_sequences: HashMap<u8, HashSet<String>> = HashMap::new();
        let mut set_lengths: HashMap<u8, usize> = HashMap::new();

        for (line_index, line) in lines.enumerate() {
            let line_number = line_index + 2;

            let line =
                line.map_err(|e| format!("Could not read barcode file line {line_number}: {e}"))?;

            let line = line.trim_end_matches('\r');

            if line.is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();

            if fields.len() != 3 {
                return Err(format!(
                    "Barcode file line {line_number} has {} columns; expected 3",
                    fields.len()
                ));
            }

            let set_text = fields[0];
            let barcode_name = fields[1];
            let sequence = fields[2].to_ascii_uppercase();

            if set_text.len() != 1 {
                return Err(format!(
                    "Invalid barcode set '{set_text}' on line {line_number}"
                ));
            }

            let symbol = set_text.as_bytes()[0];

            if !is_barcode_symbol(symbol) {
                return Err(format!(
                    "Invalid barcode set '{}'; expected an uppercase A-Z symbol other than T",
                    symbol as char
                ));
            }

            if barcode_name.is_empty() {
                return Err(format!("Barcode ID cannot be empty on line {line_number}"));
            }

            if sequence.is_empty() {
                return Err(format!(
                    "Barcode sequence cannot be empty on line {line_number}"
                ));
            }

            let encoded = encode(sequence.as_bytes())
                .map_err(|e| format!("Invalid sequence for barcode '{barcode_name}': {e}"))?;

            if encoded.n_mask != 0 {
                return Err(format!(
                    "Barcode '{barcode_name}' contains N; whitelist barcodes \
                     must contain only A, C, G, and T"
                ));
            }

            if let Some(&expected_length) = expected_lengths.get(&symbol)
                && sequence.len() != expected_length
            {
                return Err(format!(
                    "Barcode '{}' in set {} has length {}, but the \
                     read structure requires {}",
                    barcode_name,
                    symbol as char,
                    sequence.len(),
                    expected_length
                ));
            }

            if let Some(&existing_length) = set_lengths.get(&symbol) {
                if sequence.len() != existing_length {
                    return Err(format!(
                        "All barcodes in set {} must have the same length",
                        symbol as char
                    ));
                }
            } else {
                set_lengths.insert(symbol, sequence.len());
            }

            let set = sets.entry(symbol).or_insert_with(|| BarcodeSet {
                symbol,
                barcodes: Vec::new(),
                id_lookup: HashMap::new(),
            });

            if set.id_lookup.contains_key(barcode_name) {
                return Err(format!(
                    "Duplicate barcode ID '{}' in set {}",
                    barcode_name, symbol as char
                ));
            }

            // The root has no namespace in which duplicate sequences could be
            // distinguished. At child levels, duplicate physical sequences
            // are allowed globally and are validated only when they occur as
            // siblings in the same compiled routing node.
            if Some(symbol) == root_symbol {
                let sequences = seen_sequences.entry(symbol).or_default();
                if !sequences.insert(sequence.clone()) {
                    return Err(format!(
                        "Duplicate barcode sequence '{}' in root set {}",
                        sequence, symbol as char
                    ));
                }
            }

            let internal_id = u32::try_from(set.barcodes.len())
                .map_err(|_| "Too many barcodes in barcode set")?;

            set.id_lookup.insert(barcode_name.to_string(), internal_id);

            set.barcodes.push(Barcode {
                id: barcode_name.to_string(),
                sequence,
            });
        }

        for symbol in layout.barcode_symbols() {
            if !sets.contains_key(&symbol) {
                return Err(format!(
                    "Read structure requires barcode set {}, \
                     but it is missing from the barcode file",
                    symbol as char
                ));
            }
        }

        Ok(Self { sets })
    }

    pub fn resolve_id(&self, symbol: u8, barcode_name: &str) -> Option<u32> {
        self.sets.get(&symbol)?.id_lookup.get(barcode_name).copied()
    }

    pub fn set(&self, symbol: u8) -> Option<&BarcodeSet> {
        self.sets.get(&symbol)
    }
}

fn barcode_lengths(layout: &ReadLayout) -> Result<HashMap<u8, usize>, String> {
    Ok(layout
        .compile_extraction_plans()?
        .plans()
        .iter()
        .map(|plan| (plan.symbol, plan.logical_length))
        .collect())
}

impl BarcodeSet {
    pub fn barcode(&self, barcode_id: u32) -> Option<&Barcode> {
        self.barcodes.get(usize::try_from(barcode_id).ok()?)
    }
}
