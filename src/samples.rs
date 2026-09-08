use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::barcodes::BarcodeCatalog;
use crate::structure::ReadLayout;

#[derive(Debug)]
pub struct Sample {
    pub name: String,
    pub barcode_ids: Vec<u32>,
}

#[derive(Debug)]
pub struct SampleSheet {
    pub samples: Vec<Sample>,
}

impl SampleSheet {
    pub fn load(
        path: &Path,
        layout: &ReadLayout,
        barcodes: &BarcodeCatalog,
    ) -> Result<Self, String> {
        let symbols = layout.barcode_symbols();

        if symbols.is_empty() {
            return Err("Read structure contains no sample barcode segments".into());
        }

        let file = File::open(path).map_err(|e| format!("Could not open sample sheet: {e}"))?;

        let mut lines = BufReader::new(file).lines();

        let header = lines
            .next()
            .ok_or("Sample sheet is empty")?
            .map_err(|e| format!("Could not read sample sheet: {e}"))?;

        let header = header.trim_end_matches('\r');
        let columns: Vec<&str> = header.split('\t').collect();

        if columns.is_empty() || columns[0] != "Sample" {
            return Err("First sample-sheet column must be 'Sample'".into());
        }

        if columns.len() != symbols.len() + 1 {
            return Err(format!(
                "Sample sheet must contain exactly {} columns",
                symbols.len() + 1
            ));
        }

        for (index, symbol) in symbols.iter().enumerate() {
            let expected = (*symbol as char).to_string();

            if columns[index + 1] != expected {
                return Err(format!(
                    "Expected column '{expected}' at position {}",
                    index + 2
                ));
            }
        }

        let mut sample_names = HashSet::new();

        let mut combinations: HashMap<Vec<u32>, String> = HashMap::new();

        let mut samples = Vec::new();

        for (line_index, line) in lines.enumerate() {
            let line_number = line_index + 2;

            let line =
                line.map_err(|e| format!("Could not read sample sheet line {line_number}: {e}"))?;

            let line = line.trim_end_matches('\r');

            if line.is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();

            if fields.len() != columns.len() {
                return Err(format!(
                    "Line {line_number} has {} columns; expected {}",
                    fields.len(),
                    columns.len()
                ));
            }

            let sample_name = fields[0];

            if sample_name.is_empty() {
                return Err(format!("Line {line_number} has an empty sample name"));
            }

            if !sample_names.insert(sample_name.to_string()) {
                return Err(format!("Duplicate sample name '{sample_name}'"));
            }

            let mut barcode_ids = Vec::with_capacity(symbols.len());

            for (index, symbol) in symbols.iter().enumerate() {
                let barcode_name = fields[index + 1];

                if barcode_name.is_empty() {
                    return Err(format!(
                        "Sample '{sample_name}' has an empty {} barcode ID",
                        *symbol as char
                    ));
                }

                let barcode_id = barcodes.resolve_id(*symbol, barcode_name).ok_or_else(|| {
                    format!(
                        "Sample '{sample_name}' references unknown \
                             barcode '{barcode_name}' in set {}",
                        *symbol as char
                    )
                })?;

                barcode_ids.push(barcode_id);
            }

            if let Some(existing_sample) = combinations.get(&barcode_ids) {
                return Err(format!(
                    "Samples '{existing_sample}' and '{sample_name}' \
                     use the same barcode combination"
                ));
            }

            combinations.insert(barcode_ids.clone(), sample_name.to_string());

            samples.push(Sample {
                name: sample_name.to_string(),
                barcode_ids,
            });
        }

        if samples.is_empty() {
            return Err("Sample sheet contains no samples".into());
        }

        Ok(Self { samples })
    }
}
