//! Mate-specific proximal technical sequences, loaded and validated once.
use crate::structure::{ReadLayout, ReadMate, ReadStructure};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug, Default)]
pub(crate) struct TechnicalSequences(pub(crate) [Option<Vec<u8>>; 2]);

pub(crate) fn dna(sequence: &str) -> Result<Vec<u8>, String> {
    if sequence.is_empty() || !sequence.bytes().all(|b| b"ACGTacgt".contains(&b)) {
        return Err("Sequence must contain only A/C/G/T and must not be empty".into());
    }
    Ok(sequence.to_ascii_uppercase().into_bytes())
}

pub(crate) fn structures(layout: &ReadLayout) -> [Option<&ReadStructure>; 2] {
    match layout {
        ReadLayout::Single { r1 } => [r1.as_ref(), None],
        ReadLayout::Paired { r1, r2 } => [r1.as_ref(), r2.as_ref()],
    }
}

impl TechnicalSequences {
    pub(crate) fn load(path: Option<&Path>, layout: &ReadLayout) -> Result<Self, String> {
        let Some(path) = path else {
            return Ok(Self::default());
        };
        let file =
            File::open(path).map_err(|e| format!("Could not open technical sequences: {e}"))?;
        let mut lines = BufReader::new(file).lines();
        if lines
            .next()
            .transpose()
            .map_err(|e| e.to_string())?
            .as_deref()
            != Some("Mate\tSequence")
        {
            return Err("Technical sequence header must be Mate<TAB>Sequence".into());
        }
        let mut result = Self::default();
        let structures = structures(layout);
        for (index, line) in lines.enumerate() {
            let line = line.map_err(|e| e.to_string())?;
            let fields: Vec<_> = line.split('\t').collect();
            if fields.len() != 2 {
                return Err(format!("Malformed technical sequence line {}", index + 2));
            }
            let mate = match fields[0] {
                "R1" => 0,
                "R2" => 1,
                _ => return Err("Technical sequence mate must be R1 or R2".into()),
            };
            if result.0[mate].is_some() {
                return Err(format!("Duplicate technical sequence for R{}", mate + 1));
            }
            let sequence = dna(fields[1])?;
            let length: usize = structures[mate]
                .into_iter()
                .flat_map(|s| &s.segments)
                .filter(|s| s.symbol == b'T')
                .map(|s| s.end - s.start)
                .sum();
            if sequence.len() != length {
                return Err(format!(
                    "R{} technical sequence has {} bases; proximal T segments require {length}",
                    mate + 1,
                    sequence.len()
                ));
            }
            result.0[mate] = Some(sequence);
        }
        if result.0.iter().all(Option::is_none) {
            return Err("Technical sequence file has no entries".into());
        }
        Ok(result)
    }
    pub(crate) fn mate(&self, mate: ReadMate) -> Option<&[u8]> {
        self.0[usize::from(mate == ReadMate::R2)].as_deref()
    }
}
