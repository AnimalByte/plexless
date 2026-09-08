#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadMate {
    R1,
    R2,
}

#[derive(Debug, Clone, Copy)]
pub struct CompiledSegment {
    pub start: usize,
    pub end: usize,
    pub symbol: u8,
}

#[derive(Debug)]
pub struct ReadStructure {
    pub mate: ReadMate,
    pub segments: Vec<CompiledSegment>,
    pub prefix_len: usize,
}

#[derive(Debug)]
pub enum ReadLayout {
    Single {
        r1: Option<ReadStructure>,
    },
    Paired {
        r1: Option<ReadStructure>,
        r2: Option<ReadStructure>,
    },
}

impl ReadStructure {
    pub fn parse(input: &str) -> Result<Self, String> {
        if !input.is_ascii() {
            return Err("Read structure must contain only ASCII characters".into());
        }

        let (mate, structure) = input
            .split_once('_')
            .ok_or("Read structure must contain '_'")?;

        let mate = match mate {
            "R1" => ReadMate::R1,
            "R2" => ReadMate::R2,
            _ => return Err("Read structure must start with R1 or R2".into()),
        };

        if structure.is_empty() {
            return Err("Read structure cannot be empty".into());
        }

        let bytes = structure.as_bytes();

        let mut segments = Vec::new();
        let mut index = 0;
        let mut start = 0usize;

        while index < bytes.len() {
            let length_start = index;
            let mut length = 0usize;

            while index < bytes.len() && bytes[index].is_ascii_digit() {
                let digit = (bytes[index] - b'0') as usize;

                length = length
                    .checked_mul(10)
                    .and_then(|n| n.checked_add(digit))
                    .ok_or("Segment length overflow")?;

                index += 1;
            }

            if index == length_start {
                return Err("Expected segment length".into());
            }

            if length == 0 {
                return Err("Segment length must be greater than 0".into());
            }

            if index >= bytes.len() {
                return Err("Missing segment symbol".into());
            }

            let symbol = bytes[index];

            if !matches!(symbol, b'A' | b'B' | b'C' | b'T') {
                return Err(format!(
                    "Invalid segment symbol '{}'; expected A, B, C, or T",
                    symbol as char
                ));
            }

            index += 1;

            let end = start
                .checked_add(length)
                .ok_or("Segment coordinate overflow")?;

            segments.push(CompiledSegment { start, end, symbol });

            start = end;
        }

        Ok(Self {
            mate,
            segments,
            prefix_len: start,
        })
    }

    pub fn barcode_symbols(&self) -> Vec<u8> {
        let mut symbols = Vec::new();

        for segment in &self.segments {
            if segment.symbol != b'T' && !symbols.contains(&segment.symbol) {
                symbols.push(segment.symbol);
            }
        }

        symbols.sort_unstable();
        symbols
    }
}

impl ReadLayout {
    pub fn single(r1: Option<&str>) -> Result<Self, String> {
        let r1 = match r1 {
            Some(value) => {
                let structure = ReadStructure::parse(value)?;

                if structure.mate != ReadMate::R1 {
                    return Err("Single-end structure must describe R1".into());
                }

                Some(structure)
            }

            None => None,
        };

        Ok(Self::Single { r1 })
    }

    pub fn paired(r1: Option<&str>, r2: Option<&str>) -> Result<Self, String> {
        let r1 = match r1 {
            Some(value) => {
                let structure = ReadStructure::parse(value)?;

                if structure.mate != ReadMate::R1 {
                    return Err("R1 structure must describe R1".into());
                }

                Some(structure)
            }

            None => None,
        };

        let r2 = match r2 {
            Some(value) => {
                let structure = ReadStructure::parse(value)?;

                if structure.mate != ReadMate::R2 {
                    return Err("R2 structure must describe R2".into());
                }

                Some(structure)
            }

            None => None,
        };

        Ok(Self::Paired { r1, r2 })
    }

    pub fn barcode_symbols(&self) -> Vec<u8> {
        let mut symbols = Vec::new();

        match self {
            Self::Single { r1 } => {
                if let Some(r1) = r1 {
                    symbols.extend(r1.barcode_symbols());
                }
            }

            Self::Paired { r1, r2 } => {
                if let Some(r1) = r1 {
                    symbols.extend(r1.barcode_symbols());
                }

                if let Some(r2) = r2 {
                    symbols.extend(r2.barcode_symbols());
                }
            }
        }

        symbols.sort_unstable();
        symbols.dedup();

        symbols
    }
}
