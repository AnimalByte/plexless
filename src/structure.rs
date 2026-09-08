use crate::encoder::MAX_BARCODE_LENGTH;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadMate {
    R1,
    R2,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    Forward,
    ReverseComplement,
}

#[derive(Debug, Clone, Copy)]
pub struct CompiledSegment {
    pub start: usize,
    pub end: usize,
    pub symbol: u8,
    pub orientation: Orientation,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BarcodePiece {
    pub mate: ReadMate,
    pub start: usize,
    pub end: usize,
    pub orientation: Orientation,
}

#[derive(Debug)]
pub struct BarcodeExtractionPlan {
    pub symbol: u8,
    pub logical_length: usize,
    pub pieces: Vec<BarcodePiece>,
}

/// Startup-compiled locations for every logical barcode in a read layout.
///
/// Plans are ordered by symbol, and pieces within a plan are ordered as all R1
/// pieces followed by all R2 pieces. The fixed symbol index covers the generic
/// ASCII A-Z syntax; it does not limit the number of routing levels.
#[derive(Debug)]
pub struct BarcodeExtractionPlans {
    plans: Vec<BarcodeExtractionPlan>,
    symbol_to_plan: [usize; 26],
    r1_prefix_len: usize,
    r2_prefix_len: usize,
    paired: bool,
}

pub fn is_barcode_symbol(symbol: u8) -> bool {
    symbol.is_ascii_uppercase() && symbol != b'T'
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
            if symbol != b'T' && !is_barcode_symbol(symbol) {
                return Err(format!(
                    "Invalid segment symbol '{}'; expected an uppercase A-Z barcode symbol or T",
                    symbol as char
                ));
            }
            index += 1;

            let mut orientation = Orientation::Forward;
            if index < bytes.len() && bytes[index] == b'(' {
                if !bytes[index..].starts_with(b"(rc)") {
                    return Err(format!(
                        "Malformed orientation modifier after segment {}; expected '(rc)'",
                        symbol as char
                    ));
                }
                if symbol == b'T' {
                    return Err("Orientation '(rc)' is only valid for barcode segments".into());
                }
                orientation = Orientation::ReverseComplement;
                index += 4;

                if index < bytes.len() && bytes[index] == b'(' {
                    return Err(format!(
                        "Malformed orientation modifier after segment {}",
                        symbol as char
                    ));
                }
            }

            let end = start
                .checked_add(length)
                .ok_or("Segment coordinate overflow")?;

            segments.push(CompiledSegment {
                start,
                end,
                symbol,
                orientation,
            });
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
            if is_barcode_symbol(segment.symbol) && !symbols.contains(&segment.symbol) {
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

    pub fn compile_extraction_plans(&self) -> Result<BarcodeExtractionPlans, String> {
        BarcodeExtractionPlans::compile(self)
    }
}

impl BarcodeExtractionPlans {
    pub fn compile(layout: &ReadLayout) -> Result<Self, String> {
        let symbols = layout.barcode_symbols();
        let mut plans = Vec::with_capacity(symbols.len());
        let mut symbol_to_plan = [usize::MAX; 26];

        for symbol in symbols {
            let mut pieces = Vec::new();
            append_plan_pieces(layout_r1(layout), symbol, &mut pieces);
            append_plan_pieces(layout_r2(layout), symbol, &mut pieces);

            let logical_length = pieces.iter().try_fold(0usize, |total, piece| {
                total
                    .checked_add(piece.end - piece.start)
                    .ok_or("Barcode length overflow")
            })?;

            if logical_length > MAX_BARCODE_LENGTH {
                return Err(format!(
                    "Logical barcode {} is {} bases long; maximum supported length is {}",
                    symbol as char, logical_length, MAX_BARCODE_LENGTH
                ));
            }

            let plan_index = plans.len();
            symbol_to_plan[usize::from(symbol - b'A')] = plan_index;
            plans.push(BarcodeExtractionPlan {
                symbol,
                logical_length,
                pieces,
            });
        }

        let (r1_prefix_len, r2_prefix_len, paired) = match layout {
            ReadLayout::Single { r1 } => (r1.as_ref().map_or(0, |s| s.prefix_len), 0, false),
            ReadLayout::Paired { r1, r2 } => (
                r1.as_ref().map_or(0, |s| s.prefix_len),
                r2.as_ref().map_or(0, |s| s.prefix_len),
                true,
            ),
        };

        Ok(Self {
            plans,
            symbol_to_plan,
            r1_prefix_len,
            r2_prefix_len,
            paired,
        })
    }

    pub fn plans(&self) -> &[BarcodeExtractionPlan] {
        &self.plans
    }

    pub fn plan_index(&self, symbol: u8) -> Option<usize> {
        if !is_barcode_symbol(symbol) {
            return None;
        }
        let index = self.symbol_to_plan[usize::from(symbol - b'A')];
        (index != usize::MAX).then_some(index)
    }

    pub fn plan(&self, index: usize) -> Option<&BarcodeExtractionPlan> {
        self.plans.get(index)
    }

    pub fn r1_prefix_len(&self) -> usize {
        self.r1_prefix_len
    }

    pub fn r2_prefix_len(&self) -> usize {
        self.r2_prefix_len
    }

    pub fn reads_are_long_enough(&self, r1: &[u8], r2: Option<&[u8]>) -> Result<bool, String> {
        if r1.len() < self.r1_prefix_len {
            return Ok(false);
        }
        if self.paired {
            let r2 = r2.ok_or("Missing R2 sequence")?;
            if r2.len() < self.r2_prefix_len {
                return Ok(false);
            }
        }
        Ok(true)
    }

    pub fn assemble(
        &self,
        plan_index: usize,
        r1: &[u8],
        r2: Option<&[u8]>,
        buffer: &mut [u8; MAX_BARCODE_LENGTH],
    ) -> Result<Option<usize>, String> {
        let plan = self
            .plan(plan_index)
            .ok_or("Internal error: invalid barcode extraction plan")?;
        let mut output_index = 0usize;

        for piece in &plan.pieces {
            let read = match piece.mate {
                ReadMate::R1 => r1,
                ReadMate::R2 => r2.ok_or("Missing R2 sequence")?,
            };
            if piece.end > read.len() {
                return Ok(None);
            }

            let segment = &read[piece.start..piece.end];
            match piece.orientation {
                Orientation::Forward => {
                    buffer[output_index..output_index + segment.len()].copy_from_slice(segment);
                }
                Orientation::ReverseComplement => {
                    for (&base, output) in segment
                        .iter()
                        .rev()
                        .zip(buffer[output_index..output_index + segment.len()].iter_mut())
                    {
                        *output = complement(base);
                    }
                }
            }
            output_index += segment.len();
        }

        debug_assert_eq!(output_index, plan.logical_length);
        Ok(Some(output_index))
    }
}

fn layout_r1(layout: &ReadLayout) -> Option<&ReadStructure> {
    match layout {
        ReadLayout::Single { r1 } | ReadLayout::Paired { r1, .. } => r1.as_ref(),
    }
}

fn layout_r2(layout: &ReadLayout) -> Option<&ReadStructure> {
    match layout {
        ReadLayout::Single { .. } => None,
        ReadLayout::Paired { r2, .. } => r2.as_ref(),
    }
}

fn append_plan_pieces(
    structure: Option<&ReadStructure>,
    symbol: u8,
    pieces: &mut Vec<BarcodePiece>,
) {
    let Some(structure) = structure else {
        return;
    };

    pieces.extend(
        structure
            .segments
            .iter()
            .filter(|segment| segment.symbol == symbol)
            .map(|segment| BarcodePiece {
                mate: structure.mate,
                start: segment.start,
                end: segment.end,
                orientation: segment.orientation,
            }),
    );
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'a' => b't',
        b'C' => b'G',
        b'c' => b'g',
        b'G' => b'C',
        b'g' => b'c',
        b'T' => b'A',
        b't' => b'a',
        b'N' | b'n' => base,
        _ => base,
    }
}
