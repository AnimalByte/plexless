pub const MAX_BARCODE_LENGTH: usize = 32;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct EncodedBarcode {
    pub value: u64,
    pub n_mask: u32,
    pub length: u8,
}

pub fn encode(seq: &[u8]) -> Result<EncodedBarcode, String> {
    if seq.is_empty() {
        return Err("Barcode cannot be empty".into());
    }

    if seq.len() > MAX_BARCODE_LENGTH {
        return Err(format!("Barcode cannot exceed {MAX_BARCODE_LENGTH} bases"));
    }

    let length = u8::try_from(seq.len()).map_err(|_| "Barcode length exceeds u8 capacity")?;

    let mut value = 0u64;
    let mut n_mask = 0u32;

    for &base in seq {
        value <<= 2;
        n_mask <<= 1;

        match base {
            b'A' | b'a' => {}

            b'C' | b'c' => {
                value |= 0b01;
            }

            b'G' | b'g' => {
                value |= 0b10;
            }

            b'T' | b't' => {
                value |= 0b11;
            }

            b'N' | b'n' => {
                n_mask |= 1;
            }

            _ => {
                return Err(format!("Invalid DNA base: '{}'", base as char));
            }
        }
    }

    Ok(EncodedBarcode {
        value,
        n_mask,
        length,
    })
}
