use needletail::errors::ParseError;
use needletail::parser::SequenceRecord;
use needletail::{FastxReader, parse_fastx_file, parse_fastx_reader};
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

pub struct InputReader {
    reader: Box<dyn FastxReader>,
}

impl InputReader {
    pub fn open(path: &Path) -> Self {
        Self::try_open(path).expect("Could not open FASTQ file")
    }

    pub(crate) fn try_open(path: &Path) -> Result<Self, String> {
        let parse_error = match parse_fastx_file(path) {
            Ok(reader) => return Ok(Self { reader }),
            Err(error) => error,
        };

        // Needletail's XZ feature cannot coexist with HTSlib's LZMA linkage
        // because the two crates currently expose different native `links`
        // packages. Keep the qualified common-format path above unchanged and
        // use xz2 only as the compatibility fallback for an actual XZ stream.
        let mut file = File::open(path)
            .map_err(|e| format!("Could not open FASTQ input '{}': {e}", path.display()))?;
        let mut magic = [0u8; 6];
        let bytes_read = file
            .read(&mut magic)
            .map_err(|e| format!("Could not inspect FASTQ input '{}': {e}", path.display()))?;
        file.seek(SeekFrom::Start(0))
            .map_err(|e| format!("Could not rewind FASTQ input '{}': {e}", path.display()))?;
        if bytes_read == magic.len() && magic == [0xfd, b'7', b'z', b'X', b'Z', 0x00] {
            return Self::from_reader(xz2::read::XzDecoder::new_multi_decoder(file), path);
        }
        Err(format!(
            "Could not open FASTQ input '{}': {parse_error}",
            path.display()
        ))
    }

    pub(crate) fn from_reader<R>(reader: R, source: &Path) -> Result<Self, String>
    where
        R: Read + Send + 'static,
    {
        let reader = parse_fastx_reader(reader).map_err(|e| {
            format!(
                "Could not initialize FASTQ parser for '{}': {e}",
                source.display()
            )
        })?;

        Ok(Self { reader })
    }

    pub fn next_record(&mut self) -> Option<Result<SequenceRecord<'_>, ParseError>> {
        self.reader.next()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn xz_fastq_remains_supported_with_htslib_lzma_enabled() {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "plexless-xz-fastq-{}-{nonce}.fastq.xz",
            std::process::id()
        ));
        let file = File::create(&path).unwrap();
        let mut writer = xz2::write::XzEncoder::new(file, 2);
        writer.write_all(b"@read\nACGT\n+\n!I~?\n").unwrap();
        writer.finish().unwrap();

        let mut reader = InputReader::try_open(&path).unwrap();
        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.id(), b"read");
        assert_eq!(record.seq().as_ref(), b"ACGT");
        assert_eq!(record.qual(), Some(b"!I~?".as_slice()));
        assert!(reader.next_record().is_none());
        std::fs::remove_file(path).unwrap();
    }
}
