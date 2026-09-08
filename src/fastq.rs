use needletail::errors::ParseError;
use needletail::parser::SequenceRecord;
use needletail::{FastxReader, parse_fastx_file, parse_fastx_reader};
use std::io::Read;
use std::path::Path;

pub struct InputReader {
    reader: Box<dyn FastxReader>,
}

impl InputReader {
    pub fn open(path: &Path) -> Self {
        let reader = parse_fastx_file(path).expect("Could not open FASTQ file");

        Self { reader }
    }

    pub(crate) fn try_open(path: &Path) -> Result<Self, String> {
        let reader = parse_fastx_file(path)
            .map_err(|e| format!("Could not open FASTQ input '{}': {e}", path.display()))?;

        Ok(Self { reader })
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
