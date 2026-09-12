use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};

use crate::cli::{DemuxArgs, OutputFormat, ReadMode};

#[derive(Debug)]
pub(crate) enum InputFiles {
    Single(PathBuf),
    Paired { r1: PathBuf, r2: PathBuf },
}

#[derive(Debug)]
pub(crate) enum InputSource {
    Fastq(InputFiles),
    Cram {
        path: PathBuf,
        mode: ReadMode,
        output_format: OutputFormat,
    },
}

impl InputSource {
    pub(crate) fn is_paired(&self) -> bool {
        match self {
            Self::Fastq(files) => files.is_paired(),
            Self::Cram { mode, .. } => *mode == ReadMode::Paired,
        }
    }

    pub(crate) fn output_format(&self) -> OutputFormat {
        match self {
            Self::Fastq(_) => OutputFormat::Fastq,
            Self::Cram { output_format, .. } => *output_format,
        }
    }

    pub(crate) fn has_gzip(&self) -> Result<bool, String> {
        match self {
            Self::Fastq(files) => files.has_gzip(),
            Self::Cram { .. } => Ok(false),
        }
    }
}

impl InputFiles {
    pub(crate) fn is_paired(&self) -> bool {
        matches!(self, Self::Paired { .. })
    }

    pub(crate) fn has_gzip(&self) -> Result<bool, String> {
        match self {
            Self::Single(path) => is_gzip_file(path),
            Self::Paired { r1, r2 } => Ok(is_gzip_file(r1)? || is_gzip_file(r2)?),
        }
    }
}

pub(crate) fn resolve_inputs(args: &DemuxArgs) -> Result<InputFiles, String> {
    if let Some(r1) = args.r1.as_deref() {
        let r2 = args.r2.as_deref().ok_or("Missing R2 input")?;
        validate_input_path(r1)?;
        validate_input_path(r2)?;
        Ok(InputFiles::Paired {
            r1: r1.to_path_buf(),
            r2: r2.to_path_buf(),
        })
    } else {
        let reads = args
            .reads
            .as_deref()
            .ok_or("Missing single-end FASTQ input")?;
        validate_input_path(reads)?;
        Ok(InputFiles::Single(reads.to_path_buf()))
    }
}

pub(crate) fn resolve_input_source(args: &DemuxArgs) -> Result<InputSource, String> {
    if let Some(path) = args.cram.as_deref() {
        validate_input_path(path)?;
        return Ok(InputSource::Cram {
            path: path.to_path_buf(),
            mode: args
                .read_mode
                .ok_or("CRAM input requires --read-mode single or paired")?,
            output_format: args
                .output_format
                .ok_or("CRAM input requires --output-format fastq or cram")?,
        });
    }
    resolve_inputs(args).map(InputSource::Fastq)
}

fn validate_input_path(path: &Path) -> Result<(), String> {
    fs::metadata(path)
        .map(|_| ())
        .map_err(|e| format!("Could not access input '{}': {e}", path.display()))
}

pub(crate) fn is_gzip_file(path: &Path) -> Result<bool, String> {
    let mut file = fs::File::open(path)
        .map_err(|e| format!("Could not inspect input '{}': {e}", path.display()))?;
    let mut magic = [0u8; 2];
    let read = file
        .read(&mut magic)
        .map_err(|e| format!("Could not inspect input '{}': {e}", path.display()))?;

    Ok(read == magic.len() && magic == [0x1f, 0x8b])
}
