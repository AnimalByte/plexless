use clap::{ArgGroup, Args, Parser, Subcommand};
use std::path::PathBuf;

pub use crate::output::OutputMode;

fn parse_positive_usize(value: &str) -> Result<usize, String> {
    let parsed = value
        .parse::<usize>()
        .map_err(|_| format!("invalid positive integer '{value}'"))?;

    if parsed == 0 {
        return Err("value must be at least 1".into());
    }

    Ok(parsed)
}

fn parse_fraction(value: &str) -> Result<f64, String> {
    let fraction = value
        .parse::<f64>()
        .map_err(|_| format!("invalid fraction '{value}'"))?;
    if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
        return Err("fraction must be between 0 and 1".into());
    }
    Ok(fraction)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ByteSizeSetting {
    Auto,
    Bytes(usize),
}

impl ByteSizeSetting {
    pub(crate) fn explicit_bytes(self) -> Option<usize> {
        match self {
            Self::Auto => None,
            Self::Bytes(bytes) => Some(bytes),
        }
    }
}

fn parse_byte_size(value: &str) -> Result<ByteSizeSetting, String> {
    let value = value.trim();
    if value.eq_ignore_ascii_case("auto") {
        return Ok(ByteSizeSetting::Auto);
    }
    if value.is_empty() {
        return Err("byte size cannot be empty".into());
    }

    let split = value
        .find(|character: char| !character.is_ascii_digit())
        .unwrap_or(value.len());
    let (number, suffix) = value.split_at(split);
    let amount = number
        .parse::<usize>()
        .map_err(|_| format!("invalid byte size '{value}'"))?;
    if amount == 0 {
        return Err("byte size must be greater than 0".into());
    }
    let multiplier = match suffix.to_ascii_lowercase().as_str() {
        "" | "b" => 1,
        "k" | "kb" | "kib" => 1024,
        "m" | "mb" | "mib" => 1024 * 1024,
        "g" | "gb" | "gib" => 1024 * 1024 * 1024,
        _ => return Err(format!("invalid byte-size suffix in '{value}'")),
    };
    amount
        .checked_mul(multiplier)
        .map(ByteSizeSetting::Bytes)
        .ok_or_else(|| format!("byte size '{value}' is too large"))
}

#[derive(Parser, Debug)]
#[command(name = "plexless")]
#[command(about = "FASTQ nested demultiplexing tool")]
pub struct Cli {
    /// Total CPU-work budget. plexless accounts for FASTQ parsing and
    /// dynamically shares remaining capacity between gzip decompression and
    /// demultiplexing/output compression. A value of 1 preserves the serial
    /// reference pipeline.
    #[arg(
        long,
        global = true,
        default_value_t = 1,
        value_parser = parse_positive_usize
    )]
    pub threads: usize,

    #[command(subcommand)]
    pub command: Command,
}

#[derive(Subcommand, Debug)]
pub enum Command {
    Demux(DemuxArgs),
}

#[derive(Args, Debug)]
#[command(
    group(
        ArgGroup::new("input")
            .required(true)
            .multiple(false)
            .args(["reads", "r1"])
    )
)]
pub struct DemuxArgs {
    /// Single-end FASTQ input
    #[arg(long, conflicts_with_all = ["r1", "r2"])]
    pub reads: Option<PathBuf>,

    /// Paired-end R1 FASTQ input
    #[arg(long, requires = "r2")]
    pub r1: Option<PathBuf>,

    /// Paired-end R2 FASTQ input
    #[arg(long, requires = "r1")]
    pub r2: Option<PathBuf>,

    /// Read structure for single-end input
    ///
    /// Example: `R1_10A11B4T`
    #[arg(
        long,
        requires = "reads",
        conflicts_with_all = ["r1_structure", "r2_structure"]
    )]
    pub structure: Option<String>,

    /// R1 read structure for paired-end input
    ///
    /// Example: `R1_10A11B4T`
    #[arg(long, requires = "r1", conflicts_with = "structure")]
    pub r1_structure: Option<String>,

    /// R2 read structure for paired-end input
    ///
    /// Example: `R2_8C`
    #[arg(long, requires = "r2", conflicts_with = "structure")]
    pub r2_structure: Option<String>,

    /// Barcode whitelist/catalog TSV
    ///
    /// Required columns: Set, ID, Sequence
    #[arg(long)]
    pub barcodes: PathBuf,

    /// Sample routing TSV
    ///
    /// First column must be Sample, followed by barcode sets
    /// in alphabetical order
    #[arg(long)]
    pub samples: PathBuf,

    /// Output directory
    #[arg(short, long)]
    pub output: PathBuf,

    /// Gzip compression level (0-9)
    #[arg(
        long,
        default_value_t = 2,
        value_parser = clap::value_parser!(u32).range(0..=9)
    )]
    pub compression_level: u32,

    /// Output strategy: auto, direct, or buffered
    #[arg(long, value_enum, default_value_t = OutputMode::Auto)]
    pub output_mode: OutputMode,

    /// Uncompressed bytes per gzip member (auto, or e.g. 256K, 1M)
    #[arg(long, default_value = "auto", value_parser = parse_byte_size)]
    pub output_chunk_size: ByteSizeSetting,

    /// Total memory budget for active output accumulators (auto, or e.g. 512M)
    #[arg(long, default_value = "auto", value_parser = parse_byte_size)]
    pub output_buffer_memory: ByteSizeSetting,

    /// Maximum simultaneously open output files (adaptive by default)
    #[arg(long, value_parser = parse_positive_usize)]
    pub max_open_files: Option<usize>,

    /// Maximum barcode mismatches allowed during correction
    #[arg(long, default_value_t = 1)]
    pub max_mismatches: u8,

    /// Write biological-read and barcode-region FASTQ statistics
    #[arg(long)]
    pub fastq_stats: bool,

    /// Write reads that cannot be assigned to a sample
    #[arg(long)]
    pub write_unassigned: bool,

    /// Warn when a populated sample falls below this fraction of the nonzero median
    #[arg(long, default_value_t = 0.05, value_parser = parse_fraction)]
    pub low_sample_fraction: f64,
}

impl DemuxArgs {
    pub fn validate(&self) -> Result<(), String> {
        if self.reads.is_some() && self.structure.is_none() {
            return Err("Single-end input requires --structure".into());
        }

        if self.r1.is_some() && self.r1_structure.is_none() && self.r2_structure.is_none() {
            return Err("Paired-end input requires at least one of \
                 --r1-structure or --r2-structure"
                .into());
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn byte_sizes_accept_auto_and_binary_suffixes() {
        assert_eq!(parse_byte_size("auto"), Ok(ByteSizeSetting::Auto));
        assert_eq!(
            parse_byte_size("256K"),
            Ok(ByteSizeSetting::Bytes(256 * 1024))
        );
        assert_eq!(
            parse_byte_size("1MiB"),
            Ok(ByteSizeSetting::Bytes(1024 * 1024))
        );
        assert_eq!(
            parse_byte_size("2G"),
            Ok(ByteSizeSetting::Bytes(2 * 1024 * 1024 * 1024))
        );
        assert!(parse_byte_size("12XB").is_err());
        assert!(parse_byte_size("0").is_err());
    }

    #[test]
    fn low_sample_fraction_is_bounded() {
        assert_eq!(parse_fraction("0.05"), Ok(0.05));
        assert!(parse_fraction("-0.1").is_err());
        assert!(parse_fraction("1.1").is_err());
    }
}
