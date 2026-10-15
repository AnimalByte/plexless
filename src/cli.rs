use clap::{ArgGroup, Args, Parser, Subcommand};
use std::path::PathBuf;

pub use crate::output::OutputMode;

#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum ReadMode {
    Single,
    Paired,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum OutputFormat {
    Fastq,
    Cram,
}

pub(crate) fn parse_positive_usize(value: &str) -> Result<usize, String> {
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
#[command(version)]
#[command(name = "plexless")]
#[command(about = "Demultiplex sequencing reads or clean already-demultiplexed FASTQ")]
#[command(after_help = "Start here:
  plexless demux -h   Assign reads to samples using barcode/sample TSVs
  plexless trim -h    Clean FASTQ without barcode or sample files

Use -h for common options, or --help for all options and advanced settings.")]
pub struct Cli {
    /// Total CPU budget shared by input, processing, and output
    #[arg(
        long,
        global = true,
        default_value_t = 1,
        value_parser = parse_positive_usize
    )]
    #[arg(value_name = "N")]
    pub threads: usize,

    #[command(subcommand)]
    pub command: Command,
}

#[derive(Subcommand, Debug)]
pub enum Command {
    /// Demultiplex FASTQ or unmapped CRAM using barcode and sample TSVs
    #[command(after_help = r#"Example:
  plexless demux --reads reads.fq.gz --structure R1_10A11B4T \
    --barcodes barcodes.tsv --samples samples.tsv --output demux --threads 8

Defaults: declared prefixes are removed, adapter matching is on, and filters are off.
Use -h for common options, or --help for all options."#)]
    Demux(DemuxArgs),
    /// Trim and optionally filter already-demultiplexed FASTQ
    #[command(after_help = r#"Examples:
  plexless trim --r1 R1.fq.gz --r2 R2.fq.gz --fastp-defaults --output cleaned
  plexless trim --reads t7.fq.gz --r1-distal-sequence "$ADAPTER" \
    --partial-pass-through --output trimmed

Defaults: paired overlap is on; proximal clipping and quality filters are off.
--fastp-defaults enables quality/length filters and JSON/HTML reports.
Use -h for common options, or --help for all options."#)]
    Trim(TrimArgs),
}

#[derive(Args, Debug)]
#[command(
    group(
        ArgGroup::new("input")
            .required(true)
            .multiple(false)
            .args(["reads", "r1", "cram", "stdin"])
    )
)]
pub struct DemuxArgs {
    /// Single-end FASTQ input
    #[arg(long, conflicts_with_all = ["r1", "r2"])]
    #[arg(help_heading = "Input", display_order = 100, value_name = "FASTQ")]
    pub reads: Option<PathBuf>,

    /// Paired-end R1 FASTQ input
    #[arg(long, requires = "r2")]
    #[arg(help_heading = "Input", display_order = 101, value_name = "FASTQ")]
    pub r1: Option<PathBuf>,

    /// Paired-end R2 FASTQ input
    #[arg(long, requires = "r1")]
    #[arg(help_heading = "Input", display_order = 102, value_name = "FASTQ")]
    pub r2: Option<PathBuf>,

    /// Unmapped CRAM input (requires --read-mode and --output-format)
    #[arg(long, conflicts_with_all = ["reads", "r1", "r2"])]
    #[arg(
        help_heading = "Unmapped CRAM",
        display_order = 100,
        hide_short_help = true,
        value_name = "CRAM"
    )]
    pub cram: Option<PathBuf>,

    /// Record interpretation for CRAM input; never inferred from records
    #[arg(long, value_enum, requires = "cram")]
    #[arg(
        help_heading = "Unmapped CRAM",
        display_order = 101,
        hide_short_help = true
    )]
    pub read_mode: Option<ReadMode>,

    /// Output format for CRAM input; FASTQ input retains its existing output
    #[arg(long, value_enum, requires = "cram")]
    #[arg(
        help_heading = "Unmapped CRAM",
        display_order = 102,
        hide_short_help = true
    )]
    pub output_format: Option<OutputFormat>,

    /// Barcode/technical layout at the beginning of a single-end read
    ///
    /// Example: `R1_10A11B4T`
    #[arg(
        long,
        conflicts_with_all = ["r1_structure", "r2_structure"]
    )]
    #[arg(
        help_heading = "Barcode layout",
        display_order = 100,
        value_name = "LAYOUT"
    )]
    pub structure: Option<String>,

    /// Barcode/technical layout at the beginning of R1
    ///
    /// Example: `R1_10A11B4T`
    #[arg(long, conflicts_with = "structure")]
    #[arg(
        help_heading = "Barcode layout",
        display_order = 101,
        value_name = "LAYOUT"
    )]
    pub r1_structure: Option<String>,

    /// Barcode/technical layout at the beginning of R2
    ///
    /// Example: `R2_8C`
    #[arg(long, conflicts_with = "structure")]
    #[arg(
        help_heading = "Barcode layout",
        display_order = 102,
        value_name = "LAYOUT"
    )]
    pub r2_structure: Option<String>,

    /// Barcode whitelist/catalog TSV
    ///
    /// Required columns: Set, ID, Sequence
    #[arg(long)]
    #[arg(
        help_heading = "Barcode layout",
        display_order = 103,
        value_name = "TSV"
    )]
    pub barcodes: PathBuf,

    /// Sample routing TSV
    ///
    /// First column must be Sample, followed by barcode sets
    /// in alphabetical order
    #[arg(long)]
    #[arg(
        help_heading = "Barcode layout",
        display_order = 104,
        value_name = "TSV"
    )]
    pub samples: PathBuf,

    /// Output directory
    #[arg(short, long)]
    #[arg(help_heading = "Output", display_order = 100, value_name = "DIR")]
    pub output: PathBuf,

    /// Gzip compression level (0-9)
    #[arg(
        long,
        default_value_t = 2,
        value_parser = clap::value_parser!(u32).range(0..=9)
    )]
    #[arg(help_heading = "Output", display_order = 101, hide_short_help = true)]
    pub compression_level: u32,

    /// Output strategy: auto, direct, or buffered
    #[arg(long, value_enum, default_value_t = OutputMode::Auto)]
    #[arg(
        help_heading = "Advanced performance",
        display_order = 100,
        hide_short_help = true
    )]
    pub output_mode: OutputMode,

    /// Uncompressed bytes per gzip member (auto uses 32 KiB for FASTQ input)
    #[arg(long, default_value = "auto", value_parser = parse_byte_size)]
    #[arg(
        help_heading = "Advanced performance",
        display_order = 101,
        hide_short_help = true
    )]
    pub output_chunk_size: ByteSizeSetting,

    /// Total memory budget for active output accumulators (auto, or e.g. 512M)
    #[arg(long, default_value = "auto", value_parser = parse_byte_size)]
    #[arg(
        help_heading = "Advanced performance",
        display_order = 102,
        hide_short_help = true
    )]
    pub output_buffer_memory: ByteSizeSetting,

    /// FASTQ writer-cache limit; CRAM simultaneous-writer minimum (adaptive by default)
    #[arg(long, value_parser = parse_positive_usize)]
    #[arg(
        help_heading = "Advanced performance",
        display_order = 103,
        hide_short_help = true
    )]
    pub max_open_files: Option<usize>,

    /// Barcode mismatch budget (0 through barcode length); strict route-local safety checks apply
    #[arg(long, default_value_t = 1)]
    #[arg(
        help_heading = "Barcode layout",
        display_order = 105,
        hide_short_help = true
    )]
    #[arg(value_name = "N")]
    pub max_mismatches: u8,

    #[command(flatten)]
    pub preprocessing: PreprocessingArgs,

    /// Disable adapter trimming; declared proximal clipping still applies
    #[arg(long)]
    #[arg(
        help_heading = "Adapters and overlap",
        display_order = 104,
        hide_short_help = true
    )]
    pub no_pass_through: bool,

    /// Distal Hamming mismatch budget outside the exact anchor (independent of barcode correction)
    #[arg(long, default_value_t = 1)]
    #[arg(
        help_heading = "Adapters and overlap",
        display_order = 105,
        hide_short_help = true
    )]
    #[arg(value_name = "N")]
    pub pass_through_mismatches: usize,

    /// Override R1 distal structure in observed orientation, e.g. R1_12A(rc)4T
    #[arg(long, conflicts_with = "no_pass_through")]
    #[arg(
        help_heading = "Adapters and overlap",
        display_order = 106,
        hide_short_help = true,
        value_name = "LAYOUT"
    )]
    pub r1_distal_structure: Option<String>,

    /// Override R2 distal structure in observed orientation; requires paired input
    #[arg(long, conflicts_with = "no_pass_through")]
    #[arg(
        help_heading = "Adapters and overlap",
        display_order = 107,
        hide_short_help = true,
        value_name = "LAYOUT"
    )]
    pub r2_distal_structure: Option<String>,

    /// Write input, barcode-region, and retained-read QC tables
    #[arg(long)]
    #[arg(help_heading = "Reports", display_order = 100)]
    pub fastq_stats: bool,

    /// Write reads that cannot be assigned to a sample
    #[arg(long)]
    #[arg(help_heading = "Stream and auxiliary output", display_order = 115)]
    pub write_unassigned: bool,

    /// Warn when a populated sample falls below this fraction of the nonzero median
    #[arg(long, default_value_t = 0.05, value_parser = parse_fraction)]
    #[arg(help_heading = "Reports", display_order = 108, hide_short_help = true)]
    pub low_sample_fraction: f64,
}

impl DemuxArgs {
    pub fn validate(&self) -> Result<(), String> {
        if self.cram.is_some() {
            let mode = self
                .read_mode
                .ok_or("CRAM input requires --read-mode single or --read-mode paired")?;
            self.output_format
                .ok_or("CRAM input requires --output-format fastq or --output-format cram")?;
            match mode {
                ReadMode::Single if self.structure.is_none() => {
                    return Err("Single-end CRAM input requires --structure".into());
                }
                ReadMode::Paired if self.r1_structure.is_none() && self.r2_structure.is_none() => {
                    return Err("Paired-end CRAM input requires at least one of \
                         --r1-structure or --r2-structure"
                        .into());
                }
                _ => {}
            }
        } else {
            if self.read_mode.is_some() || self.output_format.is_some() {
                return Err("--read-mode and --output-format are only valid with --cram; FASTQ to CRAM is not supported".into());
            }
            if self.reads.is_some()
                && !self.preprocessing.fastp.interleaved_in
                && self.structure.is_none()
            {
                return Err("Single-end input requires --structure".into());
            }

            if (self.r1.is_some() || self.preprocessing.fastp.interleaved_in)
                && self.r1_structure.is_none()
                && self.r2_structure.is_none()
            {
                return Err("Paired-end input requires at least one of \
                     --r1-structure or --r2-structure"
                    .into());
            }
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
        assert!(parse_byte_size("18446744073709551615G").is_err());
        assert!(parse_byte_size("12XB").is_err());
        assert!(parse_byte_size("0").is_err());
    }

    #[test]
    fn low_sample_fraction_is_bounded() {
        assert_eq!(parse_fraction("0.05"), Ok(0.05));
        assert!(parse_fraction("-0.1").is_err());
        assert!(parse_fraction("1.1").is_err());
    }

    fn demux_from(arguments: &[&str]) -> DemuxArgs {
        let cli = Cli::try_parse_from(arguments).unwrap();
        match cli.command {
            Command::Demux(args) => args,
            Command::Trim(_) => panic!("expected demux"),
        }
    }

    #[test]
    fn fastq_cli_does_not_require_new_cram_options() {
        let args = demux_from(&[
            "plexless",
            "demux",
            "--reads",
            "reads.fastq.gz",
            "--structure",
            "R1_4A",
            "--barcodes",
            "barcodes.tsv",
            "--samples",
            "samples.tsv",
            "--output",
            "out",
        ]);
        assert_eq!(args.output_format, None);
        assert_eq!(args.read_mode, None);
        assert!(args.validate().is_ok());
    }

    #[test]
    fn cram_cli_requires_explicit_mode_and_output_format() {
        let args = demux_from(&[
            "plexless",
            "demux",
            "--cram",
            "reads.cram",
            "--structure",
            "R1_4A",
            "--barcodes",
            "barcodes.tsv",
            "--samples",
            "samples.tsv",
            "--output",
            "out",
        ]);
        assert!(args.validate().unwrap_err().contains("--read-mode"));

        let args = demux_from(&[
            "plexless",
            "demux",
            "--cram",
            "reads.cram",
            "--read-mode",
            "single",
            "--structure",
            "R1_4A",
            "--barcodes",
            "barcodes.tsv",
            "--samples",
            "samples.tsv",
            "--output",
            "out",
        ]);
        assert!(args.validate().unwrap_err().contains("--output-format"));
    }

    #[test]
    fn output_format_cannot_enable_fastq_to_cram() {
        let args = demux_from(&[
            "plexless",
            "demux",
            "--reads",
            "reads.fastq.gz",
            "--structure",
            "R1_4A",
            "--output-format",
            "cram",
            "--barcodes",
            "barcodes.tsv",
            "--samples",
            "samples.tsv",
            "--output",
            "out",
        ]);
        assert!(args.validate().unwrap_err().contains("FASTQ to CRAM"));
    }
}

// Shared, opt-in preprocessing controls. Quality scores are Phred+33.
#[derive(Args, Debug, Clone)]
pub struct PreprocessingArgs {
    #[command(flatten)]
    pub fastp: crate::fastp_options::FastpOptions,
    /// TSV with Mate and Sequence columns; observed proximal T bases per mate
    #[arg(long)]
    #[arg(
        help_heading = "Adapters and overlap",
        display_order = 103,
        hide_short_help = true,
        value_name = "TSV"
    )]
    pub technical_sequences: Option<PathBuf>,
    /// R1 3-prime adapter sequence, as it appears in the read (ACGT only)
    #[arg(long, conflicts_with_all = ["r1_distal_structure", "no_pass_through"])]
    #[arg(
        help_heading = "Adapters and overlap",
        display_order = 100,
        value_name = "SEQUENCE"
    )]
    pub r1_distal_sequence: Option<String>,
    /// R2 3-prime adapter sequence, as it appears in R2 (paired input only)
    #[arg(long, conflicts_with_all = ["r2_distal_structure", "no_pass_through"])]
    #[arg(
        help_heading = "Adapters and overlap",
        display_order = 101,
        value_name = "SEQUENCE"
    )]
    pub r2_distal_sequence: Option<String>,
    /// Allow exact adapter-prefix matches at the read end (20 known bases by default)
    #[arg(long, conflicts_with = "no_pass_through")]
    #[arg(help_heading = "Adapters and overlap", display_order = 102)]
    pub partial_pass_through: bool,
    /// Minimum known bases for partial matching (masked T bases do not count)
    #[arg(long, default_value_t = 20, value_parser = parse_positive_usize)]
    #[arg(
        help_heading = "Adapters and overlap",
        display_order = 108,
        hide_short_help = true
    )]
    #[arg(value_name = "BASES")]
    pub min_partial_adapter_overlap: usize,
    /// Minimum informative known positions for full adapter verification (masked T does not count)
    #[arg(long, default_value_t = 20, value_parser = parse_positive_usize)]
    #[arg(
        help_heading = "Adapters and overlap",
        display_order = 109,
        hide_short_help = true
    )]
    #[arg(value_name = "BASES")]
    pub min_full_adapter_known_bases: usize,
    /// Disable paired overlap adapter clipping; explicit merge/correction remain independent
    #[arg(long)]
    #[arg(
        help_heading = "Adapters and overlap",
        display_order = 110,
        hide_short_help = true
    )]
    pub no_overlap: bool,
    /// Minimum biological mate overlap for adapter clipping, merge, and correction
    #[arg(long, alias = "overlap_len_require", default_value_t = 30, value_parser = parse_positive_usize)]
    #[arg(
        help_heading = "Adapters and overlap",
        display_order = 111,
        hide_short_help = true
    )]
    #[arg(value_name = "BASES")]
    pub overlap_len_require: usize,
    /// Maximum mismatches across the entire mate overlap
    #[arg(long, alias = "overlap_diff_limit", default_value_t = 5)]
    #[arg(
        help_heading = "Adapters and overlap",
        display_order = 112,
        hide_short_help = true
    )]
    #[arg(value_name = "N")]
    pub overlap_diff_limit: usize,
    /// Maximum mismatch percentage across the entire mate overlap
    #[arg(long, alias = "overlap_diff_percent_limit", default_value_t = 20.0, value_parser = parse_percent)]
    #[arg(
        help_heading = "Adapters and overlap",
        display_order = 113,
        hide_short_help = true
    )]
    #[arg(value_name = "PERCENT")]
    pub overlap_diff_percent_limit: f64,
    /// Phred threshold used by --max-low-quality-percent (default: 15)
    #[arg(long, value_parser = clap::value_parser!(u8).range(0..=93))]
    #[arg(
        help_heading = "Filtering",
        display_order = 105,
        hide_short_help = true
    )]
    #[arg(value_name = "Q")]
    pub qualified_quality: Option<u8>,
    /// Reject if percent of bases below qualified quality exceeds this value
    #[arg(long, value_parser = parse_percent)]
    #[arg(help_heading = "Filtering", display_order = 104)]
    #[arg(value_name = "PERCENT")]
    pub max_low_quality_percent: Option<f64>,
    /// Reject if retained N count exceeds this value
    #[arg(long)]
    #[arg(help_heading = "Filtering", display_order = 102)]
    #[arg(value_name = "COUNT")]
    pub max_n: Option<usize>,
    /// Reject if retained length is below this value
    #[arg(long)]
    #[arg(help_heading = "Filtering", display_order = 101)]
    #[arg(value_name = "BASES")]
    pub min_length: Option<usize>,
    /// Reject if retained mean Phred is below this value (empty reads have mean 0)
    #[arg(long, value_parser = parse_quality)]
    #[arg(help_heading = "Filtering", display_order = 103)]
    #[arg(value_name = "Q")]
    pub min_average_quality: Option<f64>,
}

impl Default for PreprocessingArgs {
    fn default() -> Self {
        Self {
            fastp: Default::default(),
            technical_sequences: None,
            r1_distal_sequence: None,
            r2_distal_sequence: None,
            partial_pass_through: false,
            min_partial_adapter_overlap: 20,
            min_full_adapter_known_bases: 20,
            no_overlap: false,
            overlap_len_require: 30,
            overlap_diff_limit: 5,
            overlap_diff_percent_limit: 20.0,
            qualified_quality: None,
            max_low_quality_percent: None,
            max_n: None,
            min_length: None,
            min_average_quality: None,
        }
    }
}

pub(crate) fn parse_percent(value: &str) -> Result<f64, String> {
    let value: f64 = value.parse().map_err(|_| "invalid percentage")?;
    if value.is_finite() && (0.0..=100.0).contains(&value) {
        Ok(value)
    } else {
        Err("percentage must be between 0 and 100".into())
    }
}
fn parse_quality(value: &str) -> Result<f64, String> {
    let value: f64 = value.parse().map_err(|_| "invalid quality")?;
    if value.is_finite() && (0.0..=93.0).contains(&value) {
        Ok(value)
    } else {
        Err("quality must be between 0 and 93".into())
    }
}

#[derive(Args, Debug)]
#[command(group(ArgGroup::new("input").required(true).multiple(false).args(["reads", "r1", "stdin"])))]
pub struct TrimArgs {
    /// Single-end FASTQ input
    #[arg(long, conflicts_with_all = ["r1", "r2"])]
    #[arg(help_heading = "Input", display_order = 100, value_name = "FASTQ")]
    pub reads: Option<PathBuf>,
    /// Paired-end R1 FASTQ input
    #[arg(long, requires = "r2")]
    #[arg(help_heading = "Input", display_order = 101, value_name = "FASTQ")]
    pub r1: Option<PathBuf>,
    /// Paired-end R2 FASTQ input
    #[arg(long, requires = "r1")]
    #[arg(help_heading = "Input", display_order = 102, value_name = "FASTQ")]
    pub r2: Option<PathBuf>,
    /// Explicit proximal clipping for SE input; omitted by default
    #[arg(long, conflicts_with_all = ["r1_structure", "r2_structure"])]
    #[arg(
        help_heading = "Proximal clipping",
        display_order = 100,
        hide_short_help = true,
        value_name = "LAYOUT"
    )]
    pub structure: Option<String>,
    /// Explicit proximal clipping for R1; omitted by default
    #[arg(long, conflicts_with = "structure")]
    #[arg(
        help_heading = "Proximal clipping",
        display_order = 101,
        hide_short_help = true,
        value_name = "LAYOUT"
    )]
    pub r1_structure: Option<String>,
    /// Explicit proximal clipping for R2; omitted by default
    #[arg(long, conflicts_with = "structure")]
    #[arg(
        help_heading = "Proximal clipping",
        display_order = 102,
        hide_short_help = true,
        value_name = "LAYOUT"
    )]
    pub r2_structure: Option<String>,
    /// Output directory: trimmed.fastq.gz or trimmed_R1/R2.fastq.gz
    #[arg(short, long)]
    #[arg(help_heading = "Output", display_order = 100, value_name = "DIR")]
    pub output: PathBuf,
    /// Gzip compression level (0-9)
    #[arg(long, default_value_t = 2, value_parser = clap::value_parser!(u32).range(0..=9))]
    #[arg(help_heading = "Output", display_order = 101, hide_short_help = true)]
    pub compression_level: u32,
    /// Full adapter mismatch budget; the 12-16-base anchor must remain exact
    #[arg(long, default_value_t = 1)]
    #[arg(
        help_heading = "Adapters and overlap",
        display_order = 105,
        hide_short_help = true
    )]
    #[arg(value_name = "N")]
    pub pass_through_mismatches: usize,
    /// Disable adapter trimming; explicit proximal clipping still applies
    #[arg(long)]
    #[arg(
        help_heading = "Adapters and overlap",
        display_order = 104,
        hide_short_help = true
    )]
    pub no_pass_through: bool,
    /// Distal structure with T segments only; barcode-dependent structures require demux
    #[arg(long, conflicts_with = "no_pass_through")]
    #[arg(
        help_heading = "Adapters and overlap",
        display_order = 106,
        hide_short_help = true,
        value_name = "LAYOUT"
    )]
    pub r1_distal_structure: Option<String>,
    /// T-only R2 distal structure; requires paired input and known technical bases
    #[arg(long, conflicts_with = "no_pass_through")]
    #[arg(
        help_heading = "Adapters and overlap",
        display_order = 107,
        hide_short_help = true,
        value_name = "LAYOUT"
    )]
    pub r2_distal_structure: Option<String>,
    #[command(flatten)]
    pub preprocessing: PreprocessingArgs,
    /// Write QC tables before and after trimming
    #[arg(long)]
    #[arg(help_heading = "Reports", display_order = 100)]
    pub fastq_stats: bool,
}

impl TrimArgs {
    pub(crate) fn into_pipeline_args(self) -> DemuxArgs {
        DemuxArgs {
            reads: self.reads,
            r1: self.r1,
            r2: self.r2,
            cram: None,
            read_mode: None,
            output_format: None,
            structure: self.structure,
            r1_structure: self.r1_structure,
            r2_structure: self.r2_structure,
            barcodes: PathBuf::new(),
            samples: PathBuf::new(),
            output: self.output,
            compression_level: self.compression_level,
            output_mode: OutputMode::Auto,
            output_chunk_size: ByteSizeSetting::Auto,
            output_buffer_memory: ByteSizeSetting::Auto,
            max_open_files: None,
            max_mismatches: 1,
            preprocessing: self.preprocessing,
            no_pass_through: self.no_pass_through,
            pass_through_mismatches: self.pass_through_mismatches,
            r1_distal_structure: self.r1_distal_structure,
            r2_distal_structure: self.r2_distal_structure,
            fastq_stats: self.fastq_stats,
            write_unassigned: false,
            low_sample_fraction: 0.05,
        }
    }
}
