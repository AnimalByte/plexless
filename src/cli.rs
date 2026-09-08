use clap::{ArgGroup, Args, Parser, Subcommand};
use std::path::PathBuf;

fn parse_positive_usize(value: &str) -> Result<usize, String> {
    let threads = value
        .parse::<usize>()
        .map_err(|_| format!("invalid thread count '{value}'"))?;

    if threads == 0 {
        return Err("thread count must be at least 1".into());
    }

    Ok(threads)
}

#[derive(Parser, Debug)]
#[command(name = "simplex")]
#[command(about = "FASTQ nested demultiplexing tool")]
pub struct Cli {
    /// Total CPU-work budget. simplex accounts for FASTQ parsing and
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

    /// Maximum barcode mismatches allowed during correction
    #[arg(long, default_value_t = 1)]
    pub max_mismatches: u8,

    /// Calculate FASTQ quality statistics
    #[arg(long)]
    pub fastq_stats: bool,

    /// Write reads that cannot be assigned to a sample
    #[arg(long)]
    pub write_unassigned: bool,
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
