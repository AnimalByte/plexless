//! Configurable fastp-derived features. Original demux retention stays opt-in.
use crate::cli::{parse_percent, parse_positive_usize};
use clap::{Args, ValueEnum};
use std::path::PathBuf;

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum UmiLocation {
    Index1,
    Index2,
    Read1,
    Read2,
    #[value(alias = "per_index")]
    PerIndex,
    #[value(alias = "per_read")]
    PerRead,
}

#[derive(Debug, Clone, Args)]
pub struct FastpOptions {
    /// Enable quality/length filters, adapter discovery, and JSON/HTML reports
    ///
    /// Uses Q15, at most 40% low-quality bases, at most 5 Ns, and length >=15.
    /// Enables SE adapter discovery, instrument-based polyG detection for files,
    /// and duplicate estimation. Quality-window clipping, merging, correction,
    /// UMI extraction and duplicate removal still require their own flags.
    #[arg(long)]
    #[arg(help_heading = "Filtering", display_order = 100)]
    pub fastp_defaults: bool,
    /// Require sequence evidence as well as a unique paired overlap before adapter clipping
    #[arg(long)]
    #[arg(
        help_heading = "Adapters and overlap",
        display_order = 115,
        hide_short_help = true
    )]
    pub require_adapter_for_overlap: bool,
    /// Count adapter-supported fragments at or below this retained length as dimers
    #[arg(long, alias = "dimer_max_len", default_value_t = 2)]
    #[arg(help_heading = "Reports", display_order = 110, hide_short_help = true)]
    pub dimer_max_len: usize,
    /// Remove low-quality front windows until a passing window is reached
    #[arg(long, alias = "cut_front")]
    #[arg(
        help_heading = "Quality trimming",
        display_order = 100,
        hide_short_help = true
    )]
    pub cut_front: bool,
    /// Remove low-quality tail windows until a passing window is reached
    #[arg(long, alias = "cut_tail")]
    #[arg(
        help_heading = "Quality trimming",
        display_order = 101,
        hide_short_help = true
    )]
    pub cut_tail: bool,
    /// Clip from the first failing quality window to the read end
    #[arg(long, alias = "cut_right")]
    #[arg(
        help_heading = "Quality trimming",
        display_order = 102,
        hide_short_help = true
    )]
    pub cut_right: bool,
    /// Shared quality-window size; overridden by per-direction settings
    #[arg(long,alias="cut_window_size",default_value_t=4,value_parser=clap::value_parser!(u16).range(1..=1000))]
    #[arg(
        help_heading = "Quality trimming",
        display_order = 103,
        hide_short_help = true
    )]
    pub cut_window_size: u16,
    /// Shared minimum mean Phred quality for a passing window
    #[arg(long,alias="cut_mean_quality",default_value_t=20,value_parser=clap::value_parser!(u8).range(0..=93))]
    #[arg(
        help_heading = "Quality trimming",
        display_order = 104,
        hide_short_help = true
    )]
    pub cut_mean_quality: u8,
    /// Override front quality-window size
    #[arg(long,alias="cut_front_window_size",value_parser=clap::value_parser!(u16).range(1..=1000))]
    #[arg(
        help_heading = "Quality trimming",
        display_order = 105,
        hide_short_help = true
    )]
    pub cut_front_window_size: Option<u16>,
    /// Override tail quality-window size
    #[arg(long,alias="cut_tail_window_size",value_parser=clap::value_parser!(u16).range(1..=1000))]
    #[arg(
        help_heading = "Quality trimming",
        display_order = 106,
        hide_short_help = true
    )]
    pub cut_tail_window_size: Option<u16>,
    /// Override right quality-window size
    #[arg(long,alias="cut_right_window_size",value_parser=clap::value_parser!(u16).range(1..=1000))]
    #[arg(
        help_heading = "Quality trimming",
        display_order = 107,
        hide_short_help = true
    )]
    pub cut_right_window_size: Option<u16>,
    /// Override front window mean Phred threshold
    #[arg(long,alias="cut_front_mean_quality",value_parser=clap::value_parser!(u8).range(0..=93))]
    #[arg(
        help_heading = "Quality trimming",
        display_order = 108,
        hide_short_help = true
    )]
    pub cut_front_mean_quality: Option<u8>,
    /// Override tail window mean Phred threshold
    #[arg(long,alias="cut_tail_mean_quality",value_parser=clap::value_parser!(u8).range(0..=93))]
    #[arg(
        help_heading = "Quality trimming",
        display_order = 109,
        hide_short_help = true
    )]
    pub cut_tail_mean_quality: Option<u8>,
    /// Override right window mean Phred threshold
    #[arg(long,alias="cut_right_mean_quality",value_parser=clap::value_parser!(u8).range(0..=93))]
    #[arg(
        help_heading = "Quality trimming",
        display_order = 110,
        hide_short_help = true
    )]
    pub cut_right_mean_quality: Option<u8>,
    /// Fixed R1 front clipping after proximal structure and read UMI
    #[arg(long, alias = "trim_front1", default_value_t = 0)]
    #[arg(
        help_heading = "Quality trimming",
        display_order = 111,
        hide_short_help = true
    )]
    pub trim_front1: usize,
    /// Fixed R2 front clipping (inherits R1 when omitted)
    #[arg(long, alias = "trim_front2")]
    #[arg(
        help_heading = "Quality trimming",
        display_order = 112,
        hide_short_help = true
    )]
    pub trim_front2: Option<usize>,
    /// Fixed R1 tail clipping
    #[arg(long, alias = "trim_tail1", default_value_t = 0)]
    #[arg(
        help_heading = "Quality trimming",
        display_order = 113,
        hide_short_help = true
    )]
    pub trim_tail1: usize,
    /// Fixed R2 tail clipping (inherits R1)
    #[arg(long, alias = "trim_tail2")]
    #[arg(
        help_heading = "Quality trimming",
        display_order = 114,
        hide_short_help = true
    )]
    pub trim_tail2: Option<usize>,
    /// Clip retained R1 to this maximum length
    #[arg(long,alias="max_len1",value_parser=parse_positive_usize)]
    #[arg(
        help_heading = "Quality trimming",
        display_order = 115,
        hide_short_help = true
    )]
    pub max_len1: Option<usize>,
    /// Clip retained R2 to this maximum length (inherits R1)
    #[arg(long,alias="max_len2",value_parser=parse_positive_usize)]
    #[arg(
        help_heading = "Quality trimming",
        display_order = 116,
        hide_short_help = true
    )]
    pub max_len2: Option<usize>,
    /// Enable bounded-mismatch polyG tail removal
    #[arg(long, alias = "trim_poly_g", conflicts_with = "disable_trim_poly_g")]
    #[arg(
        help_heading = "Quality trimming",
        display_order = 117,
        hide_short_help = true
    )]
    pub trim_poly_g: bool,
    /// Disable automatic profile polyG detection
    #[arg(long, alias = "disable_trim_poly_g")]
    #[arg(
        help_heading = "Quality trimming",
        display_order = 118,
        hide_short_help = true
    )]
    pub disable_trim_poly_g: bool,
    /// Minimum polyG tail length
    #[arg(long,alias="poly_g_min_len",default_value_t=10,value_parser=parse_positive_usize)]
    #[arg(
        help_heading = "Quality trimming",
        display_order = 119,
        hide_short_help = true
    )]
    pub poly_g_min_len: usize,
    /// Enable bounded-mismatch homopolymer tail removal
    #[arg(long, alias = "trim_poly_x")]
    #[arg(
        help_heading = "Quality trimming",
        display_order = 120,
        hide_short_help = true
    )]
    pub trim_poly_x: bool,
    /// Minimum polyX tail length
    #[arg(long,alias="poly_x_min_len",default_value_t=10,value_parser=parse_positive_usize)]
    #[arg(
        help_heading = "Quality trimming",
        display_order = 121,
        hide_short_help = true
    )]
    pub poly_x_min_len: usize,
    /// Reject low adjacent-base diversity after trimming
    #[arg(long, alias = "low_complexity_filter")]
    #[arg(
        help_heading = "Filtering",
        display_order = 107,
        hide_short_help = true
    )]
    pub low_complexity_filter: bool,
    /// Minimum percent of adjacent retained bases that differ
    #[arg(long,alias="complexity_threshold",default_value_t=30.0,value_parser=parse_percent)]
    #[arg(
        help_heading = "Filtering",
        display_order = 108,
        hide_short_help = true
    )]
    pub complexity_threshold: f64,
    /// Reject retained reads longer than this limit
    #[arg(long,alias="length_limit",value_parser=parse_positive_usize)]
    #[arg(
        help_heading = "Filtering",
        display_order = 106,
        hide_short_help = true
    )]
    pub length_limit: Option<usize>,
    /// Disable quality/N filters, including profile defaults
    #[arg(long, alias = "disable_quality_filtering")]
    #[arg(
        help_heading = "Filtering",
        display_order = 109,
        hide_short_help = true
    )]
    pub disable_quality_filtering: bool,
    /// Disable minimum/maximum length filters
    #[arg(long, alias = "disable_length_filtering")]
    #[arg(
        help_heading = "Filtering",
        display_order = 110,
        hide_short_help = true
    )]
    pub disable_length_filtering: bool,
    /// Validate Phred+64 input and convert output qualities to Phred+33
    #[arg(long, alias = "phred64")]
    #[arg(
        help_heading = "Stream and auxiliary output",
        display_order = 114,
        hide_short_help = true
    )]
    pub phred64: bool,
    /// Correct an overlap mismatch only with Q30 support against a base of Q14 or lower
    #[arg(long)]
    #[arg(
        help_heading = "Merging and correction",
        display_order = 103,
        hide_short_help = true
    )]
    pub correction: bool,
    /// Extract UMI into the header before filtering
    #[arg(long, requires = "umi_loc")]
    #[arg(
        help_heading = "UMI extraction",
        display_order = 100,
        hide_short_help = true
    )]
    pub umi: bool,
    /// Read or CASAVA index location of the UMI
    #[arg(long, alias = "umi_loc", requires = "umi", value_enum)]
    #[arg(
        help_heading = "UMI extraction",
        display_order = 101,
        hide_short_help = true
    )]
    pub umi_loc: Option<UmiLocation>,
    /// Read-based UMI length
    #[arg(long, alias = "umi_len", default_value_t = 0)]
    #[arg(
        help_heading = "UMI extraction",
        display_order = 102,
        hide_short_help = true
    )]
    pub umi_len: usize,
    /// Additional bases to skip after a read-based UMI
    #[arg(long, alias = "umi_skip", default_value_t = 0)]
    #[arg(
        help_heading = "UMI extraction",
        display_order = 103,
        hide_short_help = true
    )]
    pub umi_skip: usize,
    /// Optional prefix before the extracted UMI
    #[arg(long, alias = "umi_prefix", default_value = "")]
    #[arg(
        help_heading = "UMI extraction",
        display_order = 104,
        hide_short_help = true
    )]
    pub umi_prefix: String,
    /// Delimiter between the read ID and UMI tag
    #[arg(long, alias = "umi_delim", default_value = ":")]
    #[arg(
        help_heading = "UMI extraction",
        display_order = 105,
        hide_short_help = true
    )]
    pub umi_delim: String,
    /// Remove exact duplicate original sequences/pairs; retain first input occurrence
    #[arg(long)]
    #[arg(
        help_heading = "Duplicate removal",
        display_order = 100,
        hide_short_help = true
    )]
    pub dedup: bool,
    /// Estimate original-fragment duplication without discarding reads
    #[arg(long)]
    #[arg(help_heading = "Reports", display_order = 104, hide_short_help = true)]
    pub evaluate_duplication: bool,
    /// Disable profile duplicate estimation (does not disable explicit dedup)
    #[arg(long, alias = "dont_eval_duplication")]
    #[arg(help_heading = "Reports", display_order = 105, hide_short_help = true)]
    pub dont_eval_duplication: bool,
    /// Duplicate catalog memory/cache budget in MiB; exact dedup spills to temporary SQLite
    #[arg(long,default_value_t=256,value_parser=parse_positive_usize)]
    #[arg(
        help_heading = "Duplicate removal",
        display_order = 101,
        hide_short_help = true
    )]
    pub duplication_memory_mib: usize,
    /// Report bounded sampled overrepresented substrings before/after filtering
    #[arg(long, alias = "overrepresentation_analysis")]
    #[arg(help_heading = "Reports", display_order = 106, hide_short_help = true)]
    pub overrepresentation_analysis: bool,
    /// Analyze every Nth record for overrepresentation
    #[arg(long,alias="overrepresentation_sampling",default_value_t=20,value_parser=clap::value_parser!(u32).range(1..=10000))]
    #[arg(help_heading = "Reports", display_order = 107, hide_short_help = true)]
    pub overrepresentation_sampling: u32,
    /// Write a standalone JSON QC report (Plexless schema)
    #[arg(long)]
    #[arg(help_heading = "Reports", display_order = 101)]
    #[arg(value_name = "PATH")]
    pub json: Option<PathBuf>,
    /// Write a standalone interactive HTML QC report
    #[arg(long)]
    #[arg(help_heading = "Reports", display_order = 102)]
    #[arg(value_name = "PATH")]
    pub html: Option<PathBuf>,
    /// Title of the QC reports
    #[arg(
        long,
        alias = "report_title",
        default_value = "Plexless preprocessing report"
    )]
    #[arg(help_heading = "Reports", display_order = 103, hide_short_help = true)]
    pub report_title: String,
    /// Additional observed-orientation ACGT adapters; apply to both mates
    #[arg(long, alias = "adapter_fasta")]
    #[arg(
        help_heading = "Adapters and overlap",
        display_order = 116,
        hide_short_help = true
    )]
    pub adapter_fasta: Option<PathBuf>,
    /// Infer mate-specific PE adapters from bounded startup samples
    #[arg(long, alias = "detect_adapter_for_pe")]
    #[arg(
        help_heading = "Adapters and overlap",
        display_order = 117,
        hide_short_help = true
    )]
    pub detect_adapter_for_pe: bool,
    /// Infer an adapter from a bounded startup sample (no built-in manufacturer sequences)
    #[arg(long)]
    #[arg(
        help_heading = "Adapters and overlap",
        display_order = 118,
        hide_short_help = true
    )]
    pub detect_adapter: bool,
    /// Allow one internal gap in unique short-insert adapter overlap; merge/correction stay ungapped
    #[arg(long, alias = "allow_gap_overlap_trimming")]
    #[arg(
        help_heading = "Adapters and overlap",
        display_order = 119,
        hide_short_help = true
    )]
    pub allow_gap_overlap_trimming: bool,
    /// Blacklist of CASAVA index1 sequences, one per line
    #[arg(long, alias = "filter_by_index1")]
    #[arg(
        help_heading = "Filtering",
        display_order = 111,
        hide_short_help = true
    )]
    pub filter_by_index1: Option<PathBuf>,
    /// Blacklist of CASAVA index2 sequences, one per line
    #[arg(long, alias = "filter_by_index2")]
    #[arg(
        help_heading = "Filtering",
        display_order = 112,
        hide_short_help = true
    )]
    pub filter_by_index2: Option<PathBuf>,
    /// Maximum substitutions for equal-length index blacklist matching
    #[arg(long, alias = "filter_by_index_threshold", default_value_t = 0)]
    #[arg(
        help_heading = "Filtering",
        display_order = 113,
        hide_short_help = true
    )]
    pub filter_by_index_threshold: usize,
    /// Merge uniquely overlapping paired reads before filtering
    #[arg(long)]
    #[arg(
        help_heading = "Merging and correction",
        display_order = 100,
        hide_short_help = true
    )]
    pub merge: bool,
    /// Interleave passing unmerged pairs into merged output
    #[arg(long, alias = "include_unmerged", requires = "merge")]
    #[arg(
        help_heading = "Merging and correction",
        display_order = 102,
        hide_short_help = true
    )]
    pub include_unmerged: bool,
    /// Explicit merged FASTQ path (one output sample only)
    #[arg(long, alias = "merged_out", requires = "merge")]
    #[arg(
        help_heading = "Merging and correction",
        display_order = 101,
        hide_short_help = true
    )]
    pub merged_out: Option<PathBuf>,
    /// Explicit SE/R1 output path; .gz selects gzip (one output sample only)
    #[arg(long)]
    #[arg(
        help_heading = "Stream and auxiliary output",
        display_order = 100,
        hide_short_help = true
    )]
    pub out1: Option<PathBuf>,
    /// Explicit R2 output path; .gz selects gzip (one output sample only)
    #[arg(long)]
    #[arg(
        help_heading = "Stream and auxiliary output",
        display_order = 101,
        hide_short_help = true
    )]
    pub out2: Option<PathBuf>,
    /// Write passing R1 mates from pairs rejected by a filter
    #[arg(long)]
    #[arg(
        help_heading = "Stream and auxiliary output",
        display_order = 106,
        hide_short_help = true
    )]
    pub unpaired1: Option<PathBuf>,
    /// Write passing R2 mates; can share the unpaired1 path
    #[arg(long)]
    #[arg(
        help_heading = "Stream and auxiliary output",
        display_order = 107,
        hide_short_help = true
    )]
    pub unpaired2: Option<PathBuf>,
    /// Write reads/consensus rejected by enabled filters
    #[arg(long, alias = "failed_out")]
    #[arg(
        help_heading = "Stream and auxiliary output",
        display_order = 108,
        hide_short_help = true
    )]
    pub failed_out: Option<PathBuf>,
    /// Write R1 overlap regions with zero disagreements
    #[arg(long, alias = "overlapped_out")]
    #[arg(
        help_heading = "Merging and correction",
        display_order = 104,
        hide_short_help = true
    )]
    pub overlapped_out: Option<PathBuf>,
    /// Read FASTQ from standard input; adapter autodetection is disabled
    #[arg(long)]
    #[arg(help_heading = "Input", display_order = 103)]
    pub stdin: bool,
    /// Also stream retained records to stdout; PE output is interleaved
    #[arg(long)]
    #[arg(
        help_heading = "Stream and auxiliary output",
        display_order = 103,
        hide_short_help = true
    )]
    pub stdout: bool,
    /// Treat --reads/--stdin as alternating R1/R2 FASTQ records
    #[arg(long, alias = "interleaved_in")]
    #[arg(
        help_heading = "Stream and auxiliary output",
        display_order = 104,
        hide_short_help = true
    )]
    pub interleaved_in: bool,
    /// Maximum reads/pairs to process; zero means all
    #[arg(long, alias = "reads_to_process", default_value_t = 0)]
    #[arg(
        help_heading = "Stream and auxiliary output",
        display_order = 105,
        hide_short_help = true
    )]
    pub reads_to_process: u64,
    /// Keep the default Plexless refusal to overwrite existing outputs
    #[arg(long, alias = "dont_overwrite")]
    #[arg(
        help_heading = "Stream and auxiliary output",
        display_order = 112,
        hide_short_help = true
    )]
    pub dont_overwrite: bool,
    /// Replace slash characters in the primary read ID with colons
    #[arg(long, alias = "fix_mgi_id")]
    #[arg(
        help_heading = "Stream and auxiliary output",
        display_order = 113,
        hide_short_help = true
    )]
    pub fix_mgi_id: bool,
    /// Log producer progress every million fragments
    #[arg(long)]
    #[arg(help_heading = "Reports", display_order = 109, hide_short_help = true)]
    pub verbose: bool,
    /// Target this many output parts using a full input counting pass
    #[arg(long,conflicts_with="split_by_lines",value_parser=clap::value_parser!(u16).range(2..=999))]
    #[arg(
        help_heading = "Stream and auxiliary output",
        display_order = 109,
        hide_short_help = true
    )]
    pub split: Option<u16>,
    /// Limit each output part to N FASTQ lines (multiple of 4)
    #[arg(long,alias="split_by_lines",value_parser=clap::value_parser!(u64).range(1000..))]
    #[arg(
        help_heading = "Stream and auxiliary output",
        display_order = 110,
        hide_short_help = true
    )]
    pub split_by_lines: Option<u64>,
    /// Number of digits in split filename prefixes; zero disables padding
    #[arg(long,alias="split_prefix_digits",default_value_t=4,value_parser=clap::value_parser!(u8).range(0..=10))]
    #[arg(
        help_heading = "Stream and auxiliary output",
        display_order = 111,
        hide_short_help = true
    )]
    pub split_prefix_digits: u8,
    #[arg(skip)]
    pub(crate) split_records: u64,
}
impl Default for FastpOptions {
    fn default() -> Self {
        Self {
            fastp_defaults: false,
            require_adapter_for_overlap: false,
            dimer_max_len: 2,
            cut_front: false,
            cut_tail: false,
            cut_right: false,
            cut_window_size: 4,
            cut_mean_quality: 20,
            cut_front_window_size: None,
            cut_tail_window_size: None,
            cut_right_window_size: None,
            cut_front_mean_quality: None,
            cut_tail_mean_quality: None,
            cut_right_mean_quality: None,
            trim_front1: 0,
            trim_front2: None,
            trim_tail1: 0,
            trim_tail2: None,
            max_len1: None,
            max_len2: None,
            trim_poly_g: false,
            disable_trim_poly_g: false,
            poly_g_min_len: 10,
            trim_poly_x: false,
            poly_x_min_len: 10,
            low_complexity_filter: false,
            complexity_threshold: 30.0,
            length_limit: None,
            disable_quality_filtering: false,
            disable_length_filtering: false,
            phred64: false,
            correction: false,
            umi: false,
            umi_loc: None,
            umi_len: 0,
            umi_skip: 0,
            umi_prefix: String::new(),
            umi_delim: ":".into(),
            dedup: false,
            evaluate_duplication: false,
            dont_eval_duplication: false,
            duplication_memory_mib: 256,
            overrepresentation_analysis: false,
            overrepresentation_sampling: 20,
            json: None,
            html: None,
            report_title: "Plexless preprocessing report".into(),
            adapter_fasta: None,
            detect_adapter_for_pe: false,
            detect_adapter: false,
            allow_gap_overlap_trimming: false,
            filter_by_index1: None,
            filter_by_index2: None,
            filter_by_index_threshold: 0,
            merge: false,
            include_unmerged: false,
            merged_out: None,
            out1: None,
            out2: None,
            unpaired1: None,
            unpaired2: None,
            failed_out: None,
            overlapped_out: None,
            stdin: false,
            stdout: false,
            interleaved_in: false,
            reads_to_process: 0,
            dont_overwrite: false,
            fix_mgi_id: false,
            verbose: false,
            split: None,
            split_by_lines: None,
            split_prefix_digits: 4,
            split_records: 0,
        }
    }
}

impl FastpOptions {
    pub(crate) fn needs_owned(&self) -> bool {
        self.umi
            || self.fix_mgi_id
            || self.correction
            || self.merge
            || self.failed_out.is_some()
            || self.unpaired1.is_some()
            || self.unpaired2.is_some()
            || self.overlapped_out.is_some()
            || self.stdout
            || self.filter_by_index1.is_some()
            || self.filter_by_index2.is_some()
    }
    pub(crate) fn needs_pipeline(&self) -> bool {
        self.needs_owned()
            || self.json.is_some()
            || self.html.is_some()
            || self.verbose
            || self.overrepresentation_analysis
            || self.phred64
            || self.dedup
            || self.evaluate_duplication
            || self.fastp_defaults
            || self.reads_to_process > 0
            || self.interleaved_in
            || self.stdin
            || self.split.is_some()
            || self.split_by_lines.is_some()
            || self.out1.is_some()
            || self.out2.is_some()
    }
}

/// Resolve profile defaults and validate cross-feature constraints once.
pub(crate) fn prepare(args: &mut crate::cli::DemuxArgs) -> Result<(), String> {
    let o = &args.preprocessing.fastp;
    let paired = args.r1.is_some()
        || o.interleaved_in
        || args.read_mode == Some(crate::cli::ReadMode::Paired);
    if !paired
        && (o.detect_adapter_for_pe
            || o.trim_front2.is_some()
            || o.trim_tail2.is_some()
            || o.max_len2.is_some())
    {
        return Err("R2 clipping and --detect-adapter-for-pe require paired input".into());
    }
    if !paired
        && (o.merge
            || o.correction
            || o.umi_loc
                .is_some_and(|u| matches!(u, UmiLocation::Read2 | UmiLocation::Index2))
            || o.out2.is_some()
            || o.unpaired1.is_some()
            || o.unpaired2.is_some()
            || o.overlapped_out.is_some())
    {
        return Err(
            "Merging, correction, R2 UMI, and paired output options require paired input".into(),
        );
    }
    if o.stdout && (o.split.is_some() || o.split_by_lines.is_some()) {
        return Err("--stdout cannot be combined with output splitting; choose stream output or numbered files".into());
    }
    if o.adapter_fasta.is_some() && args.no_pass_through {
        return Err("--adapter-fasta conflicts with --no-pass-through".into());
    }
    if o.allow_gap_overlap_trimming && (!paired || args.preprocessing.no_overlap) {
        return Err("--allow-gap-overlap-trimming requires paired overlap processing".into());
    }
    let mut paths = std::collections::HashMap::new();
    for (label, path) in [
        ("out1", o.out1.as_ref()),
        ("out2", o.out2.as_ref()),
        ("merged-out", o.merged_out.as_ref()),
        ("failed-out", o.failed_out.as_ref()),
        ("unpaired1", o.unpaired1.as_ref()),
        ("unpaired2", o.unpaired2.as_ref()),
        ("overlapped-out", o.overlapped_out.as_ref()),
        ("json", o.json.as_ref()),
        ("html", o.html.as_ref()),
    ] {
        if let Some(path) = path {
            let path = absolute_path(path)?;
            for reserved in [
                "preprocessing.tsv",
                "retained_fastq_stats.tsv",
                "fastq_stats.tsv",
                "barcode_stats.tsv",
                "sample_metrics.tsv",
                crate::writer::INCOMPLETE_RUN_MARKER,
            ] {
                if path == absolute_path(&args.output.join(reserved))? {
                    return Err(format!(
                        "--{label} conflicts with reserved report/marker {reserved}"
                    ));
                }
            }
            if let Some(previous) = paths.insert(path.clone(), label)
                && !(previous == "unpaired1" && label == "unpaired2")
            {
                return Err(format!(
                    "Output path collision between --{previous} and --{label}: {}",
                    path.display()
                ));
            }
            if path.exists() {
                return Err(format!("Refusing to overwrite output '{}'", path.display()));
            }
        }
    }
    if o.umi && o.umi_loc.is_none() {
        return Err("--umi requires --umi-loc".into());
    }
    if o.umi
        && o.umi_loc.is_some_and(|u| {
            matches!(
                u,
                UmiLocation::Read1 | UmiLocation::Read2 | UmiLocation::PerRead
            )
        })
        && o.umi_len == 0
    {
        return Err("Read-based UMI extraction requires --umi-len greater than zero".into());
    }
    if o.umi_delim.is_empty()
        || o.umi_delim
            .bytes()
            .chain(o.umi_prefix.bytes())
            .any(|b| b.is_ascii_whitespace() || b.is_ascii_control())
    {
        return Err("UMI delimiter and prefix must contain no whitespace/control bytes; delimiter cannot be empty".into());
    }
    o.umi_len
        .checked_add(o.umi_skip)
        .ok_or("UMI length plus skip overflow")?;
    o.duplication_memory_mib
        .checked_mul(1024 * 1024)
        .ok_or("Duplication memory budget overflow")?;
    if o.split.is_some() && (o.stdin || args.reads.as_deref() == Some(std::path::Path::new("-"))) {
        return Err("--split requires seekable input; use --split-by-lines for stdin".into());
    }
    if o.split_by_lines
        .is_some_and(|lines| !lines.is_multiple_of(4))
    {
        return Err("--split-by-lines must be a multiple of 4 FASTQ lines".into());
    }
    if o.interleaved_in && args.r1.is_some() {
        return Err(
            "--interleaved-in uses --reads or --stdin, not separate --r1/--r2 inputs".into(),
        );
    }
    if args.cram.is_some()
        && (o.needs_owned()
            || o.phred64
            || o.interleaved_in
            || o.stdin
            || o.out1.is_some()
            || o.out2.is_some()
            || o.split.is_some()
            || o.split_by_lines.is_some()
            || o.reads_to_process > 0)
    {
        return Err("FASTQ header/encoding/merging and stream-output features require FASTQ input/output; native CRAM retains its metadata policy".into());
    }
    if o.stdin {
        args.reads = Some("-".into());
    }
    if o.json.is_some() || o.html.is_some() || o.overrepresentation_analysis {
        args.fastq_stats = true;
    }
    if o.fastp_defaults {
        if !o.disable_quality_filtering {
            args.preprocessing
                .max_low_quality_percent
                .get_or_insert(40.0);
            args.preprocessing.max_n.get_or_insert(5);
        }
        if !o.disable_length_filtering {
            args.preprocessing.min_length.get_or_insert(15);
        }
        args.fastq_stats = true;
    }
    if o.disable_quality_filtering {
        args.preprocessing.max_low_quality_percent = None;
        args.preprocessing.max_n = None;
        args.preprocessing.min_average_quality = None;
    }
    if o.disable_length_filtering {
        args.preprocessing.min_length = None;
        args.preprocessing.fastp.length_limit = None;
    }
    if args.preprocessing.fastp.overrepresentation_analysis {
        args.preprocessing
            .fastp
            .json
            .get_or_insert_with(|| args.output.join("plexless.json"));
    }
    if args.preprocessing.fastp.fastp_defaults {
        args.preprocessing
            .fastp
            .json
            .get_or_insert_with(|| args.output.join("plexless.json"));
        args.preprocessing
            .fastp
            .html
            .get_or_insert_with(|| args.output.join("plexless.html"));
    }
    if let Some(parts) = args.preprocessing.fastp.split {
        let path = args
            .r1
            .as_deref()
            .or(args.reads.as_deref())
            .ok_or("--split requires FASTQ input")?;
        let mut reader = crate::fastq::InputReader::try_open(path)?;
        let mut records = 0u64;
        while let Some(record) = reader.next_record() {
            let record = record.map_err(|e| e.to_string())?;
            if record.qual().is_none() {
                return Err("Splitting requires FASTQ quality scores".into());
            }
            records = records.checked_add(1).ok_or("Input count overflow")?;
        }
        if args.preprocessing.fastp.interleaved_in {
            records = records.div_ceil(2);
        }
        if args.preprocessing.fastp.reads_to_process > 0 {
            records = records.min(args.preprocessing.fastp.reads_to_process);
        }
        args.preprocessing.fastp.split_records = records.div_ceil(u64::from(parts)).max(1);
        eprintln!(
            "Output split: counted {records} input fragments; up to {} records per part per output stream",
            args.preprocessing.fastp.split_records
        );
    } else if let Some(lines) = args.preprocessing.fastp.split_by_lines {
        args.preprocessing.fastp.split_records = lines / 4;
    }
    if args.preprocessing.fastp.fastp_defaults
        && !args.preprocessing.fastp.disable_trim_poly_g
        && !args.preprocessing.fastp.trim_poly_g
        && let Some(path) = args
            .r1
            .as_deref()
            .or(args.reads.as_deref())
            .filter(|p| *p != std::path::Path::new("-"))
    {
        let mut reader = crate::fastq::InputReader::try_open(path)?;
        if let Some(record) = reader.next_record() {
            let record = record.map_err(|e| e.to_string())?;
            let id = record.id().strip_prefix(b"@").unwrap_or(record.id());
            if [
                "FS", "MN", "SH", "NS", "NB", "NDX", "VL", "VH", "A", "NA", "LH",
            ]
            .iter()
            .any(|prefix| id.starts_with(prefix.as_bytes()))
            {
                args.preprocessing.fastp.trim_poly_g = true;
                eprintln!("PolyG trimming enabled from the input instrument identifier");
            }
        }
    }
    let detect1 = args.preprocessing.fastp.detect_adapter
        || (args.preprocessing.fastp.fastp_defaults && !paired);
    let detect2 = args.preprocessing.fastp.detect_adapter_for_pe;
    if (detect1 || detect2) && args.cram.is_none() && !args.no_pass_through {
        let paths = if args.preprocessing.fastp.interleaved_in {
            [args.reads.as_deref(), args.reads.as_deref()]
        } else if paired {
            [args.r1.as_deref(), args.r2.as_deref()]
        } else {
            [args.reads.as_deref(), None]
        };
        let layouts = if paired {
            [args.r1_structure.as_deref(), args.r2_structure.as_deref()]
        } else {
            [args.structure.as_deref(), None]
        };
        for mate in 0..if paired { 2 } else { 1 } {
            let explicit = if mate == 0 {
                &args.preprocessing.r1_distal_sequence
            } else {
                &args.preprocessing.r2_distal_sequence
            };
            let distal = if mate == 0 {
                &args.r1_distal_structure
            } else {
                &args.r2_distal_structure
            };
            if explicit.is_some() || distal.is_some() {
                continue;
            }
            if let Some(path) = paths[mate] {
                let prefix = layouts[mate]
                    .map(crate::structure::ReadStructure::parse)
                    .transpose()?
                    .map_or(0, |s| s.prefix_len);
                if let Some(adapter) = if args.preprocessing.fastp.interleaved_in {
                    crate::adapter_detection::infer_mate(path, prefix, Some(mate))?
                } else {
                    crate::adapter_detection::infer(path, prefix)?
                } {
                    let sequence = String::from_utf8(adapter).expect("inferred DNA is ASCII");
                    eprintln!("Inferred R{} distal adapter: {sequence}", mate + 1);
                    if mate == 0 {
                        args.preprocessing.r1_distal_sequence = Some(sequence);
                    } else {
                        args.preprocessing.r2_distal_sequence = Some(sequence);
                    }
                    args.preprocessing.partial_pass_through = true;
                } else {
                    eprintln!(
                        "Adapter inference: R{} has insufficient startup evidence; reads retain configured/inferred structure behavior",
                        mate + 1
                    );
                }
            }
        }
    }
    let options = &mut args.preprocessing.fastp;
    for path in [
        &mut options.out1,
        &mut options.out2,
        &mut options.merged_out,
        &mut options.unpaired1,
        &mut options.unpaired2,
        &mut options.failed_out,
        &mut options.overlapped_out,
        &mut options.json,
        &mut options.html,
    ]
    .into_iter()
    .flatten()
    {
        *path = absolute_path(path)?;
    }
    Ok(())
}

pub(crate) fn absolute_path(path: &std::path::Path) -> Result<PathBuf, String> {
    use std::path::Component;
    let path = if path.is_absolute() {
        path.to_owned()
    } else {
        std::env::current_dir()
            .map_err(|e| e.to_string())?
            .join(path)
    };
    let mut normalized = PathBuf::new();
    for component in path.components() {
        match component {
            Component::CurDir => {}
            Component::ParentDir => {
                normalized.pop();
            }
            value => normalized.push(value.as_os_str()),
        }
    }
    let mut ancestor = normalized.as_path();
    let mut suffix = Vec::new();
    while !ancestor.exists() {
        if let Some(name) = ancestor.file_name() {
            suffix.push(name.to_owned());
        }
        ancestor = ancestor
            .parent()
            .ok_or("Output path has no existing ancestor")?;
    }
    let mut resolved = ancestor.canonicalize().map_err(|e| e.to_string())?;
    for name in suffix.into_iter().rev() {
        resolved.push(name);
    }
    Ok(resolved)
}
