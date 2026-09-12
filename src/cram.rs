//! Unmapped CRAM boundary code.
//!
//! Routing never consults the HTS record stored here.  The decoded sequence and
//! normalized FASTQ-style qualities are the only fields exposed to Plexless's
//! existing demultiplexing core; the original record is a sidecar used solely
//! when the selected output format is CRAM.

use std::collections::{BTreeMap, HashMap};
use std::ffi::CString;
use std::fs;
use std::path::{Path, PathBuf};
use std::ptr::NonNull;
use std::time::{Duration, Instant};

use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, HeaderView, Read};
use rust_htslib::htslib;

use crate::cli::ReadMode;
use crate::output::{OutputMate, OutputTarget};
use crate::parallel::OwnedFastqRecord;
use crate::samples::SampleSheet;
use crate::writer::{build_sample_names, prepare_output_dir};

const CRAM_MAGIC: &[u8; 4] = b"CRAM";
const MAX_PHRED: u8 = 93;

#[derive(Debug)]
pub(crate) struct CramOwnedRecord {
    pub(crate) read: OwnedFastqRecord,
    pub(crate) source: Box<bam::Record>,
}

#[derive(Debug)]
pub(crate) enum CramInputItem {
    Single(CramOwnedRecord),
    Pair {
        r1: CramOwnedRecord,
        r2: CramOwnedRecord,
        r1_first: bool,
    },
    Orphan {
        mate: OutputMate,
        record: CramOwnedRecord,
    },
}

#[derive(Debug)]
pub(crate) enum CramOutputItem {
    Single {
        target: OutputTarget,
        record: CramOwnedRecord,
        trim_start: usize,
    },
    Pair {
        target: OutputTarget,
        r1: CramOwnedRecord,
        r2: CramOwnedRecord,
        r1_trim: usize,
        r2_trim: usize,
        r1_first: bool,
    },
    Orphan {
        target: OutputTarget,
        mate: OutputMate,
        record: CramOwnedRecord,
    },
}

impl CramOutputItem {
    pub(crate) fn target(&self) -> OutputTarget {
        match self {
            Self::Single { target, .. }
            | Self::Pair { target, .. }
            | Self::Orphan { target, .. } => *target,
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct CramHeader {
    header: bam::Header,
}

impl CramHeader {
    fn with_plexless_program(reader: &bam::Reader) -> Self {
        let view = reader.header();
        let mut header = bam::Header::from_template(view);
        let program_id = unique_program_id(view.as_bytes());
        let command = std::env::args()
            .map(|arg| arg.replace(['\t', '\n', '\r'], " "))
            .collect::<Vec<_>>()
            .join(" ");
        header.push_record(
            HeaderRecord::new(b"PG")
                .push_tag(b"ID", &program_id)
                .push_tag(b"PN", "plexless")
                .push_tag(b"VN", env!("CARGO_PKG_VERSION"))
                .push_tag(b"CL", command),
        );
        Self { header }
    }
}

pub(crate) fn inspect(path: &Path, mode: ReadMode) -> Result<CramHeader, String> {
    validate_cram_magic(path)?;
    let reader = bam::Reader::from_path(path)
        .map_err(|error| format!("Could not open CRAM input '{}': {error}", path.display()))?;
    if mode == ReadMode::Paired {
        validate_query_grouped_header(reader.header().as_bytes())?;
    }
    Ok(CramHeader::with_plexless_program(&reader))
}

pub(crate) fn read_items(
    path: &Path,
    mode: ReadMode,
    decode_threads: usize,
    mut emit: impl FnMut(CramInputItem) -> Result<(), String>,
) -> Result<(), String> {
    validate_cram_magic(path)?;
    let mut reader = bam::Reader::from_path(path)
        .map_err(|error| format!("Could not open CRAM input '{}': {error}", path.display()))?;
    if decode_threads > 0 {
        reader.set_threads(decode_threads).map_err(|error| {
            format!(
                "Could not start {decode_threads} HTSlib CRAM decode worker(s) for '{}': {error}",
                path.display()
            )
        })?;
    }

    if mode == ReadMode::Paired {
        validate_query_grouped_header(reader.header().as_bytes())?;
    }

    match mode {
        ReadMode::Single => {
            for (index, result) in reader.records().enumerate() {
                let record = result.map_err(|error| {
                    format!(
                        "CRAM decode error in '{}' at record {}: {error}",
                        path.display(),
                        index + 1
                    )
                })?;
                validate_record(&record, mode, index + 1)?;
                emit(CramInputItem::Single(to_owned_record(record)?))?;
            }
        }
        ReadMode::Paired => {
            let mut group_name: Option<Vec<u8>> = None;
            let mut group = Vec::with_capacity(2);
            for (index, result) in reader.records().enumerate() {
                let record = result.map_err(|error| {
                    format!(
                        "CRAM decode error in '{}' at record {}: {error}",
                        path.display(),
                        index + 1
                    )
                })?;
                validate_record(&record, mode, index + 1)?;
                let qname = record.qname().to_vec();
                if group_name.as_deref().is_some_and(|name| name != qname) {
                    emit_paired_group(
                        group_name.take().unwrap_or_default(),
                        &mut group,
                        &mut emit,
                    )?;
                }
                if group_name.is_none() {
                    group_name = Some(qname);
                }
                group.push(to_owned_record(record)?);
                if group.len() > 2 {
                    return Err(format!(
                        "CRAM QNAME '{}' contains more than two primary records",
                        String::from_utf8_lossy(group_name.as_deref().unwrap_or_default())
                    ));
                }
            }
            if let Some(name) = group_name {
                emit_paired_group(name, &mut group, &mut emit)?;
            }
        }
    }
    Ok(())
}

fn emit_paired_group(
    name: Vec<u8>,
    group: &mut Vec<CramOwnedRecord>,
    emit: &mut impl FnMut(CramInputItem) -> Result<(), String>,
) -> Result<(), String> {
    let display_name = String::from_utf8_lossy(&name);
    match group.len() {
        0 => Err("Internal error: empty CRAM QNAME group".into()),
        1 => {
            let record = group.pop().ok_or("Internal CRAM group error")?;
            let mate = if record.source.is_first_in_template() {
                OutputMate::R1
            } else {
                OutputMate::R2
            };
            emit(CramInputItem::Orphan { mate, record })
        }
        2 => {
            let second = group.pop().ok_or("Internal CRAM group error")?;
            let first = group.pop().ok_or("Internal CRAM group error")?;
            let first_is_r1 = first.source.is_first_in_template();
            let second_is_r1 = second.source.is_first_in_template();
            if first_is_r1 == second_is_r1 {
                let mate = if first_is_r1 { "R1" } else { "R2" };
                return Err(format!(
                    "CRAM QNAME '{display_name}' contains duplicate primary {mate} records"
                ));
            }
            let (r1, r2) = if first_is_r1 {
                (first, second)
            } else {
                (second, first)
            };
            emit(CramInputItem::Pair {
                r1,
                r2,
                r1_first: first_is_r1,
            })
        }
        _ => Err(format!(
            "CRAM QNAME '{display_name}' contains more than two primary records"
        )),
    }
}

fn validate_record(record: &bam::Record, mode: ReadMode, ordinal: usize) -> Result<(), String> {
    let name = String::from_utf8_lossy(record.qname());
    let context = format!("CRAM record {ordinal} ('{name}')");
    if record.qname().is_empty() {
        return Err(format!("{context} has an empty QNAME"));
    }
    if record.is_secondary() {
        return Err(format!(
            "{context} is secondary; only primary unmapped records are supported"
        ));
    }
    if record.is_supplementary() {
        return Err(format!(
            "{context} is supplementary; only primary unmapped records are supported"
        ));
    }
    if !record.is_unmapped()
        || record.tid() != -1
        || record.pos() != -1
        || record.mapq() != 0
        || record.cigar_len() != 0
    {
        return Err(format!(
            "{context} is aligned or carries alignment coordinates/CIGAR; aligned CRAM is not supported"
        ));
    }

    match mode {
        ReadMode::Single => {
            const PAIRED_SEMANTICS: u16 = 0x1 | 0x2 | 0x8 | 0x20 | 0x40 | 0x80;
            if record.flags() & PAIRED_SEMANTICS != 0
                || record.mtid() != -1
                || record.mpos() != -1
                || record.insert_size() != 0
            {
                return Err(format!(
                    "{context} has paired-end flags or mate fields but --read-mode single was declared"
                ));
            }
        }
        ReadMode::Paired => {
            if !record.is_paired() {
                return Err(format!(
                    "{context} is not flagged paired but --read-mode paired was declared"
                ));
            }
            if !record.is_mate_unmapped() {
                return Err(format!("{context} does not mark its mate unmapped"));
            }
            let first = record.is_first_in_template();
            let second = record.is_last_in_template();
            if first == second {
                return Err(format!(
                    "{context} must set exactly one of READ1 or READ2 in paired mode"
                ));
            }
            if record.flags() & 0x2 != 0
                || record.mtid() != -1
                || record.mpos() != -1
                || record.insert_size() != 0
            {
                return Err(format!(
                    "{context} has proper-pair or mapped-mate fields incompatible with raw unmapped CRAM"
                ));
            }
        }
    }

    if record.seq_len() == 0 {
        return Err(format!("{context} has no sequence"));
    }
    let qualities = record.qual();
    if qualities.len() != record.seq_len() {
        return Err(format!("{context} sequence/quality lengths do not match"));
    }
    if qualities.contains(&255) {
        return Err(format!("{context} has missing quality values"));
    }
    if let Some(quality) = qualities
        .iter()
        .copied()
        .find(|&quality| quality > MAX_PHRED)
    {
        return Err(format!(
            "{context} has Phred quality {quality}; the supported SAM/FASTQ range is 0..={MAX_PHRED}"
        ));
    }
    Ok(())
}

fn to_owned_record(record: bam::Record) -> Result<CramOwnedRecord, String> {
    let seq = record.seq().as_bytes();
    let qual = record
        .qual()
        .iter()
        .map(|quality| {
            quality
                .checked_add(33)
                .ok_or("CRAM quality conversion overflow")
        })
        .collect::<Result<Vec<_>, _>>()?;
    let read = OwnedFastqRecord::new(record.qname(), &seq, &qual);
    Ok(CramOwnedRecord {
        read,
        source: Box::new(record),
    })
}

fn validate_cram_magic(path: &Path) -> Result<(), String> {
    use std::io::Read as _;
    let mut file = fs::File::open(path)
        .map_err(|error| format!("Could not open CRAM input '{}': {error}", path.display()))?;
    let mut magic = [0u8; 4];
    file.read_exact(&mut magic)
        .map_err(|error| format!("Could not read CRAM header '{}': {error}", path.display()))?;
    if &magic != CRAM_MAGIC {
        return Err(format!(
            "Input '{}' is not a CRAM file (missing CRAM magic)",
            path.display()
        ));
    }
    Ok(())
}

fn validate_query_grouped_header(header: &[u8]) -> Result<(), String> {
    let hd = header
        .split(|byte| *byte == b'\n')
        .find(|line| line.starts_with(b"@HD\t"))
        .ok_or("Paired CRAM requires @HD SO:queryname or GO:query")?;
    let mut sort_order = None;
    let mut group_order = None;
    for field in hd.split(|byte| *byte == b'\t').skip(1) {
        if let Some(value) = field.strip_prefix(b"SO:") {
            sort_order = Some(value);
        }
        if let Some(value) = field.strip_prefix(b"GO:") {
            group_order = Some(value);
        }
    }
    if sort_order == Some(b"queryname") || group_order == Some(b"query") {
        return Ok(());
    }
    Err(format!(
        "Paired CRAM must be queryname-grouped (@HD SO:queryname or GO:query); found SO={} GO={}",
        sort_order
            .map(String::from_utf8_lossy)
            .as_deref()
            .unwrap_or("unspecified"),
        group_order
            .map(String::from_utf8_lossy)
            .as_deref()
            .unwrap_or("unspecified")
    ))
}

fn unique_program_id(header: &[u8]) -> String {
    let existing = header
        .split(|byte| *byte == b'\n')
        .filter(|line| line.starts_with(b"@PG\t"))
        .flat_map(|line| line.split(|byte| *byte == b'\t'))
        .filter_map(|field| field.strip_prefix(b"ID:"))
        .collect::<Vec<_>>();
    for suffix in 0u32.. {
        let candidate = if suffix == 0 {
            "plexless".to_string()
        } else {
            format!("plexless.{suffix}")
        };
        if !existing.contains(&candidate.as_bytes()) {
            return candidate;
        }
    }
    unreachable!()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum MetadataAction {
    Preserved,
    Rewritten,
    Removed,
}

impl MetadataAction {
    fn label(self) -> &'static str {
        match self {
            Self::Preserved => "preserved",
            Self::Rewritten => "rewritten",
            Self::Removed => "removed",
        }
    }
}

#[derive(Debug, Default)]
pub(crate) struct MetadataStats {
    counts: BTreeMap<(MetadataAction, [u8; 2]), u64>,
}

impl MetadataStats {
    fn observe(&mut self, action: MetadataAction, tag: &[u8]) -> Result<(), String> {
        let tag: [u8; 2] = tag.try_into().map_err(|_| "Invalid auxiliary tag length")?;
        let count = self.counts.entry((action, tag)).or_default();
        *count = count.checked_add(1).ok_or("Metadata count overflow")?;
        Ok(())
    }

    pub(crate) fn write_report(&self, output_dir: &Path) -> Result<(), String> {
        let mut report = String::from("Action\tTag\tCount\n");
        for ((action, tag), count) in &self.counts {
            report.push_str(action.label());
            report.push('\t');
            report.push_str(&String::from_utf8_lossy(tag));
            report.push('\t');
            report.push_str(&count.to_string());
            report.push('\n');
        }
        let path = output_dir.join("cram_metadata.tsv");
        fs::write(&path, report)
            .map_err(|error| format!("Could not write '{}': {error}", path.display()))?;

        for action in [
            MetadataAction::Preserved,
            MetadataAction::Rewritten,
            MetadataAction::Removed,
        ] {
            let summary = self
                .counts
                .iter()
                .filter(|((candidate, _), _)| *candidate == action)
                .map(|((_, tag), count)| format!("{}={count}", String::from_utf8_lossy(tag)))
                .collect::<Vec<_>>();
            eprintln!(
                "CRAM metadata {}: {}",
                action.label(),
                if summary.is_empty() {
                    "none".into()
                } else {
                    summary.join(", ")
                }
            );
        }
        Ok(())
    }
}

/// CRAM writer wrapper with an explicit, fallible close operation.
///
/// rust-htslib 1.0.1 closes `bam::Writer` only from `Drop` and discards the
/// return value from `hts_close`. CRAM finalization writes buffered containers
/// and the EOF block, so production success must observe that result directly.
struct FinalizingCramWriter {
    file: Option<NonNull<htslib::htsFile>>,
    header: HeaderView,
    path: PathBuf,
}

// HTSlib handles may move between threads provided one thread owns and uses a
// handle at a time. This wrapper is moved once into Plexless's writer thread.
unsafe impl Send for FinalizingCramWriter {}

impl FinalizingCramWriter {
    fn open(path: PathBuf, header: &bam::Header, compression_level: u32) -> Result<Self, String> {
        if compression_level > 9 {
            return Err(format!(
                "Invalid CRAM compression level {compression_level}; expected 0..=9"
            ));
        }
        let encoded_path = CString::new(path.as_os_str().as_encoded_bytes()).map_err(|_| {
            format!(
                "CRAM output path contains an embedded NUL byte: '{}'",
                path.display()
            )
        })?;
        let mode = c"wc";
        let raw = unsafe { htslib::hts_open(encoded_path.as_ptr(), mode.as_ptr()) };
        let file = NonNull::new(raw)
            .ok_or_else(|| format!("Could not create CRAM output '{}'", path.display()))?;
        let header = HeaderView::from_header(header);

        if unsafe { htslib::sam_hdr_write(file.as_ptr(), header.inner_ptr()) } < 0 {
            unsafe {
                htslib::hts_close(file.as_ptr());
            }
            return Err(format!(
                "Could not write CRAM header to '{}'",
                path.display()
            ));
        }
        if unsafe {
            htslib::hts_set_opt(
                file.as_ptr(),
                htslib::hts_fmt_option_HTS_OPT_COMPRESSION_LEVEL,
                compression_level,
            )
        } != 0
        {
            unsafe {
                htslib::hts_close(file.as_ptr());
            }
            return Err(format!(
                "Could not set CRAM compression level for '{}'",
                path.display()
            ));
        }

        Ok(Self {
            file: Some(file),
            header,
            path,
        })
    }

    fn write(&mut self, record: &bam::Record) -> Result<(), String> {
        let file = self.file.ok_or("CRAM writer is already finalized")?;
        if unsafe { htslib::sam_write1(file.as_ptr(), self.header.inner_ptr(), record.inner()) } < 0
        {
            return Err(format!(
                "Could not write CRAM record to '{}'",
                self.path.display()
            ));
        }
        Ok(())
    }

    fn finish(mut self) -> Result<(), String> {
        let file = self.file.take().ok_or("CRAM writer is already finalized")?;
        let status = unsafe { htslib::hts_close(file.as_ptr()) };
        check_close_status(&self.path, status)
    }
}

fn check_close_status(path: &Path, status: i32) -> Result<(), String> {
    if status != 0 {
        Err(format!(
            "Could not finalize CRAM output '{}': hts_close returned {status}: {}",
            path.display(),
            std::io::Error::last_os_error()
        ))
    } else {
        Ok(())
    }
}

impl Drop for FinalizingCramWriter {
    fn drop(&mut self) {
        if let Some(file) = self.file.take() {
            unsafe {
                htslib::hts_close(file.as_ptr());
            }
        }
    }
}

const CRAM_DESCRIPTOR_RESERVE: usize = 64;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct DescriptorPlan {
    required_soft: usize,
    raise_soft_to: Option<usize>,
}

fn descriptor_plan(
    writer_count: usize,
    requested_max_open_files: Option<usize>,
    soft_limit: usize,
    hard_limit: usize,
) -> Result<DescriptorPlan, String> {
    if let Some(requested) = requested_max_open_files
        && requested < writer_count
    {
        return Err(format!(
            "CRAM output requires {writer_count} simultaneously open writers, but --max-open-files is {requested}; CRAM append/reopen is not supported"
        ));
    }
    let required_soft = writer_count
        .checked_add(CRAM_DESCRIPTOR_RESERVE)
        .ok_or("CRAM file-descriptor requirement overflow")?;
    if hard_limit < required_soft {
        return Err(format!(
            "CRAM output requires a process file-descriptor limit of at least {required_soft} ({writer_count} writers plus {CRAM_DESCRIPTOR_RESERVE} input/report/safety descriptors), but RLIMIT_NOFILE is soft={soft_limit} hard={hard_limit}; raise the service/user hard limit or reduce the number of sample outputs"
        ));
    }
    Ok(DescriptorPlan {
        required_soft,
        raise_soft_to: (soft_limit < required_soft).then_some(required_soft),
    })
}

#[cfg(unix)]
fn preflight_cram_descriptors(
    writer_count: usize,
    requested_max_open_files: Option<usize>,
) -> Result<(), String> {
    let mut limits = libc::rlimit {
        rlim_cur: 0,
        rlim_max: 0,
    };
    if unsafe { libc::getrlimit(libc::RLIMIT_NOFILE, &mut limits) } != 0 {
        return Err(format!(
            "Could not inspect RLIMIT_NOFILE for CRAM output: {}",
            std::io::Error::last_os_error()
        ));
    }
    let soft = usize::try_from(limits.rlim_cur)
        .map_err(|_| "RLIMIT_NOFILE soft limit does not fit this platform")?;
    let hard = usize::try_from(limits.rlim_max)
        .map_err(|_| "RLIMIT_NOFILE hard limit does not fit this platform")?;
    let plan = descriptor_plan(writer_count, requested_max_open_files, soft, hard)?;
    if let Some(required) = plan.raise_soft_to {
        let raised = libc::rlimit {
            rlim_cur: required as libc::rlim_t,
            rlim_max: limits.rlim_max,
        };
        if unsafe { libc::setrlimit(libc::RLIMIT_NOFILE, &raised) } != 0 {
            return Err(format!(
                "Could not raise RLIMIT_NOFILE soft limit from {soft} to {required} for {writer_count} CRAM writers (hard limit {hard}): {}",
                std::io::Error::last_os_error()
            ));
        }
        eprintln!(
            "CRAM writer descriptors: writers={writer_count} reserve={CRAM_DESCRIPTOR_RESERVE} soft={soft}->{required} hard={hard}"
        );
    } else {
        eprintln!(
            "CRAM writer descriptors: writers={writer_count} reserve={CRAM_DESCRIPTOR_RESERVE} soft={soft} hard={hard}"
        );
    }
    Ok(())
}

#[cfg(not(unix))]
fn preflight_cram_descriptors(
    writer_count: usize,
    requested_max_open_files: Option<usize>,
) -> Result<(), String> {
    const CONSERVATIVE_NON_UNIX_WRITER_LIMIT: usize = 256;
    descriptor_plan(
        writer_count,
        requested_max_open_files,
        CONSERVATIVE_NON_UNIX_WRITER_LIMIT + CRAM_DESCRIPTOR_RESERVE,
        CONSERVATIVE_NON_UNIX_WRITER_LIMIT + CRAM_DESCRIPTOR_RESERVE,
    )
    .map(|_| ())
}

pub(crate) struct CramWriterManager {
    output_dir: PathBuf,
    sample_names: Vec<String>,
    header: bam::Header,
    writers: HashMap<OutputTarget, FinalizingCramWriter>,
    output_paths: HashMap<OutputTarget, PathBuf>,
    write_unassigned: bool,
    compression_level: u32,
    metadata: MetadataStats,
    records_written: u64,
}

impl CramWriterManager {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        output_dir: PathBuf,
        samples: &SampleSheet,
        header: CramHeader,
        write_unassigned: bool,
        max_open_files: Option<usize>,
        compression_level: u32,
    ) -> Result<Self, String> {
        let target_count = samples
            .samples
            .len()
            .checked_add(usize::from(write_unassigned))
            .ok_or("CRAM output count overflow")?;
        let sample_names = build_sample_names(samples, write_unassigned)?;
        preflight_cram_descriptors(target_count, max_open_files)?;
        prepare_output_dir(&output_dir)?;
        Ok(Self {
            output_dir,
            sample_names,
            header: header.header,
            writers: HashMap::new(),
            output_paths: HashMap::new(),
            write_unassigned,
            compression_level,
            metadata: MetadataStats::default(),
            records_written: 0,
        })
    }

    pub(crate) fn write(&mut self, item: CramOutputItem) -> Result<(), String> {
        let target = item.target();
        match item {
            CramOutputItem::Single {
                mut record,
                trim_start,
                ..
            } => {
                transform_record(&mut record, trim_start, None, None, &mut self.metadata)?;
                self.write_record(target, &record.source)?;
            }
            CramOutputItem::Pair {
                mut r1,
                mut r2,
                r1_trim,
                r2_trim,
                r1_first,
                ..
            } => {
                let r1_mate_seq = r2.read.seq[r2_trim..].to_vec();
                let r1_mate_qual = r2.read.qual[r2_trim..].to_vec();
                let r2_mate_seq = r1.read.seq[r1_trim..].to_vec();
                let r2_mate_qual = r1.read.qual[r1_trim..].to_vec();
                transform_record(
                    &mut r1,
                    r1_trim,
                    (r2_trim > 0).then_some(r1_mate_seq.as_slice()),
                    (r2_trim > 0).then_some(r1_mate_qual.as_slice()),
                    &mut self.metadata,
                )?;
                transform_record(
                    &mut r2,
                    r2_trim,
                    (r1_trim > 0).then_some(r2_mate_seq.as_slice()),
                    (r1_trim > 0).then_some(r2_mate_qual.as_slice()),
                    &mut self.metadata,
                )?;
                if r1_first {
                    self.write_record(target, &r1.source)?;
                    self.write_record(target, &r2.source)?;
                } else {
                    self.write_record(target, &r2.source)?;
                    self.write_record(target, &r1.source)?;
                }
            }
            CramOutputItem::Orphan { mut record, .. } => {
                transform_record(&mut record, 0, None, None, &mut self.metadata)?;
                self.write_record(target, &record.source)?;
            }
        }
        Ok(())
    }

    pub(crate) fn finish(mut self) -> Result<(MetadataStats, u64, u64, Duration), String> {
        let finalization_started = Instant::now();
        let mut first_close_error = None;
        for (_, writer) in self.writers.drain() {
            if let Err(error) = writer.finish()
                && first_close_error.is_none()
            {
                first_close_error = Some(error);
            }
        }
        if let Some(error) = first_close_error {
            return Err(error);
        }
        let finalization_time = finalization_started.elapsed();
        let records_written = self.records_written;
        let paths = self.output_paths.into_values().collect::<Vec<_>>();
        let output_bytes = paths.iter().try_fold(0u64, |total, path| {
            let bytes = fs::metadata(path)
                .map_err(|error| format!("Could not inspect '{}': {error}", path.display()))?
                .len();
            total
                .checked_add(bytes)
                .ok_or_else(|| "CRAM output size overflow".to_string())
        })?;
        Ok((
            self.metadata,
            records_written,
            output_bytes,
            finalization_time,
        ))
    }

    fn write_record(&mut self, target: OutputTarget, record: &bam::Record) -> Result<(), String> {
        if target == OutputTarget::Unassigned && !self.write_unassigned {
            return Ok(());
        }
        if !self.writers.contains_key(&target) {
            self.open_writer(target)?;
        }
        self.writers
            .get_mut(&target)
            .ok_or("CRAM writer cache failure")?
            .write(record)
            .map_err(|error| format!("Could not write CRAM record: {error}"))?;
        self.records_written = self
            .records_written
            .checked_add(1)
            .ok_or("CRAM output record count overflow")?;
        Ok(())
    }

    fn open_writer(&mut self, target: OutputTarget) -> Result<(), String> {
        let name = match target {
            OutputTarget::Sample { sample_id } => self
                .sample_names
                .get(usize::try_from(sample_id).map_err(|_| "Invalid sample ID")?)
                .ok_or("Unknown sample ID")?
                .as_str(),
            OutputTarget::Unassigned => "unassigned",
        };
        let path = self.output_dir.join(format!("{name}.cram"));
        let writer =
            FinalizingCramWriter::open(path.clone(), &self.header, self.compression_level)?;
        self.output_paths.insert(target, path);
        self.writers.insert(target, writer);
        Ok(())
    }
}

// Policy for an assigned record whose prefix is trimmed:
// - OQ, BQ, E2, and U2 are position-for-position strings and are trimmed identically.
// - R2 and Q2 are rewritten from the transformed mate when both mates exist.
// - MM/ML/MN and alignment/position-dependent standard tags are removed because
//   Plexless cannot safely repair their coordinates/meaning in v1.
// - all other standard and unknown/custom tags are preserved byte-for-byte by
//   rust-htslib Record::set.
const INVALIDATED_TAGS: [&[u8; 2]; 26] = [
    b"MM", b"ML", b"MN", b"MD", b"NM", b"SA", b"MC", b"AS", b"XS", b"UQ", b"MQ", b"AM", b"SM",
    b"CM", b"NH", b"HI", b"IH", b"PT", b"CC", b"CP", b"CG", b"H0", b"H1", b"H2", b"PQ", b"TS",
];

fn transform_record(
    record: &mut CramOwnedRecord,
    trim_start: usize,
    mate_seq: Option<&[u8]>,
    mate_qual: Option<&[u8]>,
    stats: &mut MetadataStats,
) -> Result<(), String> {
    if trim_start > record.read.seq.len() || trim_start > record.read.qual.len() {
        return Err("Internal error: CRAM trim prefix exceeds record length".into());
    }
    let tags = record
        .source
        .aux_iter()
        .map(|entry| {
            entry
                .map(|(tag, _)| <[u8; 2]>::try_from(tag).expect("HTS auxiliary tags are two bytes"))
                .map_err(|error| format!("Could not inspect CRAM metadata: {error}"))
        })
        .collect::<Result<Vec<_>, _>>()?;

    let mut rewritten = Vec::<[u8; 2]>::new();
    let mut removed = Vec::<[u8; 2]>::new();
    if trim_start > 0 {
        for tag in [b"OQ", b"BQ", b"E2", b"U2"] {
            if rewrite_trimmed_string(&mut record.source, tag, trim_start, record.read.seq.len())? {
                rewritten.push(*tag);
            }
        }
        for tag in INVALIDATED_TAGS {
            if remove_if_present(&mut record.source, tag)? {
                removed.push(*tag);
            }
        }
    }

    if let (Some(mate_seq), Some(mate_qual)) = (mate_seq, mate_qual) {
        if rewrite_string_if_present(&mut record.source, b"R2", mate_seq)? {
            rewritten.push(*b"R2");
        }
        if rewrite_string_if_present(&mut record.source, b"Q2", mate_qual)? {
            rewritten.push(*b"Q2");
        }
    }

    for tag in tags {
        let action = if removed.contains(&tag) {
            MetadataAction::Removed
        } else if rewritten.contains(&tag) {
            MetadataAction::Rewritten
        } else {
            MetadataAction::Preserved
        };
        stats.observe(action, &tag)?;
    }

    let seq = &record.read.seq[trim_start..];
    let normalized_qual = &record.read.qual[trim_start..];
    let numeric_qual = normalized_qual
        .iter()
        .map(|quality| {
            quality
                .checked_sub(33)
                .ok_or("Internal error: normalized CRAM quality below Phred+33")
        })
        .collect::<Result<Vec<_>, _>>()?;
    let qname = record.source.qname().to_vec();
    record.source.set(&qname, None, seq, &numeric_qual);
    Ok(())
}

fn rewrite_trimmed_string(
    record: &mut bam::Record,
    tag: &[u8; 2],
    trim_start: usize,
    original_len: usize,
) -> Result<bool, String> {
    let Some(value) = aux_string(record, tag)? else {
        return Ok(false);
    };
    if value.len() != original_len {
        return Err(format!(
            "CRAM tag {} has length {}, expected sequence length {original_len}",
            String::from_utf8_lossy(tag),
            value.len()
        ));
    }
    rewrite_string(record, tag, &value.as_bytes()[trim_start..])?;
    Ok(true)
}

fn rewrite_string_if_present(
    record: &mut bam::Record,
    tag: &[u8; 2],
    value: &[u8],
) -> Result<bool, String> {
    if aux_string(record, tag)?.is_none() {
        return Ok(false);
    }
    rewrite_string(record, tag, value)?;
    Ok(true)
}

fn rewrite_string(record: &mut bam::Record, tag: &[u8; 2], value: &[u8]) -> Result<(), String> {
    let value = std::str::from_utf8(value).map_err(|_| {
        format!(
            "CRAM tag {} rewrite produced non-UTF-8 text",
            String::from_utf8_lossy(tag)
        )
    })?;
    record.update_aux(tag, Aux::String(value)).map_err(|error| {
        format!(
            "Could not rewrite CRAM tag {}: {error}",
            String::from_utf8_lossy(tag)
        )
    })
}

fn aux_string(record: &bam::Record, tag: &[u8; 2]) -> Result<Option<String>, String> {
    match record.aux(tag) {
        Ok(Aux::String(value)) => Ok(Some(value.to_string())),
        Ok(_) => Err(format!(
            "CRAM tag {} must have SAM type Z",
            String::from_utf8_lossy(tag)
        )),
        Err(rust_htslib::errors::Error::BamAuxTagNotFound) => Ok(None),
        Err(error) => Err(format!(
            "Could not read CRAM tag {}: {error}",
            String::from_utf8_lossy(tag)
        )),
    }
}

fn remove_if_present(record: &mut bam::Record, tag: &[u8; 2]) -> Result<bool, String> {
    match record.remove_aux(tag) {
        Ok(()) => Ok(true),
        Err(rust_htslib::errors::Error::BamAuxTagNotFound) => Ok(false),
        Err(error) => Err(format!(
            "Could not remove invalidated CRAM tag {}: {error}",
            String::from_utf8_lossy(tag)
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn descriptor_preflight_raises_only_the_soft_limit() {
        let plan = descriptor_plan(1_500, None, 256, 4_096).unwrap();
        assert_eq!(plan.required_soft, 1_564);
        assert_eq!(plan.raise_soft_to, Some(1_564));

        let adequate = descriptor_plan(1_500, None, 2_048, 4_096).unwrap();
        assert_eq!(adequate.required_soft, 1_564);
        assert_eq!(adequate.raise_soft_to, None);
    }

    #[test]
    fn descriptor_preflight_rejects_insufficient_hard_limit() {
        let error = descriptor_plan(1_500, None, 256, 1_563).unwrap_err();
        assert!(error.contains("requires a process file-descriptor limit of at least 1564"));
        assert!(error.contains("soft=256 hard=1563"));
    }

    #[test]
    fn descriptor_preflight_honors_explicit_max_open_files() {
        let error = descriptor_plan(1_500, Some(750), 4_096, 4_096).unwrap_err();
        assert!(error.contains("--max-open-files is 750"));
        assert!(error.contains("append/reopen is not supported"));
    }

    #[test]
    fn explicit_finalization_status_is_fallible() {
        assert!(check_close_status(Path::new("success.cram"), 0).is_ok());
        let error = check_close_status(Path::new("failed.cram"), -1).unwrap_err();
        assert!(error.contains("Could not finalize CRAM output 'failed.cram'"));
        assert!(error.contains("hts_close returned -1"));
    }

    #[test]
    fn header_grouping_accepts_standard_queryname_declarations() {
        assert!(validate_query_grouped_header(b"@HD\tVN:1.6\tSO:queryname\n").is_ok());
        assert!(validate_query_grouped_header(b"@HD\tVN:1.6\tSO:unsorted\tGO:query\n").is_ok());
        assert!(validate_query_grouped_header(b"@HD\tVN:1.6\tSO:coordinate\n").is_err());
        assert!(validate_query_grouped_header(b"@RG\tID:one\n").is_err());
    }

    #[test]
    fn program_id_does_not_collide() {
        assert_eq!(unique_program_id(b"@HD\tVN:1.6\n"), "plexless");
        assert_eq!(
            unique_program_id(b"@PG\tID:plexless\n@PG\tID:plexless.1\n"),
            "plexless.2"
        );
    }
}
