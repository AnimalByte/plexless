use std::collections::{HashMap, HashSet};
use std::fs::{self, File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use flate2::Compression;
use flate2::write::GzEncoder;

pub use crate::output::OutputMate;
use crate::output::{CompressionJob, OutputAccumulator, OutputKey, OutputLayout, OutputTarget};
use crate::qc::SampleQc;
use crate::samples::SampleSheet;

const DEFAULT_MAX_OPEN_FILES: usize = 256;
const RESERVED_FILE_DESCRIPTORS: usize = 64;
pub(crate) const INCOMPLETE_RUN_MARKER: &str = "PLEXLESS_INCOMPLETE";

/// Durable evidence that an output directory must not be treated as complete.
///
/// The marker is deliberately not removed by `Drop`: abandoning an active run
/// leaves it behind. The run orchestrator removes it explicitly only after
/// output, reconciliation, and final reports succeed.
pub(crate) struct OutputRunMarker {
    path: PathBuf,
}

impl OutputRunMarker {
    pub(crate) fn begin(output_dir: &Path) -> Result<Self, String> {
        let path = output_dir.join(INCOMPLETE_RUN_MARKER);
        let mut file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .map_err(|error| {
                format!(
                    "Could not create incomplete-run marker '{}': {error}",
                    path.display()
                )
            })?;
        writeln!(
            file,
            "This plexless output is incomplete. Do not use it as a completed run."
        )
        .map_err(|error| format!("Could not write '{}': {error}", path.display()))?;
        file.sync_all()
            .map_err(|error| format!("Could not flush '{}': {error}", path.display()))?;
        Ok(Self { path })
    }

    pub(crate) fn complete(self) -> Result<(), String> {
        fs::remove_file(&self.path).map_err(|error| {
            format!(
                "Could not remove incomplete-run marker '{}': {error}",
                self.path.display()
            )
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum WriterKey {
    Sample { sample_id: u32, mate: OutputMate },

    Unassigned { mate: OutputMate },
}

impl WriterKey {
    fn mate(self) -> OutputMate {
        match self {
            Self::Sample { mate, .. } | Self::Unassigned { mate } => mate,
        }
    }
}

struct CachedRawWriter {
    writer: BufWriter<File>,
    last_used: u64,
}

struct CachedDirectWriter {
    writer: GzEncoder<BufWriter<File>>,
    last_used: u64,
}

pub struct WriterManager {
    accumulator: OutputAccumulator,
    writer: CompressedWriterManager,
    write_unassigned: bool,
    compression_level: u32,
    sample_qc: SampleQc,
    unassigned_records: u64,
    chunks: u64,
    uncompressed_bytes: u64,
    compressed_bytes: u64,
}

pub(crate) struct WriterCompletion {
    pub(crate) sample_qc: SampleQc,
    pub(crate) unassigned_records: u64,
    pub(crate) chunks: u64,
    pub(crate) uncompressed_bytes: u64,
    pub(crate) compressed_bytes: u64,
}

pub struct CompressedWriterManager {
    output_dir: PathBuf,
    sample_names: Vec<String>,
    cache: HashMap<WriterKey, CachedRawWriter>,
    max_open_writers: usize,
    write_unassigned: bool,
    paired: bool,
    tick: u64,
}

pub(crate) struct DirectWriterManager {
    output_dir: PathBuf,
    sample_names: Vec<String>,
    cache: HashMap<WriterKey, CachedDirectWriter>,
    max_open_writers: usize,
    write_unassigned: bool,
    paired: bool,
    compression_level: u32,
    tick: u64,
    sample_qc: SampleQc,
    unassigned_records: u64,
    chunks: u64,
    uncompressed_bytes: u64,
    compressed_bytes: u64,
}

pub(crate) fn resolve_max_open_files(
    expected_streams: usize,
    requested: Option<usize>,
) -> Result<usize, String> {
    resolve_max_open_files_with_limit(expected_streams, requested, open_file_soft_limit())
}

fn resolve_max_open_files_with_limit(
    expected_streams: usize,
    requested: Option<usize>,
    soft_limit: Option<usize>,
) -> Result<usize, String> {
    if expected_streams == 0 {
        return Err("Expected output stream count must be greater than 0".into());
    }
    let safe_limit = soft_limit
        .map(|limit| limit.saturating_sub(RESERVED_FILE_DESCRIPTORS).max(1))
        .unwrap_or(DEFAULT_MAX_OPEN_FILES);
    if let Some(requested) = requested {
        if requested == 0 {
            return Err("Maximum open files must be greater than 0".into());
        }
        if requested > safe_limit {
            return Err(format!(
                "Requested --max-open-files {requested} exceeds the safe process limit {safe_limit}"
            ));
        }
        return Ok(requested.min(expected_streams));
    }
    Ok(expected_streams
        .min(DEFAULT_MAX_OPEN_FILES)
        .min(safe_limit)
        .max(1))
}

fn open_file_soft_limit() -> Option<usize> {
    let limits = fs::read_to_string("/proc/self/limits").ok()?;
    let line = limits
        .lines()
        .find(|line| line.starts_with("Max open files"))?;
    line["Max open files".len()..]
        .split_ascii_whitespace()
        .next()?
        .parse()
        .ok()
}

impl WriterManager {
    pub fn new(
        output_dir: PathBuf,
        samples: &SampleSheet,
        paired: bool,
        write_unassigned: bool,
        max_open_writers: usize,
        compression_level: u32,
    ) -> Result<Self, String> {
        Self::new_with_policy(
            output_dir,
            samples,
            paired,
            write_unassigned,
            max_open_writers,
            compression_level,
            1024 * 1024,
            64 * 1024 * 1024,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new_with_policy(
        output_dir: PathBuf,
        samples: &SampleSheet,
        paired: bool,
        write_unassigned: bool,
        max_open_writers: usize,
        compression_level: u32,
        chunk_size: usize,
        memory_budget: usize,
    ) -> Result<Self, String> {
        let layout = OutputLayout::new(samples.samples.len(), paired, write_unassigned);
        let accumulator = OutputAccumulator::new(layout, chunk_size, memory_budget)?;
        let writer = CompressedWriterManager::new(
            output_dir,
            samples,
            paired,
            write_unassigned,
            max_open_writers,
            compression_level,
        )?;
        Ok(Self {
            accumulator,
            writer,
            write_unassigned,
            compression_level,
            sample_qc: SampleQc::new(samples.samples.len(), paired),
            unassigned_records: 0,
            chunks: 0,
            uncompressed_bytes: 0,
            compressed_bytes: 0,
        })
    }

    pub fn write_sample(
        &mut self,
        sample_id: u32,
        mate: OutputMate,
        id: &[u8],
        seq: &[u8],
        qual: &[u8],
    ) -> Result<(), String> {
        let key = OutputKey {
            target: OutputTarget::Sample { sample_id },
            mate,
        };

        self.write_record(key, id, seq, qual)?;
        self.sample_qc.observe_output(
            sample_id,
            mate,
            1,
            u64::try_from(seq.len()).map_err(|_| "Read length overflow")?,
        )
    }

    pub fn write_unassigned(
        &mut self,
        mate: OutputMate,
        id: &[u8],
        seq: &[u8],
        qual: &[u8],
    ) -> Result<(), String> {
        if !self.write_unassigned {
            return Ok(());
        }

        let key = OutputKey {
            target: OutputTarget::Unassigned,
            mate,
        };

        self.write_record(key, id, seq, qual)?;
        self.unassigned_records = self
            .unassigned_records
            .checked_add(1)
            .ok_or("Unassigned output record count overflow")?;
        Ok(())
    }

    pub fn finish(self) -> Result<(), String> {
        self.finish_with_qc().map(|_| ())
    }

    pub(crate) fn finish_with_qc(mut self) -> Result<WriterCompletion, String> {
        let jobs = self.accumulator.finish()?;
        self.write_jobs(jobs)?;
        self.writer.finish()?;
        Ok(WriterCompletion {
            sample_qc: self.sample_qc,
            unassigned_records: self.unassigned_records,
            chunks: self.chunks,
            uncompressed_bytes: self.uncompressed_bytes,
            compressed_bytes: self.compressed_bytes,
        })
    }

    fn write_record(
        &mut self,
        key: OutputKey,
        id: &[u8],
        seq: &[u8],
        qual: &[u8],
    ) -> Result<(), String> {
        let jobs = self.accumulator.append_record(key, id, seq, qual)?;
        self.write_jobs(jobs)
    }

    fn write_jobs(&mut self, jobs: Vec<CompressionJob>) -> Result<(), String> {
        for job in jobs {
            let uncompressed_len = job.fastq_bytes.len();
            let mut encoder = GzEncoder::new(Vec::new(), Compression::new(self.compression_level));
            encoder
                .write_all(&job.fastq_bytes)
                .map_err(|error| format!("Could not compress FASTQ chunk: {error}"))?;
            let member = encoder
                .finish()
                .map_err(|error| format!("Could not finish FASTQ gzip member: {error}"))?;
            self.chunks = self.chunks.checked_add(1).ok_or("Chunk count overflow")?;
            self.uncompressed_bytes = self
                .uncompressed_bytes
                .checked_add(
                    u64::try_from(uncompressed_len)
                        .map_err(|_| "Uncompressed output size overflow")?,
                )
                .ok_or("Uncompressed output size overflow")?;
            self.compressed_bytes = self
                .compressed_bytes
                .checked_add(
                    u64::try_from(member.len()).map_err(|_| "Compressed output size overflow")?,
                )
                .ok_or("Compressed output size overflow")?;
            match job.key.target {
                OutputTarget::Sample { sample_id } => {
                    self.writer
                        .append_sample_member(sample_id, job.key.mate, &member)?;
                }
                OutputTarget::Unassigned => {
                    self.writer
                        .append_unassigned_member(job.key.mate, &member)?;
                }
            }
        }
        Ok(())
    }
}

impl DirectWriterManager {
    pub(crate) fn new(
        output_dir: PathBuf,
        samples: &SampleSheet,
        paired: bool,
        write_unassigned: bool,
        max_open_writers: usize,
        compression_level: u32,
    ) -> Result<Self, String> {
        if max_open_writers == 0 {
            return Err("Maximum open writers must be greater than 0".into());
        }
        if compression_level > 9 {
            return Err("Compression level must be between 0 and 9".into());
        }
        let sample_names = build_sample_names(samples, write_unassigned)?;
        prepare_output_dir(&output_dir)?;
        Ok(Self {
            output_dir,
            sample_names,
            cache: HashMap::new(),
            max_open_writers,
            write_unassigned,
            paired,
            compression_level,
            tick: 0,
            sample_qc: SampleQc::new(samples.samples.len(), paired),
            unassigned_records: 0,
            chunks: 0,
            uncompressed_bytes: 0,
            compressed_bytes: 0,
        })
    }

    pub(crate) fn write_sample(
        &mut self,
        sample_id: u32,
        mate: OutputMate,
        id: &[u8],
        seq: &[u8],
        qual: &[u8],
    ) -> Result<(), String> {
        let key = WriterKey::Sample { sample_id, mate };
        self.write_record(key, id, seq, qual)?;
        self.sample_qc.observe_output(
            sample_id,
            mate,
            1,
            u64::try_from(seq.len()).map_err(|_| "Read length overflow")?,
        )
    }

    pub(crate) fn write_unassigned(
        &mut self,
        mate: OutputMate,
        id: &[u8],
        seq: &[u8],
        qual: &[u8],
    ) -> Result<(), String> {
        if !self.write_unassigned {
            return Ok(());
        }
        let key = WriterKey::Unassigned { mate };
        self.write_record(key, id, seq, qual)?;
        self.unassigned_records = self
            .unassigned_records
            .checked_add(1)
            .ok_or("Unassigned output record count overflow")?;
        Ok(())
    }

    pub(crate) fn finish_with_qc(mut self) -> Result<WriterCompletion, String> {
        for (_, cached) in self.cache.drain() {
            finish_direct_writer(cached.writer)?;
        }
        self.compressed_bytes = compressed_output_size(&self.output_dir)?;
        Ok(WriterCompletion {
            sample_qc: self.sample_qc,
            unassigned_records: self.unassigned_records,
            chunks: self.chunks,
            uncompressed_bytes: self.uncompressed_bytes,
            compressed_bytes: self.compressed_bytes,
        })
    }

    fn write_record(
        &mut self,
        key: WriterKey,
        id: &[u8],
        seq: &[u8],
        qual: &[u8],
    ) -> Result<(), String> {
        validate_mate(self.paired, key.mate())?;
        if let WriterKey::Sample { sample_id, .. } = key {
            validate_sample_id(&self.sample_names, sample_id)?;
        }
        if seq.len() != qual.len() {
            return Err("Sequence and quality lengths do not match".into());
        }

        self.tick = self.tick.wrapping_add(1);
        let tick = self.tick;
        if !self.cache.contains_key(&key) {
            self.open_writer(key)?;
        }
        let cached = self.cache.get_mut(&key).ok_or("Writer cache failure")?;
        cached.last_used = tick;
        let id = id.strip_prefix(b"@").unwrap_or(id);
        let record_len = serialized_record_len(id, seq, qual)?;
        write_fastq_record(&mut cached.writer, id, seq, qual)?;
        self.uncompressed_bytes = self
            .uncompressed_bytes
            .checked_add(u64::try_from(record_len).map_err(|_| "FASTQ record size overflow")?)
            .ok_or("Uncompressed output size overflow")?;
        Ok(())
    }

    fn open_writer(&mut self, key: WriterKey) -> Result<(), String> {
        if self.cache.len() >= self.max_open_writers {
            self.evict_oldest()?;
        }
        let path = output_path(&self.output_dir, &self.sample_names, key)?;
        let file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&path)
            .map_err(|error| format!("Could not open output '{}': {error}", path.display()))?;
        self.cache.insert(
            key,
            CachedDirectWriter {
                writer: GzEncoder::new(
                    BufWriter::new(file),
                    Compression::new(self.compression_level),
                ),
                last_used: self.tick,
            },
        );
        self.chunks = self.chunks.checked_add(1).ok_or("Chunk count overflow")?;
        Ok(())
    }

    fn evict_oldest(&mut self) -> Result<(), String> {
        let oldest_key = self
            .cache
            .iter()
            .min_by_key(|(_, cached)| cached.last_used)
            .map(|(&key, _)| key)
            .ok_or("Writer cache is empty")?;
        let cached = self
            .cache
            .remove(&oldest_key)
            .ok_or("Writer cache failure")?;
        finish_direct_writer(cached.writer)
    }
}

impl CompressedWriterManager {
    pub fn new(
        output_dir: PathBuf,
        samples: &SampleSheet,
        paired: bool,
        write_unassigned: bool,
        max_open_writers: usize,
        compression_level: u32,
    ) -> Result<Self, String> {
        if max_open_writers == 0 {
            return Err("Maximum open writers must be greater than 0".into());
        }

        if compression_level > 9 {
            return Err("Compression level must be between 0 and 9".into());
        }

        let sample_names = build_sample_names(samples, write_unassigned)?;
        prepare_output_dir(&output_dir)?;

        Ok(Self {
            output_dir,
            sample_names,
            cache: HashMap::new(),
            max_open_writers,
            write_unassigned,
            paired,
            tick: 0,
        })
    }

    pub fn append_sample_member(
        &mut self,
        sample_id: u32,
        mate: OutputMate,
        member: &[u8],
    ) -> Result<(), String> {
        validate_mate(self.paired, mate)?;
        validate_sample_id(&self.sample_names, sample_id)?;

        let key = WriterKey::Sample { sample_id, mate };
        self.append_member(key, member)
    }

    pub fn append_unassigned_member(
        &mut self,
        mate: OutputMate,
        member: &[u8],
    ) -> Result<(), String> {
        if !self.write_unassigned {
            return Ok(());
        }

        validate_mate(self.paired, mate)?;

        let key = WriterKey::Unassigned { mate };
        self.append_member(key, member)
    }

    pub fn finish(mut self) -> Result<(), String> {
        for (_, cached) in self.cache.drain() {
            finish_raw_writer(cached.writer)?;
        }

        Ok(())
    }

    fn append_member(&mut self, key: WriterKey, member: &[u8]) -> Result<(), String> {
        if member.is_empty() {
            return Ok(());
        }

        self.tick = self.tick.wrapping_add(1);
        let tick = self.tick;

        if !self.cache.contains_key(&key) {
            self.open_writer(key)?;
        }

        let cached = self.cache.get_mut(&key).ok_or("Writer cache failure")?;
        cached.last_used = tick;
        cached
            .writer
            .write_all(member)
            .map_err(|e| format!("Could not append compressed FASTQ member: {e}"))
    }

    fn open_writer(&mut self, key: WriterKey) -> Result<(), String> {
        if self.cache.len() >= self.max_open_writers {
            self.evict_oldest()?;
        }

        let path = output_path(&self.output_dir, &self.sample_names, key)?;

        let file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&path)
            .map_err(|e| format!("Could not open output '{}': {e}", path.display()))?;

        self.cache.insert(
            key,
            CachedRawWriter {
                writer: BufWriter::new(file),
                last_used: self.tick,
            },
        );

        Ok(())
    }

    fn evict_oldest(&mut self) -> Result<(), String> {
        let oldest_key = self
            .cache
            .iter()
            .min_by_key(|(_, cached)| cached.last_used)
            .map(|(&key, _)| key)
            .ok_or("Writer cache is empty")?;

        let cached = self
            .cache
            .remove(&oldest_key)
            .ok_or("Writer cache failure")?;

        finish_raw_writer(cached.writer)
    }
}

pub(crate) fn build_sample_names(
    samples: &SampleSheet,
    write_unassigned: bool,
) -> Result<Vec<String>, String> {
    let mut sample_names = Vec::with_capacity(samples.samples.len());
    let mut seen_names = HashSet::new();

    for sample in &samples.samples {
        let name = sanitize_name(&sample.name);

        if name.is_empty() {
            return Err(format!(
                "Sample '{}' cannot be converted to a valid filename",
                sample.name
            ));
        }

        if write_unassigned && name == "unassigned" {
            return Err(format!(
                "Sample '{}' conflicts with the reserved unassigned output filename",
                sample.name
            ));
        }

        if !seen_names.insert(name.clone()) {
            return Err(format!(
                "Sample filename collision after sanitizing '{}'",
                sample.name
            ));
        }

        sample_names.push(name);
    }

    Ok(sample_names)
}

fn validate_sample_id(sample_names: &[String], sample_id: u32) -> Result<(), String> {
    let sample_index = usize::try_from(sample_id).map_err(|_| "Invalid sample ID")?;

    if sample_index >= sample_names.len() {
        return Err(format!("Unknown sample ID {sample_id}"));
    }

    Ok(())
}

fn validate_mate(paired: bool, mate: OutputMate) -> Result<(), String> {
    match (paired, mate) {
        (false, OutputMate::Single) | (true, OutputMate::R1 | OutputMate::R2) => Ok(()),
        (false, _) => Err("Single-end output requires OutputMate::Single".into()),
        (true, OutputMate::Single) => Err("Paired-end output requires R1 or R2".into()),
    }
}

fn output_path(
    output_dir: &Path,
    sample_names: &[String],
    key: WriterKey,
) -> Result<PathBuf, String> {
    let (name, mate) = match key {
        WriterKey::Sample { sample_id, mate } => {
            let index = usize::try_from(sample_id).map_err(|_| "Invalid sample ID")?;
            let name = sample_names.get(index).ok_or("Unknown sample ID")?;
            (name.as_str(), mate)
        }
        WriterKey::Unassigned { mate } => ("unassigned", mate),
    };

    let filename = match mate {
        OutputMate::Single => format!("{name}.fastq.gz"),
        OutputMate::R1 => format!("{name}_R1.fastq.gz"),
        OutputMate::R2 => format!("{name}_R2.fastq.gz"),
    };

    Ok(output_dir.join(filename))
}

fn finish_raw_writer(mut writer: BufWriter<File>) -> Result<(), String> {
    writer
        .flush()
        .map_err(|e| format!("Could not flush output: {e}"))
}

fn finish_direct_writer(writer: GzEncoder<BufWriter<File>>) -> Result<(), String> {
    let mut buffered = writer
        .finish()
        .map_err(|error| format!("Could not finish gzip output: {error}"))?;
    buffered
        .flush()
        .map_err(|error| format!("Could not flush output: {error}"))
}

fn compressed_output_size(output_dir: &Path) -> Result<u64, String> {
    fs::read_dir(output_dir)
        .map_err(|error| {
            format!(
                "Could not read output directory '{}': {error}",
                output_dir.display()
            )
        })?
        .try_fold(0u64, |total, entry| {
            let entry = entry.map_err(|error| format!("Could not read output entry: {error}"))?;
            let is_fastq = entry
                .file_name()
                .to_str()
                .is_some_and(|name| name.ends_with(".fastq.gz"));
            if !is_fastq {
                return Ok(total);
            }
            let bytes = entry
                .metadata()
                .map_err(|error| format!("Could not inspect compressed output: {error}"))?
                .len();
            total
                .checked_add(bytes)
                .ok_or_else(|| "Compressed output size overflow".into())
        })
}

fn serialized_record_len(id: &[u8], seq: &[u8], qual: &[u8]) -> Result<usize, String> {
    1usize
        .checked_add(id.len())
        .and_then(|value| value.checked_add(1))
        .and_then(|value| value.checked_add(seq.len()))
        .and_then(|value| value.checked_add(3))
        .and_then(|value| value.checked_add(qual.len()))
        .and_then(|value| value.checked_add(1))
        .ok_or_else(|| "Serialized FASTQ record size overflow".into())
}

fn write_fastq_record(
    writer: &mut impl Write,
    id: &[u8],
    seq: &[u8],
    qual: &[u8],
) -> Result<(), String> {
    for part in [b"@".as_slice(), id, b"\n", seq, b"\n+\n", qual, b"\n"] {
        writer
            .write_all(part)
            .map_err(|error| format!("Could not write FASTQ: {error}"))?;
    }
    Ok(())
}

pub(crate) fn prepare_output_dir(path: &Path) -> Result<(), String> {
    if path.exists() {
        if !path.is_dir() {
            return Err(format!(
                "Output path '{}' exists but is not a directory",
                path.display()
            ));
        }

        let incomplete_marker = path.join(INCOMPLETE_RUN_MARKER);
        if incomplete_marker.exists() {
            return Err(format!(
                "Output directory '{}' contains {INCOMPLETE_RUN_MARKER} from an incomplete run",
                path.display()
            ));
        }

        if fs::read_dir(path)
            .map_err(|e| format!("Could not read output directory '{}': {e}", path.display()))?
            .next()
            .is_some()
        {
            return Err(format!(
                "Output directory '{}' is not empty",
                path.display()
            ));
        }
    } else {
        fs::create_dir_all(path).map_err(|e| {
            format!(
                "Could not create output directory '{}': {e}",
                path.display()
            )
        })?;
    }

    Ok(())
}

fn sanitize_name(name: &str) -> String {
    name.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, '-' | '_' | '.') {
                c
            } else {
                '_'
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn test_output_dir(name: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("plexless-{name}-{}-{nonce}", std::process::id()))
    }

    #[test]
    fn adaptive_writer_limit_reserves_descriptors_and_caps_default() {
        assert_eq!(
            resolve_max_open_files_with_limit(10, None, Some(128)),
            Ok(10)
        );
        assert_eq!(
            resolve_max_open_files_with_limit(1_000, None, Some(128)),
            Ok(64)
        );
        assert_eq!(
            resolve_max_open_files_with_limit(1_000, None, Some(10_000)),
            Ok(DEFAULT_MAX_OPEN_FILES)
        );
        assert!(resolve_max_open_files_with_limit(1_000, Some(65), Some(128)).is_err());
        assert_eq!(
            resolve_max_open_files_with_limit(1_000, Some(32), Some(128)),
            Ok(32)
        );
    }

    #[test]
    fn incomplete_marker_survives_drop_and_blocks_output_reuse() {
        let output = test_output_dir("incomplete-marker");
        fs::create_dir_all(&output).unwrap();
        let marker = OutputRunMarker::begin(&output).unwrap();
        assert!(output.join(INCOMPLETE_RUN_MARKER).is_file());
        drop(marker);
        let error = prepare_output_dir(&output).unwrap_err();
        assert!(error.contains(INCOMPLETE_RUN_MARKER));
        fs::remove_dir_all(output).unwrap();
    }

    #[test]
    fn completing_run_removes_incomplete_marker() {
        let output = test_output_dir("complete-marker");
        fs::create_dir_all(&output).unwrap();
        OutputRunMarker::begin(&output).unwrap().complete().unwrap();
        assert!(!output.join(INCOMPLETE_RUN_MARKER).exists());
        assert!(prepare_output_dir(&output).is_ok());
        fs::remove_dir_all(output).unwrap();
    }
}
