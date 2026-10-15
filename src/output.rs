use std::fs;

use clap::ValueEnum;

pub(crate) const FASTQ_BUFFERED_CHUNK_TARGET: usize = 32 * 1024;
pub(crate) const AUTO_MIN_CHUNK_SIZE: usize = 128 * 1024;
pub(crate) const AUTO_MAX_CHUNK_SIZE: usize = 1024 * 1024;
const AUTO_MIN_BUFFER_BUDGET: usize = 16 * 1024 * 1024;
const AUTO_MAX_BUFFER_BUDGET: usize = 2 * 1024 * 1024 * 1024;
const MIN_EXPLICIT_CHUNK_SIZE: usize = 32 * 1024;
const MAX_EXPLICIT_CHUNK_SIZE: usize = 64 * 1024 * 1024;
const MIN_EXPLICIT_BUFFER_BUDGET: usize = 1024 * 1024;

/// Use one stable output-stream threshold so paired and unassigned outputs are
/// accounted for without coupling selection to sample count, input layout, or
/// machine-specific timing.
pub(crate) const AUTO_BUFFERED_STREAM_THRESHOLD: usize = 384;

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum OutputMode {
    Auto,
    Direct,
    Buffered,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ResolvedOutputMode {
    Direct,
    Buffered,
}

impl ResolvedOutputMode {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::Direct => "direct",
            Self::Buffered => "buffered",
        }
    }
}

pub(crate) fn report_output_mode(
    requested: OutputMode,
    resolved: ResolvedOutputMode,
    expected_streams: usize,
) {
    match requested {
        OutputMode::Auto => eprintln!(
            "Output mode: auto -> {} ({expected_streams} expected streams)",
            resolved.as_str()
        ),
        OutputMode::Direct | OutputMode::Buffered => {
            eprintln!("Output mode: {}", resolved.as_str());
        }
    }
}

pub(crate) fn resolve_output_mode(
    requested: OutputMode,
    expected_streams: usize,
) -> Result<ResolvedOutputMode, String> {
    if expected_streams == 0 {
        return Err("Expected output stream count must be greater than 0".into());
    }
    Ok(match requested {
        OutputMode::Auto if expected_streams < AUTO_BUFFERED_STREAM_THRESHOLD => {
            ResolvedOutputMode::Direct
        }
        OutputMode::Auto => ResolvedOutputMode::Buffered,
        OutputMode::Direct => ResolvedOutputMode::Direct,
        OutputMode::Buffered => ResolvedOutputMode::Buffered,
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum OutputMate {
    Single,
    R1,
    R2,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub(crate) enum AuxiliaryStream {
    Failed = 0,
    Unpaired1 = 1,
    Unpaired2 = 2,
    Overlap = 3,
    Stdout = 4,
}
impl AuxiliaryStream {
    pub(crate) const ALL: [Self; 5] = [
        Self::Failed,
        Self::Unpaired1,
        Self::Unpaired2,
        Self::Overlap,
        Self::Stdout,
    ];
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Failed => "failed",
            Self::Unpaired1 => "unpaired_R1",
            Self::Unpaired2 => "unpaired_R2",
            Self::Overlap => "overlapped",
            Self::Stdout => "stdout",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(crate) enum OutputTarget {
    Merged { sample_id: u32 },
    Auxiliary { stream: AuxiliaryStream },
    Sample { sample_id: u32 },
    Unassigned,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(crate) struct OutputKey {
    pub(crate) target: OutputTarget,
    pub(crate) mate: OutputMate,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct OutputLayout {
    sample_count: usize,
    paired: bool,
    write_unassigned: bool,
    merged: bool,
    auxiliary: u8,
    plain_primary: u8,
    plain_merged: bool,
    plain_auxiliary: u8,
    split_records: u64,
    shared_unpaired: bool,
}

impl OutputLayout {
    pub(crate) fn new(sample_count: usize, paired: bool, write_unassigned: bool) -> Self {
        Self {
            sample_count,
            paired,
            write_unassigned,
            merged: false,
            auxiliary: 0,
            plain_primary: 0,
            plain_merged: false,
            plain_auxiliary: 0,
            split_records: 0,
            shared_unpaired: false,
        }
    }

    pub(crate) fn canonical_target(self, target: OutputTarget) -> OutputTarget {
        if self.shared_unpaired
            && target
                == (OutputTarget::Auxiliary {
                    stream: AuxiliaryStream::Unpaired2,
                })
        {
            OutputTarget::Auxiliary {
                stream: AuxiliaryStream::Unpaired1,
            }
        } else {
            target
        }
    }
    pub(crate) fn with_fastp(mut self, o: &crate::fastp_options::FastpOptions) -> Self {
        let plain = |path: Option<&std::path::Path>| {
            path.is_some_and(|p| {
                p.extension()
                    .is_none_or(|ext| !ext.eq_ignore_ascii_case("gz"))
            })
        };
        if plain(o.out1.as_deref()) {
            self.plain_primary |= 1;
        }
        if plain(o.out2.as_deref()) {
            self.plain_primary |= 2;
        }
        self.plain_merged = plain(o.merged_out.as_deref());
        for (stream, path) in [
            (AuxiliaryStream::Failed, o.failed_out.as_deref()),
            (AuxiliaryStream::Unpaired1, o.unpaired1.as_deref()),
            (AuxiliaryStream::Unpaired2, o.unpaired2.as_deref()),
            (AuxiliaryStream::Overlap, o.overlapped_out.as_deref()),
        ] {
            if plain(path) {
                self.plain_auxiliary |= 1 << stream as u8;
            }
        }
        if o.stdout {
            self.plain_auxiliary |= 1 << AuxiliaryStream::Stdout as u8;
        }
        self.shared_unpaired = o.unpaired1.is_some() && o.unpaired1 == o.unpaired2;
        self.split_records = o.split_records;
        self.merged = o.merge;
        for (stream, enabled) in [
            (AuxiliaryStream::Failed, o.failed_out.is_some()),
            (AuxiliaryStream::Unpaired1, o.unpaired1.is_some()),
            (AuxiliaryStream::Unpaired2, o.unpaired2.is_some()),
            (AuxiliaryStream::Overlap, o.overlapped_out.is_some()),
            (AuxiliaryStream::Stdout, o.stdout),
        ] {
            if enabled {
                self.auxiliary |= 1 << stream as u8;
            }
        }
        self
    }
    pub(crate) fn gzip(self, key: OutputKey) -> bool {
        match key.target {
            OutputTarget::Merged { .. } => !self.plain_merged,
            OutputTarget::Auxiliary { stream } => self.plain_auxiliary & (1 << stream as u8) == 0,
            OutputTarget::Sample { .. } => {
                self.plain_primary & if key.mate == OutputMate::R2 { 2 } else { 1 } == 0
            }
            OutputTarget::Unassigned => true,
        }
    }
    fn primary_count(self) -> Result<usize, String> {
        self.sample_count
            .checked_add(usize::from(self.write_unassigned))
            .and_then(|n| n.checked_mul(if self.paired { 2 } else { 1 }))
            .ok_or_else(|| "Output stream count overflow".into())
    }
    pub(crate) fn stream_count(self) -> Result<usize, String> {
        self.primary_count()?
            .checked_add(if self.merged { self.sample_count } else { 0 })
            .and_then(|n| n.checked_add(self.auxiliary.count_ones() as usize))
            .ok_or_else(|| "Output stream count overflow".into())
    }
    pub(crate) fn stream_index(self, key: OutputKey) -> Result<usize, String> {
        let target = match key.target {
            OutputTarget::Merged { sample_id } => {
                if !self.merged
                    || sample_id as usize >= self.sample_count
                    || key.mate != OutputMate::Single
                {
                    return Err("Invalid merged output stream".into());
                }
                return Ok(self.primary_count()? + sample_id as usize);
            }
            OutputTarget::Auxiliary { stream } => {
                let bit = 1u8 << stream as u8;
                if self.auxiliary & bit == 0 || key.mate != OutputMate::Single {
                    return Err("Auxiliary output is disabled or has invalid mate".into());
                }
                return Ok(self.primary_count()?
                    + if self.merged { self.sample_count } else { 0 }
                    + (self.auxiliary & (bit - 1)).count_ones() as usize);
            }
            OutputTarget::Sample { sample_id } => {
                if sample_id as usize >= self.sample_count {
                    return Err(format!("Unknown sample ID {sample_id}"));
                }
                sample_id as usize
            }
            OutputTarget::Unassigned if self.write_unassigned => self.sample_count,
            OutputTarget::Unassigned => return Err("Unassigned output is disabled".into()),
        };
        if self.paired {
            let mate = match key.mate {
                OutputMate::R1 => 0,
                OutputMate::R2 => 1,
                OutputMate::Single => return Err("Paired-end output requires R1 or R2".into()),
            };
            Ok(target * 2 + mate)
        } else if key.mate == OutputMate::Single {
            Ok(target)
        } else {
            Err("Single-end output requires OutputMate::Single".into())
        }
    }
    pub(crate) fn key_at(self, index: usize) -> Result<OutputKey, String> {
        if index >= self.stream_count()? {
            return Err("Output stream index is out of range".into());
        }
        let primary = self.primary_count()?;
        if index >= primary {
            let index = index - primary;
            if self.merged && index < self.sample_count {
                return Ok(OutputKey {
                    target: OutputTarget::Merged {
                        sample_id: u32::try_from(index).map_err(|_| "Too many samples")?,
                    },
                    mate: OutputMate::Single,
                });
            }
            let aux_index = index - if self.merged { self.sample_count } else { 0 };
            let stream = AuxiliaryStream::ALL
                .into_iter()
                .filter(|s| self.auxiliary & (1 << *s as u8) != 0)
                .nth(aux_index)
                .ok_or("Invalid auxiliary stream")?;
            return Ok(OutputKey {
                target: OutputTarget::Auxiliary { stream },
                mate: OutputMate::Single,
            });
        }
        let (target_index, mate) = if self.paired {
            (
                index / 2,
                if index.is_multiple_of(2) {
                    OutputMate::R1
                } else {
                    OutputMate::R2
                },
            )
        } else {
            (index, OutputMate::Single)
        };
        let target = if target_index < self.sample_count {
            OutputTarget::Sample {
                sample_id: u32::try_from(target_index).map_err(|_| "Too many samples")?,
            }
        } else {
            OutputTarget::Unassigned
        };
        Ok(OutputKey { target, mate })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct OutputBufferPolicy {
    pub(crate) chunk_size: usize,
    pub(crate) memory_budget: usize,
    pub(crate) expected_streams: usize,
    pub(crate) available_memory: usize,
}

pub(crate) fn available_memory_bytes() -> Option<usize> {
    let contents = fs::read_to_string("/proc/meminfo").ok()?;
    let kib = contents.lines().find_map(|line| {
        let mut fields = line.split_ascii_whitespace();
        (fields.next()? == "MemAvailable:")
            .then(|| fields.next()?.parse::<usize>().ok())
            .flatten()
    })?;
    let host_available = kib.checked_mul(1024)?;
    Some(
        cgroup_available_memory()
            .map(|available| available.min(host_available))
            .unwrap_or(host_available),
    )
}

fn cgroup_available_memory() -> Option<usize> {
    cgroup_available_memory_at("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory.current")
        .or_else(|| {
            cgroup_available_memory_at(
                "/sys/fs/cgroup/memory/memory.limit_in_bytes",
                "/sys/fs/cgroup/memory/memory.usage_in_bytes",
            )
        })
}

fn cgroup_available_memory_at(limit_path: &str, used_path: &str) -> Option<usize> {
    let limit_text = fs::read_to_string(limit_path).ok()?;
    if limit_text.trim() == "max" {
        return None;
    }
    let limit = limit_text.trim().parse::<usize>().ok()?;
    let used = fs::read_to_string(used_path)
        .ok()?
        .trim()
        .parse::<usize>()
        .ok()?;
    limit.checked_sub(used)
}

/// Reserves most available memory for input, queues, compression, writers, and
/// process/OS overhead. The accumulator itself is also pressure-flushed at the
/// returned budget, so `streams * chunk_size` is not preallocated or required.
pub(crate) fn calculate_output_buffer_policy(
    expected_streams: usize,
    available_memory: usize,
    chunk_override: Option<usize>,
    memory_override: Option<usize>,
) -> Result<OutputBufferPolicy, String> {
    if expected_streams == 0 {
        return Err("Expected output stream count must be greater than 0".into());
    }
    if available_memory == 0 {
        return Err("Available memory must be greater than 0".into());
    }
    let automatic_budget = (available_memory / 4)
        .clamp(
            AUTO_MIN_BUFFER_BUDGET.min(available_memory),
            AUTO_MAX_BUFFER_BUDGET,
        )
        .min((available_memory / 2).max(1));
    if let Some(memory) = memory_override {
        if memory < MIN_EXPLICIT_BUFFER_BUDGET {
            return Err("Explicit output buffer memory must be at least 1 MiB".into());
        }
        if memory > (available_memory / 2).max(1) {
            return Err(
                "Explicit output buffer memory exceeds half of currently available memory".into(),
            );
        }
    }
    let memory_budget = memory_override.unwrap_or(automatic_budget);

    let automatic_chunk =
        (memory_budget / expected_streams).clamp(AUTO_MIN_CHUNK_SIZE, AUTO_MAX_CHUNK_SIZE);
    if let Some(chunk) = chunk_override {
        if !(MIN_EXPLICIT_CHUNK_SIZE..=MAX_EXPLICIT_CHUNK_SIZE).contains(&chunk) {
            return Err("Explicit output chunk size must be between 32 KiB and 64 MiB".into());
        }
        if chunk > memory_budget {
            return Err("Output chunk size cannot exceed output buffer memory".into());
        }
    }
    let chunk_size = chunk_override.unwrap_or(automatic_chunk);

    Ok(OutputBufferPolicy {
        chunk_size,
        memory_budget,
        expected_streams,
        available_memory,
    })
}

/// FASTQ input uses the qualified fixed buffered chunk target unless the user
/// explicitly overrides it. Keeping this policy in one helper ensures serial
/// and parallel buffered FASTQ paths cannot drift apart.
pub(crate) fn calculate_fastq_output_buffer_policy(
    expected_streams: usize,
    available_memory: usize,
    chunk_override: Option<usize>,
    memory_override: Option<usize>,
) -> Result<OutputBufferPolicy, String> {
    calculate_output_buffer_policy(
        expected_streams,
        available_memory,
        Some(chunk_override.unwrap_or(FASTQ_BUFFERED_CHUNK_TARGET)),
        memory_override,
    )
}

#[derive(Debug)]
pub(crate) struct CompressionJob {
    pub(crate) part: u64,
    pub(crate) gzip: bool,
    pub(crate) key: OutputKey,
    pub(crate) chunk_id: u64,
    pub(crate) fastq_bytes: Vec<u8>,
    pub(crate) flush_reason: CompressionFlushReason,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CompressionFlushReason {
    ChunkThreshold,
    MemoryPressure,
    FinalDrain,
}

#[derive(Debug, Default)]
struct Accumulator {
    bytes: Vec<u8>,
    next_chunk_id: u64,
    records: u64,
    part: u64,
}

/// Persistent buffers are indexed by output identity and span input work
/// batches. `memory_budget` bounds their combined logical byte length; bounded
/// pipeline channels separately cap queued jobs and compressed members.
pub(crate) struct OutputAccumulator {
    layout: OutputLayout,
    buffers: Vec<Accumulator>,
    chunk_size: usize,
    memory_budget: usize,
    buffered_bytes: usize,
    active_buffers: usize,
}

impl OutputAccumulator {
    pub(crate) fn new(
        layout: OutputLayout,
        chunk_size: usize,
        memory_budget: usize,
    ) -> Result<Self, String> {
        if chunk_size == 0 || memory_budget == 0 {
            return Err("Output chunk and memory sizes must be greater than 0".into());
        }
        let stream_count = layout.stream_count()?;
        Ok(Self {
            layout,
            buffers: (0..stream_count).map(|_| Accumulator::default()).collect(),
            chunk_size,
            memory_budget,
            buffered_bytes: 0,
            active_buffers: 0,
        })
    }

    pub(crate) fn append(
        &mut self,
        key: OutputKey,
        bytes: Vec<u8>,
    ) -> Result<Vec<CompressionJob>, String> {
        if bytes.is_empty() {
            return Ok(Vec::new());
        }
        if self.layout.split_records == 0 {
            return self.append_bytes(key, bytes);
        }
        let index = self.layout.stream_index(key)?;
        let mut jobs = Vec::new();
        self.advance_part(index, &mut jobs)?;
        let mut begin = 0;
        let mut lines = 0usize;
        for (i, &byte) in bytes.iter().enumerate() {
            if byte != b'\n' {
                continue;
            }
            lines += 1;
            if lines.is_multiple_of(4) {
                self.buffers[index].records += 1;
                if self.buffers[index].records == self.layout.split_records && i + 1 < bytes.len() {
                    jobs.extend(self.append_bytes(key, bytes[begin..=i].to_vec())?);
                    begin = i + 1;
                    self.advance_part(index, &mut jobs)?;
                }
            }
        }
        if !lines.is_multiple_of(4) || bytes.last() != Some(&b'\n') {
            return Err("Split output received an incomplete FASTQ record".into());
        }
        if begin == 0 {
            jobs.extend(self.append_bytes(key, bytes)?);
        } else {
            jobs.extend(self.append_bytes(key, bytes[begin..].to_vec())?);
        }
        Ok(jobs)
    }
    fn advance_part(&mut self, index: usize, jobs: &mut Vec<CompressionJob>) -> Result<(), String> {
        if self.layout.split_records > 0 && self.buffers[index].records == self.layout.split_records
        {
            if !self.buffers[index].bytes.is_empty() {
                jobs.push(self.flush_index(index, CompressionFlushReason::ChunkThreshold)?);
            }
            self.buffers[index].part = self.buffers[index]
                .part
                .checked_add(1)
                .ok_or("Split part overflow")?;
            self.buffers[index].records = 0;
        }
        Ok(())
    }
    fn append_bytes(
        &mut self,
        key: OutputKey,
        mut bytes: Vec<u8>,
    ) -> Result<Vec<CompressionJob>, String> {
        if bytes.is_empty() {
            return Ok(Vec::new());
        }
        let index = self.layout.stream_index(key)?;
        self.buffered_bytes = self
            .buffered_bytes
            .checked_add(bytes.len())
            .ok_or("Output buffer byte count overflow")?;
        let buffer = self
            .buffers
            .get_mut(index)
            .ok_or("Output accumulator index is missing")?;
        if buffer.bytes.is_empty() {
            self.active_buffers = self
                .active_buffers
                .checked_add(1)
                .ok_or("Active output buffer count overflow")?;
            std::mem::swap(&mut buffer.bytes, &mut bytes);
        } else {
            buffer.bytes.extend_from_slice(&bytes);
        }

        self.flush_ready(index)
    }

    pub(crate) fn append_record(
        &mut self,
        key: OutputKey,
        id: &[u8],
        seq: &[u8],
        qual: &[u8],
    ) -> Result<Vec<CompressionJob>, String> {
        if seq.len() != qual.len() {
            return Err("Sequence and quality lengths do not match".into());
        }
        let id = id.strip_prefix(b"@").unwrap_or(id);
        let record_len = 1usize
            .checked_add(id.len())
            .and_then(|value| value.checked_add(1))
            .and_then(|value| value.checked_add(seq.len()))
            .and_then(|value| value.checked_add(3))
            .and_then(|value| value.checked_add(qual.len()))
            .and_then(|value| value.checked_add(1))
            .ok_or("Serialized FASTQ record size overflow")?;
        let index = self.layout.stream_index(key)?;
        let mut jobs = Vec::new();
        if self.layout.split_records > 0 {
            self.advance_part(index, &mut jobs)?;
            self.buffers[index].records += 1;
        }
        let accumulator = self
            .buffers
            .get_mut(index)
            .ok_or("Output accumulator index is missing")?;
        if accumulator.bytes.is_empty() {
            self.active_buffers = self
                .active_buffers
                .checked_add(1)
                .ok_or("Active output buffer count overflow")?;
        }
        accumulator.bytes.reserve(record_len);
        accumulator.bytes.extend_from_slice(b"@");
        accumulator.bytes.extend_from_slice(id);
        accumulator.bytes.extend_from_slice(b"\n");
        accumulator.bytes.extend_from_slice(seq);
        // The FASTQ parser does not expose source '+' metadata. Plexless
        // intentionally emits the standards-compliant normalized separator.
        accumulator.bytes.extend_from_slice(b"\n+\n");
        accumulator.bytes.extend_from_slice(qual);
        accumulator.bytes.extend_from_slice(b"\n");
        self.buffered_bytes = self
            .buffered_bytes
            .checked_add(record_len)
            .ok_or("Output buffer byte count overflow")?;
        jobs.extend(self.flush_ready(index)?);
        Ok(jobs)
    }

    fn flush_ready(&mut self, appended_index: usize) -> Result<Vec<CompressionJob>, String> {
        let mut jobs = Vec::new();
        if self.buffers[appended_index].bytes.len() >= self.chunk_size {
            jobs.push(self.flush_index(appended_index, CompressionFlushReason::ChunkThreshold)?);
        }
        while self.buffered_bytes > self.memory_budget {
            let largest = self
                .buffers
                .iter()
                .enumerate()
                .filter(|(_, buffer)| !buffer.bytes.is_empty())
                .max_by_key(|(_, buffer)| buffer.bytes.len())
                .map(|(index, _)| index)
                .ok_or("Output memory budget exceeded without buffered output")?;
            jobs.push(self.flush_index(largest, CompressionFlushReason::MemoryPressure)?);
        }
        Ok(jobs)
    }

    pub(crate) fn buffered_bytes(&self) -> usize {
        self.buffered_bytes
    }

    pub(crate) fn active_buffers(&self) -> usize {
        self.active_buffers
    }

    pub(crate) fn layout(&self) -> OutputLayout {
        self.layout
    }

    pub(crate) fn finish(&mut self) -> Result<Vec<CompressionJob>, String> {
        let mut jobs = Vec::new();
        for index in 0..self.buffers.len() {
            if !self.buffers[index].bytes.is_empty() {
                jobs.push(self.flush_index(index, CompressionFlushReason::FinalDrain)?);
            }
        }
        debug_assert_eq!(self.buffered_bytes, 0);
        debug_assert_eq!(self.active_buffers, 0);
        Ok(jobs)
    }

    fn flush_index(
        &mut self,
        index: usize,
        flush_reason: CompressionFlushReason,
    ) -> Result<CompressionJob, String> {
        let key = self.layout.key_at(index)?;
        let accumulator = self
            .buffers
            .get_mut(index)
            .ok_or("Output accumulator index is missing")?;
        let fastq_bytes = std::mem::take(&mut accumulator.bytes);
        self.buffered_bytes = self
            .buffered_bytes
            .checked_sub(fastq_bytes.len())
            .ok_or("Output buffer byte count underflow")?;
        self.active_buffers = self
            .active_buffers
            .checked_sub(1)
            .ok_or("Active output buffer count underflow")?;
        let chunk_id = accumulator.next_chunk_id;
        accumulator.next_chunk_id = accumulator
            .next_chunk_id
            .checked_add(1)
            .ok_or("Output chunk ID overflow")?;
        Ok(CompressionJob {
            part: accumulator.part,
            gzip: self.layout.gzip(key),
            key,
            chunk_id,
            fastq_bytes,
            flush_reason,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn auto_policy_scales_and_clamps_by_stream_count() {
        let gib = 1024 * 1024 * 1024;
        let small = calculate_output_buffer_policy(2, 16 * gib, None, None).unwrap();
        let medium = calculate_output_buffer_policy(768, 16 * gib, None, None).unwrap();
        let many = calculate_output_buffer_policy(6_000, 512 * 1024 * 1024, None, None).unwrap();
        assert_eq!(small.chunk_size, AUTO_MAX_CHUNK_SIZE);
        assert!((AUTO_MIN_CHUNK_SIZE..=AUTO_MAX_CHUNK_SIZE).contains(&medium.chunk_size));
        assert_eq!(many.chunk_size, AUTO_MIN_CHUNK_SIZE);
    }

    #[test]
    fn fastq_auto_policy_uses_qualified_chunk_target() {
        let mib = 1024 * 1024;
        let policy = calculate_fastq_output_buffer_policy(768, 8 * 1024 * mib, None, None).unwrap();
        assert_eq!(policy.chunk_size, FASTQ_BUFFERED_CHUNK_TARGET);

        let overridden = calculate_fastq_output_buffer_policy(
            768,
            8 * 1024 * mib,
            Some(64 * 1024),
            Some(512 * mib),
        )
        .unwrap();
        assert_eq!(overridden.chunk_size, 64 * 1024);
        assert_eq!(overridden.memory_budget, 512 * mib);
    }

    #[test]
    fn explicit_policy_honors_valid_values_and_rejects_unsafe_ones() {
        let mib = 1024 * 1024;
        let policy = calculate_output_buffer_policy(
            1_500,
            8 * 1024 * mib,
            Some(256 * 1024),
            Some(512 * mib),
        )
        .unwrap();
        assert_eq!(policy.chunk_size, 256 * 1024);
        assert_eq!(policy.memory_budget, 512 * mib);
        assert!(calculate_output_buffer_policy(1, 64 * mib, None, Some(40 * mib)).is_err());
        assert!(calculate_output_buffer_policy(1, 1024 * mib, Some(16 * 1024), None).is_err());
    }

    #[test]
    fn explicit_budget_may_exceed_two_gib_when_memory_allows() {
        let gib = 1024 * 1024 * 1024;
        let policy =
            calculate_output_buffer_policy(768, 32 * gib, Some(4 * 1024 * 1024), Some(8 * gib))
                .unwrap();
        assert_eq!(policy.memory_budget, 8 * gib);
        assert_eq!(policy.chunk_size, 4 * 1024 * 1024);
        assert!(calculate_output_buffer_policy(768, 16 * gib, None, Some(9 * gib)).is_err());
    }

    #[test]
    fn layout_counts_se_pe_and_unassigned_streams() {
        assert_eq!(OutputLayout::new(16, false, false).stream_count(), Ok(16));
        assert_eq!(OutputLayout::new(16, true, false).stream_count(), Ok(32));
        assert_eq!(OutputLayout::new(16, false, true).stream_count(), Ok(17));
        assert_eq!(OutputLayout::new(16, true, true).stream_count(), Ok(34));
        assert_eq!(OutputLayout::new(384, false, false).stream_count(), Ok(384));
        assert_eq!(OutputLayout::new(384, true, false).stream_count(), Ok(768));
        assert_eq!(OutputLayout::new(384, true, true).stream_count(), Ok(770));
    }

    #[test]
    fn auto_mode_uses_the_exact_stream_threshold() {
        assert_eq!(
            resolve_output_mode(OutputMode::Auto, AUTO_BUFFERED_STREAM_THRESHOLD - 1),
            Ok(ResolvedOutputMode::Direct)
        );
        assert_eq!(
            resolve_output_mode(OutputMode::Auto, AUTO_BUFFERED_STREAM_THRESHOLD),
            Ok(ResolvedOutputMode::Buffered)
        );
        assert_eq!(
            resolve_output_mode(OutputMode::Auto, AUTO_BUFFERED_STREAM_THRESHOLD + 1),
            Ok(ResolvedOutputMode::Buffered)
        );
    }

    #[test]
    fn forced_output_modes_bypass_auto_selection() {
        assert_eq!(
            resolve_output_mode(OutputMode::Direct, usize::MAX),
            Ok(ResolvedOutputMode::Direct)
        );
        assert_eq!(
            resolve_output_mode(OutputMode::Buffered, 1),
            Ok(ResolvedOutputMode::Buffered)
        );
    }

    #[test]
    fn memory_pressure_flushes_without_waiting_for_chunk_target() {
        let layout = OutputLayout::new(2, false, false);
        let mut accumulator = OutputAccumulator::new(layout, 1024, 10).unwrap();
        let key = OutputKey {
            target: OutputTarget::Sample { sample_id: 0 },
            mate: OutputMate::Single,
        };
        let jobs = accumulator.append(key, vec![b'A'; 11]).unwrap();
        assert_eq!(jobs.len(), 1);
        assert_eq!(jobs[0].chunk_id, 0);
        assert_eq!(jobs[0].fastq_bytes.len(), 11);
        assert_eq!(jobs[0].flush_reason, CompressionFlushReason::MemoryPressure);
    }

    #[test]
    fn accumulates_across_appends_and_sequences_chunks_per_output() {
        let layout = OutputLayout::new(1, false, false);
        let mut accumulator = OutputAccumulator::new(layout, 10, 100).unwrap();
        let key = OutputKey {
            target: OutputTarget::Sample { sample_id: 0 },
            mate: OutputMate::Single,
        };
        assert!(accumulator.append(key, vec![b'A'; 6]).unwrap().is_empty());
        let first = accumulator.append(key, vec![b'B'; 5]).unwrap();
        assert_eq!(first[0].chunk_id, 0);
        assert_eq!(first[0].fastq_bytes, b"AAAAAABBBBB");
        assert_eq!(
            first[0].flush_reason,
            CompressionFlushReason::ChunkThreshold
        );
        assert!(accumulator.append(key, vec![b'C'; 3]).unwrap().is_empty());
        let final_jobs = accumulator.finish().unwrap();
        assert_eq!(final_jobs[0].chunk_id, 1);
        assert_eq!(final_jobs[0].fastq_bytes, b"CCC");
        assert_eq!(
            final_jobs[0].flush_reason,
            CompressionFlushReason::FinalDrain
        );
    }
}
