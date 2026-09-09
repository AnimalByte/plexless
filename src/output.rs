use std::fs;

use clap::ValueEnum;

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
pub(crate) enum OutputTarget {
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
}

impl OutputLayout {
    pub(crate) fn new(sample_count: usize, paired: bool, write_unassigned: bool) -> Self {
        Self {
            sample_count,
            paired,
            write_unassigned,
        }
    }

    pub(crate) fn stream_count(self) -> Result<usize, String> {
        let mates = if self.paired { 2 } else { 1 };
        self.sample_count
            .checked_add(usize::from(self.write_unassigned))
            .and_then(|targets| targets.checked_mul(mates))
            .ok_or_else(|| "Output stream count overflow".to_string())
    }

    pub(crate) fn stream_index(self, key: OutputKey) -> Result<usize, String> {
        let target = match key.target {
            OutputTarget::Sample { sample_id } => {
                let index = usize::try_from(sample_id).map_err(|_| "Invalid sample ID")?;
                if index >= self.sample_count {
                    return Err(format!("Unknown sample ID {sample_id}"));
                }
                index
            }
            OutputTarget::Unassigned if self.write_unassigned => self.sample_count,
            OutputTarget::Unassigned => {
                return Err("Unassigned output is disabled".into());
            }
        };

        if self.paired {
            let mate = match key.mate {
                OutputMate::R1 => 0,
                OutputMate::R2 => 1,
                OutputMate::Single => {
                    return Err("Paired-end output requires R1 or R2".into());
                }
            };
            target
                .checked_mul(2)
                .and_then(|value| value.checked_add(mate))
                .ok_or_else(|| "Output stream index overflow".to_string())
        } else {
            if key.mate != OutputMate::Single {
                return Err("Single-end output requires OutputMate::Single".into());
            }
            Ok(target)
        }
    }

    pub(crate) fn key_at(self, index: usize) -> Result<OutputKey, String> {
        if index >= self.stream_count()? {
            return Err("Output stream index is out of range".into());
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

#[derive(Debug)]
pub(crate) struct CompressionJob {
    pub(crate) key: OutputKey,
    pub(crate) chunk_id: u64,
    pub(crate) fastq_bytes: Vec<u8>,
}

#[derive(Debug, Default)]
struct Accumulator {
    bytes: Vec<u8>,
    next_chunk_id: u64,
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
        })
    }

    pub(crate) fn append(
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
        let accumulator = self
            .buffers
            .get_mut(index)
            .ok_or("Output accumulator index is missing")?;
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
        self.flush_ready(index)
    }

    fn flush_ready(&mut self, appended_index: usize) -> Result<Vec<CompressionJob>, String> {
        let mut jobs = Vec::new();
        if self.buffers[appended_index].bytes.len() >= self.chunk_size {
            jobs.push(self.flush_index(appended_index)?);
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
            jobs.push(self.flush_index(largest)?);
        }
        Ok(jobs)
    }

    pub(crate) fn finish(&mut self) -> Result<Vec<CompressionJob>, String> {
        let mut jobs = Vec::new();
        for index in 0..self.buffers.len() {
            if !self.buffers[index].bytes.is_empty() {
                jobs.push(self.flush_index(index)?);
            }
        }
        debug_assert_eq!(self.buffered_bytes, 0);
        Ok(jobs)
    }

    fn flush_index(&mut self, index: usize) -> Result<CompressionJob, String> {
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
        let chunk_id = accumulator.next_chunk_id;
        accumulator.next_chunk_id = accumulator
            .next_chunk_id
            .checked_add(1)
            .ok_or("Output chunk ID overflow")?;
        Ok(CompressionJob {
            key,
            chunk_id,
            fastq_bytes,
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
        assert!(accumulator.append(key, vec![b'C'; 3]).unwrap().is_empty());
        let final_jobs = accumulator.finish().unwrap();
        assert_eq!(final_jobs[0].chunk_id, 1);
        assert_eq!(final_jobs[0].fastq_bytes, b"CCC");
    }
}
