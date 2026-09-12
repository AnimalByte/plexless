use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io::Write;
use std::path::Path;

use crate::barcodes::BarcodeCatalog;
use crate::cli::DemuxArgs;
use crate::fastq::InputReader;
use crate::input::{InputFiles, InputSource, resolve_input_source, resolve_inputs};
use crate::output::{
    OutputLayout, ResolvedOutputMode, available_memory_bytes, calculate_output_buffer_policy,
    report_output_mode, resolve_output_mode,
};
use crate::routing::{RouteResult, RoutingTree};
use crate::samples::SampleSheet;
use crate::stats::{FastqStats, MateQcStats};
use crate::structure::{Orientation, ReadLayout};
use crate::writer::{
    DirectWriterManager, OutputMate, OutputRunMarker, WriterCompletion, WriterManager,
    resolve_max_open_files,
};

const PAIRED_RESYNC_WINDOW: usize = 1024;

#[derive(Debug, Default)]
pub(crate) struct DemuxCounts {
    pub(crate) assigned: u64,
    pub(crate) unmatched: u64,
    pub(crate) ambiguous: u64,
    pub(crate) unrouted: u64,
    pub(crate) short_reads: u64,
    pub(crate) orphan_r1: u64,
    pub(crate) orphan_r2: u64,
}

impl DemuxCounts {
    pub(crate) fn merge(&mut self, other: &Self) {
        self.assigned += other.assigned;
        self.unmatched += other.unmatched;
        self.ambiguous += other.ambiguous;
        self.unrouted += other.unrouted;
        self.short_reads += other.short_reads;
        self.orphan_r1 += other.orphan_r1;
        self.orphan_r2 += other.orphan_r2;
    }
}

enum SerialWriter {
    Direct(DirectWriterManager),
    Buffered(WriterManager),
}

impl SerialWriter {
    fn write_sample(
        &mut self,
        sample_id: u32,
        mate: OutputMate,
        id: &[u8],
        seq: &[u8],
        qual: &[u8],
    ) -> Result<(), String> {
        match self {
            Self::Direct(writer) => writer.write_sample(sample_id, mate, id, seq, qual),
            Self::Buffered(writer) => writer.write_sample(sample_id, mate, id, seq, qual),
        }
    }

    fn write_unassigned(
        &mut self,
        mate: OutputMate,
        id: &[u8],
        seq: &[u8],
        qual: &[u8],
    ) -> Result<(), String> {
        match self {
            Self::Direct(writer) => writer.write_unassigned(mate, id, seq, qual),
            Self::Buffered(writer) => writer.write_unassigned(mate, id, seq, qual),
        }
    }

    fn finish_with_qc(self) -> Result<WriterCompletion, String> {
        match self {
            Self::Direct(writer) => writer.finish_with_qc(),
            Self::Buffered(writer) => writer.finish_with_qc(),
        }
    }
}

pub fn run_with_threads(args: DemuxArgs, threads: usize) -> Result<(), String> {
    if threads == 0 {
        return Err("Thread count must be greater than 0".into());
    }

    let input = resolve_input_source(&args)?;

    if let InputSource::Fastq(inputs) = input {
        if threads == 1 {
            return run_resolved(args, inputs);
        }
        return crate::parallel::run(args, threads, InputSource::Fastq(inputs));
    }

    crate::parallel::run(args, threads, input)
}

pub fn run(args: DemuxArgs) -> Result<(), String> {
    let inputs = resolve_inputs(&args)?;
    run_resolved(args, inputs)
}

fn run_resolved(args: DemuxArgs, inputs: InputFiles) -> Result<(), String> {
    let paired = inputs.is_paired();

    let layout = if paired {
        ReadLayout::paired(args.r1_structure.as_deref(), args.r2_structure.as_deref())?
    } else {
        ReadLayout::single(args.structure.as_deref())?
    };

    let catalog = BarcodeCatalog::load(&args.barcodes, &layout)?;

    let samples = SampleSheet::load(&args.samples, &layout, &catalog)?;

    let routing = RoutingTree::new(&layout, &catalog, &samples, args.max_mismatches)?;

    let output_dir = args.output.clone();
    let output_layout = OutputLayout::new(samples.samples.len(), paired, args.write_unassigned);
    let expected_streams = output_layout.stream_count()?;
    let output_mode = resolve_output_mode(args.output_mode, expected_streams)?;
    report_output_mode(args.output_mode, output_mode, expected_streams);
    let max_open_files = resolve_max_open_files(expected_streams, args.max_open_files)?;
    let mut writer = match output_mode {
        ResolvedOutputMode::Direct => {
            eprintln!("Output writers: max-open-files={max_open_files}");
            SerialWriter::Direct(DirectWriterManager::new(
                args.output,
                &samples,
                paired,
                args.write_unassigned,
                max_open_files,
                args.compression_level,
            )?)
        }
        ResolvedOutputMode::Buffered => {
            let output_policy = calculate_output_buffer_policy(
                expected_streams,
                available_memory_bytes().unwrap_or(512 * 1024 * 1024),
                args.output_chunk_size.explicit_bytes(),
                args.output_buffer_memory.explicit_bytes(),
            )?;
            eprintln!(
                "Output buffering: streams={} chunk={} bytes budget={} bytes max-open-files={}",
                output_policy.expected_streams,
                output_policy.chunk_size,
                output_policy.memory_budget,
                max_open_files,
            );
            SerialWriter::Buffered(WriterManager::new_with_policy(
                args.output,
                &samples,
                paired,
                args.write_unassigned,
                max_open_files,
                args.compression_level,
                output_policy.chunk_size,
                output_policy.memory_budget,
            )?)
        }
    };
    let run_marker = OutputRunMarker::begin(&output_dir)?;

    let (r1_qc, r2_qc) = MateQcStats::for_layout(&layout);
    let mut r1_stats = args.fastq_stats.then_some(r1_qc);
    let mut r2_stats = if args.fastq_stats { r2_qc } else { None };

    let mut counts = DemuxCounts::default();

    match &inputs {
        InputFiles::Single(reads) => {
            run_single(reads, &routing, &mut writer, r1_stats.as_mut(), &mut counts)?
        }
        InputFiles::Paired { r1, r2 } => run_paired(
            r1,
            r2,
            &routing,
            &mut writer,
            r1_stats.as_mut(),
            r2_stats.as_mut(),
            &mut counts,
        )?,
    }

    let completion = writer.finish_with_qc()?;
    completion.sample_qc.verify(counts.assigned)?;
    verify_serial_unassigned_counts(
        &counts,
        paired,
        args.write_unassigned,
        completion.unassigned_records,
    )?;

    if let Some(r1) = &r1_stats {
        write_fastq_stats(&output_dir, r1, r2_stats.as_ref())?;
    }

    let qc_summary =
        completion
            .sample_qc
            .write_report(&output_dir, &samples, args.low_sample_fraction)?;

    run_marker.complete()?;

    print_summary(&counts);
    eprintln!("  expected samples:  {}", qc_summary.expected);
    eprintln!("  populated samples: {}", qc_summary.populated);
    eprintln!("  missing samples:   {}", qc_summary.missing);
    eprintln!("  low samples:       {}", qc_summary.low);
    eprintln!(
        "  gzip chunks: {} ({:.1} KiB average uncompressed, {} compressed bytes)",
        completion.chunks,
        if completion.chunks == 0 {
            0.0
        } else {
            completion.uncompressed_bytes as f64 / completion.chunks as f64 / 1024.0
        },
        completion.compressed_bytes,
    );

    Ok(())
}

fn verify_serial_unassigned_counts(
    counts: &DemuxCounts,
    paired: bool,
    write_unassigned: bool,
    observed_records: u64,
) -> Result<(), String> {
    if !write_unassigned {
        return (observed_records == 0)
            .then_some(())
            .ok_or_else(|| "Unassigned records were emitted while output was disabled".into());
    }
    let classified = counts
        .unmatched
        .checked_add(counts.ambiguous)
        .and_then(|value| value.checked_add(counts.unrouted))
        .ok_or("Unassigned count overflow")?;
    let expected = if paired {
        classified
            .checked_mul(2)
            .and_then(|value| value.checked_add(counts.orphan_r1))
            .and_then(|value| value.checked_add(counts.orphan_r2))
            .ok_or("Unassigned output count overflow")?
    } else {
        classified
    };
    if expected != observed_records {
        return Err(format!(
            "Unassigned-count reconciliation failed: expected={expected}, output={observed_records}"
        ));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_single(
    path: &Path,
    routing: &RoutingTree,
    writer: &mut SerialWriter,
    mut stats: Option<&mut MateQcStats>,
    counts: &mut DemuxCounts,
) -> Result<(), String> {
    let mut reader = InputReader::try_open(path)?;

    let prefix_len = routing.r1_prefix_len();

    while let Some(record) = reader.next_record() {
        let record =
            record.map_err(|e| format!("FASTQ parse error in '{}': {e}", path.display()))?;

        let seq_cow = record.seq();
        let seq = seq_cow.as_ref();

        let qual = record.qual().ok_or_else(|| {
            format!(
                "Input '{}' is not FASTQ or has no quality scores",
                path.display()
            )
        })?;

        if let Some(stats) = stats.as_mut() {
            stats.update(seq, qual)?;
        }

        let route = if let Some(route) = routing.route_read(seq, None)? {
            route
        } else {
            counts.short_reads += 1;
            RouteResult::Unmatched
        };

        match route {
            RouteResult::Assigned { sample_id } => {
                counts.assigned += 1;

                let trimmed_seq = &seq[prefix_len..];

                let trimmed_qual = &qual[prefix_len..];

                writer.write_sample(
                    sample_id,
                    OutputMate::Single,
                    record.id(),
                    trimmed_seq,
                    trimmed_qual,
                )?;
            }

            RouteResult::Unmatched => {
                counts.unmatched += 1;

                writer.write_unassigned(OutputMate::Single, record.id(), seq, qual)?;
            }

            RouteResult::Ambiguous => {
                counts.ambiguous += 1;

                writer.write_unassigned(OutputMate::Single, record.id(), seq, qual)?;
            }

            RouteResult::Unrouted => {
                counts.unrouted += 1;

                writer.write_unassigned(OutputMate::Single, record.id(), seq, qual)?;
            }
        }
    }

    Ok(())
}

#[derive(Debug)]
struct OwnedFastqRecord {
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

impl OwnedFastqRecord {
    fn new(id: &[u8], seq: &[u8], qual: &[u8]) -> Self {
        Self {
            id: id.to_vec(),
            seq: seq.to_vec(),
            qual: qual.to_vec(),
        }
    }
}

fn read_owned_record(
    reader: &mut InputReader,
    path: &Path,
    stats: Option<&mut MateQcStats>,
) -> Result<Option<OwnedFastqRecord>, String> {
    let Some(record) = reader.next_record() else {
        return Ok(None);
    };

    let record = record.map_err(|e| format!("FASTQ parse error in '{}': {e}", path.display()))?;

    let seq_cow = record.seq();
    let seq = seq_cow.as_ref();

    let qual = record
        .qual()
        .ok_or_else(|| format!("Input '{}' has no FASTQ quality scores", path.display()))?;

    if let Some(stats) = stats {
        stats.update(seq, qual)?;
    }

    Ok(Some(OwnedFastqRecord::new(record.id(), seq, qual)))
}

fn write_orphan_parts(
    writer: &mut SerialWriter,
    mate: OutputMate,
    id: &[u8],
    seq: &[u8],
    qual: &[u8],
    counts: &mut DemuxCounts,
) -> Result<(), String> {
    match mate {
        OutputMate::R1 => counts.orphan_r1 += 1,
        OutputMate::R2 => counts.orphan_r2 += 1,
        OutputMate::Single => {
            return Err("Internal error: paired orphan cannot be single-end".into());
        }
    }

    writer.write_unassigned(mate, id, seq, qual)
}

fn write_owned_orphan(
    writer: &mut SerialWriter,
    mate: OutputMate,
    record: OwnedFastqRecord,
    counts: &mut DemuxCounts,
) -> Result<(), String> {
    write_orphan_parts(writer, mate, &record.id, &record.seq, &record.qual, counts)
}

fn drain_owned_orphans(
    queue: &mut VecDeque<OwnedFastqRecord>,
    count: usize,
    mate: OutputMate,
    writer: &mut SerialWriter,
    counts: &mut DemuxCounts,
) -> Result<(), String> {
    for _ in 0..count {
        let record = queue
            .pop_front()
            .ok_or("Internal error: missing buffered orphan record")?;

        write_owned_orphan(writer, mate, record, counts)?;
    }

    Ok(())
}

fn drain_reader_as_orphans(
    reader: &mut InputReader,
    path: &Path,
    mate: OutputMate,
    mut stats: Option<&mut MateQcStats>,
    writer: &mut SerialWriter,
    counts: &mut DemuxCounts,
) -> Result<(), String> {
    while let Some(record) = reader.next_record() {
        let record =
            record.map_err(|e| format!("FASTQ parse error in '{}': {e}", path.display()))?;

        let seq_cow = record.seq();
        let seq = seq_cow.as_ref();

        let qual = record
            .qual()
            .ok_or_else(|| format!("Input '{}' has no FASTQ quality scores", path.display()))?;

        if let Some(stats) = stats.as_mut() {
            stats.update(seq, qual)?;
        }

        write_orphan_parts(writer, mate, record.id(), seq, qual, counts)?;
    }

    Ok(())
}

fn find_resync_offsets(
    r1_queue: &VecDeque<OwnedFastqRecord>,
    r2_queue: &VecDeque<OwnedFastqRecord>,
) -> Option<(usize, usize)> {
    if r1_queue.is_empty() || r2_queue.is_empty() {
        return None;
    }

    let mut r2_positions: HashMap<&[u8], usize> = HashMap::with_capacity(r2_queue.len());

    for (index, record) in r2_queue.iter().enumerate() {
        r2_positions
            .entry(core_read_id(&record.id))
            .or_insert(index);
    }

    let mut best: Option<(usize, usize)> = None;

    for (r1_index, record) in r1_queue.iter().enumerate() {
        let Some(&r2_index) = r2_positions.get(core_read_id(&record.id)) else {
            continue;
        };

        let candidate = (r1_index, r2_index);

        if best.is_none_or(|current| r1_index + r2_index < current.0 + current.1) {
            best = Some(candidate);
        }
    }

    best
}

#[allow(clippy::too_many_arguments)]
fn process_pair_parts(
    r1_id: &[u8],
    r1_seq: &[u8],
    r1_qual: &[u8],
    r2_id: &[u8],
    r2_seq: &[u8],
    r2_qual: &[u8],
    routing: &RoutingTree,
    writer: &mut SerialWriter,
    r1_prefix: usize,
    r2_prefix: usize,
    counts: &mut DemuxCounts,
) -> Result<(), String> {
    let route = if let Some(route) = routing.route_read(r1_seq, Some(r2_seq))? {
        route
    } else {
        counts.short_reads += 1;
        RouteResult::Unmatched
    };

    match route {
        RouteResult::Assigned { sample_id } => {
            counts.assigned += 1;

            writer.write_sample(
                sample_id,
                OutputMate::R1,
                r1_id,
                &r1_seq[r1_prefix..],
                &r1_qual[r1_prefix..],
            )?;

            writer.write_sample(
                sample_id,
                OutputMate::R2,
                r2_id,
                &r2_seq[r2_prefix..],
                &r2_qual[r2_prefix..],
            )?;
        }

        RouteResult::Unmatched => {
            counts.unmatched += 1;
            writer.write_unassigned(OutputMate::R1, r1_id, r1_seq, r1_qual)?;
            writer.write_unassigned(OutputMate::R2, r2_id, r2_seq, r2_qual)?;
        }

        RouteResult::Ambiguous => {
            counts.ambiguous += 1;
            writer.write_unassigned(OutputMate::R1, r1_id, r1_seq, r1_qual)?;
            writer.write_unassigned(OutputMate::R2, r2_id, r2_seq, r2_qual)?;
        }

        RouteResult::Unrouted => {
            counts.unrouted += 1;
            writer.write_unassigned(OutputMate::R1, r1_id, r1_seq, r1_qual)?;
            writer.write_unassigned(OutputMate::R2, r2_id, r2_seq, r2_qual)?;
        }
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn recover_pairing(
    first_r1: OwnedFastqRecord,
    first_r2: OwnedFastqRecord,
    r1_reader: &mut InputReader,
    r2_reader: &mut InputReader,
    r1_path: &Path,
    r2_path: &Path,
    routing: &RoutingTree,
    writer: &mut SerialWriter,
    r1_prefix: usize,
    r2_prefix: usize,
    mut r1_stats: Option<&mut MateQcStats>,
    mut r2_stats: Option<&mut MateQcStats>,
    counts: &mut DemuxCounts,
) -> Result<(), String> {
    let mut r1_queue = VecDeque::from([first_r1]);
    let mut r2_queue = VecDeque::from([first_r2]);
    let mut r1_eof = false;
    let mut r2_eof = false;

    loop {
        if let Some((r1_orphans, r2_orphans)) = find_resync_offsets(&r1_queue, &r2_queue) {
            drain_owned_orphans(&mut r1_queue, r1_orphans, OutputMate::R1, writer, counts)?;

            drain_owned_orphans(&mut r2_queue, r2_orphans, OutputMate::R2, writer, counts)?;

            let r1_record = r1_queue
                .pop_front()
                .ok_or("Internal error: missing resynchronized R1 record")?;
            let r2_record = r2_queue
                .pop_front()
                .ok_or("Internal error: missing resynchronized R2 record")?;

            process_pair_parts(
                &r1_record.id,
                &r1_record.seq,
                &r1_record.qual,
                &r2_record.id,
                &r2_record.seq,
                &r2_record.qual,
                routing,
                writer,
                r1_prefix,
                r2_prefix,
                counts,
            )?;

            if r1_queue.is_empty() && r2_queue.is_empty() {
                return Ok(());
            }

            continue;
        }

        if r1_eof && r2_eof {
            let r1_orphans = r1_queue.len();
            let r2_orphans = r2_queue.len();

            drain_owned_orphans(&mut r1_queue, r1_orphans, OutputMate::R1, writer, counts)?;

            drain_owned_orphans(&mut r2_queue, r2_orphans, OutputMate::R2, writer, counts)?;

            return Ok(());
        }

        if r1_eof && r1_queue.is_empty() {
            let r2_orphans = r2_queue.len();

            drain_owned_orphans(&mut r2_queue, r2_orphans, OutputMate::R2, writer, counts)?;

            drain_reader_as_orphans(
                r2_reader,
                r2_path,
                OutputMate::R2,
                r2_stats.as_deref_mut(),
                writer,
                counts,
            )?;

            return Ok(());
        }

        if r2_eof && r2_queue.is_empty() {
            let r1_orphans = r1_queue.len();

            drain_owned_orphans(&mut r1_queue, r1_orphans, OutputMate::R1, writer, counts)?;

            drain_reader_as_orphans(
                r1_reader,
                r1_path,
                OutputMate::R1,
                r1_stats.as_deref_mut(),
                writer,
                counts,
            )?;

            return Ok(());
        }

        let mut progressed = false;

        if !r1_eof && r1_queue.len() < PAIRED_RESYNC_WINDOW {
            match read_owned_record(r1_reader, r1_path, r1_stats.as_deref_mut())? {
                Some(record) => {
                    r1_queue.push_back(record);
                    progressed = true;
                }
                None => r1_eof = true,
            }

            if find_resync_offsets(&r1_queue, &r2_queue).is_some() {
                continue;
            }
        }

        if !r2_eof && r2_queue.len() < PAIRED_RESYNC_WINDOW {
            match read_owned_record(r2_reader, r2_path, r2_stats.as_deref_mut())? {
                Some(record) => {
                    r2_queue.push_back(record);
                    progressed = true;
                }
                None => r2_eof = true,
            }

            if find_resync_offsets(&r1_queue, &r2_queue).is_some() {
                continue;
            }
        }

        if !progressed && !(r1_eof && r2_eof) {
            let r1_id = r1_queue
                .front()
                .map_or(b"<empty>".as_slice(), |record| record.id.as_slice());
            let r2_id = r2_queue
                .front()
                .map_or(b"<empty>".as_slice(), |record| record.id.as_slice());

            return Err(format!(
                "Could not resynchronize paired FASTQ records within {PAIRED_RESYNC_WINDOW} records: '{}' vs '{}'",
                String::from_utf8_lossy(r1_id),
                String::from_utf8_lossy(r2_id),
            ));
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn run_paired(
    r1_path: &Path,
    r2_path: &Path,
    routing: &RoutingTree,
    writer: &mut SerialWriter,
    mut r1_stats: Option<&mut MateQcStats>,
    mut r2_stats: Option<&mut MateQcStats>,
    counts: &mut DemuxCounts,
) -> Result<(), String> {
    let mut r1_reader = InputReader::try_open(r1_path)?;
    let mut r2_reader = InputReader::try_open(r2_path)?;

    let (r1_prefix, r2_prefix) = (routing.r1_prefix_len(), routing.r2_prefix_len());

    loop {
        let r1_next = r1_reader.next_record();
        let r2_next = r2_reader.next_record();

        match (r1_next, r2_next) {
            (None, None) => break,

            (Some(r1), Some(r2)) => {
                let r1_record =
                    r1.map_err(|e| format!("FASTQ parse error in '{}': {e}", r1_path.display()))?;

                let r2_record =
                    r2.map_err(|e| format!("FASTQ parse error in '{}': {e}", r2_path.display()))?;

                let r1_seq_cow = r1_record.seq();
                let r2_seq_cow = r2_record.seq();
                let r1_seq = r1_seq_cow.as_ref();
                let r2_seq = r2_seq_cow.as_ref();

                let r1_qual = r1_record.qual().ok_or_else(|| {
                    format!("Input '{}' has no FASTQ quality scores", r1_path.display())
                })?;

                let r2_qual = r2_record.qual().ok_or_else(|| {
                    format!("Input '{}' has no FASTQ quality scores", r2_path.display())
                })?;

                if let Some(stats) = r1_stats.as_mut() {
                    stats.update(r1_seq, r1_qual)?;
                }

                if let Some(stats) = r2_stats.as_mut() {
                    stats.update(r2_seq, r2_qual)?;
                }

                if core_read_id(r1_record.id()) == core_read_id(r2_record.id()) {
                    process_pair_parts(
                        r1_record.id(),
                        r1_seq,
                        r1_qual,
                        r2_record.id(),
                        r2_seq,
                        r2_qual,
                        routing,
                        writer,
                        r1_prefix,
                        r2_prefix,
                        counts,
                    )?;

                    continue;
                }

                let owned_r1 = OwnedFastqRecord::new(r1_record.id(), r1_seq, r1_qual);
                let owned_r2 = OwnedFastqRecord::new(r2_record.id(), r2_seq, r2_qual);

                recover_pairing(
                    owned_r1,
                    owned_r2,
                    &mut r1_reader,
                    &mut r2_reader,
                    r1_path,
                    r2_path,
                    routing,
                    writer,
                    r1_prefix,
                    r2_prefix,
                    r1_stats.as_deref_mut(),
                    r2_stats.as_deref_mut(),
                    counts,
                )?;
            }

            (Some(r1), None) => {
                let r1_record =
                    r1.map_err(|e| format!("FASTQ parse error in '{}': {e}", r1_path.display()))?;

                let seq_cow = r1_record.seq();
                let seq = seq_cow.as_ref();

                let qual = r1_record.qual().ok_or_else(|| {
                    format!("Input '{}' has no FASTQ quality scores", r1_path.display())
                })?;

                if let Some(stats) = r1_stats.as_mut() {
                    stats.update(seq, qual)?;
                }

                write_orphan_parts(writer, OutputMate::R1, r1_record.id(), seq, qual, counts)?;

                drain_reader_as_orphans(
                    &mut r1_reader,
                    r1_path,
                    OutputMate::R1,
                    r1_stats.as_deref_mut(),
                    writer,
                    counts,
                )?;

                break;
            }

            (None, Some(r2)) => {
                let r2_record =
                    r2.map_err(|e| format!("FASTQ parse error in '{}': {e}", r2_path.display()))?;

                let seq_cow = r2_record.seq();
                let seq = seq_cow.as_ref();

                let qual = r2_record.qual().ok_or_else(|| {
                    format!("Input '{}' has no FASTQ quality scores", r2_path.display())
                })?;

                if let Some(stats) = r2_stats.as_mut() {
                    stats.update(seq, qual)?;
                }

                write_orphan_parts(writer, OutputMate::R2, r2_record.id(), seq, qual, counts)?;

                drain_reader_as_orphans(
                    &mut r2_reader,
                    r2_path,
                    OutputMate::R2,
                    r2_stats.as_deref_mut(),
                    writer,
                    counts,
                )?;

                break;
            }
        }
    }

    Ok(())
}

pub(crate) fn core_read_id(id: &[u8]) -> &[u8] {
    let token = id.split(u8::is_ascii_whitespace).next().unwrap_or(id);

    if token.ends_with(b"/1") || token.ends_with(b"/2") {
        &token[..token.len() - 2]
    } else {
        token
    }
}

pub(crate) fn write_fastq_stats(
    output_dir: &Path,
    r1: &MateQcStats,
    r2: Option<&MateQcStats>,
) -> Result<(), String> {
    let path = output_dir.join("fastq_stats.tsv");

    let mut file =
        File::create(&path).map_err(|e| format!("Could not create '{}': {e}", path.display()))?;

    writeln!(
        file,
        "Mate\tReads\tBases\tMinLength\tMaxLength\tMeanLength\tGCPercent\tNPercent\tMeanQuality\tQ20Percent\tQ30Percent"
    )
    .map_err(|e| {
        format!(
            "Could not write FASTQ stats: {e}"
        )
    })?;

    write_stats_row(&mut file, "R1", r1.biological())?;

    if let Some(r2) = r2 {
        write_stats_row(&mut file, "R2", r2.biological())?;
    }

    file.flush()
        .map_err(|e| format!("Could not flush FASTQ stats: {e}"))?;

    write_barcode_stats(output_dir, r1, r2)
}

fn write_barcode_stats(
    output_dir: &Path,
    r1: &MateQcStats,
    r2: Option<&MateQcStats>,
) -> Result<(), String> {
    let path = output_dir.join("barcode_stats.tsv");
    let mut file =
        File::create(&path).map_err(|e| format!("Could not create '{}': {e}", path.display()))?;

    writeln!(
        file,
        "Mate\tBarcodeSymbol\tMatePiece\tStartCycle\tEndCycle\tOrientation\tReads\tBases\tMinLength\tMaxLength\tMeanLength\tGCPercent\tNPercent\tMeanQuality\tQ20Percent\tQ30Percent"
    )
    .map_err(|e| format!("Could not write barcode stats: {e}"))?;

    write_barcode_stats_rows(&mut file, "R1", r1)?;
    if let Some(r2) = r2 {
        write_barcode_stats_rows(&mut file, "R2", r2)?;
    }

    file.flush()
        .map_err(|e| format!("Could not flush barcode stats: {e}"))
}

fn write_barcode_stats_rows(
    file: &mut File,
    mate: &str,
    mate_stats: &MateQcStats,
) -> Result<(), String> {
    for segment in mate_stats.barcode_segments() {
        let orientation = match segment.orientation {
            Orientation::Forward => "forward",
            Orientation::ReverseComplement => "reverse-complement",
        };
        let stats = &segment.stats;
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.4}\t{:.4}\t{:.2}\t{:.4}\t{:.4}",
            mate,
            char::from(segment.symbol),
            segment.piece,
            segment.start + 1,
            segment.end,
            orientation,
            stats.reads,
            stats.bases,
            stats.min_length.unwrap_or(0),
            stats.max_length,
            stats.mean_length(),
            stats.gc_percent(),
            stats.n_percent(),
            stats.mean_quality(),
            stats.q20_percent(),
            stats.q30_percent(),
        )
        .map_err(|e| format!("Could not write barcode stats: {e}"))?;
    }
    Ok(())
}

fn write_stats_row(file: &mut File, mate: &str, stats: &FastqStats) -> Result<(), String> {
    writeln!(
        file,
        "{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.4}\t{:.4}\t{:.2}\t{:.4}\t{:.4}",
        mate,
        stats.reads,
        stats.bases,
        stats.min_length.unwrap_or(0),
        stats.max_length,
        stats.mean_length(),
        stats.gc_percent(),
        stats.n_percent(),
        stats.mean_quality(),
        stats.q20_percent(),
        stats.q30_percent(),
    )
    .map_err(|e| format!("Could not write FASTQ stats: {e}"))
}

pub(crate) fn print_summary(counts: &DemuxCounts) {
    let total = counts.assigned
        + counts.unmatched
        + counts.ambiguous
        + counts.unrouted
        + counts.orphan_r1
        + counts.orphan_r2;

    eprintln!("Demultiplexing complete");
    eprintln!("  total:      {total}");
    eprintln!("  assigned:   {}", counts.assigned);
    eprintln!("  unmatched:  {}", counts.unmatched);
    eprintln!("  ambiguous:  {}", counts.ambiguous);
    eprintln!("  unrouted:   {}", counts.unrouted);

    if counts.short_reads > 0 {
        eprintln!("  short reads: {}", counts.short_reads);
    }

    if counts.orphan_r1 > 0 {
        eprintln!("  orphan R1:   {}", counts.orphan_r1);
    }

    if counts.orphan_r2 > 0 {
        eprintln!("  orphan R2:   {}", counts.orphan_r2);
    }
}
