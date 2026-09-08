use std::collections::{BTreeMap, HashMap, VecDeque};
use std::io::Write;
use std::path::Path;
use std::sync::{Arc, Condvar, Mutex};
use std::thread;

use crossbeam_channel::{Receiver, Sender, bounded};
use flate2::Compression;
use flate2::write::GzEncoder;

use crate::barcodes::BarcodeCatalog;
use crate::cli::DemuxArgs;
use crate::demux::{DemuxCounts, core_read_id, print_summary, write_fastq_stats};
use crate::input::InputFiles;
use crate::parallel_input::ParallelInput;
use crate::routing::{RouteResult, RoutingTree};
use crate::samples::SampleSheet;
use crate::stats::FastqStats;
use crate::structure::ReadLayout;
use crate::thread_plan::{AdaptiveAllocation, AllocationSnapshot, ThreadPlan};
use crate::writer::{CompressedWriterManager, OutputMate};

const MAX_OPEN_WRITERS: usize = 64;
const PAIRED_RESYNC_WINDOW: usize = 1024;
const PAIRED_READER_QUEUE: usize = 256;
const BATCH_SIZE: usize = 1024;
const QUEUE_DEPTH_PER_WORKER: usize = 2;
const MAX_QUEUED_BATCHES: usize = 32;

struct SharedState {
    routing: RoutingTree,
    r1_prefix: usize,
    r2_prefix: usize,
    write_unassigned: bool,
    compression_level: u32,
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

#[derive(Debug)]
enum WorkItem {
    Single(OwnedFastqRecord),
    Pair {
        r1: OwnedFastqRecord,
        r2: OwnedFastqRecord,
    },
    Orphan {
        mate: OutputMate,
        record: OwnedFastqRecord,
    },
}

#[derive(Debug)]
struct WorkBatch {
    id: u64,
    items: Vec<WorkItem>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum OutputTarget {
    Sample { sample_id: u32 },
    Unassigned,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct OutputKey {
    target: OutputTarget,
    mate: OutputMate,
}

#[derive(Debug)]
struct CompressedOutput {
    target: OutputTarget,
    mate: OutputMate,
    member: Vec<u8>,
}

#[derive(Debug)]
struct ProcessedBatch {
    id: u64,
    outputs: Vec<CompressedOutput>,
    counts: DemuxCounts,
}

#[derive(Default)]
struct BatchBuffers {
    buffers: HashMap<OutputKey, Vec<u8>>,
}

impl BatchBuffers {
    fn push_record(
        &mut self,
        target: OutputTarget,
        mate: OutputMate,
        record: &OwnedFastqRecord,
        trim_start: usize,
    ) -> Result<(), String> {
        if trim_start > record.seq.len() || trim_start > record.qual.len() {
            return Err("Internal error: trim prefix exceeds FASTQ record length".into());
        }

        let seq = &record.seq[trim_start..];
        let qual = &record.qual[trim_start..];

        if seq.len() != qual.len() {
            return Err("Sequence and quality lengths do not match".into());
        }

        let key = OutputKey { target, mate };
        let buffer = self.buffers.entry(key).or_default();
        let id = record.id.strip_prefix(b"@").unwrap_or(record.id.as_slice());

        buffer.extend_from_slice(b"@");
        buffer.extend_from_slice(id);
        buffer.extend_from_slice(b"\n");
        buffer.extend_from_slice(seq);
        buffer.extend_from_slice(b"\n+\n");
        buffer.extend_from_slice(qual);
        buffer.extend_from_slice(b"\n");

        Ok(())
    }

    fn compress(self, compression_level: u32) -> Result<Vec<CompressedOutput>, String> {
        let mut outputs = Vec::with_capacity(self.buffers.len());

        for (key, fastq) in self.buffers {
            let mut encoder = GzEncoder::new(Vec::new(), Compression::new(compression_level));
            encoder
                .write_all(&fastq)
                .map_err(|e| format!("Could not compress FASTQ batch: {e}"))?;
            let member = encoder
                .finish()
                .map_err(|e| format!("Could not finish FASTQ gzip member: {e}"))?;

            outputs.push(CompressedOutput {
                target: key.target,
                mate: key.mate,
                member,
            });
        }

        Ok(outputs)
    }
}

struct BatchSender {
    sender: Sender<WorkBatch>,
    next_batch_id: u64,
    items: Vec<WorkItem>,
    allocation: Option<AllocationController>,
}

impl BatchSender {
    fn new(sender: Sender<WorkBatch>) -> Self {
        Self {
            sender,
            next_batch_id: 0,
            items: Vec::with_capacity(BATCH_SIZE),
            allocation: None,
        }
    }

    fn with_allocation_controller(
        sender: Sender<WorkBatch>,
        allocation: AllocationController,
    ) -> Self {
        Self {
            allocation: Some(allocation),
            ..Self::new(sender)
        }
    }

    fn push(&mut self, item: WorkItem) -> Result<(), String> {
        self.items.push(item);

        if self.items.len() >= BATCH_SIZE {
            self.flush()?;
        }

        Ok(())
    }

    fn finish(mut self) -> Result<u64, String> {
        self.flush()?;
        Ok(self.next_batch_id)
    }

    fn flush(&mut self) -> Result<(), String> {
        if self.items.is_empty() {
            return Ok(());
        }

        let items = std::mem::replace(&mut self.items, Vec::with_capacity(BATCH_SIZE));
        let batch = WorkBatch {
            id: self.next_batch_id,
            items,
        };

        self.sender
            .send(batch)
            .map_err(|_| "Parallel worker queue disconnected".to_string())?;

        if let Some(allocation) = &self.allocation {
            allocation.observe(self.sender.len())?;
        }

        self.next_batch_id = self
            .next_batch_id
            .checked_add(1)
            .ok_or("Batch ID overflow")?;

        Ok(())
    }
}

#[derive(Clone)]
struct WorkerGate {
    state: Arc<(Mutex<usize>, Condvar)>,
}

impl WorkerGate {
    fn new(initial_workers: usize) -> Self {
        Self {
            state: Arc::new((Mutex::new(initial_workers), Condvar::new())),
        }
    }

    fn wait_until_active(&self, worker_index: usize) -> Result<(), String> {
        let (limit, changed) = &*self.state;
        let mut limit = limit
            .lock()
            .map_err(|_| "Demultiplexing worker gate was poisoned".to_string())?;

        while worker_index >= *limit {
            limit = changed
                .wait(limit)
                .map_err(|_| "Demultiplexing worker gate was poisoned".to_string())?;
        }

        Ok(())
    }

    fn set_limit(&self, workers: usize) -> Result<(), String> {
        let (limit, changed) = &*self.state;
        let mut limit = limit
            .lock()
            .map_err(|_| "Demultiplexing worker gate was poisoned".to_string())?;
        *limit = workers;
        changed.notify_all();
        Ok(())
    }
}

#[derive(Clone)]
struct AllocationController {
    allocation: Arc<Mutex<AdaptiveAllocation>>,
    parallel_input: ParallelInput,
    worker_gate: WorkerGate,
    queue_capacity: usize,
}

impl AllocationController {
    fn new(
        allocation: AdaptiveAllocation,
        parallel_input: ParallelInput,
        worker_gate: WorkerGate,
        queue_capacity: usize,
    ) -> Self {
        Self {
            allocation: Arc::new(Mutex::new(allocation)),
            parallel_input,
            worker_gate,
            queue_capacity,
        }
    }

    fn observe(&self, queued_batches: usize) -> Result<(), String> {
        let mut allocation = self
            .allocation
            .lock()
            .map_err(|_| "Adaptive thread allocation was poisoned".to_string())?;
        let previous = allocation.current();
        let Some(next) = allocation.observe(queued_batches, self.queue_capacity) else {
            return Ok(());
        };
        drop(allocation);

        self.apply(previous, next)
    }

    fn apply(&self, previous: AllocationSnapshot, next: AllocationSnapshot) -> Result<(), String> {
        if next.input_threads < previous.input_threads {
            self.parallel_input.set_worker_limit(next.input_threads)?;
            self.worker_gate.set_limit(next.worker_threads)?;
        } else {
            self.worker_gate.set_limit(next.worker_threads)?;
            self.parallel_input.set_worker_limit(next.input_threads)?;
        }

        Ok(())
    }
}

pub(crate) fn run(args: DemuxArgs, threads: usize, inputs: InputFiles) -> Result<(), String> {
    if threads == 0 {
        return Err("Parallel pipeline requires at least 1 worker thread".into());
    }

    let paired = inputs.is_paired();
    let has_gzip = inputs.has_gzip()?;
    let thread_plan = ThreadPlan::new(threads, paired, has_gzip)?;
    let parallel_input = ParallelInput::new(
        thread_plan.max_input_threads,
        thread_plan.initial_input_threads,
        thread_plan.parallel_gzip,
    )?;
    let worker_threads = thread_plan.initial_worker_threads;
    let worker_headroom = thread_plan.worker_headroom;

    if threads >= 2 {
        eprintln!(
            "Thread allocation: budget={} fastq-parsing={} input-decompression={} demux/compression={} overcommit={}",
            thread_plan.requested_threads,
            thread_plan.parser_threads,
            thread_plan.initial_input_threads,
            thread_plan.initial_worker_threads,
            thread_plan.budget_overcommit,
        );
    }

    let layout = if paired {
        ReadLayout::paired(args.r1_structure.as_deref(), args.r2_structure.as_deref())?
    } else {
        ReadLayout::single(args.structure.as_deref())?
    };

    let catalog = BarcodeCatalog::load(&args.barcodes, &layout)?;
    let samples = SampleSheet::load(&args.samples, &layout, &catalog)?;
    let routing = RoutingTree::new(&layout, &catalog, &samples, args.max_mismatches)?;
    let (r1_prefix, r2_prefix) = (routing.r1_prefix_len(), routing.r2_prefix_len());

    let state = SharedState {
        routing,
        r1_prefix,
        r2_prefix,
        write_unassigned: args.write_unassigned,
        compression_level: args.compression_level,
    };

    let output_dir = args.output.clone();
    let writer = CompressedWriterManager::new(
        args.output.clone(),
        &samples,
        paired,
        args.write_unassigned,
        MAX_OPEN_WRITERS,
        args.compression_level,
    )?;

    let mut r1_stats = if args.fastq_stats {
        Some(FastqStats::default())
    } else {
        None
    };

    let mut r2_stats = if paired && args.fastq_stats {
        Some(FastqStats::default())
    } else {
        None
    };

    let queue_capacity = worker_threads
        .checked_mul(QUEUE_DEPTH_PER_WORKER)
        .ok_or("Worker queue capacity overflow")?
        .clamp(1, MAX_QUEUED_BATCHES);

    let (work_tx, work_rx) = bounded::<WorkBatch>(queue_capacity);
    let (result_tx, result_rx) = bounded::<Result<ProcessedBatch, String>>(queue_capacity);
    let worker_gate = WorkerGate::new(worker_threads);
    let allocation_controller = thread_plan.adaptive_allocation().map(|allocation| {
        AllocationController::new(
            allocation,
            parallel_input.clone(),
            worker_gate.clone(),
            queue_capacity,
        )
    });

    let pipeline_result = thread::scope(|scope| {
        let writer_handle = scope.spawn(move || writer_loop(writer, result_rx));

        let mut worker_handles = Vec::with_capacity(worker_headroom);

        for worker_index in 0..worker_headroom {
            let worker_rx = work_rx.clone();
            let worker_tx = result_tx.clone();
            let worker_state = &state;
            let gate = worker_gate.clone();

            worker_handles.push(scope.spawn(move || {
                worker_loop(worker_index, worker_rx, worker_tx, worker_state, gate)
            }));
        }

        drop(work_rx);
        drop(result_tx);

        let mut batcher = match allocation_controller {
            Some(controller) => BatchSender::with_allocation_controller(work_tx, controller),
            None => BatchSender::new(work_tx),
        };

        let producer_result = match &inputs {
            InputFiles::Single(path) => {
                produce_single(path, &mut batcher, r1_stats.as_mut(), &parallel_input)
            }
            InputFiles::Paired { r1, r2 } => produce_paired(
                r1,
                r2,
                &mut batcher,
                r1_stats.as_mut(),
                r2_stats.as_mut(),
                &parallel_input,
            ),
        };

        let batches_sent = match producer_result {
            Ok(()) => batcher.finish(),
            Err(error) => {
                drop(batcher);
                Err(error)
            }
        };

        let gate_result = worker_gate.set_limit(worker_headroom);

        let mut worker_error = None;

        for handle in worker_handles {
            match handle.join() {
                Ok(Ok(())) => {}
                Ok(Err(error)) => {
                    if worker_error.is_none() {
                        worker_error = Some(error);
                    }
                }
                Err(_) => {
                    if worker_error.is_none() {
                        worker_error = Some("Parallel worker thread panicked".into());
                    }
                }
            }
        }

        let writer_result = match writer_handle.join() {
            Ok(result) => result,
            Err(_) => Err("Writer thread panicked".into()),
        };

        let (counts, written_batches) = writer_result?;

        if let Some(error) = worker_error {
            return Err(error);
        }

        gate_result?;

        let expected_batches = batches_sent?;

        if written_batches != expected_batches {
            return Err(format!(
                "Parallel pipeline lost work: sent {expected_batches} batches but wrote {written_batches}"
            ));
        }

        Ok(counts)
    });

    let counts = pipeline_result?;

    if let Some(r1) = &r1_stats {
        write_fastq_stats(&output_dir, r1, r2_stats.as_ref())?;
    }

    print_summary(&counts);

    Ok(())
}

fn worker_loop(
    worker_index: usize,
    receiver: Receiver<WorkBatch>,
    sender: Sender<Result<ProcessedBatch, String>>,
    state: &SharedState,
    gate: WorkerGate,
) -> Result<(), String> {
    loop {
        gate.wait_until_active(worker_index)?;
        let Ok(batch) = receiver.recv() else {
            return Ok(());
        };

        match process_batch(batch, state) {
            Ok(processed) => {
                if sender.send(Ok(processed)).is_err() {
                    return Ok(());
                }
            }
            Err(error) => {
                let _ = sender.send(Err(error.clone()));
                return Err(error);
            }
        }
    }
}

fn process_batch(batch: WorkBatch, state: &SharedState) -> Result<ProcessedBatch, String> {
    let mut buffers = BatchBuffers::default();
    let mut counts = DemuxCounts::default();

    for item in batch.items {
        match item {
            WorkItem::Single(record) => {
                process_single_record(record, state, &mut buffers, &mut counts)?;
            }
            WorkItem::Pair { r1, r2 } => {
                process_pair(r1, r2, state, &mut buffers, &mut counts)?;
            }
            WorkItem::Orphan { mate, record } => {
                match mate {
                    OutputMate::R1 => counts.orphan_r1 += 1,
                    OutputMate::R2 => counts.orphan_r2 += 1,
                    OutputMate::Single => {
                        return Err("Internal error: paired orphan cannot be single-end".into());
                    }
                }

                if state.write_unassigned {
                    buffers.push_record(OutputTarget::Unassigned, mate, &record, 0)?;
                }
            }
        }
    }

    let outputs = buffers.compress(state.compression_level)?;

    Ok(ProcessedBatch {
        id: batch.id,
        outputs,
        counts,
    })
}

fn process_single_record(
    record: OwnedFastqRecord,
    state: &SharedState,
    buffers: &mut BatchBuffers,
    counts: &mut DemuxCounts,
) -> Result<(), String> {
    let route = if let Some(route) = state.routing.route_read(&record.seq, None)? {
        route
    } else {
        counts.short_reads += 1;
        RouteResult::Unmatched
    };

    let (target, trim_start) = route_to_output(route, state.r1_prefix, counts);

    if should_emit(target, state.write_unassigned) {
        buffers.push_record(target, OutputMate::Single, &record, trim_start)?;
    }

    Ok(())
}

fn process_pair(
    r1: OwnedFastqRecord,
    r2: OwnedFastqRecord,
    state: &SharedState,
    buffers: &mut BatchBuffers,
    counts: &mut DemuxCounts,
) -> Result<(), String> {
    let route = if let Some(route) = state.routing.route_read(&r1.seq, Some(&r2.seq))? {
        route
    } else {
        counts.short_reads += 1;
        RouteResult::Unmatched
    };

    let assigned = matches!(route, RouteResult::Assigned { .. });
    let target = route_to_target(route, counts);

    if should_emit(target, state.write_unassigned) {
        buffers.push_record(
            target,
            OutputMate::R1,
            &r1,
            if assigned { state.r1_prefix } else { 0 },
        )?;
        buffers.push_record(
            target,
            OutputMate::R2,
            &r2,
            if assigned { state.r2_prefix } else { 0 },
        )?;
    }

    Ok(())
}

fn should_emit(target: OutputTarget, write_unassigned: bool) -> bool {
    matches!(target, OutputTarget::Sample { .. }) || write_unassigned
}

fn route_to_output(
    route: RouteResult,
    assigned_prefix: usize,
    counts: &mut DemuxCounts,
) -> (OutputTarget, usize) {
    let assigned = matches!(route, RouteResult::Assigned { .. });
    let target = route_to_target(route, counts);
    (target, if assigned { assigned_prefix } else { 0 })
}

fn route_to_target(route: RouteResult, counts: &mut DemuxCounts) -> OutputTarget {
    match route {
        RouteResult::Assigned { sample_id } => {
            counts.assigned += 1;
            OutputTarget::Sample { sample_id }
        }
        RouteResult::Unmatched => {
            counts.unmatched += 1;
            OutputTarget::Unassigned
        }
        RouteResult::Ambiguous => {
            counts.ambiguous += 1;
            OutputTarget::Unassigned
        }
        RouteResult::Unrouted => {
            counts.unrouted += 1;
            OutputTarget::Unassigned
        }
    }
}

fn writer_loop(
    mut writer: CompressedWriterManager,
    receiver: Receiver<Result<ProcessedBatch, String>>,
) -> Result<(DemuxCounts, u64), String> {
    let mut pending = BTreeMap::<u64, ProcessedBatch>::new();
    let mut next_batch_id = 0u64;
    let mut counts = DemuxCounts::default();

    while let Ok(message) = receiver.recv() {
        let batch = message?;

        if batch.id < next_batch_id || pending.insert(batch.id, batch).is_some() {
            return Err("Duplicate or stale processed batch ID".into());
        }

        while let Some(batch) = pending.remove(&next_batch_id) {
            write_processed_batch(&mut writer, batch, &mut counts)?;
            next_batch_id = next_batch_id.checked_add(1).ok_or("Batch ID overflow")?;
        }
    }

    if !pending.is_empty() {
        return Err(format!(
            "Parallel output ended with a missing batch before batch {next_batch_id}"
        ));
    }

    writer.finish()?;

    Ok((counts, next_batch_id))
}

fn write_processed_batch(
    writer: &mut CompressedWriterManager,
    batch: ProcessedBatch,
    counts: &mut DemuxCounts,
) -> Result<(), String> {
    for output in batch.outputs {
        match output.target {
            OutputTarget::Sample { sample_id } => {
                writer.append_sample_member(sample_id, output.mate, &output.member)?;
            }
            OutputTarget::Unassigned => {
                writer.append_unassigned_member(output.mate, &output.member)?;
            }
        }
    }

    counts.merge(&batch.counts);
    Ok(())
}

fn produce_single(
    path: &Path,
    batcher: &mut BatchSender,
    mut stats: Option<&mut FastqStats>,
    parallel_input: &ParallelInput,
) -> Result<(), String> {
    let mut reader = parallel_input.open(path)?;

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

        batcher.push(WorkItem::Single(OwnedFastqRecord::new(
            record.id(),
            seq,
            qual,
        )))?;
    }

    Ok(())
}

fn stream_owned_records(
    path: &Path,
    sender: Sender<Result<OwnedFastqRecord, String>>,
    parallel_input: ParallelInput,
) {
    let mut reader = match parallel_input.open(path) {
        Ok(reader) => reader,
        Err(error) => {
            let _ = sender.send(Err(error));
            return;
        }
    };

    while let Some(record) = reader.next_record() {
        let record = match record {
            Ok(record) => record,
            Err(error) => {
                let message = format!("FASTQ parse error in '{}': {error}", path.display());
                let _ = sender.send(Err(message));
                return;
            }
        };

        let seq_cow = record.seq();
        let seq = seq_cow.as_ref();
        let Some(qual) = record.qual() else {
            let message = format!("Input '{}' has no FASTQ quality scores", path.display());
            let _ = sender.send(Err(message));
            return;
        };

        let owned = OwnedFastqRecord::new(record.id(), seq, qual);

        if sender.send(Ok(owned)).is_err() {
            return;
        }
    }
}

fn recv_owned_record(
    receiver: &Receiver<Result<OwnedFastqRecord, String>>,
    stats: Option<&mut FastqStats>,
) -> Result<Option<OwnedFastqRecord>, String> {
    match receiver.recv() {
        Ok(Ok(record)) => {
            if let Some(stats) = stats {
                stats.update(&record.seq, &record.qual)?;
            }
            Ok(Some(record))
        }
        Ok(Err(error)) => Err(error),
        Err(_) => Ok(None),
    }
}

fn push_orphan(
    batcher: &mut BatchSender,
    mate: OutputMate,
    record: OwnedFastqRecord,
) -> Result<(), String> {
    batcher.push(WorkItem::Orphan { mate, record })
}

fn drain_owned_orphans(
    queue: &mut VecDeque<OwnedFastqRecord>,
    count: usize,
    mate: OutputMate,
    batcher: &mut BatchSender,
) -> Result<(), String> {
    for _ in 0..count {
        let record = queue
            .pop_front()
            .ok_or("Internal error: missing buffered orphan record")?;
        push_orphan(batcher, mate, record)?;
    }
    Ok(())
}

fn drain_receiver_as_orphans(
    receiver: &Receiver<Result<OwnedFastqRecord, String>>,
    mate: OutputMate,
    mut stats: Option<&mut FastqStats>,
    batcher: &mut BatchSender,
) -> Result<(), String> {
    while let Some(record) = recv_owned_record(receiver, stats.as_deref_mut())? {
        push_orphan(batcher, mate, record)?;
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
fn recover_pairing(
    first_r1: OwnedFastqRecord,
    first_r2: OwnedFastqRecord,
    r1_receiver: &Receiver<Result<OwnedFastqRecord, String>>,
    r2_receiver: &Receiver<Result<OwnedFastqRecord, String>>,
    batcher: &mut BatchSender,
    mut r1_stats: Option<&mut FastqStats>,
    mut r2_stats: Option<&mut FastqStats>,
) -> Result<(), String> {
    let mut r1_queue = VecDeque::from([first_r1]);
    let mut r2_queue = VecDeque::from([first_r2]);
    let mut r1_eof = false;
    let mut r2_eof = false;

    loop {
        if let Some((r1_orphans, r2_orphans)) = find_resync_offsets(&r1_queue, &r2_queue) {
            drain_owned_orphans(&mut r1_queue, r1_orphans, OutputMate::R1, batcher)?;
            drain_owned_orphans(&mut r2_queue, r2_orphans, OutputMate::R2, batcher)?;

            let r1 = r1_queue
                .pop_front()
                .ok_or("Internal error: missing resynchronized R1 record")?;
            let r2 = r2_queue
                .pop_front()
                .ok_or("Internal error: missing resynchronized R2 record")?;

            batcher.push(WorkItem::Pair { r1, r2 })?;

            if r1_queue.is_empty() && r2_queue.is_empty() {
                return Ok(());
            }
            continue;
        }

        if r1_eof && r2_eof {
            let r1_orphans = r1_queue.len();
            let r2_orphans = r2_queue.len();
            drain_owned_orphans(&mut r1_queue, r1_orphans, OutputMate::R1, batcher)?;
            drain_owned_orphans(&mut r2_queue, r2_orphans, OutputMate::R2, batcher)?;
            return Ok(());
        }

        if r1_eof && r1_queue.is_empty() {
            let r2_orphans = r2_queue.len();
            drain_owned_orphans(&mut r2_queue, r2_orphans, OutputMate::R2, batcher)?;
            drain_receiver_as_orphans(
                r2_receiver,
                OutputMate::R2,
                r2_stats.as_deref_mut(),
                batcher,
            )?;
            return Ok(());
        }

        if r2_eof && r2_queue.is_empty() {
            let r1_orphans = r1_queue.len();
            drain_owned_orphans(&mut r1_queue, r1_orphans, OutputMate::R1, batcher)?;
            drain_receiver_as_orphans(
                r1_receiver,
                OutputMate::R1,
                r1_stats.as_deref_mut(),
                batcher,
            )?;
            return Ok(());
        }

        let mut progressed = false;

        if !r1_eof && r1_queue.len() < PAIRED_RESYNC_WINDOW {
            match recv_owned_record(r1_receiver, r1_stats.as_deref_mut())? {
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
            match recv_owned_record(r2_receiver, r2_stats.as_deref_mut())? {
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

fn produce_paired_from_receivers(
    r1_receiver: &Receiver<Result<OwnedFastqRecord, String>>,
    r2_receiver: &Receiver<Result<OwnedFastqRecord, String>>,
    batcher: &mut BatchSender,
    mut r1_stats: Option<&mut FastqStats>,
    mut r2_stats: Option<&mut FastqStats>,
) -> Result<(), String> {
    loop {
        let r1_next = recv_owned_record(r1_receiver, r1_stats.as_deref_mut())?;
        let r2_next = recv_owned_record(r2_receiver, r2_stats.as_deref_mut())?;

        match (r1_next, r2_next) {
            (None, None) => break,
            (Some(r1), Some(r2)) => {
                if core_read_id(&r1.id) == core_read_id(&r2.id) {
                    batcher.push(WorkItem::Pair { r1, r2 })?;
                    continue;
                }

                recover_pairing(
                    r1,
                    r2,
                    r1_receiver,
                    r2_receiver,
                    batcher,
                    r1_stats.as_deref_mut(),
                    r2_stats.as_deref_mut(),
                )?;
            }
            (Some(r1), None) => {
                push_orphan(batcher, OutputMate::R1, r1)?;
                drain_receiver_as_orphans(
                    r1_receiver,
                    OutputMate::R1,
                    r1_stats.as_deref_mut(),
                    batcher,
                )?;
                break;
            }
            (None, Some(r2)) => {
                push_orphan(batcher, OutputMate::R2, r2)?;
                drain_receiver_as_orphans(
                    r2_receiver,
                    OutputMate::R2,
                    r2_stats.as_deref_mut(),
                    batcher,
                )?;
                break;
            }
        }
    }

    Ok(())
}

fn produce_paired(
    r1_path: &Path,
    r2_path: &Path,
    batcher: &mut BatchSender,
    r1_stats: Option<&mut FastqStats>,
    r2_stats: Option<&mut FastqStats>,
    parallel_input: &ParallelInput,
) -> Result<(), String> {
    thread::scope(|scope| {
        let (r1_tx, r1_rx) = bounded::<Result<OwnedFastqRecord, String>>(PAIRED_READER_QUEUE);
        let (r2_tx, r2_rx) = bounded::<Result<OwnedFastqRecord, String>>(PAIRED_READER_QUEUE);

        let r1_input = parallel_input.clone();
        let r2_input = parallel_input.clone();

        let r1_handle = scope.spawn(move || stream_owned_records(r1_path, r1_tx, r1_input));
        let r2_handle = scope.spawn(move || stream_owned_records(r2_path, r2_tx, r2_input));

        let producer_result =
            produce_paired_from_receivers(&r1_rx, &r2_rx, batcher, r1_stats, r2_stats);

        drop(r1_rx);
        drop(r2_rx);

        let r1_result = match r1_handle.join() {
            Ok(()) => Ok(()),
            Err(_) => Err("R1 reader thread panicked".to_string()),
        };
        let r2_result = match r2_handle.join() {
            Ok(()) => Ok(()),
            Err(_) => Err("R2 reader thread panicked".to_string()),
        };

        producer_result?;
        r1_result?;
        r2_result?;

        Ok(())
    })
}
