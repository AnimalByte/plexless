use std::any::Any;
use std::collections::{BTreeMap, HashMap, VecDeque};
use std::io::Write;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::path::Path;
#[cfg(test)]
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Condvar, Mutex, Weak};
use std::thread;

use crossbeam_channel::{Receiver, Sender, TrySendError, bounded};
use flate2::Compression;
use flate2::write::GzEncoder;

use crate::barcodes::BarcodeCatalog;
use crate::cli::DemuxArgs;
use crate::demux::{DemuxCounts, core_read_id, print_summary, write_fastq_stats};
use crate::input::InputFiles;
use crate::output::{
    CompressionJob, OutputAccumulator, OutputKey, OutputLayout, OutputMate, OutputTarget,
    ResolvedOutputMode, available_memory_bytes, calculate_output_buffer_policy, report_output_mode,
    resolve_output_mode,
};
use crate::parallel_input::ParallelInput;
use crate::qc::{SampleQc, SampleQcSummary};
use crate::routing::{RouteResult, RoutingTree};
use crate::samples::SampleSheet;
use crate::stats::MateQcStats;
use crate::structure::ReadLayout;
use crate::thread_plan::{AdaptiveAllocation, AllocationSnapshot, ThreadPlan};
use crate::writer::{CompressedWriterManager, OutputRunMarker, resolve_max_open_files};

const PAIRED_RESYNC_WINDOW: usize = 1024;
const PAIRED_READER_QUEUE: usize = 256;
const BATCH_SIZE: usize = 1024;
const QUEUE_DEPTH_PER_WORKER: usize = 2;
const MAX_QUEUED_BATCHES: usize = 32;

#[derive(Clone)]
struct PipelineControl {
    inner: Arc<PipelineControlInner>,
}

struct PipelineControlInner {
    cancelled: AtomicBool,
    first_error: Mutex<Option<String>>,
    cancel_sender: Mutex<Option<Sender<()>>>,
    cancel_receiver: Receiver<()>,
    gate_wakers: Mutex<Vec<Weak<WorkerGateState>>>,
}

impl PipelineControl {
    fn new() -> Self {
        let (cancel_sender, cancel_receiver) = bounded::<()>(0);
        Self {
            inner: Arc::new(PipelineControlInner {
                cancelled: AtomicBool::new(false),
                first_error: Mutex::new(None),
                cancel_sender: Mutex::new(Some(cancel_sender)),
                cancel_receiver,
                gate_wakers: Mutex::new(Vec::new()),
            }),
        }
    }

    fn is_cancelled(&self) -> bool {
        self.inner.cancelled.load(Ordering::Acquire)
    }

    fn check(&self) -> Result<(), String> {
        if self.is_cancelled() {
            Err(self.root_error())
        } else {
            Ok(())
        }
    }

    fn cancel_receiver(&self) -> &Receiver<()> {
        &self.inner.cancel_receiver
    }

    fn fail(&self, error: impl Into<String>) -> String {
        let error = error.into();
        let first = {
            let mut first_error = self
                .inner
                .first_error
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            if first_error.is_none() {
                *first_error = Some(error);
                true
            } else {
                false
            }
        };

        if first {
            self.inner.cancelled.store(true, Ordering::Release);
            let sender = self
                .inner
                .cancel_sender
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner)
                .take();
            drop(sender);
            self.wake_worker_gates();
        }

        self.root_error()
    }

    fn first_error(&self) -> Option<String> {
        self.inner
            .first_error
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .clone()
    }

    fn root_error(&self) -> String {
        self.first_error()
            .unwrap_or_else(|| "Parallel pipeline cancelled".into())
    }

    fn register_worker_gate(&self, gate: &Arc<WorkerGateState>) {
        self.inner
            .gate_wakers
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .push(Arc::downgrade(gate));
    }

    fn wake_worker_gates(&self) {
        let gates: Vec<_> = {
            let mut registered = self
                .inner
                .gate_wakers
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            let gates = registered.iter().filter_map(Weak::upgrade).collect();
            registered.retain(|gate| gate.strong_count() > 0);
            gates
        };
        for gate in gates {
            // Taking the gate lock before notification prevents a cancellation
            // signal from racing between the wait predicate and Condvar::wait.
            let _limit = gate
                .limit
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            gate.changed.notify_all();
        }
    }
}

fn panic_payload(payload: Box<dyn Any + Send>) -> String {
    if let Some(message) = payload.downcast_ref::<&str>() {
        (*message).to_string()
    } else if let Some(message) = payload.downcast_ref::<String>() {
        message.clone()
    } else {
        "non-string panic payload".into()
    }
}

fn describe_output_key(key: OutputKey) -> String {
    let target = match key.target {
        OutputTarget::Sample { sample_id } => format!("sample {sample_id}"),
        OutputTarget::Unassigned => "unassigned output".into(),
    };
    let mate = match key.mate {
        OutputMate::Single => "single-end",
        OutputMate::R1 => "R1",
        OutputMate::R2 => "R2",
    };
    format!("{target} {mate}")
}

#[cfg(not(test))]
#[inline]
fn inject_test_failure(_stage: &str, _unit_id: u64) -> Result<(), String> {
    Ok(())
}

#[cfg(test)]
fn inject_test_failure(stage: &str, unit_id: u64) -> Result<(), String> {
    test_failure::checkpoint(stage, unit_id)
}

fn stage_boundary<T>(
    control: &PipelineControl,
    panic_context: &str,
    operation: impl FnOnce() -> Result<T, String>,
) -> Result<T, String> {
    match catch_unwind(AssertUnwindSafe(operation)) {
        Ok(Ok(value)) => Ok(value),
        Ok(Err(error)) => Err(control.fail(format!("{panic_context} failed: {error}"))),
        Err(payload) => Err(control.fail(format!(
            "{panic_context} panicked: {}",
            panic_payload(payload)
        ))),
    }
}

fn recv_or_cancel<T>(
    receiver: &Receiver<T>,
    control: &PipelineControl,
) -> Result<Option<T>, String> {
    control.check()?;
    crossbeam_channel::select_biased! {
        recv(control.cancel_receiver()) -> _ => Err(control.root_error()),
        recv(receiver) -> message => match message {
            Ok(value) => {
                control.check()?;
                Ok(Some(value))
            }
            Err(_) => {
                control.check()?;
                Ok(None)
            }
        }
    }
}

fn send_or_cancel<T>(
    sender: &Sender<T>,
    value: T,
    control: &PipelineControl,
    disconnected_error: &'static str,
) -> Result<(), String> {
    control.check()?;
    crossbeam_channel::select_biased! {
        recv(control.cancel_receiver()) -> _ => Err(control.root_error()),
        send(sender, value) -> result => match result {
            Ok(()) => control.check(),
            Err(_) if control.is_cancelled() => Err(control.root_error()),
            Err(_) => Err(disconnected_error.into()),
        }
    }
}

fn acquire_permit(
    receiver: &Receiver<()>,
    return_sender: &Sender<()>,
    control: &PipelineControl,
    name: &'static str,
) -> Result<PermitGuard, String> {
    control.check()?;
    let received = crossbeam_channel::select_biased! {
        recv(control.cancel_receiver()) -> _ => return Err(control.root_error()),
        recv(receiver) -> message => message,
    };
    received.map_err(|_| {
        if control.is_cancelled() {
            control.root_error()
        } else {
            format!("{name} permit channel disconnected")
        }
    })?;
    let permit = PermitGuard::new(return_sender.clone(), control.clone(), name);
    control.check()?;
    Ok(permit)
}

struct PermitGuard {
    return_sender: Option<Sender<()>>,
    control: Option<PipelineControl>,
    name: &'static str,
}

impl PermitGuard {
    fn new(return_sender: Sender<()>, control: PipelineControl, name: &'static str) -> Self {
        Self {
            return_sender: Some(return_sender),
            control: Some(control),
            name,
        }
    }

    #[cfg(test)]
    fn detached(name: &'static str) -> Self {
        Self {
            return_sender: None,
            control: None,
            name,
        }
    }
}

impl std::fmt::Debug for PermitGuard {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("PermitGuard")
            .field("name", &self.name)
            .finish_non_exhaustive()
    }
}

impl Drop for PermitGuard {
    fn drop(&mut self) {
        let Some(sender) = self.return_sender.take() else {
            return;
        };
        match sender.try_send(()) {
            Ok(()) | Err(TrySendError::Disconnected(())) => {}
            Err(TrySendError::Full(())) => {
                if !thread::panicking()
                    && let Some(control) = &self.control
                {
                    control.fail(format!(
                        "Internal error: {} permit was released more than once",
                        self.name
                    ));
                }
            }
        }
    }
}

fn compression_outstanding_limit(
    queue_capacity: usize,
    chunk_size: usize,
    memory_budget: usize,
) -> Result<usize, String> {
    if queue_capacity == 0 || chunk_size == 0 || memory_budget == 0 {
        return Err("Compression queue sizes must be greater than zero".into());
    }
    let count_limit = queue_capacity
        .checked_mul(2)
        .ok_or("Outstanding compression limit overflow")?;
    let byte_limit = (memory_budget / chunk_size).max(1);
    Ok(count_limit.min(byte_limit).max(1))
}

struct SharedState {
    routing: RoutingTree,
    output_layout: OutputLayout,
    r1_prefix: usize,
    r2_prefix: usize,
    write_unassigned: bool,
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
    permit: PermitGuard,
}

#[derive(Debug)]
struct BufferedOutput {
    key: OutputKey,
    fastq: Vec<u8>,
    records: u64,
    bases: u64,
}

#[derive(Debug)]
struct ProcessedBatch {
    id: u64,
    outputs: Vec<BufferedOutput>,
    counts: DemuxCounts,
    permit: PermitGuard,
}

#[derive(Debug)]
struct DirectCompressedOutput {
    key: OutputKey,
    member: Vec<u8>,
    records: u64,
    bases: u64,
    uncompressed_len: usize,
}

#[derive(Debug)]
struct DirectProcessedBatch {
    id: u64,
    outputs: Vec<DirectCompressedOutput>,
    counts: DemuxCounts,
    permit: PermitGuard,
}

#[derive(Debug)]
struct QueuedCompressionJob {
    job: CompressionJob,
    permit: PermitGuard,
}

#[derive(Debug)]
struct CompressedChunk {
    key: OutputKey,
    chunk_id: u64,
    member: Vec<u8>,
    uncompressed_len: usize,
    _permit: PermitGuard,
}

struct AggregationResult {
    counts: DemuxCounts,
    sample_qc: SampleQc,
    batches: u64,
    chunks: u64,
    uncompressed_bytes: u64,
    unassigned_records: u64,
}

#[derive(Debug, Default)]
struct WriterResult {
    chunks: u64,
    compressed_bytes: u64,
}

#[derive(Debug, Default)]
struct OutputBuffer {
    fastq: Vec<u8>,
    records: u64,
    bases: u64,
}

struct BatchBuffers {
    layout: OutputLayout,
    buffers: Vec<Option<OutputBuffer>>,
}

impl BatchBuffers {
    fn new(layout: OutputLayout) -> Result<Self, String> {
        Ok(Self {
            layout,
            buffers: (0..layout.stream_count()?).map(|_| None).collect(),
        })
    }

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
        let index = self.layout.stream_index(key)?;
        let buffer = self
            .buffers
            .get_mut(index)
            .ok_or("Output batch buffer index is missing")?
            .get_or_insert_with(OutputBuffer::default);
        let id = record.id.strip_prefix(b"@").unwrap_or(record.id.as_slice());

        buffer.fastq.extend_from_slice(b"@");
        buffer.fastq.extend_from_slice(id);
        buffer.fastq.extend_from_slice(b"\n");
        buffer.fastq.extend_from_slice(seq);
        buffer.fastq.extend_from_slice(b"\n+\n");
        buffer.fastq.extend_from_slice(qual);
        buffer.fastq.extend_from_slice(b"\n");
        buffer.records = buffer
            .records
            .checked_add(1)
            .ok_or("Output record count overflow")?;
        buffer.bases = buffer
            .bases
            .checked_add(u64::try_from(seq.len()).map_err(|_| "Read length overflow")?)
            .ok_or("Output base count overflow")?;

        Ok(())
    }

    fn finish(self) -> Result<Vec<BufferedOutput>, String> {
        self.buffers
            .into_iter()
            .enumerate()
            .filter_map(|(index, buffer)| buffer.map(|buffer| (index, buffer)))
            .map(|(index, buffer)| {
                Ok(BufferedOutput {
                    key: self.layout.key_at(index)?,
                    fastq: buffer.fastq,
                    records: buffer.records,
                    bases: buffer.bases,
                })
            })
            .collect()
    }
}

struct BatchSender {
    sender: Sender<WorkBatch>,
    permits: Receiver<()>,
    permit_returns: Sender<()>,
    next_batch_id: u64,
    items: Vec<WorkItem>,
    allocation: Option<AllocationController>,
    control: PipelineControl,
}

impl BatchSender {
    fn new(
        sender: Sender<WorkBatch>,
        permits: Receiver<()>,
        permit_returns: Sender<()>,
        control: PipelineControl,
    ) -> Self {
        Self {
            sender,
            permits,
            permit_returns,
            next_batch_id: 0,
            items: Vec::with_capacity(BATCH_SIZE),
            allocation: None,
            control,
        }
    }

    fn with_allocation_controller(
        sender: Sender<WorkBatch>,
        permits: Receiver<()>,
        permit_returns: Sender<()>,
        control: PipelineControl,
        allocation: AllocationController,
    ) -> Self {
        Self {
            allocation: Some(allocation),
            ..Self::new(sender, permits, permit_returns, control)
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
        let permit = acquire_permit(&self.permits, &self.permit_returns, &self.control, "batch")?;
        let batch = WorkBatch {
            id: self.next_batch_id,
            items,
            permit,
        };

        send_or_cancel(
            &self.sender,
            batch,
            &self.control,
            "Parallel worker queue disconnected",
        )?;

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
    state: Arc<WorkerGateState>,
    control: PipelineControl,
}

struct WorkerGateState {
    limit: Mutex<usize>,
    changed: Condvar,
    #[cfg(test)]
    waiting: AtomicUsize,
}

impl WorkerGate {
    fn new(initial_workers: usize, control: PipelineControl) -> Self {
        let state = Arc::new(WorkerGateState {
            limit: Mutex::new(initial_workers),
            changed: Condvar::new(),
            #[cfg(test)]
            waiting: AtomicUsize::new(0),
        });
        control.register_worker_gate(&state);
        Self { state, control }
    }

    fn wait_until_active(&self, worker_index: usize) -> Result<(), String> {
        let mut limit = self
            .state
            .limit
            .lock()
            .map_err(|_| "Demultiplexing worker gate was poisoned".to_string())?;

        while worker_index >= *limit {
            self.control.check()?;
            #[cfg(test)]
            self.state.waiting.fetch_add(1, Ordering::Release);
            limit = self
                .state
                .changed
                .wait(limit)
                .map_err(|_| "Demultiplexing worker gate was poisoned".to_string())?;
            #[cfg(test)]
            self.state.waiting.fetch_sub(1, Ordering::Release);
        }

        self.control.check()
    }

    fn set_limit(&self, workers: usize) -> Result<(), String> {
        self.control.check()?;
        let mut limit = self
            .state
            .limit
            .lock()
            .map_err(|_| "Demultiplexing worker gate was poisoned".to_string())?;
        *limit = workers;
        self.state.changed.notify_all();
        Ok(())
    }
}

#[derive(Clone)]
struct AllocationController {
    allocation: Arc<Mutex<AdaptiveAllocation>>,
    parallel_input: ParallelInput,
    worker_gates: AllocationGates,
    queue_capacity: usize,
}

#[derive(Clone)]
enum AllocationGates {
    Split(StageGates),
    Combined(WorkerGate),
}

#[derive(Clone)]
struct StageGates {
    demux: WorkerGate,
    compression: WorkerGate,
}

fn split_worker_slots(total: usize) -> (usize, usize) {
    if total <= 1 {
        return (1, 1);
    }
    let demux = total.div_ceil(3).max(1);
    (demux, total - demux)
}

fn report_thread_allocation(thread_plan: ThreadPlan, output_mode: ResolvedOutputMode) {
    if thread_plan.requested_threads < 2 {
        return;
    }
    match output_mode {
        ResolvedOutputMode::Direct => eprintln!(
            "Thread allocation: budget={} fastq-parsing={} input-decompression={} direct-workers={} overcommit={}",
            thread_plan.requested_threads,
            thread_plan.parser_threads,
            thread_plan.initial_input_threads,
            thread_plan.initial_worker_threads,
            thread_plan.budget_overcommit,
        ),
        ResolvedOutputMode::Buffered => {
            let (demux, compression) = split_worker_slots(thread_plan.initial_worker_threads);
            let overcommit = thread_plan
                .parser_threads
                .saturating_add(thread_plan.initial_input_threads)
                .saturating_add(demux)
                .saturating_add(compression)
                .saturating_sub(thread_plan.requested_threads)
                .max(thread_plan.budget_overcommit);
            eprintln!(
                "Thread allocation: budget={} fastq-parsing={} input-decompression={} demux={} compression={} overcommit={}",
                thread_plan.requested_threads,
                thread_plan.parser_threads,
                thread_plan.initial_input_threads,
                demux,
                compression,
                overcommit,
            );
        }
    }
}

impl StageGates {
    fn new(total: usize, control: PipelineControl) -> Self {
        let (demux, compression) = split_worker_slots(total);
        Self {
            demux: WorkerGate::new(demux, control.clone()),
            compression: WorkerGate::new(compression, control),
        }
    }

    fn set_total(&self, total: usize) -> Result<(), String> {
        let (demux, compression) = split_worker_slots(total);
        self.demux.set_limit(demux)?;
        self.compression.set_limit(compression)
    }
}

impl AllocationController {
    fn new_split(
        allocation: AdaptiveAllocation,
        parallel_input: ParallelInput,
        stage_gates: StageGates,
        queue_capacity: usize,
    ) -> Self {
        Self {
            allocation: Arc::new(Mutex::new(allocation)),
            parallel_input,
            worker_gates: AllocationGates::Split(stage_gates),
            queue_capacity,
        }
    }

    fn new_combined(
        allocation: AdaptiveAllocation,
        parallel_input: ParallelInput,
        worker_gate: WorkerGate,
        queue_capacity: usize,
    ) -> Self {
        Self {
            allocation: Arc::new(Mutex::new(allocation)),
            parallel_input,
            worker_gates: AllocationGates::Combined(worker_gate),
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
            self.set_worker_total(next.worker_threads)?;
        } else {
            self.set_worker_total(next.worker_threads)?;
            self.parallel_input.set_worker_limit(next.input_threads)?;
        }

        Ok(())
    }

    fn set_worker_total(&self, workers: usize) -> Result<(), String> {
        match &self.worker_gates {
            AllocationGates::Split(gates) => gates.set_total(workers),
            AllocationGates::Combined(gate) => gate.set_limit(workers),
        }
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

    let layout = if paired {
        ReadLayout::paired(args.r1_structure.as_deref(), args.r2_structure.as_deref())?
    } else {
        ReadLayout::single(args.structure.as_deref())?
    };

    let catalog = BarcodeCatalog::load(&args.barcodes, &layout)?;
    let samples = SampleSheet::load(&args.samples, &layout, &catalog)?;
    let routing = RoutingTree::new(&layout, &catalog, &samples, args.max_mismatches)?;
    let (r1_qc, r2_qc) = MateQcStats::for_layout(&layout);
    let mut r1_stats = args.fastq_stats.then_some(r1_qc);
    let mut r2_stats = if args.fastq_stats { r2_qc } else { None };
    let (r1_prefix, r2_prefix) = (routing.r1_prefix_len(), routing.r2_prefix_len());
    let output_layout = OutputLayout::new(samples.samples.len(), paired, args.write_unassigned);
    let expected_streams = output_layout.stream_count()?;
    let output_mode = resolve_output_mode(args.output_mode, expected_streams)?;
    report_output_mode(args.output_mode, output_mode, expected_streams);
    report_thread_allocation(thread_plan, output_mode);
    let max_open_files = resolve_max_open_files(expected_streams, args.max_open_files)?;

    let state = SharedState {
        routing,
        output_layout,
        r1_prefix,
        r2_prefix,
        write_unassigned: args.write_unassigned,
    };
    let control = PipelineControl::new();

    if output_mode == ResolvedOutputMode::Direct {
        return run_direct_parallel(
            args,
            inputs,
            paired,
            thread_plan,
            parallel_input,
            state,
            samples,
            max_open_files,
            control,
            r1_stats,
            r2_stats,
        );
    }

    let available_memory = available_memory_bytes().unwrap_or(512 * 1024 * 1024);
    let output_policy = calculate_output_buffer_policy(
        expected_streams,
        available_memory,
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

    let output_dir = args.output.clone();
    let writer = CompressedWriterManager::new(
        args.output.clone(),
        &samples,
        paired,
        args.write_unassigned,
        max_open_files,
        args.compression_level,
    )?;
    let accumulator = OutputAccumulator::new(
        output_layout,
        output_policy.chunk_size,
        output_policy.memory_budget,
    )?;
    let run_marker = OutputRunMarker::begin(&output_dir)?;

    let queue_capacity = worker_threads
        .checked_mul(QUEUE_DEPTH_PER_WORKER)
        .ok_or("Worker queue capacity overflow")?
        .clamp(1, MAX_QUEUED_BATCHES);

    let (work_tx, work_rx) = bounded::<WorkBatch>(queue_capacity);
    let (result_tx, result_rx) = bounded::<ProcessedBatch>(queue_capacity);
    let batch_outstanding_limit = queue_capacity
        .checked_mul(2)
        .ok_or("Outstanding batch limit overflow")?;
    let (batch_permit_tx, batch_permit_rx) = bounded::<()>(batch_outstanding_limit);
    for _ in 0..batch_outstanding_limit {
        batch_permit_tx
            .send(())
            .map_err(|_| "Could not initialize batch permits")?;
    }
    let (compression_tx, compression_rx) = bounded::<QueuedCompressionJob>(queue_capacity);
    let (compressed_tx, compressed_rx) = bounded::<CompressedChunk>(queue_capacity);
    let outstanding_limit = compression_outstanding_limit(
        queue_capacity,
        output_policy.chunk_size,
        output_policy.memory_budget,
    )?;
    let (permit_tx, permit_rx) = bounded::<()>(outstanding_limit);
    for _ in 0..outstanding_limit {
        permit_tx
            .send(())
            .map_err(|_| "Could not initialize compression permits")?;
    }
    let stage_gates = StageGates::new(worker_threads, control.clone());
    let allocation_controller = thread_plan.adaptive_allocation().map(|allocation| {
        AllocationController::new_split(
            allocation,
            parallel_input.clone(),
            stage_gates.clone(),
            queue_capacity,
        )
    });
    let (demux_headroom, compression_headroom) = split_worker_slots(worker_headroom);
    let compression_level = args.compression_level;
    let sample_count = samples.samples.len();

    let pipeline_result = thread::scope(|scope| {
        stage_boundary(&control, "Buffered pipeline coordinator", || {
            let writer_control = control.clone();
            let writer_handle = scope.spawn(move || {
                stage_boundary(&writer_control, "Output writer thread", || {
                    writer_loop(writer, &compressed_rx, output_layout, &writer_control)
                })
            });

            let mut compression_handles = Vec::with_capacity(compression_headroom);
            for worker_index in 0..compression_headroom {
                let receiver = compression_rx.clone();
                let sender = compressed_tx.clone();
                let gate = stage_gates.compression.clone();
                let worker_control = control.clone();
                compression_handles.push(scope.spawn(move || {
                    stage_boundary(
                        &worker_control,
                        &format!("Compression worker {worker_index}"),
                        || {
                            compression_loop(
                                worker_index,
                                &receiver,
                                &sender,
                                compression_level,
                                gate,
                                &worker_control,
                            )
                        },
                    )
                }));
            }
            drop(compression_rx);
            drop(compressed_tx);

            let aggregator_control = control.clone();
            let compression_permit_returns = permit_tx.clone();
            let aggregator_handle = scope.spawn(move || {
                stage_boundary(&aggregator_control, "Output aggregator thread", || {
                    aggregation_loop(
                        &result_rx,
                        &compression_tx,
                        &permit_rx,
                        &compression_permit_returns,
                        accumulator,
                        sample_count,
                        paired,
                        &aggregator_control,
                    )
                })
            });

            let mut worker_handles = Vec::with_capacity(demux_headroom);
            for worker_index in 0..demux_headroom {
                let worker_rx = work_rx.clone();
                let worker_tx = result_tx.clone();
                let worker_state = &state;
                let gate = stage_gates.demux.clone();
                let worker_control = control.clone();
                worker_handles.push(scope.spawn(move || {
                    stage_boundary(
                        &worker_control,
                        &format!("Routing worker {worker_index}"),
                        || {
                            worker_loop(
                                worker_index,
                                &worker_rx,
                                &worker_tx,
                                worker_state,
                                gate,
                                &worker_control,
                            )
                        },
                    )
                }));
            }
            drop(work_rx);
            drop(result_tx);

            let mut batcher = match allocation_controller {
                Some(controller) => BatchSender::with_allocation_controller(
                    work_tx,
                    batch_permit_rx,
                    batch_permit_tx,
                    control.clone(),
                    controller,
                ),
                None => {
                    BatchSender::new(work_tx, batch_permit_rx, batch_permit_tx, control.clone())
                }
            };

            let producer_result = stage_boundary(&control, "Input producer", || match &inputs {
                InputFiles::Single(path) => produce_single(
                    path,
                    &mut batcher,
                    r1_stats.as_mut(),
                    &parallel_input,
                    &control,
                ),
                InputFiles::Paired { r1, r2 } => produce_paired(
                    r1,
                    r2,
                    &mut batcher,
                    r1_stats.as_mut(),
                    r2_stats.as_mut(),
                    &parallel_input,
                    &control,
                ),
            });

            let batches_sent = match producer_result {
                Ok(()) => stage_boundary(&control, "Input batch finalization", || batcher.finish()),
                Err(error) => {
                    drop(batcher);
                    Err(error)
                }
            };

            let gate_result = if control.is_cancelled() {
                Err(control.root_error())
            } else {
                stage_boundary(&control, "Worker gate finalization", || {
                    stage_gates.set_total(worker_headroom)
                })
            };

            let mut worker_results = Vec::with_capacity(worker_handles.len());
            for handle in worker_handles {
                worker_results.push(match handle.join() {
                    Ok(result) => result,
                    Err(payload) => Err(control.fail(format!(
                        "Routing worker boundary panicked: {}",
                        panic_payload(payload)
                    ))),
                });
            }
            let aggregation_result = match aggregator_handle.join() {
                Ok(result) => result,
                Err(payload) => Err(control.fail(format!(
                    "Output aggregator boundary panicked: {}",
                    panic_payload(payload)
                ))),
            };
            let mut compression_results = Vec::with_capacity(compression_handles.len());
            for handle in compression_handles {
                compression_results.push(match handle.join() {
                    Ok(result) => result,
                    Err(payload) => Err(control.fail(format!(
                        "Compression worker boundary panicked: {}",
                        panic_payload(payload)
                    ))),
                });
            }
            let writer_result = match writer_handle.join() {
                Ok(result) => result,
                Err(payload) => Err(control.fail(format!(
                    "Output writer boundary panicked: {}",
                    panic_payload(payload)
                ))),
            };

            if let Some(error) = control.first_error() {
                return Err(error);
            }
            for result in worker_results {
                result?;
            }
            for result in compression_results {
                result?;
            }
            gate_result?;
            let aggregation = aggregation_result?;
            let writer_result = writer_result?;
            let expected_batches = batches_sent?;
            if aggregation.batches != expected_batches {
                return Err(format!(
                    "Parallel pipeline lost work: sent {expected_batches} batches but aggregated {}",
                    aggregation.batches
                ));
            }
            if writer_result.chunks != aggregation.chunks {
                return Err(format!(
                    "Parallel pipeline lost compressed output: emitted {} chunks but wrote {}",
                    aggregation.chunks, writer_result.chunks
                ));
            }
            Ok((aggregation, writer_result))
        })
    });

    let (aggregation, writer_result) = pipeline_result?;
    let counts = aggregation.counts;
    aggregation.sample_qc.verify(counts.assigned)?;
    verify_unassigned_output_counts(
        &counts,
        paired,
        args.write_unassigned,
        aggregation.unassigned_records,
    )?;

    if let Some(r1) = &r1_stats {
        write_fastq_stats(&output_dir, r1, r2_stats.as_ref())?;
    }

    let qc_summary =
        aggregation
            .sample_qc
            .write_report(&output_dir, &samples, args.low_sample_fraction)?;

    run_marker.complete()?;

    print_summary(&counts);
    print_sample_summary(qc_summary);
    print_chunk_summary(aggregation.uncompressed_bytes, &writer_result);

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_direct_parallel(
    args: DemuxArgs,
    inputs: InputFiles,
    paired: bool,
    thread_plan: ThreadPlan,
    parallel_input: ParallelInput,
    state: SharedState,
    samples: SampleSheet,
    max_open_files: usize,
    control: PipelineControl,
    mut r1_stats: Option<MateQcStats>,
    mut r2_stats: Option<MateQcStats>,
) -> Result<(), String> {
    let output_dir = args.output.clone();
    eprintln!("Output writers: max-open-files={max_open_files}");
    let writer = CompressedWriterManager::new(
        args.output.clone(),
        &samples,
        paired,
        args.write_unassigned,
        max_open_files,
        args.compression_level,
    )?;
    let run_marker = OutputRunMarker::begin(&output_dir)?;
    let worker_threads = thread_plan.initial_worker_threads;
    let worker_headroom = thread_plan.worker_headroom;
    let queue_capacity = worker_threads
        .checked_mul(QUEUE_DEPTH_PER_WORKER)
        .ok_or("Worker queue capacity overflow")?
        .clamp(1, MAX_QUEUED_BATCHES);
    let (work_tx, work_rx) = bounded::<WorkBatch>(queue_capacity);
    let (result_tx, result_rx) = bounded::<DirectProcessedBatch>(queue_capacity);
    let batch_outstanding_limit = queue_capacity
        .checked_mul(2)
        .ok_or("Outstanding batch limit overflow")?;
    let (batch_permit_tx, batch_permit_rx) = bounded::<()>(batch_outstanding_limit);
    for _ in 0..batch_outstanding_limit {
        batch_permit_tx
            .send(())
            .map_err(|_| "Could not initialize batch permits")?;
    }
    let worker_gate = WorkerGate::new(worker_threads, control.clone());
    let allocation_controller = thread_plan.adaptive_allocation().map(|allocation| {
        AllocationController::new_combined(
            allocation,
            parallel_input.clone(),
            worker_gate.clone(),
            queue_capacity,
        )
    });
    let compression_level = args.compression_level;
    let sample_count = samples.samples.len();

    let pipeline_result = thread::scope(|scope| {
        stage_boundary(&control, "Direct pipeline coordinator", || {
            let writer_control = control.clone();
            let writer_handle = scope.spawn(move || {
                stage_boundary(&writer_control, "Direct output writer thread", || {
                    direct_writer_loop(writer, &result_rx, sample_count, paired, &writer_control)
                })
            });
            let mut worker_handles = Vec::with_capacity(worker_headroom);
            for worker_index in 0..worker_headroom {
                let receiver = work_rx.clone();
                let sender = result_tx.clone();
                let worker_state = &state;
                let gate = worker_gate.clone();
                let worker_control = control.clone();
                worker_handles.push(scope.spawn(move || {
                    stage_boundary(
                        &worker_control,
                        &format!("Direct routing/compression worker {worker_index}"),
                        || {
                            direct_worker_loop(
                                worker_index,
                                &receiver,
                                &sender,
                                worker_state,
                                compression_level,
                                gate,
                                &worker_control,
                            )
                        },
                    )
                }));
            }
            drop(work_rx);
            drop(result_tx);

            let mut batcher = match allocation_controller {
                Some(controller) => BatchSender::with_allocation_controller(
                    work_tx,
                    batch_permit_rx,
                    batch_permit_tx,
                    control.clone(),
                    controller,
                ),
                None => {
                    BatchSender::new(work_tx, batch_permit_rx, batch_permit_tx, control.clone())
                }
            };
            let producer_result = stage_boundary(&control, "Input producer", || match &inputs {
                InputFiles::Single(path) => produce_single(
                    path,
                    &mut batcher,
                    r1_stats.as_mut(),
                    &parallel_input,
                    &control,
                ),
                InputFiles::Paired { r1, r2 } => produce_paired(
                    r1,
                    r2,
                    &mut batcher,
                    r1_stats.as_mut(),
                    r2_stats.as_mut(),
                    &parallel_input,
                    &control,
                ),
            });
            let batches_sent = match producer_result {
                Ok(()) => stage_boundary(&control, "Input batch finalization", || batcher.finish()),
                Err(error) => {
                    drop(batcher);
                    Err(error)
                }
            };
            let gate_result = if control.is_cancelled() {
                Err(control.root_error())
            } else {
                stage_boundary(&control, "Worker gate finalization", || {
                    worker_gate.set_limit(worker_headroom)
                })
            };

            let mut worker_results = Vec::with_capacity(worker_handles.len());
            for handle in worker_handles {
                worker_results.push(match handle.join() {
                    Ok(result) => result,
                    Err(payload) => Err(control.fail(format!(
                        "Direct worker boundary panicked: {}",
                        panic_payload(payload)
                    ))),
                });
            }
            let writer_result = match writer_handle.join() {
                Ok(result) => result,
                Err(payload) => Err(control.fail(format!(
                    "Direct output writer boundary panicked: {}",
                    panic_payload(payload)
                ))),
            };

            if let Some(error) = control.first_error() {
                return Err(error);
            }
            for result in worker_results {
                result?;
            }
            gate_result?;
            let (aggregation, writer_result) = writer_result?;
            let expected_batches = batches_sent?;
            if aggregation.batches != expected_batches {
                return Err(format!(
                    "Direct pipeline lost work: sent {expected_batches} batches but wrote {}",
                    aggregation.batches
                ));
            }
            if aggregation.chunks != writer_result.chunks {
                return Err(format!(
                    "Direct pipeline lost compressed output: emitted {} chunks but wrote {}",
                    aggregation.chunks, writer_result.chunks
                ));
            }
            Ok((aggregation, writer_result))
        })
    });

    let (aggregation, writer_result) = pipeline_result?;
    let counts = aggregation.counts;
    aggregation.sample_qc.verify(counts.assigned)?;
    verify_unassigned_output_counts(
        &counts,
        paired,
        args.write_unassigned,
        aggregation.unassigned_records,
    )?;
    if let Some(r1) = &r1_stats {
        write_fastq_stats(&output_dir, r1, r2_stats.as_ref())?;
    }
    let qc_summary =
        aggregation
            .sample_qc
            .write_report(&output_dir, &samples, args.low_sample_fraction)?;
    run_marker.complete()?;
    print_summary(&counts);
    print_sample_summary(qc_summary);
    print_chunk_summary(aggregation.uncompressed_bytes, &writer_result);
    Ok(())
}

fn direct_worker_loop(
    worker_index: usize,
    receiver: &Receiver<WorkBatch>,
    sender: &Sender<DirectProcessedBatch>,
    state: &SharedState,
    compression_level: u32,
    gate: WorkerGate,
    control: &PipelineControl,
) -> Result<(), String> {
    loop {
        gate.wait_until_active(worker_index)?;
        let Some(batch) = recv_or_cancel(receiver, control)? else {
            return Ok(());
        };
        let batch_id = batch.id;
        let processed = stage_boundary(
            control,
            &format!(
                "Direct routing/compression worker {worker_index} while processing batch {batch_id}"
            ),
            || {
                inject_test_failure("direct-worker", batch_id)?;
                process_direct_batch(batch, state, compression_level)
            },
        )?;
        send_or_cancel(
            sender,
            processed,
            control,
            "Direct output result queue disconnected",
        )?;
    }
}

fn process_direct_batch(
    batch: WorkBatch,
    state: &SharedState,
    compression_level: u32,
) -> Result<DirectProcessedBatch, String> {
    let processed = process_batch(batch, state)?;
    let ProcessedBatch {
        id,
        outputs,
        counts,
        permit,
    } = processed;
    let outputs = outputs
        .into_iter()
        .map(|output| {
            let uncompressed_len = output.fastq.len();
            let mut encoder = GzEncoder::new(Vec::new(), Compression::new(compression_level));
            encoder
                .write_all(&output.fastq)
                .map_err(|error| format!("Could not compress direct FASTQ batch: {error}"))?;
            let member = encoder
                .finish()
                .map_err(|error| format!("Could not finish direct FASTQ member: {error}"))?;
            Ok(DirectCompressedOutput {
                key: output.key,
                member,
                records: output.records,
                bases: output.bases,
                uncompressed_len,
            })
        })
        .collect::<Result<_, String>>()?;
    Ok(DirectProcessedBatch {
        id,
        outputs,
        counts,
        permit,
    })
}

fn direct_writer_loop(
    mut writer: CompressedWriterManager,
    receiver: &Receiver<DirectProcessedBatch>,
    sample_count: usize,
    paired: bool,
    control: &PipelineControl,
) -> Result<(AggregationResult, WriterResult), String> {
    let mut pending = BTreeMap::<u64, DirectProcessedBatch>::new();
    let mut next_batch_id = 0u64;
    let mut counts = DemuxCounts::default();
    let mut sample_qc = SampleQc::new(sample_count, paired);
    let mut writer_result = WriterResult::default();
    let mut uncompressed_bytes = 0u64;
    let mut unassigned_records = 0u64;

    while let Some(batch) = recv_or_cancel(receiver, control)? {
        let received_batch_id = batch.id;
        stage_boundary(
            control,
            &format!("Direct output writer while handling batch {received_batch_id}"),
            || {
                inject_test_failure("direct-writer", received_batch_id)?;
                if received_batch_id < next_batch_id
                    || pending.insert(received_batch_id, batch).is_some()
                {
                    return Err("Duplicate or stale direct batch ID".into());
                }
                while let Some(batch) = pending.remove(&next_batch_id) {
                    let DirectProcessedBatch {
                        outputs,
                        counts: batch_counts,
                        permit,
                        ..
                    } = batch;
                    for output in outputs {
                        control.check()?;
                        observe_output_metrics(
                            output.key,
                            output.records,
                            output.bases,
                            &mut sample_qc,
                            &mut unassigned_records,
                        )?;
                        control.check()?;
                        append_member(&mut writer, output.key, &output.member)?;
                        writer_result.chunks = writer_result
                            .chunks
                            .checked_add(1)
                            .ok_or("Chunk count overflow")?;
                        writer_result.compressed_bytes = writer_result
                            .compressed_bytes
                            .checked_add(
                                u64::try_from(output.member.len())
                                    .map_err(|_| "Compressed output size overflow")?,
                            )
                            .ok_or("Compressed output size overflow")?;
                        uncompressed_bytes = uncompressed_bytes
                            .checked_add(
                                u64::try_from(output.uncompressed_len)
                                    .map_err(|_| "Uncompressed output size overflow")?,
                            )
                            .ok_or("Uncompressed output size overflow")?;
                    }
                    counts.merge(&batch_counts);
                    drop(permit);
                    next_batch_id = next_batch_id.checked_add(1).ok_or("Batch ID overflow")?;
                }
                Ok(())
            },
        )?;
    }
    if !pending.is_empty() {
        return Err(format!(
            "Direct output ended with a missing batch before batch {next_batch_id}"
        ));
    }
    writer.finish()?;
    Ok((
        AggregationResult {
            counts,
            sample_qc,
            batches: next_batch_id,
            chunks: writer_result.chunks,
            uncompressed_bytes,
            unassigned_records,
        },
        writer_result,
    ))
}

fn worker_loop(
    worker_index: usize,
    receiver: &Receiver<WorkBatch>,
    sender: &Sender<ProcessedBatch>,
    state: &SharedState,
    gate: WorkerGate,
    control: &PipelineControl,
) -> Result<(), String> {
    loop {
        gate.wait_until_active(worker_index)?;
        let Some(batch) = recv_or_cancel(receiver, control)? else {
            return Ok(());
        };
        let batch_id = batch.id;
        let processed = stage_boundary(
            control,
            &format!("Routing worker {worker_index} while processing batch {batch_id}"),
            || {
                inject_test_failure("routing-worker", batch_id)?;
                process_batch(batch, state)
            },
        )?;
        send_or_cancel(
            sender,
            processed,
            control,
            "Processed batch queue disconnected",
        )?;
    }
}

fn process_batch(batch: WorkBatch, state: &SharedState) -> Result<ProcessedBatch, String> {
    let WorkBatch { id, items, permit } = batch;
    let mut buffers = BatchBuffers::new(state.output_layout)?;
    let mut counts = DemuxCounts::default();

    for item in items {
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

    let outputs = buffers.finish()?;

    Ok(ProcessedBatch {
        id,
        outputs,
        counts,
        permit,
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

#[allow(clippy::too_many_arguments)]
fn aggregation_loop(
    receiver: &Receiver<ProcessedBatch>,
    sender: &Sender<QueuedCompressionJob>,
    permits: &Receiver<()>,
    permit_returns: &Sender<()>,
    mut accumulator: OutputAccumulator,
    sample_count: usize,
    paired: bool,
    control: &PipelineControl,
) -> Result<AggregationResult, String> {
    let mut pending = BTreeMap::<u64, ProcessedBatch>::new();
    let mut next_batch_id = 0u64;
    let mut counts = DemuxCounts::default();
    let mut sample_qc = SampleQc::new(sample_count, paired);
    let mut chunks = 0u64;
    let mut uncompressed_bytes = 0u64;
    let mut unassigned_records = 0u64;

    while let Some(batch) = recv_or_cancel(receiver, control)? {
        let received_batch_id = batch.id;
        stage_boundary(
            control,
            &format!("Output aggregator while handling batch {received_batch_id}"),
            || {
                inject_test_failure("aggregation", received_batch_id)?;
                if received_batch_id < next_batch_id
                    || pending.insert(received_batch_id, batch).is_some()
                {
                    return Err("Duplicate or stale processed batch ID".into());
                }

                while let Some(batch) = pending.remove(&next_batch_id) {
                    aggregate_processed_batch(
                        batch,
                        &mut accumulator,
                        sender,
                        permits,
                        permit_returns,
                        &mut counts,
                        &mut sample_qc,
                        &mut chunks,
                        &mut uncompressed_bytes,
                        &mut unassigned_records,
                        control,
                    )?;
                    next_batch_id = next_batch_id.checked_add(1).ok_or("Batch ID overflow")?;
                }
                Ok(())
            },
        )?;
    }

    if !pending.is_empty() {
        return Err(format!(
            "Parallel output ended with a missing batch before batch {next_batch_id}"
        ));
    }

    send_compression_jobs(
        sender,
        permits,
        permit_returns,
        accumulator.finish()?,
        &mut chunks,
        &mut uncompressed_bytes,
        control,
    )?;

    Ok(AggregationResult {
        counts,
        sample_qc,
        batches: next_batch_id,
        chunks,
        uncompressed_bytes,
        unassigned_records,
    })
}

#[allow(clippy::too_many_arguments)]
fn aggregate_processed_batch(
    batch: ProcessedBatch,
    accumulator: &mut OutputAccumulator,
    sender: &Sender<QueuedCompressionJob>,
    permits: &Receiver<()>,
    permit_returns: &Sender<()>,
    counts: &mut DemuxCounts,
    sample_qc: &mut SampleQc,
    chunks: &mut u64,
    uncompressed_bytes: &mut u64,
    unassigned_records: &mut u64,
    control: &PipelineControl,
) -> Result<(), String> {
    let ProcessedBatch {
        outputs,
        counts: batch_counts,
        permit,
        ..
    } = batch;
    for output in outputs {
        control.check()?;
        observe_output_metrics(
            output.key,
            output.records,
            output.bases,
            sample_qc,
            unassigned_records,
        )?;
        let jobs = accumulator.append(output.key, output.fastq)?;
        send_compression_jobs(
            sender,
            permits,
            permit_returns,
            jobs,
            chunks,
            uncompressed_bytes,
            control,
        )?;
    }

    counts.merge(&batch_counts);
    drop(permit);
    Ok(())
}

fn send_compression_jobs(
    sender: &Sender<QueuedCompressionJob>,
    permits: &Receiver<()>,
    permit_returns: &Sender<()>,
    jobs: Vec<CompressionJob>,
    chunks: &mut u64,
    uncompressed_bytes: &mut u64,
    control: &PipelineControl,
) -> Result<(), String> {
    for job in jobs {
        // A credit is returned only after the writer appends this chunk. This
        // bounds queued raw jobs, compressed results, and reorder buffers as a
        // single end-to-end pool even if one compression task is delayed.
        let permit = acquire_permit(permits, permit_returns, control, "compression")?;
        *chunks = chunks
            .checked_add(1)
            .ok_or("Compression job count overflow")?;
        *uncompressed_bytes = uncompressed_bytes
            .checked_add(
                u64::try_from(job.fastq_bytes.len())
                    .map_err(|_| "Uncompressed output size overflow")?,
            )
            .ok_or("Uncompressed output size overflow")?;
        send_or_cancel(
            sender,
            QueuedCompressionJob { job, permit },
            control,
            "Compression job queue disconnected",
        )?;
    }
    Ok(())
}

fn compression_loop(
    worker_index: usize,
    receiver: &Receiver<QueuedCompressionJob>,
    sender: &Sender<CompressedChunk>,
    compression_level: u32,
    gate: WorkerGate,
    control: &PipelineControl,
) -> Result<(), String> {
    loop {
        gate.wait_until_active(worker_index)?;
        let Some(queued) = recv_or_cancel(receiver, control)? else {
            return Ok(());
        };
        let key = queued.job.key;
        let chunk_id = queued.job.chunk_id;
        let context = format!(
            "Compression worker {worker_index} for {} chunk {chunk_id}",
            describe_output_key(key)
        );
        let chunk = stage_boundary(control, &context, || {
            inject_test_failure("compression-worker", chunk_id)?;
            compress_job(queued.job, compression_level, queued.permit)
        })?;
        send_or_cancel(
            sender,
            chunk,
            control,
            "Compressed chunk queue disconnected",
        )?;
    }
}

fn compress_job(
    job: CompressionJob,
    compression_level: u32,
    permit: PermitGuard,
) -> Result<CompressedChunk, String> {
    let uncompressed_len = job.fastq_bytes.len();
    let mut encoder = GzEncoder::new(Vec::new(), Compression::new(compression_level));
    encoder
        .write_all(&job.fastq_bytes)
        .map_err(|error| format!("Could not compress FASTQ chunk: {error}"))?;
    let member = encoder
        .finish()
        .map_err(|error| format!("Could not finish FASTQ gzip member: {error}"))?;
    Ok(CompressedChunk {
        key: job.key,
        chunk_id: job.chunk_id,
        member,
        uncompressed_len,
        _permit: permit,
    })
}

#[derive(Default)]
struct OutputOrderState {
    next_chunk_id: u64,
    pending: BTreeMap<u64, CompressedChunk>,
}

fn writer_loop(
    mut writer: CompressedWriterManager,
    receiver: &Receiver<CompressedChunk>,
    layout: OutputLayout,
    control: &PipelineControl,
) -> Result<WriterResult, String> {
    let mut ordering: Vec<OutputOrderState> = (0..layout.stream_count()?)
        .map(|_| OutputOrderState::default())
        .collect();
    let mut result = WriterResult::default();

    while let Some(chunk) = recv_or_cancel(receiver, control)? {
        let key = chunk.key;
        let chunk_id = chunk.chunk_id;
        stage_boundary(
            control,
            &format!(
                "Output writer while handling {} chunk {chunk_id}",
                describe_output_key(key)
            ),
            || {
                inject_test_failure("buffered-writer", chunk_id)?;
                let index = layout.stream_index(key)?;
                let state = ordering
                    .get_mut(index)
                    .ok_or("Output ordering state is missing")?;
                if chunk_id < state.next_chunk_id || state.pending.insert(chunk_id, chunk).is_some()
                {
                    return Err("Duplicate or stale output chunk ID".into());
                }
                while let Some(chunk) = state.pending.remove(&state.next_chunk_id) {
                    control.check()?;
                    append_compressed_chunk(&mut writer, &chunk)?;
                    result.chunks = result.chunks.checked_add(1).ok_or("Chunk count overflow")?;
                    result.compressed_bytes = result
                        .compressed_bytes
                        .checked_add(
                            u64::try_from(chunk.member.len())
                                .map_err(|_| "Compressed output size overflow")?,
                        )
                        .ok_or("Compressed output size overflow")?;
                    state.next_chunk_id = state
                        .next_chunk_id
                        .checked_add(1)
                        .ok_or("Output chunk ID overflow")?;
                    drop(chunk);
                }
                Ok(())
            },
        )?;
    }

    if let Some((index, state)) = ordering
        .iter()
        .enumerate()
        .find(|(_, state)| !state.pending.is_empty())
    {
        return Err(format!(
            "Parallel output ended before chunk {} for output stream {index}",
            state.next_chunk_id
        ));
    }
    writer.finish()?;
    Ok(result)
}

fn append_compressed_chunk(
    writer: &mut CompressedWriterManager,
    chunk: &CompressedChunk,
) -> Result<(), String> {
    debug_assert!(chunk.uncompressed_len > 0);
    append_member(writer, chunk.key, &chunk.member)
}

fn append_member(
    writer: &mut CompressedWriterManager,
    key: OutputKey,
    member: &[u8],
) -> Result<(), String> {
    match key.target {
        OutputTarget::Sample { sample_id } => {
            writer.append_sample_member(sample_id, key.mate, member)
        }
        OutputTarget::Unassigned => writer.append_unassigned_member(key.mate, member),
    }
}

fn observe_output_metrics(
    key: OutputKey,
    records: u64,
    bases: u64,
    sample_qc: &mut SampleQc,
    unassigned_records: &mut u64,
) -> Result<(), String> {
    match key.target {
        OutputTarget::Sample { sample_id } => {
            sample_qc.observe_output(sample_id, key.mate, records, bases)
        }
        OutputTarget::Unassigned => {
            *unassigned_records = unassigned_records
                .checked_add(records)
                .ok_or("Unassigned output record count overflow")?;
            Ok(())
        }
    }
}

fn verify_unassigned_output_counts(
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

fn print_sample_summary(summary: SampleQcSummary) {
    eprintln!("  expected samples:  {}", summary.expected);
    eprintln!("  populated samples: {}", summary.populated);
    eprintln!("  missing samples:   {}", summary.missing);
    eprintln!("  low samples:       {}", summary.low);
}

fn print_chunk_summary(uncompressed_bytes: u64, writer: &WriterResult) {
    eprintln!(
        "  gzip chunks: {} ({:.1} KiB average uncompressed, {} compressed bytes)",
        writer.chunks,
        if writer.chunks == 0 {
            0.0
        } else {
            uncompressed_bytes as f64 / writer.chunks as f64 / 1024.0
        },
        writer.compressed_bytes,
    );
}

fn produce_single(
    path: &Path,
    batcher: &mut BatchSender,
    mut stats: Option<&mut MateQcStats>,
    parallel_input: &ParallelInput,
    control: &PipelineControl,
) -> Result<(), String> {
    let mut reader = parallel_input.open(path)?;

    loop {
        control.check()?;
        let Some(record) = reader.next_record() else {
            break;
        };
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
    sender: &Sender<OwnedFastqRecord>,
    parallel_input: ParallelInput,
    control: &PipelineControl,
) -> Result<(), String> {
    let mut reader = parallel_input.open(path)?;

    loop {
        control.check()?;
        let Some(record) = reader.next_record() else {
            break;
        };
        let record = record
            .map_err(|error| format!("FASTQ parse error in '{}': {error}", path.display()))?;

        let seq_cow = record.seq();
        let seq = seq_cow.as_ref();
        let qual = record
            .qual()
            .ok_or_else(|| format!("Input '{}' has no FASTQ quality scores", path.display()))?;

        let owned = OwnedFastqRecord::new(record.id(), seq, qual);
        send_or_cancel(
            sender,
            owned,
            control,
            "Paired FASTQ reader queue disconnected",
        )?;
    }

    Ok(())
}

fn recv_owned_record(
    receiver: &Receiver<OwnedFastqRecord>,
    stats: Option<&mut MateQcStats>,
    control: &PipelineControl,
) -> Result<Option<OwnedFastqRecord>, String> {
    match recv_or_cancel(receiver, control)? {
        Some(record) => {
            if let Some(stats) = stats {
                stats.update(&record.seq, &record.qual)?;
            }
            Ok(Some(record))
        }
        None => Ok(None),
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
    receiver: &Receiver<OwnedFastqRecord>,
    mate: OutputMate,
    mut stats: Option<&mut MateQcStats>,
    batcher: &mut BatchSender,
    control: &PipelineControl,
) -> Result<(), String> {
    while let Some(record) = recv_owned_record(receiver, stats.as_deref_mut(), control)? {
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
    r1_receiver: &Receiver<OwnedFastqRecord>,
    r2_receiver: &Receiver<OwnedFastqRecord>,
    batcher: &mut BatchSender,
    mut r1_stats: Option<&mut MateQcStats>,
    mut r2_stats: Option<&mut MateQcStats>,
    control: &PipelineControl,
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
                control,
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
                control,
            )?;
            return Ok(());
        }

        let mut progressed = false;

        if !r1_eof && r1_queue.len() < PAIRED_RESYNC_WINDOW {
            match recv_owned_record(r1_receiver, r1_stats.as_deref_mut(), control)? {
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
            match recv_owned_record(r2_receiver, r2_stats.as_deref_mut(), control)? {
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
    r1_receiver: &Receiver<OwnedFastqRecord>,
    r2_receiver: &Receiver<OwnedFastqRecord>,
    batcher: &mut BatchSender,
    mut r1_stats: Option<&mut MateQcStats>,
    mut r2_stats: Option<&mut MateQcStats>,
    control: &PipelineControl,
) -> Result<(), String> {
    loop {
        let r1_next = recv_owned_record(r1_receiver, r1_stats.as_deref_mut(), control)?;
        let r2_next = recv_owned_record(r2_receiver, r2_stats.as_deref_mut(), control)?;

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
                    control,
                )?;
            }
            (Some(r1), None) => {
                push_orphan(batcher, OutputMate::R1, r1)?;
                drain_receiver_as_orphans(
                    r1_receiver,
                    OutputMate::R1,
                    r1_stats.as_deref_mut(),
                    batcher,
                    control,
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
                    control,
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
    r1_stats: Option<&mut MateQcStats>,
    r2_stats: Option<&mut MateQcStats>,
    parallel_input: &ParallelInput,
    control: &PipelineControl,
) -> Result<(), String> {
    thread::scope(|scope| {
        stage_boundary(control, "Paired input scope coordinator", || {
            let (r1_tx, r1_rx) = bounded::<OwnedFastqRecord>(PAIRED_READER_QUEUE);
            let (r2_tx, r2_rx) = bounded::<OwnedFastqRecord>(PAIRED_READER_QUEUE);

            let r1_input = parallel_input.clone();
            let r2_input = parallel_input.clone();
            let r1_control = control.clone();
            let r2_control = control.clone();

            let r1_handle = scope.spawn(move || {
                stage_boundary(&r1_control, "R1 input reader", || {
                    stream_owned_records(r1_path, &r1_tx, r1_input, &r1_control)
                })
            });
            let r2_handle = scope.spawn(move || {
                stage_boundary(&r2_control, "R2 input reader", || {
                    stream_owned_records(r2_path, &r2_tx, r2_input, &r2_control)
                })
            });

            let producer_result = stage_boundary(control, "Paired input coordinator", || {
                produce_paired_from_receivers(&r1_rx, &r2_rx, batcher, r1_stats, r2_stats, control)
            });

            drop(r1_rx);
            drop(r2_rx);

            let r1_result = match r1_handle.join() {
                Ok(result) => result,
                Err(payload) => Err(control.fail(format!(
                    "R1 reader boundary panicked: {}",
                    panic_payload(payload)
                ))),
            };
            let r2_result = match r2_handle.join() {
                Ok(result) => result,
                Err(payload) => Err(control.fail(format!(
                    "R2 reader boundary panicked: {}",
                    panic_payload(payload)
                ))),
            };

            if let Some(error) = control.first_error() {
                return Err(error);
            }
            producer_result?;
            r1_result?;
            r2_result?;

            Ok(())
        })
    })
}

#[cfg(test)]
mod test_failure {
    use std::collections::HashSet;
    use std::sync::{Condvar, Mutex, OnceLock};
    use std::time::{Duration, Instant};

    struct FailurePoint {
        stage: String,
        target_unit: u64,
        wait_for_later_units: usize,
        return_error: bool,
        later_units: Mutex<HashSet<u64>>,
        changed: Condvar,
    }

    static FAILURE: OnceLock<Option<FailurePoint>> = OnceLock::new();

    fn configured_failure() -> Option<&'static FailurePoint> {
        FAILURE
            .get_or_init(|| {
                let specification = std::env::var("PLEXLESS_TEST_FAILURE").ok()?;
                let mut parts = specification.split(':');
                let stage = parts.next()?.to_string();
                let target_unit = parts.next()?.parse().ok()?;
                let wait_for_later_units = parts.next()?.parse().ok()?;
                if stage.is_empty() || parts.next().is_some() {
                    return None;
                }
                Some(FailurePoint {
                    stage,
                    target_unit,
                    wait_for_later_units,
                    return_error: std::env::var("PLEXLESS_TEST_FAILURE_KIND").as_deref()
                        == Ok("error"),
                    later_units: Mutex::new(HashSet::new()),
                    changed: Condvar::new(),
                })
            })
            .as_ref()
    }

    pub(super) fn checkpoint(stage: &str, unit_id: u64) -> Result<(), String> {
        let Some(failure) = configured_failure() else {
            return Ok(());
        };
        if failure.stage != stage {
            return Ok(());
        }

        if unit_id != failure.target_unit {
            if unit_id > failure.target_unit {
                let mut units = failure
                    .later_units
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
                if units.insert(unit_id) {
                    failure.changed.notify_all();
                }
            }
            return Ok(());
        }

        let deadline = Instant::now() + Duration::from_secs(10);
        let mut units = failure
            .later_units
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        while units.len() < failure.wait_for_later_units {
            let now = Instant::now();
            assert!(
                now < deadline,
                "failure injection for {stage} could not observe {} later units; saw {}",
                failure.wait_for_later_units,
                units.len()
            );
            let timeout = deadline.saturating_duration_since(now);
            let (next_units, wait) = failure
                .changed
                .wait_timeout(units, timeout)
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            units = next_units;
            assert!(
                !wait.timed_out() || units.len() >= failure.wait_for_later_units,
                "failure injection for {stage} could not observe {} later units; saw {}",
                failure.wait_for_later_units,
                units.len()
            );
        }
        if failure.return_error {
            Err(format!(
                "injected {stage} fatal error for unit {unit_id} after observing {} later units",
                units.len()
            ))
        } else {
            panic!(
                "injected {stage} panic for unit {unit_id} after observing {} later units",
                units.len()
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::{Read, Write};
    use std::process::{Command, Stdio};
    use std::time::{Duration, Instant};
    use std::time::{SystemTime, UNIX_EPOCH};

    use crate::cli::{ByteSizeSetting, OutputMode};
    use crate::samples::Sample;
    use crate::writer::INCOMPLETE_RUN_MARKER;
    use flate2::read::MultiGzDecoder;

    fn test_output_dir(name: &str) -> std::path::PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("plexless-{name}-{}-{nonce}", std::process::id()))
    }

    fn chunk(sample_id: u32, chunk_id: u64, text: &str) -> CompressedChunk {
        compress_job(
            CompressionJob {
                key: OutputKey {
                    target: OutputTarget::Sample { sample_id },
                    mate: OutputMate::Single,
                },
                chunk_id,
                fastq_bytes: text.as_bytes().to_vec(),
            },
            2,
            PermitGuard::detached("test compression"),
        )
        .unwrap()
    }

    fn read_gzip(path: &Path) -> String {
        let file = fs::File::open(path).unwrap();
        let mut decoder = MultiGzDecoder::new(file);
        let mut text = String::new();
        decoder.read_to_string(&mut text).unwrap();
        text
    }

    #[test]
    fn per_output_ordering_allows_independent_sample_progress() {
        let output = test_output_dir("per-output-ordering");
        let samples = SampleSheet {
            samples: vec![
                Sample {
                    name: "sample_0".into(),
                    barcode_ids: Vec::new(),
                },
                Sample {
                    name: "sample_1".into(),
                    barcode_ids: Vec::new(),
                },
                Sample {
                    name: "sample_2".into(),
                    barcode_ids: Vec::new(),
                },
            ],
        };
        let writer =
            CompressedWriterManager::new(output.clone(), &samples, false, false, 1, 2).unwrap();
        let layout = OutputLayout::new(3, false, false);
        let (sender, receiver) = bounded(8);
        let control = PipelineControl::new();

        let result = thread::scope(|scope| {
            let writer_control = control.clone();
            let handle =
                scope.spawn(move || writer_loop(writer, &receiver, layout, &writer_control));
            sender.send(chunk(0, 1, "@a1\nTGCA\n+\nIIII\n")).unwrap();
            sender.send(chunk(1, 0, "@b0\nCCCC\n+\nIIII\n")).unwrap();
            // Opening sample 2 evicts and flushes sample 1 while sample 0 is
            // still waiting for chunk 0, proving independent output progress.
            sender.send(chunk(2, 0, "@c0\nGGGG\n+\nIIII\n")).unwrap();

            let sample_1_path = output.join("sample_1.fastq.gz");
            let deadline = Instant::now() + Duration::from_secs(2);
            while fs::metadata(&sample_1_path).map_or(true, |metadata| metadata.len() == 0) {
                assert!(
                    Instant::now() < deadline,
                    "sample 1 did not progress independently"
                );
                thread::sleep(Duration::from_millis(5));
            }
            assert_eq!(read_gzip(&sample_1_path), "@b0\nCCCC\n+\nIIII\n");

            sender.send(chunk(0, 0, "@a0\nACGT\n+\nIIII\n")).unwrap();
            drop(sender);
            handle.join().unwrap().unwrap()
        });
        assert_eq!(result.chunks, 4);
        assert_eq!(
            read_gzip(&output.join("sample_0.fastq.gz")),
            "@a0\nACGT\n+\nIIII\n@a1\nTGCA\n+\nIIII\n"
        );
        assert_eq!(
            read_gzip(&output.join("sample_1.fastq.gz")),
            "@b0\nCCCC\n+\nIIII\n"
        );
        assert_eq!(
            read_gzip(&output.join("sample_2.fastq.gz")),
            "@c0\nGGGG\n+\nIIII\n"
        );
        fs::remove_dir_all(output).unwrap();
    }

    #[test]
    fn ordered_writers_fail_instead_of_skipping_missing_units() {
        let samples = SampleSheet {
            samples: vec![Sample {
                name: "sample_0".into(),
                barcode_ids: Vec::new(),
            }],
        };

        let direct_output = test_output_dir("direct-missing-batch");
        let direct_writer =
            CompressedWriterManager::new(direct_output.clone(), &samples, false, false, 1, 2)
                .unwrap();
        let (direct_sender, direct_receiver) = bounded(1);
        direct_sender
            .send(DirectProcessedBatch {
                id: 1,
                outputs: Vec::new(),
                counts: DemuxCounts::default(),
                permit: PermitGuard::detached("test batch"),
            })
            .unwrap();
        drop(direct_sender);
        let direct_error = match direct_writer_loop(
            direct_writer,
            &direct_receiver,
            samples.samples.len(),
            false,
            &PipelineControl::new(),
        ) {
            Err(error) => error,
            Ok(_) => panic!("direct writer skipped a missing batch"),
        };
        assert!(direct_error.contains("missing batch before batch 0"));

        let buffered_output = test_output_dir("buffered-missing-chunk");
        let buffered_writer =
            CompressedWriterManager::new(buffered_output.clone(), &samples, false, false, 1, 2)
                .unwrap();
        let (buffered_sender, buffered_receiver) = bounded(1);
        buffered_sender
            .send(chunk(0, 1, "@read\nACGT\n+\nIIII\n"))
            .unwrap();
        drop(buffered_sender);
        let buffered_error = writer_loop(
            buffered_writer,
            &buffered_receiver,
            OutputLayout::new(1, false, false),
            &PipelineControl::new(),
        )
        .unwrap_err();
        assert!(buffered_error.contains("before chunk 0"));

        fs::remove_dir_all(direct_output).unwrap();
        fs::remove_dir_all(buffered_output).unwrap();
    }

    #[test]
    fn worker_slot_split_keeps_routing_and_parallel_compression_live() {
        assert_eq!(split_worker_slots(1), (1, 1));
        assert_eq!(split_worker_slots(2), (1, 1));
        assert_eq!(split_worker_slots(8), (3, 5));
    }

    #[test]
    fn compression_credits_respect_count_and_byte_limits() {
        let mib = 1024 * 1024;
        assert_eq!(compression_outstanding_limit(32, mib, 512 * mib), Ok(64));
        assert_eq!(compression_outstanding_limit(32, mib, 16 * mib), Ok(16));
        assert_eq!(compression_outstanding_limit(32, 64 * mib, 64 * mib), Ok(1));
    }

    #[test]
    fn first_pipeline_failure_wins_and_wakes_worker_gates() {
        let control = PipelineControl::new();
        let gate = WorkerGate::new(0, control.clone());
        let gate_probe = gate.clone();
        let handle = thread::spawn(move || gate.wait_until_active(0));
        let deadline = Instant::now() + Duration::from_secs(2);
        while gate_probe.state.waiting.load(Ordering::Acquire) == 0 {
            assert!(Instant::now() < deadline, "worker did not enter gate wait");
            thread::yield_now();
        }

        assert_eq!(control.fail("root panic"), "root panic");
        assert_eq!(control.fail("secondary disconnect"), "root panic");
        assert_eq!(handle.join().unwrap(), Err("root panic".into()));
        assert_eq!(control.first_error().as_deref(), Some("root panic"));
    }

    #[test]
    fn permit_guard_returns_credit_during_unwind() {
        let control = PipelineControl::new();
        let (returns, permits) = bounded(1);
        returns.send(()).unwrap();
        let guard = acquire_permit(&permits, &returns, &control, "test").unwrap();
        assert!(permits.is_empty());

        let unwind = catch_unwind(AssertUnwindSafe(|| {
            let _guard = guard;
            panic!("test unwind");
        }));
        assert!(unwind.is_err());
        assert_eq!(permits.try_recv(), Ok(()));
        assert!(!control.is_cancelled());
    }

    #[test]
    fn injected_pipeline_failure_child() {
        let Ok(failure) = std::env::var("PLEXLESS_TEST_FAILURE") else {
            return;
        };
        let stage = failure.split(':').next().unwrap();
        let failure_kind =
            std::env::var("PLEXLESS_TEST_FAILURE_KIND").unwrap_or_else(|_| "panic".into());
        let paired = std::env::var("PLEXLESS_TEST_PAIRED").as_deref() == Ok("1");
        let mode = if stage.starts_with("direct-") {
            OutputMode::Direct
        } else {
            OutputMode::Buffered
        };

        let root = test_output_dir(&format!("panic-child-{stage}"));
        fs::create_dir_all(&root).unwrap();
        let barcodes = root.join("barcodes.tsv");
        let samples = root.join("samples.tsv");
        let reads = root.join("reads.fastq");
        let output = root.join("output");
        fs::write(&barcodes, "Set\tID\tSequence\nA\tA0\tACGT\n").unwrap();
        fs::write(&samples, "Sample\tA\nsample_0\tA0\n").unwrap();
        let mut fastq = fs::File::create(&reads).unwrap();
        let sequence = format!("ACGT{}", "T".repeat(196));
        let quality = "I".repeat(sequence.len());
        for index in 0..36_000 {
            writeln!(fastq, "@read{index} metadata").unwrap();
            writeln!(fastq, "{sequence}").unwrap();
            writeln!(fastq, "+source metadata").unwrap();
            writeln!(fastq, "{quality}").unwrap();
        }
        fastq.flush().unwrap();

        let r2 = paired.then(|| {
            let path = root.join("R2.fastq");
            let mut fastq = fs::File::create(&path).unwrap();
            let r2_sequence = "G".repeat(200);
            for index in 0..36_000 {
                writeln!(fastq, "@read{index} metadata").unwrap();
                writeln!(fastq, "{r2_sequence}").unwrap();
                writeln!(fastq, "+source metadata").unwrap();
                writeln!(fastq, "{quality}").unwrap();
            }
            fastq.flush().unwrap();
            path
        });
        let inputs = if let Some(r2) = &r2 {
            InputFiles::Paired {
                r1: reads.clone(),
                r2: r2.clone(),
            }
        } else {
            InputFiles::Single(reads.clone())
        };

        let args = DemuxArgs {
            reads: (!paired).then_some(reads.clone()),
            r1: paired.then_some(reads),
            r2,
            structure: (!paired).then(|| "R1_4A2T".into()),
            r1_structure: paired.then(|| "R1_4A2T".into()),
            r2_structure: paired.then(|| "R2_2T".into()),
            barcodes,
            samples,
            output: output.clone(),
            compression_level: 2,
            output_mode: mode,
            output_chunk_size: ByteSizeSetting::Bytes(128 * 1024),
            output_buffer_memory: ByteSizeSetting::Bytes(4 * 1024 * 1024),
            max_open_files: Some(2),
            max_mismatches: 0,
            fastq_stats: true,
            write_unassigned: false,
            low_sample_fraction: 0.05,
        };

        let error = run(args, 8, inputs).unwrap_err();
        let expected_root = if failure_kind == "error" {
            assert!(error.contains("failed"), "unexpected error: {error}");
            format!("injected {stage} fatal error")
        } else {
            assert!(error.contains("panicked"), "unexpected error: {error}");
            format!("injected {stage} panic")
        };
        assert!(
            error.contains(&expected_root),
            "root failure was not preserved: {error}"
        );
        match stage {
            "direct-worker" | "routing-worker" | "direct-writer" | "aggregation" => {
                assert!(error.contains("batch 0"), "batch context missing: {error}");
            }
            "compression-worker" | "buffered-writer" => {
                assert!(
                    error.contains("sample 0"),
                    "sample context missing: {error}"
                );
                assert!(error.contains("chunk 0"), "chunk context missing: {error}");
                assert!(
                    error.contains("single-end"),
                    "mate context missing: {error}"
                );
            }
            _ => panic!("unknown injected stage {stage}"),
        }
        assert!(
            !error.contains("queue disconnected") && !error.contains("missing batch"),
            "secondary shutdown error replaced the root cause: {error}"
        );
        assert!(output.join(INCOMPLETE_RUN_MARKER).is_file());
        assert!(!output.join("sample_metrics.tsv").exists());
        assert!(!output.join("fastq_stats.tsv").exists());
        assert!(!output.join("barcode_stats.tsv").exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn injected_failures_cancel_saturated_pipelines_without_deadlock() {
        let cases = [
            // Batch zero stays owned while every other batch credit reaches an
            // ordered pending map, so the producer is blocked on credit 29.
            ("direct-worker:0:27", false, "panic", 3),
            ("direct-writer:0:0", false, "panic", 1),
            ("routing-worker:0:27", false, "panic", 3),
            ("routing-worker:0:23", true, "panic", 3),
            // Chunk zero stays active while all 27 remaining end-to-end
            // compression credits are held by later ordered chunks.
            ("compression-worker:0:27", false, "panic", 3),
            ("aggregation:0:0", false, "panic", 1),
            ("buffered-writer:0:0", false, "panic", 1),
            // Ordinary fatal stage errors use the same first-error shutdown
            // path without going through panic capture.
            ("compression-worker:0:27", false, "error", 1),
        ];

        for (failure, paired, failure_kind, repetitions) in cases {
            for repetition in 0..repetitions {
                let log_path = test_output_dir(&format!(
                    "failure-parent-{}-{failure_kind}-{repetition}",
                    failure.replace(':', "-")
                ));
                let log = fs::File::create(&log_path).unwrap();
                let mut child = Command::new(std::env::current_exe().unwrap())
                    .args([
                        "--exact",
                        "parallel::tests::injected_pipeline_failure_child",
                        "--nocapture",
                        "--test-threads=1",
                    ])
                    .env("PLEXLESS_TEST_FAILURE", failure)
                    .env("PLEXLESS_TEST_FAILURE_KIND", failure_kind)
                    .env("PLEXLESS_TEST_PAIRED", if paired { "1" } else { "0" })
                    .stdout(Stdio::null())
                    .stderr(Stdio::from(log))
                    .spawn()
                    .unwrap();
                let deadline = Instant::now() + Duration::from_secs(20);
                let status = loop {
                    if let Some(status) = child.try_wait().unwrap() {
                        break status;
                    }
                    if Instant::now() >= deadline {
                        let _ = child.kill();
                        let _ = child.wait();
                        let diagnostics = fs::read_to_string(&log_path).unwrap_or_default();
                        panic!(
                            "injected {failure_kind} failure {failure} deadlocked:\n{diagnostics}"
                        );
                    }
                    thread::sleep(Duration::from_millis(10));
                };
                let diagnostics = fs::read_to_string(&log_path).unwrap();
                fs::remove_file(log_path).unwrap();
                assert!(
                    status.success(),
                    "injected failure {failure} failed its child assertions:\n{diagnostics}"
                );
            }
        }
    }
}
