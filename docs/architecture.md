# Processing architecture

## Streaming data flow

Plexless uses a bounded streaming pipeline:

```text
FASTQ input
  -> Needletail parsing / parallel gzip decoding
  -> bounded record batches
  -> demultiplexing workers
  -> selected direct or buffered output strategy
  -> bounded output-writer cache
```

Paired input validates normalized record IDs and uses bounded lookahead to
resynchronize, emitting records present in only one input as orphans. Assigned
reads have each mate's declared structured prefix trimmed; unassigned and
orphan reads remain untrimmed.

## Compiled hierarchical model

`ReadLayout::compile_extraction_plans` converts parsed read structures into one
`BarcodeExtractionPlan` per logical symbol. Each plan contains its logical
length and an ordered list of `(mate, start, end, orientation)` pieces. The
piece order is all R1 segments in structure order followed by all R2 segments
in structure order.

`RoutingTree::new` compiles the ordinary tabular sample sheet into immutable
contiguous `CompiledNode` values. Each node stores:

- its generic barcode symbol;
- the index of its extraction plan;
- a `Decoder` built only from barcode candidates reachable at that parent
  path;
- a target array indexed directly by the decoder-local barcode ID.

A target is another node, a sample ID, or an unrouted terminal. The root uses
the complete root whitelist so a known but unused root barcode remains
unrouted. Every deeper node contains only sample-sheet children beneath its
specific parent path.

When the sample sheet is a complete Cartesian product, startup compilation
selects a dense representation instead. Each decoder still maps packed bases
and the N mask to a dense whitelist-local index; row-major index arithmetic
then selects a precompiled sample slot. Sparse and irregular sheets continue
to use hierarchical nodes, including their parent-local correction rules.

## Read hot path

```text
check R1/R2 structured-prefix lengths
  -> assemble root symbol into a 32-byte stack buffer
  -> 2-bit encode and decode in the root node
  -> index the resolved target
  -> if it is a node, lazily assemble that node's symbol
  -> repeat until sample, unmatched, ambiguous, or unrouted
```

Traversal allocates no strings, vectors, or maps per read. An exact/corrected
decoder result is a node-local integer used directly to index the target
array. Invalid observed bases classify the read as unmatched. Observed `N`
bases retain the existing mismatch-budget behavior, but the N-aware scan sees
only the current node's candidates.

## Output mode selection

Before input processing, `OutputLayout` counts sample/mate streams as
`(samples + enabled unassigned target) * mates`. `--output-mode auto` selects
the direct pipeline below 384 streams and the buffered pipeline at 384 or more.
Explicit `direct` and `buffered` values bypass this stable stream-count
threshold. The startup report records the decision.

Parallel direct mode routes, serializes, and compresses each populated output
fragment inside its work-batch worker. Its bounded writer reorders completed
batches by input batch ID before appending their ordinary gzip members. This
path has low coordination and memory overhead, but creates many small members
as multiplexing increases. With one thread, direct mode instead writes records
to long-lived `GzEncoder`s held by the bounded LRU cache; eviction finishes a
member and reopening appends a new standards-compliant member.

## Buffered output aggregation and ordering

The buffered pipeline separates scheduling batches from compression:

```text
bounded 1,024-record work batches
  -> parallel route/serialize into per-output batch fragments
  -> input-ordered persistent sample/mate accumulation
  -> bounded compression jobs with per-output chunk IDs
  -> parallel independent gzip members
  -> per-output ordered append through an adaptive LRU writer cache
```

The accumulator briefly orders completed input batches so fragments enter each
sample in source order. Compression completion is not globally ordered. The
writer keeps a separate next-chunk ID for every `(sample_id, mate)` output, so
one ready sample can progress while another sample's earlier chunk is still
compressing. Paired R1 and R2 fragments are accumulated from the same ordered
batch stream and therefore retain correspondence.

Each gzip member is independently valid. Appending members creates an ordinary
concatenated `.fastq.gz` readable by standard multi-member gzip tools.

## Bounded-memory model

Buffered auto sizing reads Linux `MemAvailable` and, when present, tighter cgroup memory
headroom. One quarter is assigned to active output bytes, clamped to 16 MiB–2
GiB and never above half of available memory. Dividing that budget by expected
sample/mate/unassigned streams produces a 128 KiB–1 MiB target chunk. A global
pressure check flushes the largest accumulator when the logical budget is
exceeded; vector growth can transiently reserve more capacity, which is why the
policy leaves substantial headroom.

Work, processed-batch, compression-job, and compressed-result channels are all
bounded. One end-to-end credit pool bounds active and completed batches plus
the pre-aggregation reorder map. A second bounds queued jobs, active
compression, completed members, and writer reorder state. Its credit count is
the smaller of twice the worker-queue depth and the accumulator budget divided
by the target chunk size. A chunk can exceed its target by at most one
serialized work-batch fragment, and an active encoder briefly holds compressed
and uncompressed forms; these bounded transients are included in the process
headroom rather than the accumulator's logical-byte counter. Paired
resynchronization is separately bounded to 1,024 records per mate, plus
256-record reader queues.

The writer cache defaults to the smaller of expected streams, 256, and the
process descriptor soft limit minus a 64-descriptor reserve. Reopening a file
uses append mode and does not invalidate its concatenated gzip stream.
Linux reads the limit from `/proc/self/limits`; other platforms use the safe
256-file fallback unless `--max-open-files` selects a lower value.

Direct mode does not allocate persistent accumulators or separate compression
queues. Its work/result queues and end-to-end batch credits remain bounded, and
each in-flight batch contains at most one serialized fragment per populated
output. Serial direct caches at most its configured number of gzip encoders;
parallel direct caches that many raw writers while each active worker owns one
transient encoder.

Buffered parallel runs reserve at least one routing worker and one compression worker.
When the requested thread budget leaves only one shared CPU slot, this creates
a one-thread minimum pipeline overcommit; the startup allocation line reports
it explicitly. Larger budgets split the existing shared worker allocation and
remain within it. `--threads 1` continues to use the serial reference path.
Direct parallel runs use the shared worker pool for combined routing and
compression rather than splitting it into two stages.

## Completion and QC invariants

Success requires every submitted input batch to reach aggregation, every
compression job to reach the writer, all partial accumulators to flush, and all
writers to finish. Per-sample fragment totals must equal the global assigned
count; paired samples must have equal fragment, R1, and R2 counts. When
unassigned output is enabled, emitted record counts are reconciled with
unmatched, ambiguous, unrouted, and orphan counts before `sample_metrics.tsv`
and the terminal QC summary are written.

With `--fastq-stats`, input records are projected through the compiled read
layout before aggregation. `fastq_stats.tsv` covers only the biological suffix
after the complete structured prefix. `barcode_stats.tsv` covers each physical
barcode segment by mate, symbol, mate-local repeated-piece number, and cycle
range; technical `T` segments are not included. Barcode observations remain in
raw sequencer orientation because GC, N, and aggregate quality statistics are
orientation invariant. Cycle ranges are one-based and inclusive, and symbols
refer to the read layout rather than decoded whitelist IDs. Collection stays
in the input producer so assigned, unassigned, and orphan records have
identical QC semantics in serial, direct, and buffered modes.

The run creates `PLEXLESS_INCOMPLETE` after output setup and never removes it
through `Drop`. Successful orchestration explicitly removes it only after
writer finalization, count reconciliation, and flushed statistics and sample
metrics. Thus any Rust error, caught worker panic, process interruption, or
final-report failure before that completion point leaves an unmistakable
partial-run marker. Existing marked directories are rejected rather than
appended to.

## Failure and cancellation model

All parallel stages share a cancellation object containing an atomic fast path,
a first-error slot, and a disconnected zero-capacity channel used as a
broadcast wakeup. The first fatal error wins; later shutdown symptoms cannot
replace its diagnostic. Work/result sends and receives, paired-reader queues,
batch credits, and compression credits select between their normal operation
and cancellation. Adaptive worker gates are registered with the same object
and notified while holding their predicate mutex, avoiding lost wakeups.

Batch and compression credits are RAII guards. Ownership moves with a batch or
chunk and `Drop` returns its token on success, ordinary error, or unwind. This
bounds memory during normal operation and prevents a panic from leaking the
credit needed by shutdown. A missing ordered batch or per-output chunk is never
skipped: the run cancels, drops queued work, joins all scoped threads, and
returns the preserved root cause.

Panic capture is limited to clear boundaries: routing/compression jobs,
aggregation and writing units, input-reader threads, and stage/coordinator
threads. Expected parse, configuration, compression, and filesystem failures
remain ordinary `Result` errors. Cancellation can wake application-level
queues and condition variables; as with conventional synchronous I/O, it
cannot forcibly interrupt a kernel call or third-party decompressor that never
returns.

## Node-local correction safety

Each node independently validates minimum sibling Hamming distance
`2 * max_mismatches + 1` and independently precomputes exact and correction
lookups. An error identifies the parent path plus both conflicting IDs.

Consequently, close or identical child sequences are allowed in independent
parent namespaces. They conflict only when both are candidates in the same
node. This keeps correction validation scoped to barcodes that can compete at
the same routing position.

## Extension point

`Decoder` still exposes exact, corrected, ambiguous, and unmatched outcomes,
while `RoutingTree` owns traversal and node targets. A future path-aware rescue
implementation can maintain multiple `(node, score)` states over the same
compiled representation instead of replacing the extraction or tree model.
