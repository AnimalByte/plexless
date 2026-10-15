# Performance and tuning

Plexless is designed to scale with highly multiplexed sequencing runs. In
`--output-mode auto`, it uses the low-coordination direct path below 384
expected FASTQ output streams and the batch-owned buffered pipeline at 384
streams or more.

The table below reports median wall time from three repetitions of a synthetic
10-million-fragment demultiplexing benchmark on an Intel Core i7-13700K desktop
with 24 logical CPUs and DDR5 memory. Single-end (SE) runs contain 10 million
reads; paired-end (PE) runs contain 10 million read pairs. Both tools used gzip
compression level 2. Plexless used the portable release build
(`cargo build --release`), without `-C target-cpu=native`, and was compared
against fqtk 0.4.0.

These historical measurements predate the current preprocessing features,
including default paired overlap. They do not measure the new feature bundle.
The older v0.6.0 distal-scanning comparison is linked below.
See the [v0.6.0 benchmark report](pass-through-validation.md) for the current
baseline/disabled/enabled comparison; do not assume the historical rates include
the new scan.

| Mode | Samples | Plexless | fqtk | Plexless speed vs fqtk |
| --- | ---: | ---: | ---: | ---: |
| SE | 384 | 3.01 s | 6.43 s | 2.14x |
| SE | 700 | 3.06 s | 6.29 s | 2.06x |
| SE | 1,500 | 3.34 s | 6.44 s | 1.93x |
| PE | 384 | 6.16 s | 6.55 s | 1.06x |
| PE | 700 | 6.55 s | 6.63 s | 1.01x |
| PE | 1,500 | 6.92 s | 7.11 s | 1.03x |

These results are workload- and hardware-specific. fqtk was faster in this
benchmark at lower fan-out and at lower thread counts. Plexless's strongest
scaling appeared in the high-multiplexing buffered path as additional CPU
parallelism became available.

For the 1,500-sample fqtk comparisons, the process open-file soft limit was
raised from 1,024 to 8,192 after fqtk exhausted the default limit. Plexless
FASTQ output uses a bounded writer cache, so FASTQ descriptor use is capped
instead of requiring every output file to remain open simultaneously.


## Parallelism and compression

`--threads` is a total work budget, not simply a worker count. Plexless accounts
for FASTQ parser threads, parallel gzip decompression, routing, buffered
accumulation, and output compression. Buffered FASTQ fixes that allocation at
startup; it does not dynamically remove compression capacity while running.
The compiled routing representation is immutable and shared without hot-path
locks.

`--output-mode auto` selects `direct` below 384 expected output streams and
`buffered` at 384 or more. A stream is one sample/mate FASTQ, so 192 paired-end
samples are 384 streams; enabled unassigned outputs are included. Startup logs
the requested and selected mode. In parallel runs, `direct` compresses each
batch's populated sample buffers with low coordination overhead; one-thread
runs write through bounded, long-lived gzip writers instead. `buffered`
accumulates records across input batches, then parallel-compresses larger
per-output chunks; this avoids pathological gzip fragmentation at high
multiplexing. Either mode can be forced for diagnostics or unusual workloads.

For FASTQ input, buffered mode copies parser records once into immutable,
batch-owned byte blocks and routes compact checked descriptors into one ordered
accumulator. Trimming is represented by offsets, and records are serialized
directly into their destination buffers without per-record output allocation.
The default target is 32 KiB of uncompressed FASTQ per gzip member. One bounded
ready-job reservoir feeds shared zlib-rs compression workers; completed members
may finish out of order, but per-output sequence IDs preserve original record
order without globally ordering unrelated outputs. When the worker budget is
too small to reserve a dedicated accumulator while retaining routing and
compression, ordering and accumulation remain combined.

The accumulator memory budget is a pressure limit, not a preallocation. Auto
sizing reserves a conservative part of currently available memory, capped at
2 GiB. Bounded queues and compression credits apply backpressure throughout the
pipeline. CRAM-to-FASTQ retains its separately qualified chunk and scheduling
policy.

Buffering overrides apply to buffered mode and accept binary byte suffixes such
as `256K`, `512MiB`, and `2G`.
Explicit chunks must be 32 KiB–64 MiB and cannot exceed the buffer budget;
explicit buffer budgets must be at least 1 MiB and no more than half of memory
currently available at startup.

Output is standard concatenated gzip. The adaptive LRU writer cache uses at
most 256 files by default, further limited by the process file-descriptor limit
with 64 descriptors reserved on Linux. Other platforms use the conservative
256-file cap. Evicted outputs are safely reopened for append.

CRAM input uses rust-htslib/HTSlib background decode workers within the same
global budget. For FASTQ output, budgets of four or more assign approximately
one quarter to CRAM decode (capped at eight). CRAM output keeps HTSlib decoding
inline because controlled measurements found that background decoding reduced
whole-pipeline throughput. At budgets of eight or more, an ordered coordinator
dispatches to two destination-owner writer threads; each output CRAM belongs to
exactly one writer for the entire run. Smaller staged budgets use one combined
order/writer thread. The remaining budget goes to Plexless routing workers and
the exact allocation is logged. At a budget of one, the caller thread reads,
routes, and writes inline. A budget of two reports its one-thread minimum stage
overcommit.

Unlike concatenated gzip, CRAM output is not reopened for append. Plexless
preflights all possible output writers plus a 64-descriptor reserve before
creating the output directory. On Unix it raises a low soft `RLIMIT_NOFILE`
within the existing hard limit; if the hard limit is insufficient, it fails
before processing with exact required/soft/hard values. Normal execution
explicitly finalizes writers, reconciles successful writes, writes reports,
and removes `PLEXLESS_INCOMPLETE`; it does not decode all output files again.
The per-owner queues hold at most two ordered batch chunks each, so a hot
destination applies bounded backpressure rather than consuming unbounded
memory. Full rereading and `samtools` checks are qualification operations.

Parallel stages share one first-error cancellation state. Worker panics are
caught at batch or chunk boundaries with stage and unit context; cancellation
wakes bounded channel operations and worker gates. In-flight batch
and compression credits are ownership guards, so early returns and panic
unwinding release them. The process joins every spawned stage and returns the
original failure instead of reporting a disconnected queue or normal success.
