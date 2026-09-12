# plexless

`plexless` is a high-throughput hierarchical FASTQ and unmapped-CRAM
demultiplexer for single-end and paired-end sequencing data. Read structures describe physical
barcode pieces, while an ordinary tabular sample sheet defines parent-to-child
barcode namespaces and sample routes.

Plexless resolves one routing level at a time. After barcode A resolves, B is
matched only against B barcodes reachable beneath that A; C is then matched
only beneath the resolved A/B path, and so on. Unrelated child namespaces never
participate in matching or mismatch-safety validation.

Plexless is under active development. Validate demultiplexing results against
an established workflow before using it in production or clinical pipelines.

## Highlights

- True hierarchical routing compiled once at startup
- Plain or gzip-compressed FASTQ input
- Single-end and paired-end demultiplexing
- Generic barcode symbols `A-Z`, with `T` reserved for technical sequence
- Repeated barcode pieces within and across mates
- Per-segment forward or reverse-complement orientation
- Parent-local exact matching and correction by 0, 1, or 2 mismatches
- Parent-local handling of observed `N` bases
- Child sequence reuse across independent parent namespaces
- Structured-prefix trimming for assigned reads
- Optional unassigned output and separate biological/barcode FASTQ statistics
- Paired-read ID validation, bounded resynchronization, and orphan handling
- Bounded multithreaded parsing, demultiplexing, compression, and ordered output
- Unmapped CRAM input with FASTQ or metadata-preserving unmapped CRAM output

## Installation

Plexless links against the system zlib library. Source builds also compile
HTSlib for CRAM support and require a C toolchain, `make`, `pkg-config`, Clang,
and libclang in addition to a recent stable Rust toolchain and Cargo:

```bash
# Ubuntu/Debian
sudo apt install build-essential clang libclang-dev make pkg-config zlib1g-dev

# Fedora/RHEL
sudo dnf install clang clang-devel gcc make pkgconf-pkg-config zlib-devel

# Arch Linux
sudo pacman -S base-devel clang pkgconf zlib
```

### Install from crates.io

Install Plexless and verify the CLI:

```bash
cargo install plexless
plexless --help
```

`cargo install plexless` installs only the `plexless` executable into Cargo's
binary directory, normally `$HOME/.cargo/bin`.

### Download a prebuilt binary

Linux x86-64 archives are available from the public
[GitHub releases page](https://github.com/AnimalByte/plexless/releases). After
downloading the archive for the desired version:

```bash
tar -xzf plexless-vX.Y.Z-x86_64-unknown-linux-gnu.tar.gz
mkdir -p "$HOME/.local/bin"
install -m 0755 \
  plexless-vX.Y.Z-x86_64-unknown-linux-gnu/plexless \
  "$HOME/.local/bin/plexless"
plexless --help
```

Prebuilt Linux binaries require the zlib runtime library (`zlib1g` on
Ubuntu/Debian and `zlib` on Fedora/Arch).

### Build from source

Clone the public repository and build an optimized binary:

```bash
git clone https://github.com/AnimalByte/plexless.git
cd plexless
cargo build --release
```

The resulting executable is `target/release/plexless`. Alternatively, install
the current checkout through Cargo:

```bash
cargo install --path .
```

## Command-line reference

Plexless has one subcommand, `demux`:

```text
plexless [--threads <N>] demux [OPTIONS] \
  --barcodes <TSV> --samples <TSV> --output <DIR> \
  <--reads <FASTQ> | --r1 <FASTQ> --r2 <FASTQ> | --cram <CRAM>>
```

`--threads` is global and may also appear after `demux`. Show the authoritative
help for the installed version with:

```bash
plexless --help
plexless demux --help
```

### Global options

| Flag | Default | Description |
| --- | ---: | --- |
| `--threads <N>` | `1` | Global work budget shared by input decoding, routing, and output. FASTQ `1` uses the serial reference; CRAM uses no HTSlib background workers at `1`. Values must be positive integers. |
| `-h`, `--help` | — | Print help and exit. |

### Input and routing options

| Flag | Required | Description |
| --- | :---: | --- |
| `--reads <FASTQ>` | SE | Single-end FASTQ or FASTQ.gz input. Conflicts with `--r1` and `--r2`. |
| `--r1 <FASTQ>` | PE | Paired-end R1 FASTQ or FASTQ.gz input. Requires `--r2`. |
| `--r2 <FASTQ>` | PE | Paired-end R2 FASTQ or FASTQ.gz input. Requires `--r1`. |
| `--cram <CRAM>` | CRAM | One unmapped CRAM input. Conflicts with all FASTQ inputs. |
| `--read-mode <single\|paired>` | CRAM | Required declaration of CRAM record interpretation; never inferred. |
| `--output-format <fastq\|cram>` | CRAM | Required output choice. Existing FASTQ input needs no new option and remains FASTQ output. |
| `--structure <LAYOUT>` | SE | Single-end read structure, such as `R1_10A11B4T`. |
| `--r1-structure <LAYOUT>` | PE¹ | R1 read structure, such as `R1_10A11B4T`. |
| `--r2-structure <LAYOUT>` | PE¹ | R2 read structure, such as `R2_8C`. |
| `--barcodes <TSV>` | Always | Barcode whitelist with `Set`, `ID`, and `Sequence` columns. |
| `--samples <TSV>` | Always | Routing table with `Sample` first, followed by barcode-set columns. |
| `--max-mismatches <0-2>` | No (`1`) | Maximum substitutions allowed at each routing node. Observed `N` positions consume this budget. |

¹ Paired input requires one or both mate structures. A mate without a structure
is retained as biological sequence but supplies no barcode pieces.

Exactly one input mode is allowed. Single-end runs require `--reads` and
`--structure`. Paired-end runs require `--r1`, `--r2`, and at least one of
`--r1-structure` or `--r2-structure`.

CRAM input requires both `--read-mode` and `--output-format`. Only raw,
unmapped CRAM is accepted. Paired CRAM must declare queryname ordering/grouping
in `@HD` with `SO:queryname` or `GO:query`. See the complete
[CRAM contract](docs/cram.md), including validation and metadata policy.

### Output and performance options

| Flag | Default | Description |
| --- | ---: | --- |
| `-o`, `--output <DIR>` | Required | Output directory. Plexless creates it if absent; an existing directory must be empty. |
| `--compression-level <0-9>` | `2` | Compression level for output FASTQ gzip streams or CRAM files. Lower levels generally favor throughput. |
| `--output-mode <MODE>` | `auto` | `auto`, `direct`, or `buffered`. Auto selects direct below 384 expected streams and buffered at 384 or more. |
| `--output-chunk-size <SIZE>` | `auto` | Buffered-mode target for uncompressed bytes per gzip member. Explicit values must be 32 KiB–64 MiB and no larger than the buffer budget. |
| `--output-buffer-memory <SIZE>` | `auto` | Buffered-mode memory budget for active output accumulators. Explicit values must be at least 1 MiB and no more than half of currently available memory. |
| `--max-open-files <N>` | adaptive | FASTQ: positive override for the bounded writer cache. CRAM: must be at least the simultaneous destination count; CRAM writers cannot be reopened for append. |

`SIZE` is a positive integer followed optionally by `B`, `K`/`KB`/`KiB`,
`M`/`MB`/`MiB`, or `G`/`GB`/`GiB`; suffixes are case-insensitive and use
binary multiples. Examples: `262144`, `256K`, `512MiB`, and `2G`. Decimal
quantities such as `1.5G` are not accepted. `auto` is valid for both size
options.

The defaults are recommended. Buffer-size overrides affect only buffered mode.
Expected streams count sample/mate files, plus enabled unassigned files: 16 SE
samples produce 16 streams, while 16 PE samples produce 32. See
[Parallelism and compression](#parallelism-and-compression) for memory sizing,
ordering, and file-cache details.

### QC and unassigned-read options

| Flag | Default | Description |
| --- | ---: | --- |
| `--fastq-stats` | Off | Write biological-read statistics to `fastq_stats.tsv` and physical barcode-region statistics to `barcode_stats.tsv`. |
| `--write-unassigned` | Off | Write unmatched, ambiguous, unrouted, and orphan reads untrimmed to `unassigned*.fastq.gz`. |
| `--low-sample-fraction <0-1>` | `0.05` | Mark a populated sample `LOW_REPRESENTATION` below this fraction of the nonzero sample median. `0` disables low-representation classifications; missing-sample warnings remain mandatory. |
| `-h`, `--help` | — | Print `demux` help and exit. |

Every successful run writes `sample_metrics.tsv`; `--fastq-stats` adds the
biological and barcode reports. Output filenames and report fields are described in
[Output behavior](#output-behavior).

### Single-end example

```bash
plexless --threads 8 demux \
  --reads reads.fastq.gz \
  --structure R1_10A11B4T \
  --barcodes barcodes.tsv \
  --samples samples.tsv \
  --output demux_out \
  --max-mismatches 1 \
  --compression-level 2 \
  --write-unassigned
```

### Paired-end example with orientation normalization

```bash
plexless --threads 16 demux \
  --r1 reads_R1.fastq.gz \
  --r2 reads_R2.fastq.gz \
  --r1-structure R1_10A11B4T \
  --r2-structure 'R2_10A(rc)11B(rc)' \
  --barcodes barcodes.tsv \
  --samples samples.tsv \
  --output demux_out \
  --fastq-stats
```

Quote structures containing `(rc)` so shells do not interpret the
parentheses.

### Unmapped CRAM examples

Single-end CRAM to FASTQ requires an explicit lossy format choice:

```bash
plexless --threads 8 demux \
  --cram reads.cram --read-mode single --output-format fastq \
  --structure R1_10A11B4T \
  --barcodes barcodes.tsv --samples samples.tsv \
  --output demux_fastq --write-unassigned
```

Plexless reports that SAM/CRAM metadata not representable in FASTQ will not be
present in this output. It does not encode arbitrary tags in FASTQ headers.

Paired, queryname-grouped CRAM to metadata-preserving CRAM:

```bash
plexless --threads 16 demux \
  --cram paired.queryname.cram --read-mode paired --output-format cram \
  --r1-structure R1_10A11B4T \
  --r2-structure 'R2_10A(rc)11B(rc)' \
  --barcodes barcodes.tsv --samples samples.tsv \
  --output demux_cram --write-unassigned
```

FASTQ to CRAM and aligned CRAM input/output are not supported.

### Explicit output tuning

Auto mode needs no tuning. To force the low-coordination path for a small run:

```bash
plexless --threads 8 demux [INPUT AND ROUTING OPTIONS] \
  --output demux_out --output-mode direct
```

To force buffered output with explicit advanced limits:

```bash
plexless --threads 16 demux [INPUT AND ROUTING OPTIONS] \
  --output demux_out \
  --output-mode buffered \
  --output-chunk-size 512KiB \
  --output-buffer-memory 512MiB \
  --max-open-files 128
```

`[INPUT AND ROUTING OPTIONS]` is illustrative shell notation, not literal
syntax; supply one of the complete input forms shown above.

## Input configuration

### Read structures

A read structure describes the structured prefix as repeated `<length><symbol>`
segments. For example:

```text
R1_10A11B4T
```

means 10 bases of barcode A, 11 bases of barcode B, and 4 technical bases.
The remaining sequence is biological insert. The entire declared prefix is
removed from assigned output.

| Symbol | Meaning |
| --- | --- |
| `A-S`, `U-Z` | Generic logical barcode/routing levels. |
| `T` | Technical sequence that is trimmed but not decoded. |

Append `(rc)` to normalize only that barcode segment by reverse complement:

```text
R1_10A11B
R2_10A(rc)11B(rc)
```

Without `(rc)`, orientation is forward. Entire reads and quality strings are
never reverse-complemented. Whitelist sequences are always canonical.

Repeated occurrences of a symbol form one logical barcode. Pieces are
concatenated in this order:

1. All matching R1 pieces in R1 structure order.
2. All matching R2 pieces in R2 structure order.

For example, with `R1_2A2B` and `R2_2A2B`, R1 A sequence `AC`, and normalized
R2 A sequence `GT`, logical A is `ACGT`. Asymmetric structures such as
`R1_10A11B` plus `R2_10A(rc)` are supported.

### Barcode whitelist

The barcode file is tab-separated with this exact header:

```tsv
Set	ID	Sequence
A	A01	ACGT
A	A02	TGCA
B	B01	GATC
B	B02	CTAG
```

Whitelist sequences must contain only `A`, `C`, `G`, and `T`. IDs must be
unique within each set, all sequences in a set must have one length, and that
length must equal the logical length from the read structure.

Root sequences must be unique. A child sequence may be reused under
independent parent paths, including with different IDs. Reuse is rejected only
when both sequences become candidates in the same routing node.

### Sample sheet

The sample sheet remains a simple tab-separated table:

```tsv
Sample	A	B	C
sample_1	A01	B01	C01
sample_2	A01	B01	C02
sample_3	A01	B02	C01
sample_4	A02	B01	C03
```

The first column must be `Sample`. Remaining columns are the read-structure
barcode symbols in alphabetical order, skipping reserved `T`. Every row must
specify every level and every barcode ID must exist in the whitelist.

The columns define routing order. The table above compiles to:

```text
A decoder
  -> resolved A child
     -> parent-specific B decoder
        -> resolved A/B child
           -> parent-specific C decoder
              -> sample
```

## Hierarchical matching and routing

At runtime Plexless:

1. Executes the compiled extraction plan for the root symbol.
2. Encodes the logical barcode using the packed 2-bit representation.
3. Decodes against that node's exact/correction index.
4. Uses the decoder-local integer to select a child node or sample.
5. Lazily repeats only if another routing level is reached.

No barcode from an unrelated parent namespace is scanned. Correction safety is
also node-local: mismatch tolerance `e` requires sibling Hamming distance
`2e + 1`, but imposes no distance requirement across independent parents.

Observed `N` bases consume mismatch budget rather than acting as free
wildcards. N-aware matching scans only the current node's candidates.

| Outcome | Meaning |
| --- | --- |
| Assigned | Every level resolves through the hierarchy to a sample. |
| Unmatched | A barcode cannot resolve against candidates below its resolved parent. |
| Ambiguous | The current node does not have one deterministic match. |
| Unrouted | A valid root barcode has no compiled sample route. |

Plexless never silently selects the first ambiguous candidate. Cross-level
candidate-path rescue is not currently performed.

## Output behavior

Assigned single-end reads are written as `<sample>.fastq.gz`. Paired output is
written as `<sample>_R1.fastq.gz` and `<sample>_R2.fastq.gz`.

For CRAM output, each sample is one `<sample>.cram`; paired mates remain in the
same file. Enabled unassigned output is `unassigned.cram`. The preserved input
header gains one uniquely named Plexless `@PG` entry. No CRAM index is created.
See [CRAM metadata preservation](docs/cram.md#metadata-preservation-policy).

Assigned reads have their declared structured prefixes removed independently
from each mate. With `--write-unassigned`, unassigned and orphan records are
written untrimmed to corresponding `unassigned*.fastq.gz` files.
Because those filenames are reserved, a worklist sample that sanitizes to
`unassigned` is rejected when unassigned output is enabled.

`--fastq-stats` writes two input-level QC reports. `fastq_stats.tsv` summarizes
only the biological suffix after each mate's complete declared prefix.
`barcode_stats.tsv` summarizes every physical barcode segment separately and
identifies its mate, symbol, mate-local repeated-piece number, sequencer cycles,
and declared orientation. Technical `T` segments are excluded from both
reports. Both reports include read/base counts, length statistics, GC and N
percentages, mean quality, and Q20/Q30 percentages.

Statistics cover every parsed input record, regardless of whether it is
assigned, unmatched, ambiguous, unrouted, or orphaned. A mate without a read
structure is entirely biological. Barcode statistics use the raw observed
bases and qualities before mismatch correction or reverse-complement
normalization; orientation is retained as report metadata. `BarcodeSymbol`
refers to the read-structure symbol, not a decoded whitelist ID, and each row
aggregates all input observations. `StartCycle` and `EndCycle` are one-based and
inclusive. For reads shorter than the declared prefix, available barcode cycles
are counted and the biological observation has length zero. `Reads` in a
barcode row is the number of parsed records for that mate, including zero-length
observations when the record ends before that segment.

Every successful run writes `sample_metrics.tsv`, including all worklist
samples—even those with zero reads. It reports assigned fragments, mate-level
record and base counts, fractions of the run total and nonzero median, and an
`OK`, `LOW_REPRESENTATION`, or `MISSING` status. For paired data, one assigned
fragment means one R1/R2 pair. Missing samples always produce a warning;
low-representation warnings default to below 5% of the nonzero median.

When output processing starts, Plexless creates `PLEXLESS_INCOMPLETE` in the
output directory. It removes this marker only after all output writers finish,
counts reconcile, and optional statistics plus sample metrics are flushed. Any
failure before that completion point leaves the marker behind, and a later run
refuses to reuse that nonempty directory. Treat every directory that contains
the marker as partial output.

Plexless preserves the complete parser-provided read header and trims sequence
and quality by identical prefix lengths. Output uses a normalized `+` separator
line; metadata from the input `+` line is not retained.

## Paired-read synchronization

R1 and R2 IDs are normalized and checked before demultiplexing. When one mate
is missing, Plexless searches ahead up to 1,024 records for the next shared ID,
emits intervening records as R1 or R2 orphans, and resumes in order. It returns
an error instead of guessing if synchronization cannot be recovered within the
window. Orphans are never assigned to a sample or supplied with missing barcode
pieces.

## Parallelism and compression

`--threads` is a total work budget, not simply a worker count. Plexless accounts
for FASTQ parser threads and dynamically shares remaining capacity between
parallel gzip decompression and demultiplexing/output compression. The compiled
routing representation is immutable and shared without hot-path locks.

`--output-mode auto` selects `direct` below 384 expected output streams and
`buffered` at 384 or more. A stream is one sample/mate FASTQ, so 192 paired-end
samples are 384 streams; enabled unassigned outputs are included. Startup logs
the requested and selected mode. In parallel runs, `direct` compresses each
batch's populated sample buffers with low coordination overhead; one-thread
runs write through bounded, long-lived gzip writers instead. `buffered`
accumulates records across input batches, then parallel-compresses larger
per-output chunks; this avoids pathological gzip fragmentation at high
multiplexing. Either mode can be forced for diagnostics or unusual workloads.

In buffered mode, input scheduling batches remain small while auto sizing
budgets a conservative part of currently available memory and normally chooses
128 KiB–1 MiB uncompressed chunks. Bounded queues apply backpressure through
compression and writing. Chunks may finish out of order, but per-output
sequence IDs preserve original record order without globally ordering
unrelated outputs.

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
wakes bounded channel operations and adaptive worker gates. In-flight batch
and compression credits are ownership guards, so early returns and panic
unwinding release them. The process joins every spawned stage and returns the
original failure instead of reporting a disconnected queue or normal success.

## Limits and validation

- Maximum logical barcode length: 32 bases
- Barcode symbols: `A-S` and `U-Z`; `T` is reserved
- Mismatch correction: 0, 1, or 2 substitutions
- Sample-sheet rows must contain complete paths; partial paths are unsupported
- Ambiguity stops traversal; path-aware rescue is not yet implemented
- Paired resynchronization lookahead: 1,024 records
- Automatic gzip chunk target: 128 KiB–1 MiB uncompressed
- Automatic output buffer budget: 25% of available memory, capped at 2 GiB
- Automatic output mode: direct below 384 streams; buffered at 384 or more
- Automatic output writer cache: expected streams, capped at 256 and the safe descriptor limit
- CRAM input: unmapped records only; missing sequence or quality is rejected
- Paired CRAM: authoritative queryname grouping and at most one primary R1/R2 per QNAME
- CRAM output: all potentially populated files remain open; descriptor limits are preflighted and writer memory scales approximately linearly with output count

Startup validation covers read-structure and `(rc)` syntax, required barcode
sets, barcode bases and logical lengths, duplicate IDs and complete paths,
unknown sample-sheet IDs, root and sibling sequence conflicts, node-local
correction safety with parent-path errors, and the 32-base packed-encoding
limit. Short reads are safely classified as unmatched.

## Architecture and tests

- [`docs/architecture.md`](docs/architecture.md) documents compilation and the
  allocation-free read hot path.
- [`docs/cram.md`](docs/cram.md) documents the supported CRAM matrix, validation,
  metadata transformation contract, and limitations.
- [`tests/hierarchical_routing.rs`](tests/hierarchical_routing.rs) covers local
  namespaces, dense-grid selection, correction safety, reuse, and generic symbols.
- [`tests/scalable_output_qc.rs`](tests/scalable_output_qc.rs) covers high
  multiplexing, SE/PE ordering, header preservation, and sample QC.

Development checks:

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all
cargo build --release --all-targets
cargo test --release
```

## License

Plexless is dual-licensed under the Apache License, Version 2.0, or the MIT
license, at your option. See `LICENSE-APACHE` and `LICENSE-MIT`.
