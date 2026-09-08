# plexless

`plexless` is a streaming FASTQ demultiplexer written in Rust.

It is designed for barcode layouts that are described directly from the read structure, including single-end and paired-end data where a logical barcode may be split across R1 and R2.

## Features

- Single-end and paired-end FASTQ demultiplexing
- FASTQ and gzip-compressed FASTQ input
- User-defined read structures
- Hierarchical barcode routing compiled from a tabular sample sheet
- Generic logical barcode symbols `A-Z` except reserved technical symbol `T`
- Technical sequence segments with `T`
- Barcode concatenation across repeated segments and across mates
- Per-segment forward or reverse-complement orientation
- Exact barcode matching
- Configurable mismatch correction
- Explicit handling of `N` bases
- Sample routing from a separate sample sheet
- Structured-prefix trimming
- Optional raw-input FASTQ statistics
- Optional unassigned FASTQ output
- Paired-read orphan detection and bounded resynchronization
- Bounded gzip writer cache for large sample counts
- Gzip-compressed per-sample output

## Prebuilt binary

Prebuilt Linux x86-64 binaries are available from the
[GitHub releases page](https://github.com/AnimalByte/plexless/releases). The
release binary includes zlib, so zlib does not need to be installed separately.

## Build from source

### System dependencies

`plexless` uses the system zlib implementation through `flate2` for output gzip compression and for the single-thread reference input path. Parallel gzip input uses `rapidgzip-core`, which validates gzip framing/checksums while decoding ordinary gzip streams in parallel.

### Ubuntu / Debian

```bash
sudo apt update
sudo apt install zlib1g-dev
```

### Fedora / RHEL

```bash
sudo dnf install zlib-devel
```

### Arch Linux

```bash
sudo pacman -S zlib
```

A recent stable Rust toolchain is also required.

### Build

```bash
cargo build --release
```

The executable will be:

```text
target/release/plexless
```

For migration compatibility, the build also produces
`target/release/simplex` with the same command-line interface.

For development:

```bash
cargo run -- --help
cargo run -- demux --help
```

## Read structures

A read structure describes the structured prefix of a read.

For example:

```text
R1_10A11B4T
```

means:

```text
R1
├── 10 bases of barcode A
├── 11 bases of barcode B
├── 4 bases of technical sequence T
└── remaining sequence is biological insert
```

The structured prefix is removed from assigned reads before output.

### Segment symbols

| Symbol | Meaning |
| --- | --- |
| `A-S`, `U-Z` | Generic logical barcode/routing levels |
| `T` | Technical sequence; trimmed but not used for sample identity |

`T` does not appear in the barcode whitelist or sample sheet.

### Repeated barcode symbols

If the same barcode symbol appears more than once, its pieces form one logical barcode.

For paired reads, pieces are concatenated in R1 order followed by R2 order.

Example:

```text
R1_2A2B2T
R2_2A2B2T
```

If:

```text
R1 A = AC
R2 A = GT
```

then the logical A barcode is:

```text
ACGT
```

The same rule applies to every logical barcode symbol.

By default, pieces are interpreted in forward orientation. Append `(rc)` to
reverse-complement only that barcode piece before concatenation:

```text
R1_4A4B
R2_4A(rc)4B(rc)
```

Orientation is per segment and can be mixed between symbols and mates. Quality
strings are not transformed because barcode decoding does not use qualities.
Whitelist sequences are always written in canonical orientation.

## Barcode whitelist

Barcodes are supplied as a tab-separated file.

Example:

```tsv
Set	ID	Sequence
A	A01	ACGT
A	A02	TGCA
B	B01	GATC
B	B02	CTAG
```

Required columns:

- `Set`
- `ID`
- `Sequence`

The barcode sequence length must match the logical length defined by the read structure.

Whitelist barcode sequences:

- must contain valid DNA bases
- cannot contain `N`
- must have unique IDs within their set

The root barcode set must have unique sequences. A child sequence may be
reused under independent parent paths, including under different IDs. Reused
sequences are rejected only if they become siblings in the same routing node.

Logical barcodes are currently limited to 32 bases.

## Sample sheet

The sample sheet maps complete hierarchical barcode paths to samples.

Example:

```tsv
Sample	A	B
sample_1	A01	B01
sample_2	A01	B02
sample_3	A02	B01
sample_4	A02	B02
```

Only barcode sets present in the read structure belong in the sample sheet.

Columns follow deterministic symbol order (`A`, then `B`, then `C`, and so
on, skipping `T`). Every row must specify every level. At startup, the table is
compiled into a routing tree:

```text
decode A at the root
  -> choose the A parent
  -> decode B only against children beneath that A
  -> continue until the path reaches a sample
```

The barcode whitelist and sample sheet remain deliberately separate:

- the barcode whitelist defines the complete valid barcode space
- the sample sheet defines parent-local child namespaces and sample routes

An unused root barcode remains `unrouted`. At child levels, a barcode not
reachable from the resolved parent is `unmatched`; unrelated namespaces are
never scanned to reinterpret it.

## Single-end example

```bash
plexless demux \
    --reads reads.fastq.gz \
    --structure R1_10A11B4T \
    --barcodes barcodes.tsv \
    --samples samples.tsv \
    --output demux_out
```

Equivalent with the short output option:

```bash
plexless demux \
    --reads reads.fastq.gz \
    --structure R1_10A11B4T \
    --barcodes barcodes.tsv \
    --samples samples.tsv \
    -o demux_out
```

## Paired-end example

```bash
plexless demux \
    --r1 reads_R1.fastq.gz \
    --r2 reads_R2.fastq.gz \
    --r1-structure R1_2A2B2T \
    --r2-structure R2_2A2B2T \
    --barcodes barcodes.tsv \
    --samples samples.tsv \
    --output demux_out
```

At least one paired read structure must be supplied.

## Barcode correction

The mismatch threshold is configurable:

```bash
--max-mismatches 1
```

The default is:

```text
1
```

Supported correction distances are currently 0, 1, and 2 mismatches.

### Unique correction

`plexless` validates each effective routing-node candidate set before
processing reads.

For a mismatch tolerance of `e`, sibling barcodes must be separated
sufficiently to make correction unique. For example, one-mismatch correction
requires sibling distance of at least 3. Barcodes beneath different parents do
not need to be mutually distinguishable because they are never decoded
together.

An unsafe node is rejected before demultiplexing begins, with its parent path
and conflicting barcode IDs in the error.

### `N` bases

An observed `N` is not treated as a free wildcard.

It consumes mismatch budget.

For example, with:

```text
--max-mismatches 1
```

a barcode containing one `N` may still be recoverable, while a barcode containing two `N` positions cannot be rescued solely within that budget.

## Routing outcomes

A read or pair can end in one of four routing states.

### Assigned

Every routing level resolves within its current parent namespace and the leaf
maps to a sample.

### Unmatched

At least one barcode cannot be decoded within the configured mismatch limit
against the candidates reachable from its resolved parent.

### Ambiguous

Barcode evidence does not support a unique correction.

### Unrouted

A valid root barcode has no compiled sample route. Deeper barcodes outside a
parent's candidate namespace are unmatched rather than globally decoded.

## Unassigned reads

Use:

```bash
--write-unassigned
```

to write unmatched, ambiguous, unrouted, and orphan records.

Unassigned records are written in their original untrimmed form.

## Paired-end synchronization and orphans

Normal paired reads are processed directly when their normalized read IDs match.

For example:

```text
R1: read1  read2  read3
R2: read1         read3
```

`read2` is treated as an R1 orphan. `plexless` searches forward for the next shared read ID, resynchronizes at `read3`, and continues processing.

Orphan records:

- are counted separately for R1 and R2
- are never assigned to a sample
- are written raw when `--write-unassigned` is enabled
- do not cause immediate failure merely because one mate is missing

Resynchronization uses a bounded lookahead window of 1024 records per mate. If synchronization cannot be recovered within that bound, `plexless` returns an error instead of guessing.

Trailing records present in only one mate file are also treated as orphans.

## Output

Single-end assigned output:

```text
sample_1.fastq.gz
sample_2.fastq.gz
...
```

Paired-end assigned output:

```text
sample_1_R1.fastq.gz
sample_1_R2.fastq.gz
sample_2_R1.fastq.gz
sample_2_R2.fastq.gz
...
```

With unassigned output enabled, corresponding unassigned FASTQ files are also created.

Assigned reads have the full structured prefix removed.

Unassigned and orphan reads remain untrimmed.

## Compression

`plexless` uses `flate2` with the system `zlib` backend for output gzip compression. With `--threads 1`, gzip input also follows the original Needletail/system-zlib reference path. With a larger CPU budget, ordinary seekable gzip input is decoded in parallel internally and then passed to Needletail as one ordered byte stream.

The gzip compression level is configurable:

```bash
--compression-level 2
```

The default is level 2.

`plexless` maintains a bounded cache of open output writers so large sample sheets do not require every output file to remain open simultaneously.

The current maximum number of open writers is 64.

When an output is evicted and later reopened, another gzip member is appended. Concatenated gzip members are valid gzip and decompress as one continuous FASTQ stream.

## FASTQ statistics

Enable raw-input FASTQ statistics with:

```bash
--fastq-stats
```

Statistics are computed before structured-prefix trimming.

The report includes:

- read count
- base count
- minimum read length
- maximum read length
- mean read length
- GC percentage
- N percentage
- mean quality
- Q20 percentage
- Q30 percentage

For paired data, R1 and R2 are reported separately.

## Validation and safety checks

Before or during processing, `plexless` validates conditions including:

- malformed read structures
- malformed `(rc)` orientation modifiers
- invalid barcode symbols
- missing barcode sets
- duplicate barcode IDs
- duplicate root or node-local sibling sequences
- barcode length mismatch
- duplicate sample names
- duplicate sample barcode combinations
- unknown barcode IDs in the sample sheet
- logical barcodes longer than 32 bases
- unsafe node-local mismatch-correction geometry with parent-path context
- FASTQ parse failures
- missing quality scores
- reads shorter than the declared structured prefix

Short reads are classified as unmatched rather than causing an out-of-bounds trim.

## Quality gates

The project currently uses:

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test
cargo audit
cargo build --release
```

The integration suite covers:

- single-end A+B+T demultiplexing
- paired-end A+B+T demultiplexing
- barcode assembly across mates
- repeated-symbol assembly in R1-then-R2 order
- segment-local reverse-complement normalization
- asymmetric paired layouts
- hierarchical exact routing through two and three levels
- child sequence reuse across independent parent namespaces
- parent-local correction safety and true sibling collision rejection
- generic hierarchies deeper than three levels
- technical-sequence trimming
- one-mismatch correction
- `N` handling
- unmatched reads
- unrouted barcode combinations
- short-read safety
- FASTQ statistics
- configuration validation
- orphan recovery in both directions
- trailing orphan records
- writer-cache eviction and gzip reopening

## Current implementation limits

- Maximum logical barcode length: 32 bases
- Barcode symbols: uppercase `A-Z`, with `T` reserved for technical sequence
- Paired resynchronization lookahead: 1024 records
- Output writer cache: 64 open writers
- Mismatch correction: 0, 1, or 2 mismatches

## Development status

`plexless` is currently under active development. Conservative traversal stops
on an ambiguous decision at the current node; cross-level candidate-path
rescue/scoring is intentionally not implemented yet.

The parallel pipeline uses bounded queues, ordered batch emission, parallel output gzip compression, paired-read resynchronization, and automatic parallel gzip input. Performance work remains performance-tested so that correctness, deterministic routing, bounded memory use, and output integrity remain the primary constraints.

## Automatic parallel gzip input

Pass one FASTQ file for single-end data or one R1/R2 pair for paired-end data.
`plexless` handles parallel input processing internally.

For example:

```bash
plexless --threads 16 demux \
  --reads reads.fastq.gz \
  --structure R1_10A11B4T \
  --barcodes barcodes.tsv \
  --samples samples.tsv \
  --output out
```

For a normal seekable `.fastq.gz`, `plexless` internally parallelizes gzip
decompression and presents Needletail with one ordered decompressed stream. It
does not write temporary decompressed FASTQ files.

For paired-end data, R1 and R2 are decoded as independent ordered streams under
one shared decompression budget. They are then passed through the existing
read-ID pairing and bounded resynchronization logic. R1 and R2 are therefore
not independently partitioned into FASTQ-record ranges that could break orphan
handling.

### CPU budget

`--threads` is a total CPU-work budget rather than a raw demultiplexing-worker
count. Any positive integer is accepted; the value does not need to be a power
of two. `plexless` accounts for one FASTQ parser for single-end input and two
FASTQ parsers for paired-end input. It shares the remaining budget between
input decompression and demultiplexing/output compression.

For plain FASTQ input, the full post-parser budget is assigned to
demultiplexing/output compression. For gzip input with a sufficiently large
budget, the initial decompression allocation is approximately one quarter of
the shared single-end budget or one third of the shared paired-end budget. The
rest begins in demultiplexing/output compression.

For example:

| `--threads` | Single-end: parser / input / workers | Paired-end: parsers / input / workers |
|------------:|--------------------------------------:|----------------------------------------:|
| 16 | 1 / 4 / 11 | 2 / 4 / 10 |
| 32 | 1 / 8 / 23 | 2 / 10 / 20 |
| 37 | 1 / 9 / 27 | 2 / 11 / 24 |
| 48 | 1 / 12 / 35 | 2 / 15 / 31 |
| 64 | 1 / 16 / 47 | 2 / 20 / 42 |

These are initial allocations. During processing, sustained work-queue
pressure moves capacity one slot at a time between parallel decompression and
demultiplexing/output compression. Paired R1/R2 decompression shares one
aggregate input budget. The initial allocation is printed for parallel runs.

The budget controls active CPU work; it is not a hard cap on the number of
operating-system threads. Dormant demultiplexing workers provide headroom for
adaptive changes, and parallel gzip decoder workers are created lazily as work
becomes available. Large budgets should therefore be performance-tested for both
throughput and memory use on the target system.
