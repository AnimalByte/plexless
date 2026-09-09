# plexless

`plexless` is a high-throughput hierarchical FASTQ demultiplexer for
single-end and paired-end sequencing data. Read structures describe physical
barcode pieces, while an ordinary tabular sample sheet defines parent-to-child
barcode namespaces and sample routes.

Plexless resolves one routing level at a time. After barcode A resolves, B is
matched only against B barcodes reachable beneath that A; C is then matched
only beneath the resolved A/B path, and so on. Unrelated child namespaces never
participate in matching or mismatch-safety validation.

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
- Optional unassigned output and raw-input FASTQ statistics
- Paired-read ID validation, bounded resynchronization, and orphan handling
- Bounded multithreaded parsing, demultiplexing, compression, and ordered output

## Installation

### Download the prebuilt binary

Plexless is hosted in a private GitHub repository, so downloads require an
account with repository access. Install the
[GitHub CLI](https://cli.github.com/), authenticate, and download the Linux
x86-64 release:

```bash
gh auth login
gh release download v0.2.0 \
  --repo AnimalByte/plexless \
  --pattern 'plexless-v0.2.0-x86_64-unknown-linux-gnu.tar.gz'

tar -xzf plexless-v0.2.0-x86_64-unknown-linux-gnu.tar.gz
mkdir -p "$HOME/.local/bin"
install -m 0755 \
  plexless-v0.2.0-x86_64-unknown-linux-gnu/plexless \
  "$HOME/.local/bin/plexless"
```

Ensure `$HOME/.local/bin` is on `PATH`, then verify the installation:

```bash
plexless --help
```

Users with repository access can also download the archive in a browser from
the [private releases page](https://github.com/AnimalByte/plexless/releases).
The release archive includes a compatibility executable named `simplex` for
existing command lines.

The Linux binaries dynamically use the standard `libz.so.1` runtime library,
which is installed by default on most distributions (`zlib1g` on
Ubuntu/Debian and `zlib` on Fedora/Arch).

### Build from source

A recent stable Rust toolchain and zlib development package are required.

```bash
# Ubuntu/Debian
sudo apt install zlib1g-dev

# Fedora/RHEL
sudo dnf install zlib-devel

# Arch Linux
sudo pacman -S zlib
```

Authenticate to GitHub, clone the private repository, and build:

```bash
gh auth login
gh repo clone AnimalByte/plexless
cd plexless
cargo build --release
```

The primary executable is `target/release/plexless`. To install it for your
user:

```bash
mkdir -p "$HOME/.local/bin"
install -m 0755 target/release/plexless "$HOME/.local/bin/plexless"
```

## Usage

```text
plexless [GLOBAL OPTIONS] demux [DEMUX OPTIONS]
```

Show command help:

```bash
plexless --help
plexless demux --help
```

### Global options

| Flag | Default | Description |
| --- | ---: | --- |
| `--threads <N>` | `1` | Total CPU-work budget shared between FASTQ parsing, gzip input, demultiplexing, and output compression. Must be at least 1. |
| `-h`, `--help` | — | Print help. |

### `demux` options

| Flag | Required | Default | Description |
| --- | :---: | ---: | --- |
| `--reads <FASTQ>` | Single-end | — | Single-end FASTQ or FASTQ.gz input. Cannot be combined with `--r1` or `--r2`. |
| `--r1 <FASTQ>` | Paired-end | — | Paired R1 FASTQ or FASTQ.gz input. Requires `--r2`. |
| `--r2 <FASTQ>` | Paired-end | — | Paired R2 FASTQ or FASTQ.gz input. Requires `--r1`. |
| `--structure <LAYOUT>` | Single-end | — | R1 read structure, for example `R1_10A11B4T`. |
| `--r1-structure <LAYOUT>` | Paired-end¹ | — | R1 read structure. |
| `--r2-structure <LAYOUT>` | Paired-end¹ | — | R2 read structure. |
| `--barcodes <TSV>` | Yes | — | Barcode whitelist with `Set`, `ID`, and `Sequence` columns. |
| `--samples <TSV>` | Yes | — | Hierarchical sample routing table whose first column is `Sample`. |
| `-o`, `--output <DIR>` | Yes | — | New or empty output directory. |
| `--compression-level <0-9>` | No | `2` | Gzip compression level for output FASTQs. |
| `--max-mismatches <0-2>` | No | `1` | Maximum substitutions, including observed `N` positions, allowed at each routing node. |
| `--fastq-stats` | No | Off | Write raw-input FASTQ quality and composition statistics. |
| `--write-unassigned` | No | Off | Write unmatched, ambiguous, unrouted, and orphan reads. |
| `-h`, `--help` | No | — | Print demultiplexing help. |

¹ Paired input requires at least one of `--r1-structure` or `--r2-structure`.

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

Assigned reads have their declared structured prefixes removed independently
from each mate. With `--write-unassigned`, unassigned and orphan records are
written untrimmed to corresponding `unassigned*.fastq.gz` files.

`--fastq-stats` writes `fastq_stats.tsv` from raw input before trimming. It
includes read/base counts, length statistics, GC and N percentages, mean
quality, and Q20/Q30 percentages.

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
routing tree is immutable and shared without hot-path locks. Output order is
deterministic across bounded batches.

Output is always gzip-compressed. The writer keeps at most 64 outputs open;
evicted files are safely reopened as concatenated gzip members.

## Limits and validation

- Maximum logical barcode length: 32 bases
- Barcode symbols: `A-S` and `U-Z`; `T` is reserved
- Mismatch correction: 0, 1, or 2 substitutions
- Sample-sheet rows must contain complete paths; partial paths are unsupported
- Ambiguity stops traversal; path-aware rescue is not yet implemented
- Paired resynchronization lookahead: 1,024 records
- Output writer cache: 64 open writers

Startup validation covers read-structure and `(rc)` syntax, required barcode
sets, barcode bases and logical lengths, duplicate IDs and complete paths,
unknown sample-sheet IDs, root and sibling sequence conflicts, node-local
correction safety with parent-path errors, and the 32-base packed-encoding
limit. Short reads are safely classified as unmatched.

## Architecture and tests

- [`docs/architecture.md`](docs/architecture.md) documents compilation and the
  allocation-free read hot path.
- [`tests/hierarchical_routing.rs`](tests/hierarchical_routing.rs) covers local
  namespaces, correction safety, reuse, three levels, and generic symbols.

Development checks:

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test
cargo build --release --all-targets
```
