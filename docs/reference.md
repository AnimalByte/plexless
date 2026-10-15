# Detailed usage reference

Use `plexless demux` for sample assignment and `plexless trim` for
already-demultiplexed FASTQ. This guide covers demux layouts, routing files,
output contracts and advanced settings. For trim and preprocessing options, see
[preprocessing.md](preprocessing.md) and
[fastp-compatibility.md](fastp-compatibility.md).

Demux syntax:

```text
plexless [--threads <N>] demux [OPTIONS] \
  --barcodes <TSV> --samples <TSV> --output <DIR> \
  <--reads <FASTQ> | --r1 <FASTQ> --r2 <FASTQ> | --cram <CRAM>>
```

`--threads` is global and may also appear after `demux`. Show the authoritative
help for the installed version with:

```bash
plexless --help
plexless demux -h
plexless demux --help
plexless trim -h
plexless trim --help
```

### Global options

| Flag | Default | Description |
| --- | ---: | --- |
| `--threads <N>` | `1` | Global work budget shared by input decoding, routing, and output. FASTQ `1` uses the serial reference; CRAM uses no HTSlib background workers at `1`. Values must be positive integers. |
| `-h`, `--help` | — | Print help and exit. |

### Demux input and routing options

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
| `--max-mismatches <N>` | No (`1`) | Maximum substitutions allowed at each routing node. Observed `N` positions consume this budget. |
| `--no-pass-through` | No | Disable default-on distal read-through detection; retain proximal-only trimming. |
| `--pass-through-mismatches <N>` | No (`1`) | Distal substitution budget, independent of barcode correction. The anchor must remain exact. |
| `--r1-distal-structure <LAYOUT>` | No | Override inferred R1 distal geometry in observed orientation. |
| `--r2-distal-structure <LAYOUT>` | No | Override inferred R2 distal geometry in observed orientation (PE only). |

¹ Paired input requires one or both mate structures. A mate without a structure
is retained as biological sequence but supplies no barcode pieces.

The requirements in this table apply to `demux`; `trim` requires neither a
barcode/sample TSV nor a proximal structure. Exactly one input mode is allowed.
Single-end demux runs require `--reads` and
`--structure`. Paired-end runs require `--r1`, `--r2`, and at least one of
`--r1-structure` or `--r2-structure`.

CRAM input requires both `--read-mode` and `--output-format`. Only raw,
unmapped CRAM is accepted. Paired CRAM must declare queryname ordering/grouping
in `@HD` with `SO:queryname` or `GO:query`. See the complete
[CRAM contract](cram.md), including validation and metadata policy.

### Output and performance options

| Flag | Default | Description |
| --- | ---: | --- |
| `-o`, `--output <DIR>` | Required | Output directory. Plexless creates it if absent; an existing directory must be empty. |
| `--compression-level <0-9>` | `2` | Compression level for output FASTQ gzip streams or CRAM files. Lower levels generally favor throughput. |
| `--output-mode <MODE>` | `auto` | `auto`, `direct`, or `buffered`. Auto selects direct below 384 expected streams and buffered at 384 or more. |
| `--output-chunk-size <SIZE>` | `auto` | Buffered-mode target for uncompressed bytes per gzip member. For FASTQ input, auto uses the qualified 32 KiB target. Explicit values must be 32 KiB–64 MiB and no larger than the buffer budget. |
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
[Performance and tuning](performance.md) for memory sizing,
ordering, and file-cache details.

### QC and unassigned-read options

| Flag | Default | Description |
| --- | ---: | --- |
| `--fastq-stats` | Off | Write biological-read statistics to `fastq_stats.tsv` and physical barcode-region statistics to `barcode_stats.tsv`. |
| `--write-unassigned` | Off | Write unmatched, ambiguous, unrouted, and orphan reads untrimmed to `unassigned*.fastq.gz`. |
| `--low-sample-fraction <0-1>` | `0.05` | Mark a populated sample `LOW_REPRESENTATION` below this fraction of the nonzero sample median. `0` disables low-representation classifications; missing-sample warnings remain mandatory. |
| `-h`, `--help` | — | Print `demux` help and exit. |

Every successful demux run writes `sample_metrics.tsv`; `--fastq-stats` adds the
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
  --output-chunk-size 32KiB \
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
A	A01	ACGTTGCACT
A	A02	TGCAAGGTCA
B	B01	GATCCGTAGTC
B	B02	CTAGGACTTGA
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
Sample	A	B
sample_1	A01	B01
sample_2	A01	B02
sample_3	A02	B01
```

The first column must be `Sample`. Remaining columns are the read-structure
barcode symbols in alphabetical order, skipping reserved `T`. Every row must
specify every level and every barcode ID must exist in the whitelist.

These example barcode lengths match `R1_10A11B4T`. The four T bases are
technical sequence and do not appear as a barcode column.

The columns define routing order. The table above compiles to:

```text
A decoder
  -> resolved A child
     -> parent-specific B decoder
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
See [CRAM metadata preservation](cram.md#metadata-preservation-policy).

Assigned reads have their declared structured prefixes removed independently
from each mate, followed by default-on distal pass-through trimming as described
in the [preprocessing guide](preprocessing.md). With `--write-unassigned`, unassigned and orphan records are
written untrimmed to corresponding `unassigned*.fastq.gz` files.
Because those filenames are reserved, a worklist sample that sanitizes to
`unassigned` is rejected when unassigned output is enabled.

`--fastq-stats` preserves two input-level QC reports and adds
`retained_fastq_stats.tsv` for biological reads after trimming and filtering. `fastq_stats.tsv` summarizes
the input suffix after each mate's complete declared prefix, before distal
trimming (these input-level reports retain their existing meaning).
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

Every successful demux run writes `sample_metrics.tsv`, including all worklist
samples—even those with zero reads. It reports retained assigned fragments (after optional filtering), mate-level
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

## Limits and validation

- Maximum logical barcode length: 32 bases
- Barcode symbols: `A-S` and `U-Z`; `T` is reserved
- Mismatch correction: 0 through the actual barcode length, with strict distance validation
- Sample-sheet rows must contain complete paths; partial paths are unsupported
- Ambiguity stops traversal; path-aware rescue is not yet implemented
- Paired resynchronization lookahead: 1,024 records
- Buffered FASTQ-input gzip chunk target: 32 KiB uncompressed
- CRAM-to-FASTQ automatic gzip chunk target: 128 KiB–1 MiB uncompressed
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

- [`docs/architecture.md`](architecture.md) documents compilation and the
  batch-owned, allocation-reduced FASTQ hot path.
- [`docs/cram.md`](cram.md) documents the supported CRAM matrix, validation,
  metadata transformation contract, and limitations.
- [`tests/hierarchical_routing.rs`](../tests/hierarchical_routing.rs) covers local
  namespaces, dense-grid selection, correction safety, reuse, and generic symbols.
- [`tests/scalable_output_qc.rs`](../tests/scalable_output_qc.rs) covers high
  multiplexing, SE/PE ordering, header preservation, and sample QC.

Development checks:

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all
cargo build --release --all-targets
cargo test --release
```
