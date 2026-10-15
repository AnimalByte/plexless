# Plexless

Demultiplex sequencing reads into samples, or trim and filter already-demultiplexed
FASTQ. Plexless supports single-end and paired-end FASTQ, hierarchical barcode
routing, and unmapped CRAM demultiplexing.

| Your input | Use | Required configuration |
| --- | --- | --- |
| Reads containing sample barcodes | `plexless demux` | Read layout, barcode TSV, sample TSV |
| FASTQ already assigned to a sample | `plexless trim` | Input and output paths; adapter sequence for SE adapter trimming |

## Quick start

### 1. Demultiplex reads into samples

```bash
plexless demux --reads reads.fq.gz --structure R1_10A11B4T \
  --barcodes barcodes.tsv --samples samples.tsv \
  --output demux --threads 8
```

`R1_10A11B4T` means **10 barcode-A bases, 11 barcode-B bases, 4 technical bases,
then biological sequence**. Adapt the layout to your library. Plexless assigns the
sample from the original barcode regions, removes that prefix, and checks for
adapter read-through. Filtering is off unless requested.

Outputs: `demux/<sample>.fastq.gz`, `sample_metrics.tsv`, and `preprocessing.tsv`.
For paired input, use `--r1`, `--r2`, and mate-specific structures; see the
[paired demux example](docs/reference.md#paired-end-example-with-orientation-normalization).

### 2. Clean paired-end FASTQ with filtering and reports

```bash
plexless trim --r1 R1.fq.gz --r2 R2.fq.gz \
  --fastp-defaults --output cleaned --threads 8
```

Plexless uses a unique, informative mate overlap to identify read-through.
`--fastp-defaults` adds quality/N/length filters, duplicate estimation, and reports.
It also enables automatic polyG trimming for recognized instrument IDs in files.

Outputs: `cleaned/trimmed_R1.fastq.gz`, `trimmed_R2.fastq.gz`, QC tables,
`plexless.json`, and **`plexless.html`**. Open the HTML file to inspect the results.
If either mate fails a filter, the primary output rejects the pair.

Omit `--fastp-defaults` for adapter cleanup without the profile's filters or reports.
Add `--fastq-stats` for QC tables without enabling filtering.

### 3. Trim an already-demultiplexed T7 FASTQ

Put your library's complete distal adapter sequence in `r1_distal.txt`, in the
orientation it appears at the read's 3-prime end. Use only A/C/G/T bases.

```bash
R1_DISTAL=$(tr -d '\r\n' < r1_distal.txt)
plexless trim --reads t7.fq.gz --r1-distal-sequence "$R1_DISTAL" \
  --partial-pass-through --fastq-stats --output t7-trimmed --threads 8
```

This keeps the read's 5-prime bases and removes supported distal adapter matches.
Partial matching requires at least **20 informative known adapter bases** by
default. No manufacturer adapter sequences are built in.

A **structure describes the beginning of a read**. A **distal sequence describes
its 3-prime adapter**. Supply `--structure` in trim mode only when you intend to
remove a prefix. Barcode and sample files are unnecessary for trim.

## Defaults and optional features

| Behavior | Ordinary run | With `--fastp-defaults` |
| --- | --- | --- |
| Paired overlap adapter detection | On | On |
| Proximal clipping in trim | Only with an explicit structure | Same |
| Quality/N/length filters | Off | Q15; at most 40% low-quality bases; at most 5 Ns; minimum length 15 |
| Automatic adapter discovery | Off | On for SE; PE requires `--detect-adapter-for-pe` |
| Automatic polyG trimming | Off | Instrument-based detection for file input |
| Duplicate estimation | Off | On |
| JSON/HTML reports | Only when requested | On, in the output directory |
| Quality-window clipping, merging, correction, UMI extraction, duplicate removal | Only when requested | Same |

`--fastp-defaults` enables the profile settings listed above. Other transformations
require their own flags. Explicit filter thresholds override the profile defaults.

Full sequence adapter matching is enabled when a usable signature is configured
or inferred from the demux layout. Sequence-based partial matching requires
`--partial-pass-through`; successful startup adapter discovery also enables it for
that run. Paired overlap can identify a boundary independently of sequence-based
partial matching. `--no-overlap` selects sequence-only adapter detection.
See [adapter matching and orientation](docs/preprocessing.md) for the evidence
rules and [all feature defaults](docs/fastp-compatibility.md) for advanced controls.

## Common additions

Append only the options you need; both commands share these controls.

| Goal | Add |
| --- | --- |
| QC tables without filtering | `--fastq-stats` |
| Standalone reports | `--json report.json --html report.html` |
| Keep reads at least 30 bases long | `--min-length 30` |
| Reject reads with more than 5 Ns | `--max-n 5` |
| Cut a low-quality tail using a quality window | `--cut-tail` |
| Trim partial adapters with more evidence | `--partial-pass-through --min-partial-adapter-overlap 30` |
| Supply known technical T bases | `--technical-sequences technical.tsv` |
| Merge uniquely overlapping pairs | `--merge` |
| Extract an 8-base UMI from each read | `--umi --umi-loc per-read --umi-len 8` |
| Remove exact duplicate original fragments | `--dedup` |

Quality clipping changes the retained sequence. Quality filtering accepts or
rejects the resulting read/pair. Merging, UMI extraction, and duplicate removal
are separate optional steps; their [ordering and behavior](docs/fastp-compatibility.md#workflow-and-pairing)
are documented in the feature guide.

## Demux configuration files

The quick-start layout above uses these example files. They are tab-separated;
replace the synthetic sequences and sample routes with your library's values.

`barcodes.tsv`:

```tsv
Set	ID	Sequence
A	A01	ACGTTGCACT
A	A02	TGCAAGGTCA
B	B01	GATCCGTAGTC
B	B02	CTAGGACTTGA
```

`samples.tsv`:

```tsv
Sample	A	B
sample_1	A01	B01
sample_2	A01	B02
sample_3	A02	B01
```

Each barcode set's sequence length must match its logical length in the layout.
Correction defaults to one mismatch and rejects unsafe or ambiguous assignments.
A resolves first, then B is matched within that A's sample routes.
[Layouts, orientation, and routing details](docs/reference.md#input-configuration).

## Install and get help

From the current source checkout, with Rust and the native build dependencies
installed:

```bash
cargo install --path . --locked --force
plexless --version
```

See [installation](docs/installation.md) for system prerequisites, source builds,
and published packages/binaries.

```bash
plexless --help       # Choose demux or trim
plexless demux -h     # Common demux options and an example
plexless trim -h      # Common trim options and examples
plexless trim --help  # Every trim option, grouped by purpose
plexless demux --help # Every demux option, including CRAM and tuning
```

Use a new or empty output directory. Existing outputs are protected from
overwrite. A `PLEXLESS_INCOMPLETE` marker means a run did not finish successfully.
The default output is gzip-compressed FASTQ; the [output reference](docs/reference.md#output-behavior)
explains filenames, reports, and optional output modes.

## Detailed guides

- [Usage reference](docs/reference.md): paired demux, structures, barcode/sample files, routing, outputs, and limits.
- [Adapters and preprocessing](docs/preprocessing.md): known T sequences, orientation, full/partial matching, filtering, and T7 examples.
- [Feature guide](docs/fastp-compatibility.md): quality trimming, UMIs, duplicates, merging, reports, streaming, and splitting.
- [CRAM support](docs/cram.md): accepted inputs and metadata behavior.
- [Performance and tuning](docs/performance.md): historical measurements, CPU/memory settings, and compression.
- [Architecture](docs/architecture.md) and [validation report](docs/validation.md): implementation and tests. Current preprocessing changes have correctness comparisons; throughput comparisons remain to be run.

## License

MIT or Apache-2.0, at your option. See [LICENSE-MIT](LICENSE-MIT),
[LICENSE-APACHE](LICENSE-APACHE), and [third-party notices](THIRD_PARTY.md).
