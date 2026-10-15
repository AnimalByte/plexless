# Preprocessing change report and validation

The local `AnimalByte/plexless` checkout was fast-forwarded from `f96cba2` to
upstream v0.6.0, `eb622a4`, before implementation. This report covers the
v0.7.0 preprocessing release. `parallel.rs` keeps its existing
orchestration; it was not reorganized as an unrelated cleanup.

## 1. Files changed

The complete manifest is at the end of this report. Main integration points are
`cli`, `demux`, `parallel`, `output`, `writer`, `stats`, `qc`, `routing`, `decoder`,
`pass_through`, and the FASTQ/CRAM input adapters. New modules separate technical
sequence loading, shared preprocessing, trim entry, overlap, quality bounds,
mutable correction/merge operations, duplicate/UMI runtime state, adapter discovery,
detailed distributions, overrepresentation, and reports. Existing tests were
retained and updated for new optional argument defaults and explicit legacy
sequence-only matching where their fixtures require it.

## 2. CLI syntax

`plexless demux` retains its original required barcode/sample inputs.
`plexless trim --reads input.fq.gz --output trimmed` needs neither file; paired
trim accepts `--r1/--r2`. Both commands share the new preprocessing flags.
See [preprocessing.md](preprocessing.md) for original/known-T/filter/partial/SE/PE/
T7 examples and [fastp-compatibility.md](fastp-compatibility.md) for the full feature
matrix, profile defaults, and advanced commands. All three help pages were checked.

## 3. Technical sequences

Strict `Mate<TAB>Sequence` TSV, one optional entry per mate, A/C/G/T normalization,
duplicate/malformed/length errors at startup. Each sequence concatenates that
mate's proximal T segments in observed order. Missing T knowledge remains masked.
There is no whitelist format change, adapter barcode, or larger logical barcode.
Known sequences and adapter masks/chunks are compiled once.

## 4. Orientation

Inferred opposite signatures reverse physical segment order and complement
barcode/technical pieces according to the existing extraction orientation.
Explicit distal sequences are already in their destination mate's observed
orientation. Multiple T pieces, asymmetric mates, `(rc)` pieces, and long technical
strings have orientation tests. Read UMI/skip cycles are masked at the beginning
of inferred opposite distal signatures and add no known-base evidence. Explicit
complete distal configurations retain their supplied physical meaning.

## 5. Full matching

Exact diverse 12–16-base anchor, complete masked verification in <=32-base packed
chunks, independent configurable mismatch budget, conservative observed-N handling,
exact whole-distal boundary. Full signatures require 20 known positions by default,
configurable through `--min-full-adapter-known-bases`. Unknown expected positions
never count as evidence. Full matching takes precedence over partial matching.

## 6. Partial matching

Optional bounded suffix/prefix Shift-And matching uses startup masks and reused
scratch, without requiring the full anchor to appear in the sequenced prefix.
Known positions match exactly; full mismatch tolerance does not weaken partial
specificity. Candidates require base/dinucleotide diversity and reject low
complexity. Across adapters, maximize matching known evidence, then favor less
clipping. No arbitrary internal adapter fragment is eligible.

## 7. Minimum overlap and default mate evidence

Sequence-based partial matching is off by default. Its default minimum is 20
known bases: 20/30 known bases can qualify, 8 cannot, and masked positions do not
count. Default PE overlap separately supports boundaries without a complete
adapter: unique informative alignment, minimum 30 biological bases, <=5 mismatches
and <=20% differences across the entire overlap. Ambiguity preserves reads.
`--no-overlap` selects sequence-only adapter matching; `--require-adapter-for-overlap`
requires both forms of evidence. Optional one-gap fallback is restricted to adapter
bounds; merge/correction remain ungapped.

## 8. Barcode mismatch scaling

Independent `--max-mismatches` and `--pass-through-mismatches`, both default 1.
Barcode budgets allow 0 through the actual <=32-base barcode length. Route-local
strict distance `2*k+1` validation and exact precedence remain. Observed Ns consume
budget; ties never silently assign. Neighborhood estimates above one million
entries/node use a packed scan rather than exponential enumeration. High budgets
can cost more per nonexact observation; pairwise startup safety checks remain
quadratic in the route-local catalog size.

## 9. Trim architecture

Trim constructs one output descriptor, omits barcode extraction/decoding/routing,
and reuses the existing bounded FASTQ parser, batches, workers, compression and
ordered writers, including at one thread. Proximal clipping requires an explicit
structure. Distal sequence-only T7 trimming works without proximal/barcode inputs.
There is no independent slow trim reader or separate adapter implementation.

## 10. Shared filtering, QC and added features

Both commands use shared bounds/filter/QC logic. Quality/N/length filters run on
retained sequence; default runs add no filtering. Empty reads have length/mean
quality/low-quality percentage zero; enabled thresholds determine retention.
Reports include processed/retained/rejected fragments, overlapping rejection
reasons, mate-specific full/partial/overlap trimming and bases, proximal/quality/
polyG/polyX removal, dimers, correction and merging, plus before/after biological QC.
Detailed JSON and standalone interactive HTML are opt-in or profile-enabled.

Additional fastp-derived families are fixed/window quality clipping, polyG/polyX,
complexity filtering, Phred+64 conversion, UMI/header handling, index blacklists,
startup adapter inference, correction, merging, exact duplicate removal with
bounded-memory disk spill, bounded duplicate estimates and overrepresentation,
stdin/interleaving/stdout, auxiliary outputs, processing limits and record-safe
split files. The feature guide records deliberate behavior differences and costs.

## 11. Paired behavior

Original barcode regions are assigned first. Default mate overlap verifies the
whole alignment and accounts for asymmetric fronts. Mates can retain unequal
lengths. Without merging, either failed mate rejects the primary pair; requested
auxiliary streams can retain passing mates/failing reads without changing paired
counts. Merge emits a consensus fragment, or retains unmerged pairs according to
configuration. Per-sample merging and sample-QC reconciliation have tests. Asymmetric
front clipping preserves biological overhangs covered by the other mate. Sequence/
quality coordinates, gzip validity and input order remain paired and ordered.

## 12. Validation results

**200 tests passed in each of the debug and optimized release profiles, with
zero failures or ignored tests, across 24 targets per profile.**
Formatting and the formatting check passed. Clippy passed all targets/features
with `-D warnings`. The root, demux and trim help pages succeeded and were
inspected. All **nine official fastp 1.4.0 sequence/quality comparisons passed**.
No speed or production-workload equivalence is claimed.

`cargo install --path . --locked --force` completed in the optimized release
profile. The optimized executable exposes `demux` and `trim` and passed all nine
fastp comparisons again. The release workflow also runs the complete test suite
in the release profile before packaging the Linux executable, documentation,
licenses and SHA-256 checksums. See the [release notes](release-v0.7.0.md) and
[installation guide](installation.md) for the downloadable build's requirements.

```bash
cargo fmt --all
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-targets --all-features
cargo test --locked --release --all-targets --all-features
./target/debug/plexless --help
./target/debug/plexless demux --help
./target/debug/plexless trim --help
python3 scripts/compare_fastp.py --fastp /path/to/fastp \
  --plexless target/debug/plexless
```

This machine's Clang installation lacks resource headers. Compile/check/test
commands used the environment-only workaround below; no machine-specific include
path was added to the repository:

```bash
BINDGEN_EXTRA_CLANG_ARGS="-I/usr/lib/gcc/x86_64-linux-gnu/15/include" \
  cargo test --all-targets --all-features
```

The official fastp 1.4.0 comparison covers nine synthetic transformations with
matching output sequences/qualities: fixed clipping, front/tail/right windows,
polyG, polyX, UMI, merging, and short-insert overlap trimming. Their golden outputs
are checked by Cargo tests offline. Header/report formats intentionally differ.
Other tests cover technical TSV/orientation, 42/64/100/128+ adapters, >255 full
mismatch counts, genomic/random negatives, default and one-gap overlaps, short/
ambiguous/invalid/late-mismatch cases, barcode budgets 0/1/2/higher, route-local
safety and bounded fallback, demux assignment before transformations, filtering,
auxiliary and per-sample merged output, exact dedup determinism/spill, serial/
parallel FASTQ/CRAM, report escaping, stdin/interleaving, splitting, output-path
collisions/overwrite protection and cancellation of saturated pipelines.

## 13. Performance considerations

Packed barcodes and exact lookup, full-signature anchor scanning, chunked Hamming
verification, precompiled suffix masks, bounded queues, parallel gzip/processing/
compression, ordered aggregation and output caches remain. Overlap scratch is
reused, and unchanged biological bounds reuse the preceding overlap result for
merge/correction. Mutation/header/auxiliary features acquire owned records only
when needed. Split output is divided before compression without recompression.

Optional detailed reports and heavy-hitter analysis scan additional bases.
Exact dedup adds producer work and temporary disk proportional to unique original
fragments after memory fills; reporting-only Bloom estimates cannot discard reads.
`--split` adds an input counting pass. Catalog/cache settings bound their own
state, not total process RSS. No throughput improvement is claimed.

## 14. Remaining limits and deliberate differences

This is fastp feature-family coverage rather than drop-in command, report-schema,
or bit-for-bit algorithm compatibility. No manufacturer adapters are built in.
Full and partial sequence matching still reject indels; only optional mate adapter
overlap permits one internal gap. Unique/diverse whole-overlap evidence is stricter
than inspected fastp prefix/first-hit shortcuts, so very short, low-complexity or
repeated inserts can retain contamination. `--no-overlap` provides sequence-only
behavior for such libraries when appropriate. Original resynchronization orphans
retain their legacy routing rather than entering new fragment-level filters.

Duplicate estimates and overrepresentation are bounded/sampled; exact duplicate
removal stores original sequences and preserves the first occurrence even if that
occurrence later fails a filter. Temporary disk can limit exact removal. Adapter
discovery is bounded and heuristic, disabled for stdin. Custom primary paths require
one output sample; multisample output retains per-sample filenames. `--split` is
an input-based target and can yield fewer/more parts after fanout/merging; the line
limit is the strict per-file record control. Stdout mirrors retained file output.
FASTQ-only transformations are explicitly rejected for native CRAM. Library-specific
truth data and sustained production workloads remain to be evaluated.

## 15. Recommended fastp benchmarks

First establish equivalent retained IDs, sequences, qualities, pair counts and
rejection metrics. Match adapters, evidence/mismatch policies, inference, filters,
quality encoding, duplication mode, compression, storage and total CPU budget.
Compare SE; ordinary PE; short-insert PE; long known/unknown T structures; high-N/
low-quality and repetitive negatives; many-sample fanout; and duplicate-rich data
with/without disk spill. Record versions/commands, wall/CPU time, peak RSS, temporary
disk, output bytes and accuracy at true biological boundaries. Separate startup
sampling/counting from steady-state processing. The nine correctness comparisons
are not equivalent-workload throughput benchmarks.

## Usability update

The README now starts with three runnable recipes: demux, paired FASTQ cleanup,
and already-demultiplexed T7 trimming. It is 181 lines, down from roughly 750.
A defaults/profile table distinguishes routine adapter cleanup from optional
filtering, clipping, merging, UMI and duplicate operations. Installation, detailed
usage and performance material moved to their own guides without removing the
advanced controls.

`plexless demux -h` and `plexless trim -h` show common options and examples.
`--help` lists every option under task headings, including advanced settings.
Short parameter labels and plain descriptions explain read inputs, 3-prime
adapters, filters and reports. Existing option names, parsing and processing
settings are preserved.

All **200 tests**, strict Clippy and formatting checks passed after the help
changes. Each of the three README recipes was executed unchanged against
synthetic gzip FASTQ, using the README's TSVs and checking exact output sequences,
qualities and report files. Local Markdown paths and heading anchors were also
checked across all 14 Markdown files. The optimized installed executable was
updated and all three recipes passed again through the installed command. Its
common-help views contain 55 lines for demux and 47 for trim; every advanced
option remains available in grouped `--help`.

## Complete changed-file manifest

- [`Cargo.lock`](../Cargo.lock)
- [`Cargo.toml`](../Cargo.toml)
- [`.github/workflows/release.yml`](../.github/workflows/release.yml)
- [`README.md`](../README.md)
- [`THIRD_PARTY.md`](../THIRD_PARTY.md)
- [`docs/architecture.md`](../docs/architecture.md)
- [`docs/fastp-compatibility.md`](../docs/fastp-compatibility.md)
- [`docs/installation.md`](../docs/installation.md)
- [`docs/performance.md`](../docs/performance.md)
- [`docs/preprocessing.md`](../docs/preprocessing.md)
- [`docs/reference.md`](../docs/reference.md)
- [`docs/release-v0.7.0.md`](../docs/release-v0.7.0.md)
- [`docs/validation.md`](../docs/validation.md)
- [`licenses/fastp-MIT.txt`](../licenses/fastp-MIT.txt)
- [`scripts/compare_fastp.py`](../scripts/compare_fastp.py)
- [`src/adapter_detection.rs`](../src/adapter_detection.rs)
- [`src/cli.rs`](../src/cli.rs)
- [`src/cram.rs`](../src/cram.rs)
- [`src/decoder.rs`](../src/decoder.rs)
- [`src/demux.rs`](../src/demux.rs)
- [`src/distribution.rs`](../src/distribution.rs)
- [`src/fastp_options.rs`](../src/fastp_options.rs)
- [`src/fastp_runtime.rs`](../src/fastp_runtime.rs)
- [`src/fastq.rs`](../src/fastq.rs)
- [`src/input.rs`](../src/input.rs)
- [`src/lib.rs`](../src/lib.rs)
- [`src/main.rs`](../src/main.rs)
- [`src/output.rs`](../src/output.rs)
- [`src/overlap.rs`](../src/overlap.rs)
- [`src/overrepresentation.rs`](../src/overrepresentation.rs)
- [`src/parallel.rs`](../src/parallel.rs)
- [`src/parallel_input.rs`](../src/parallel_input.rs)
- [`src/pass_through.rs`](../src/pass_through.rs)
- [`src/preprocessing.rs`](../src/preprocessing.rs)
- [`src/qc.rs`](../src/qc.rs)
- [`src/quality_trim.rs`](../src/quality_trim.rs)
- [`src/read_ops.rs`](../src/read_ops.rs)
- [`src/report.rs`](../src/report.rs)
- [`src/report_template.html`](../src/report_template.html)
- [`src/routing.rs`](../src/routing.rs)
- [`src/stats.rs`](../src/stats.rs)
- [`src/technical.rs`](../src/technical.rs)
- [`src/trim.rs`](../src/trim.rs)
- [`src/writer.rs`](../src/writer.rs)
- [`tests/barcode_behavior.rs`](../tests/barcode_behavior.rs)
- [`tests/cram_support.rs`](../tests/cram_support.rs)
- [`tests/fastp_features.rs`](../tests/fastp_features.rs)
- [`tests/fixtures/README.md`](../tests/fixtures/README.md)
- [`tests/fixtures/ecoli_U00096_3_1_10000.fa`](../tests/fixtures/ecoli_U00096_3_1_10000.fa)
- [`tests/fixtures/fastp_v1_4_0.json`](../tests/fixtures/fastp_v1_4_0.json)
- [`tests/fixtures/human_NC_000022_11_20000001_20010000.fa`](../tests/fixtures/human_NC_000022_11_20000001_20010000.fa)
- [`tests/hierarchical_routing.rs`](../tests/hierarchical_routing.rs)
- [`tests/orientation_demux.rs`](../tests/orientation_demux.rs)
- [`tests/output_modes.rs`](../tests/output_modes.rs)
- [`tests/overlap_preprocessing.rs`](../tests/overlap_preprocessing.rs)
- [`tests/paired_end.rs`](../tests/paired_end.rs)
- [`tests/paired_errors.rs`](../tests/paired_errors.rs)
- [`tests/parallel_gzip_input.rs`](../tests/parallel_gzip_input.rs)
- [`tests/parallel_paired.rs`](../tests/parallel_paired.rs)
- [`tests/parallel_single_end.rs`](../tests/parallel_single_end.rs)
- [`tests/pass_through.rs`](../tests/pass_through.rs)
- [`tests/preprocessing.rs`](../tests/preprocessing.rs)
- [`tests/scalable_output_qc.rs`](../tests/scalable_output_qc.rs)
- [`tests/single_end.rs`](../tests/single_end.rs)
- [`tests/threads_cli.rs`](../tests/threads_cli.rs)
