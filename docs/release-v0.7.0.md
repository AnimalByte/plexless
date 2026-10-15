# Plexless v0.7.0

Plexless now cleans already-demultiplexed FASTQ with `plexless trim` and shares
adapter detection, filtering and QC with demultiplexing. The README starts with
three runnable recipes; `-h` shows common options and `--help` includes all controls.

## Main changes

- Known mate-specific technical sequences, explicit observed-orientation distal
  adapters, long packed signatures and optional conservative partial matching.
- Default paired-end overlap evidence for short-insert adapter detection;
  ambiguous or insufficiently informative overlaps preserve reads.
- Configurable barcode correction budgets with bounded correction indexes and a
  packed-distance fallback, preserving exact precedence and ambiguity safeguards.
- Opt-in quality/N/length filtering and a `--fastp-defaults` profile with reports.
- Quality-window and fixed clipping, polyG/polyX trimming, complexity filters,
  adapter discovery, UMI handling, index blacklists and Phred+64 conversion.
- Optional overlap correction/merging, exact duplicate removal with disk spill,
  bounded duplicate estimates and sampled overrepresentation reports.
- JSON and standalone interactive HTML, interleaved input, stdin/stdout,
  auxiliary outputs and record-safe output splitting.
- Existing routing, paired output ordering, bounded parallel processing and
  compression remain shared between workflows. Native CRAM supports shared
  bounds/filter/QC and exact dedup; FASTQ-only transformations are rejected there.

## Defaults and migration

Barcode and full-adapter mismatch budgets still default to 1. Ordinary runs do
not enable quality/N/length filtering. Sequence-based partial adapter matching is
off unless requested or startup adapter discovery succeeds. Its default evidence
threshold is 20 known bases, matched exactly; masked positions do not count.

Paired-end overlap adapter detection is now on by default. Use `--no-overlap` for
sequence-only adapter matching, or `--require-adapter-for-overlap` to require both
sequence and overlap evidence. Full sequence signatures require 20 known bases
by default; the threshold is configurable. A proximal structure still describes
the beginning of a read. Trim mode removes no proximal bases without an explicit
structure; T7 trimming can use only its real distal adapter sequence.

The fastp-derived features are implemented in Rust, without a fastp runtime
dependency. Command options, report schemas and conservative matching policies
deliberately differ; this is not a drop-in fastp replacement.

## Validation and limits

The local validation passed 200 tests in each of the debug and optimized release
profiles, formatting, all-target/all-feature Clippy with warnings denied, all
three README recipes, and nine sequence/quality
comparisons with official fastp 1.4.0. The Linux packaging workflow additionally
runs all-target/all-feature tests in the optimized release profile.

No equivalent-workload speed improvement is claimed. Library-specific truth
data and sustained production workloads remain to be evaluated. Very short or
repetitive inserts can lack sufficient trimming evidence. Adapter discovery is
heuristic; exact duplicate removal may use temporary disk. See the
[validation report](https://github.com/AnimalByte/plexless/blob/v0.7.0/docs/validation.md)
and [fastp feature guide](https://github.com/AnimalByte/plexless/blob/v0.7.0/docs/fastp-compatibility.md).

## Downloads

The Linux x86-64 GNU archive includes the optimized executable, documentation,
licenses and build provenance. A raw `plexless` executable, `SHA256SUMS`, and an
archive-specific checksum are also attached. The binary is built on Ubuntu 22.04,
requires glibc 2.35 or newer and `libz.so.1`, and is dynamically linked. For other
platforms, build from source. This release does not publish a new crates.io version.

[Full changelog](https://github.com/AnimalByte/plexless/compare/v0.6.0...v0.7.0)
