# Adapter trimming, filtering, and trim-only FASTQ

`plexless demux` assigns barcodes from the original read, clips the declared
proximal prefix, detects distal adapters, optionally filters, and writes ordered
sample outputs. `plexless trim` skips barcode extraction, correction, and sample
routing entirely. It uses the same preprocessing, batched FASTQ parser, bounded
queues, gzip input, compression, and ordered writer as demux, with one output
stream per mate. No whitelist or sample file is opened in trim mode.

## Proximal and distal configuration

A structure such as `R1_10A33T` describes the **beginning** of a read: a logical
barcode A and 33 technical bases. The rest is biological sequence. It is not a
3-prime adapter declaration. Trim mode clips no proximal bases unless you
explicitly supply `--structure` (SE) or `--r1-structure`/`--r2-structure` (PE).
An already-demultiplexed T7 FASTQ needs only its actual distal sequence.

`--r1-distal-sequence ACGT...` and `--r2-distal-sequence ACGT...` specify a complete
expected distal sequence in the **observed orientation of the destination read**.
They accept nonempty A/C/G/T strings, case-insensitively, with no barcode-length
limit. They do not trigger proximal clipping. There are no built-in manufacturer
sequences. Supply the complete sequence appropriate to your library.

For each mate, an explicit distal sequence or `--r1-distal-structure` /
`--r2-distal-structure` replaces inferred distal configuration. Sequence and
structure options for the same mate conflict; configuration for R2 requires
paired input. Distal options and partial matching conflict with
`--no-pass-through`. Distal structures must name their destination mate.

Demux retains its existing inference: R1 expects the reverse complement of the
physical R2 prefix, and R2 expects the reverse complement of the physical R1
prefix. SE inference assumes symmetric ends and reverses its own prefix.
Physical segment order reverses; each barcode piece's declared `(rc)` orientation
is inverted. Explicit distal structures are already in observed order: no
whole-structure reversal is applied, and `(rc)` applies to that barcode piece.
Each referenced logical barcode must be declared proximally and consumed in full.

Trim cannot resolve barcode-dependent distal structures; use explicit distal
sequences instead. For a T-only distal structure, known technical bases from the
same mate must have matching length and are used in observed order. For trim's
optional inference, only a wholly known T-only opposite proximal prefix can
supply a signature; proximal barcode symbols never become fake assignments.

## Known technical sequences

`--technical-sequences technical.tsv` loads a strict tab-separated file once:

```tsv
Mate	Sequence
R1	ACGT
R2	TGCA
```

These short sequences are illustrative synthetic examples, not manufacturer
adapters. Use only the mates present in your layout. Header names are exact.
Each mate can occur once. Sequences are normalized to uppercase and must be
nonempty A/C/G/T; invalid bases, extra columns, blank records, files with no entries,
duplicate mates, and incorrect lengths produce startup errors. A file may
specify just one mate; unspecified mates remain masked.

Each entry concatenates **all proximal T segments for that mate, in read order**.
For `R1_2T10A3T`, the five technical bases map first to the 2T piece and then the
3T piece. They are supplied as seen in R1, independently of barcode orientation.
When compiling an opposite-mate signature, Plexless reverses the segment order
and reverse-complements each technical piece. No barcode TSV changes are needed;
logical barcodes retain their 32-base limit. Technical sequences can exceed 100
bases and signatures use multiple packed verification chunks.

For explicit distal structure overrides, T pieces consume the destination mate's
known technical string in occurrence order, without reversing it; their total T
length must match that string. If the intended distal technical order or length
differs, specify the complete observed distal sequence instead. Without a known
technical entry, T positions remain expected unknowns and are masked.

## Default paired overlap

Paired reads use unique biological mate overlap by default, with at least 30
bases and both absolute (5) and percentage (20%) mismatch limits across the
complete alignment. A unique informative overlap can establish a boundary even
without adapter sequence evidence. Ambiguous or low-complexity alignments retain
the reads. `--require-adapter-for-overlap` requires both forms of evidence;
`--no-overlap` restores independent sequence-only matching. These policies also
apply to demux after assignment. Optional one-gap adapter overlap, correction,
merging, quality trimming, UMI, deduplication, reports, and all output extensions
are detailed in the [fastp feature guide](fastp-compatibility.md).

## Full adapter detection

Full detection always runs first. Plexless scans for a startup-selected exact
12–16-base anchor, then verifies every known position in the complete expected
sequence using packed Hamming chunks of at most 32 bases. Full verification
requires at least `--min-full-adapter-known-bases` known positions (default 20);
masked T positions do not count. The independent
`--pass-through-mismatches N` budget (default 1) applies across all chunks; it
never relaxes the exact anchor. There is no hardcoded two-substitution ceiling.
The counter uses `usize`, including for signatures longer than 255 bases.

The full sequence must fit. Failed candidates do not hide later valid full
matches; the earliest verified boundary wins in sequence-only mode. Default
paired mode requires agreement with the unique overlap boundary and can find a
later full candidate agreeing with that boundary. Trimming starts at the complete
distal boundary, including unknown technical positions before the anchor, and
removes everything following it. Unknown observed bases at known expected
positions reject a candidate even with a large mismatch budget. Unknown expected
T positions accept any observed base and add no mismatch penalty.

Anchors retain the existing complexity requirements: at least three base types,
at least six distinct dinucleotides, and no base above 75%. Ends lacking an
anchor are unavailable for full matching. With partial matching disabled,
sequence-only mode preserves incomplete adapters and reads without a verified
full signature. Default paired overlap is a separate evidence source.

## Conservative partial read-through

Enable with `--partial-pass-through`. It is **off by default**. Set
`--min-partial-adapter-overlap N` to control evidence (default **20 known bases**).
A 30- or 20-known-base adapter prefix can qualify; an eight-base prefix cannot
qualify with the default. Masks do not count toward evidence, so 30 observed
bases with 11 masked positions are insufficient.

Partial detection runs only if full detection finds no complete signature.
It compares a **3-prime read suffix** to the **prefix** of the expected adapter.
Internal adapter fragments, adapter suffixes, and adapter-like internal matches
with unrelated sequence after them do not qualify. Complete full matches always
win, even if a partial match elsewhere has a longer overlap.

All known partial positions must match **exactly**. Neither the full-signature
mismatch budget nor barcode correction weakens this rule. Observed N or invalid
bases at known positions reject that candidate. Expected unknown T positions are
masked. Candidate prefixes also require three base types, six distinct
contiguous known dinucleotides, and no base above 75% of known evidence. Masked
positions break dinucleotides. Short, homopolymeric, and repetitive evidence is
insufficient. Failure preserves the sequence-based candidate. Paired overlap
can independently supply a boundary even when partial sequence matching is off.

If multiple suffix/prefix candidates qualify, maximize the number of matching
known bases; for equal evidence choose the shortest overlap (least clipping).
This deterministic rule uses only sequenced evidence. It does not require the
complete adapter or the full-match anchor to have been sequenced.

The matcher compiles Shift-And masks and prefix evidence at startup. It scans
only the last `min(read length, adapter length - 1)` bases with a bit-parallel
prefix automaton, propagating states across 64-bit words. This avoids comparing
full strings at every read position and works with 42-, 64-, 100-, and 128-base
adapters. Up to 256 observed adapter bases use stack storage with no per-read
allocation; longer matches reuse scratch storage within a worker batch. Runtime
is bounded by the configured adapter length and observed suffix length.

These safeguards reduce false trimming but cannot distinguish biological DNA
identical to a qualifying adapter prefix. Validate against library-specific
negative reads before production rollout. The included randomized negatives and fixed human/E. coli reference windows
are regression tests, not a library-specific empirical false-positive estimate.

## Barcode mismatch budgets

`--max-mismatches N` defaults to 1 and accepts 0 through the actual logical barcode
length (at most 32), without an arbitrary ceiling of 2. Each route-local sibling
namespace must retain minimum pairwise Hamming distance `2*N + 1`. Unsafe
geometry is rejected at startup with barcode names, route, and measured distance.
Reusing sequences under unrelated parents remains supported. Exact lookup takes
precedence; observed Ns each consume one mismatch, and runtime ties never assign.

Plexless estimates neighborhood size before building a correction index. Nodes
whose estimated total exceeds 1,000,000 correction entries use packed-distance
scans over their route-local candidate list instead. Exact lookup stays O(1).
The fallback uses O(candidate count) storage and O(candidate count) work for a
nonexact observation; it does not enumerate exponential neighborhoods. Existing
pairwise geometry validation is O(candidate count squared) at startup. Large
catalogs with high budgets can therefore be slower; no throughput improvement is
claimed. The entry cap is a resource policy, not a mismatch threshold.

## Filtering and paired behavior

Filtering is opt-in for both commands and runs after proximal and distal clipping:

| Option | Rejection condition |
| --- | --- |
| `--qualified-quality Q` | Sets the low-quality definition: Phred **below** Q (default 15 when percentage filtering is enabled); alone discards nothing. |
| `--max-low-quality-percent P` | Low-quality fraction times 100 **exceeds** P. |
| `--max-n COUNT` | Retained N/n count **exceeds** COUNT. |
| `--min-length BASES` | Retained length is **below** BASES. |
| `--min-average-quality Q` | Retained mean Phred is **below** Q. |

Qualities are Phred+33. Percentages accept 0–100 and quality thresholds 0–93.
Exact thresholds pass. Empty reads have mean quality 0, N count 0, low-quality
percentage 0, and length 0. They remain unless an enabled filter rejects them.
Quality/window clipping, polyG/polyX, correction and merging are separately
opt-in; see the [fastp feature guide](fastp-compatibility.md).

Mates trim independently and can retain different lengths. If either mate fails,
the entire pair is rejected. Output order and sequence/quality correspondence are
preserved. Existing demux assignment happens first, using original barcode regions;
unassigned reads keep their original bounds and can also be filtered. Existing
orphan routing remains unchanged; resynchronization orphans do not enter
fragment-level filters or post-filter QC. Optional unpaired outputs capture the
passing mate of a pair rejected by an enabled filter.

## Reports

Every run writes `preprocessing.tsv`: processed, retained, rejected fragments
(SE reads or PE pairs), rejection reasons, and R1/R2 full/partial detection counts,
full/partial bases removed, and proximal bases removed. Reason counts can overlap,
while a rejected pair counts once. Detection/removal counts include subsequent
filter rejections. Orphans are reported separately by the existing summary.
Retained fragment counts describe filter decisions, including unassigned reads
whose output is disabled. Sample output counts describe emitted records.

`--fastq-stats` never enables filtering. Existing `fastq_stats.tsv` keeps its input
biological meaning: after the declared proximal region, before distal clipping or
filtering, across all input records. `barcode_stats.tsv` keeps raw barcode-piece
QC. New `retained_fastq_stats.tsv` describes biological reads after trimming and
filtering, excluding orphans, with reads/bases, length, GC%, N%, mean quality,
Q20%, and Q30%. Filter measurements and post-filter QC share the same base scan.
The input QC scan remains separate because it measures different coordinates.

## Commands

Original demux (existing barcode/sample TSVs are unchanged):

```bash
plexless demux --reads reads.fastq.gz --structure R1_10A11B4T \
  --barcodes barcodes.tsv --samples samples.tsv --output demux-out --threads 8

# technical.tsv contains R1's four proximal T bases in observed orientation
plexless demux --reads reads.fastq.gz --structure R1_10A11B4T \
  --barcodes barcodes.tsv --samples samples.tsv --technical-sequences technical.tsv \
  --output demux-known --threads 8

plexless demux --reads reads.fastq.gz --structure R1_10A11B4T \
  --barcodes barcodes.tsv --samples samples.tsv --output demux-filtered \
  --min-length 30 --max-n 2 --qualified-quality 15 --max-low-quality-percent 20 \
  --min-average-quality 20 --fastq-stats --threads 8

plexless demux --reads reads.fastq.gz --structure R1_10A11B4T \
  --barcodes barcodes.tsv --samples samples.tsv --technical-sequences technical.tsv \
  --partial-pass-through --min-partial-adapter-overlap 20 --output demux-partial
```

Read your library's actual complete distal sequences from text files containing
only their observed A/C/G/T bases (line endings are removed below):

```bash
R1_DISTAL=$(tr -d '\r\n' < r1_distal.txt)
R2_DISTAL=$(tr -d '\r\n' < r2_distal.txt)

# Trim-only SE: output trim-se/trimmed.fastq.gz
plexless trim --reads sample.fastq.gz --r1-distal-sequence "$R1_DISTAL" \
  --output trim-se --threads 8 --fastq-stats

# Trim-only PE: trim-pe/trimmed_R1.fastq.gz and trimmed_R2.fastq.gz
plexless trim --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz \
  --r1-distal-sequence "$R1_DISTAL" --r2-distal-sequence "$R2_DISTAL" \
  --partial-pass-through --min-length 30 --output trim-pe --threads 8

# Already-demultiplexed T7 FASTQ: no proximal clipping is implied
plexless trim --reads t7_sample.fastq.gz --r1-distal-sequence "$R1_DISTAL" \
  --output t7-full --pass-through-mismatches 1 --fastq-stats --threads 8
plexless trim --reads t7_sample.fastq.gz --r1-distal-sequence "$R1_DISTAL" \
  --partial-pass-through --min-partial-adapter-overlap 20 \
  --output t7-partial --fastq-stats --threads 8
```

## Performance validation against fastp

Use the same compressed input, library-specific explicit adapters, output
compression level, storage, and thread budget. Compare full-only and partial
read-through separately, across SE and PE, varied inserts, adapter lengths
42/64/100/128, GC/repetitiveness, unknown bases, and output fanout. Include
negative genomic controls and known truth boundaries, and evaluate false trimming,
missed adapters, exact retained sequences/qualities, paired counts, wall time,
CPU time, peak RSS, and output bytes. Warm up tools and repeat runs, retaining
versions and exact commands. Test filtering disabled first, then equivalent
filter thresholds and QC. Do not compare demux fanout with trim-only work or
claim equivalent behavior without checking overlap specificity and ambiguity policies. Different
partial-match specificity requires a correctness comparison before throughput
comparisons. No equivalent-workload fastp benchmark has been run for this change.


## Validation of this implementation

The local checkout was fast-forwarded from `f96cba2` to upstream v0.6.0
`eb622a4` before the v0.7.0 preprocessing work.

The required checks cover all targets and features; the final local validation
results are recorded in [validation.md](validation.md). Existing tests were retained; Rust callers were updated to
supply defaults for the new optional preprocessing fields and handle the new CLI
subcommand variant. The three help pages are included in validation.

```bash
cargo fmt --all
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-targets --all-features
./target/debug/plexless --help
./target/debug/plexless demux --help
./target/debug/plexless trim --help
```

This machine's Clang installation lacks resource headers, causing HTSlib bindgen
to report missing `stddef.h`. Clippy and tests were run with the existing GCC
header directory provided as a command environment override:

```bash
BINDGEN_EXTRA_CLANG_ARGS="-I/usr/lib/gcc/x86_64-linux-gnu/15/include" \
  cargo test --all-targets --all-features
```

The same prefix was used for Clippy. No machine-specific build path was added
to repository configuration. Tests cover technical TSV validation and multiple
T pieces, SE/PE physical orientation, full/partial adapter lengths 42/64/100/128,
full mismatch counts above 255, exact partial thresholds, low complexity, real
human/bacterial negative windows on both strands, bounded barcode correction,
route-local safety, pair rejection, gzip ordering across batches, before/after QC,
empty retained output, and FASTQ/CRAM regression paths.
