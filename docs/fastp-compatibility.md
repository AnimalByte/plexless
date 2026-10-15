# fastp-derived preprocessing

Plexless implements the fastp preprocessing feature families in Rust while keeping
its barcode routing, bounded pipeline, ordered output, and CRAM metadata policy.
This is feature coverage, not a drop-in fastp CLI, report schema, or bit-for-bit
algorithm replacement. Quality trimming and filtering remain opt-in under ordinary
Plexless defaults. Paired overlap analysis is enabled by default.

The implementation was reviewed against OpenGene/fastp source revision
[`ee125a52`](https://github.com/OpenGene/fastp/tree/ee125a52fddae9d8189e13cd8a33ad04988a1351).
The official fastp **1.4.0** executable was used for sequence/quality comparisons.
See [THIRD_PARTY.md](../THIRD_PARTY.md) and the retained
[MIT license](../licenses/fastp-MIT.txt). These features do not invoke fastp at runtime.

## Feature coverage and defaults

Both `demux` and `trim` accept these controls. All barcode assignment uses original
sequence, before any UMI, quality trimming, correction, merging, or filtering.

| Feature | Plexless controls and behavior |
| --- | --- |
| Standard quality/length profile | `--fastp-defaults`: Q15 qualified bases, at most 40% low quality, at most 5 Ns, minimum length 15; duplicate estimation and JSON/HTML enabled. Enables SE adapter inference and file-based instrument detection for polyG. Explicit filter thresholds take precedence. |
| Quality filtering | `--qualified-quality`, `--max-low-quality-percent`, `--max-n`, `--min-average-quality`; after trimming. |
| Length filtering | `--min-length`, `--length-limit`; zero-length reads pass only thresholds that allow zero. |
| Disable profile filters | `--disable-quality-filtering`, `--disable-length-filtering`; the latter also removes `--length-limit`. |
| Fixed clipping | `--trim-front1/2`, `--trim-tail1/2`, `--max-len1/2`; unspecified R2 settings inherit R1. Coordinates begin after explicitly declared proximal structure and read UMI. |
| Quality windows | `--cut-front`, `--cut-tail`, `--cut-right`; default window 4, mean Q20. Common `--cut-window-size/--cut-mean-quality` and per-direction overrides. Right cutting takes precedence over tail cutting. |
| PolyG | `--trim-poly-g`, `--poly-g-min-len` (10); profile detection from supported instrument prefixes for seekable FASTQ. `--disable-trim-poly-g` disables automatic detection. |
| PolyX | `--trim-poly-x`, `--poly-x-min-len` (10); permits the bounded tail mismatch policy derived from fastp. |
| Low complexity | `--low-complexity-filter`, `--complexity-threshold` (30% adjacent-base changes). Reads shorter than two bases fail this enabled filter. |
| Phred+64 input | `--phred64`; validate input encoding and convert once in the producer. Output remains Phred+33. |
| Explicit adapters | `--r1-distal-sequence`, `--r2-distal-sequence`, `--adapter-fasta`; strict A/C/G/T observed orientation, no manufacturer catalog. FASTA entries apply to both mates in their supplied orientation. |
| Automatic adapters | `--detect-adapter` (SE or both PE mates), `--detect-adapter-for-pe`; bounded discovery from at most 5,000 records/mate and 200 tail bases/record. Includes seekable interleaved input. Disabled for stdin. |
| Paired overlap | Default on, minimum 30 bases, at most 5 differences and 20% differences over the **entire** overlap. `--no-overlap` disables adapter overlap analysis. `--require-adapter-for-overlap` additionally requires explicit/compiled sequence evidence. |
| One-gap adapter overlap | `--allow-gap-overlap-trimming`; short-insert fallback only after no ungapped alignment, at least 16 configured overlap bases, one internal gap at least 8 matched positions from either edge. Not used for correction or merging. |
| Adapter dimers | `--dimer-max-len` (2) counts adapter-supported short retained reads per mate; does not itself discard reads. |
| Base correction | `--correction`: overlap disagreement is corrected only with Q30-or-higher support against Q14-or-lower quality; copy supported base and quality. Default off. |
| Merge | `--merge`, `--merged-out`, `--include-unmerged`; unique ungapped overlap, quality-supported consensus. Default off. |
| UMI | `--umi --umi-loc read1/read2/per-read/index1/index2/per-index`, `--umi-len`, `--umi-skip`, `--umi-prefix`, `--umi-delim`. Read UMI follows any declared proximal structure. Same tag on both mates, before `/1` and `/2` suffixes. Inferred opposite distal signatures mask read UMI/skip cycles at their physical position; these add no known evidence. |
| Exact duplicate removal | `--dedup`; original sequences, case-insensitive, pair-aware, independent of IDs and qualities. Deterministically retain the first fragment in input order (it still must pass subsequent filters). |
| Duplicate reporting | `--evaluate-duplication`; bounded four-hash Bloom estimate, also on in the profile unless `--dont-eval-duplication`. Reporting never rejects reads. |
| Duplicate resource control | `--duplication-memory-mib` (256); explicit MiB budget replaces fastp's accuracy-level preset. Exact removal spills to transient SQLite; estimation uses a fixed bitset. |
| Index blacklist | `--filter-by-index1/2`, `--filter-by-index-threshold`; one DNA sequence per line, exact lengths, independent mismatch policy. Uses CASAVA indexes after the final colon; index2 requires `index1+index2`. Missing/malformed requested indexes cause an error. |
| Overrepresented sequences | `--overrepresentation-analysis`, sampling every 20 records by default (`--overrepresentation-sampling`); before/after sampled substring occurrence counts with bounded candidate storage and replacement-error bounds. Creates JSON if no report path is supplied. |
| Reports | `--json`, `--html`, `--report-title`; standalone interactive HTML, per-cycle bases/mean quality, length/quality histograms, filtering/trimming/duplicate summaries, before/after biological QC. Profile creates `plexless.json` and `plexless.html`. |
| Input modes | `--reads`, separate `--r1/--r2`, `--stdin`, `--interleaved-in`; paired interleaved IDs are checked and incomplete pairs error. |
| Output modes | Default gzip or explicit `--out1/--out2`; `.gz` selects gzip, other extensions plain FASTQ. `--stdout` additionally streams retained output, interleaved for PE, preserving Plexless file output. |
| Auxiliary outputs | `--failed-out`, `--unpaired1/2`, `--overlapped-out`; passing mate from a filtered pair goes only to requested auxiliary output. Equal unpaired paths share a single ordered stream. Perfect-overlap output contains R1's overlap in its observed orientation. |
| Processing limit | `--reads-to-process N`; SE records or PE pairs, zero means all. Pair resynchronization limits each input to N records, so missing mates can yield fewer than N pairs. |
| Split files | `--split-by-lines N` (multiple of 4, at least 1000), or `--split N` (2–999); `--split-prefix-digits` (0–10). Split at record boundaries before compression. Paired primary parts have corresponding records. |
| Miscellaneous | `--fix-mgi-id`, `--verbose`, existing `--threads` and `--compression-level`. Outputs are protected from overwrite by default; `--dont-overwrite` is accepted and preserves that policy. |

The original barcode/sample TSV formats and 32-base logical barcode limit are
unchanged. Snake-case aliases are accepted for most fastp-derived long options;
UMI location values accept `per_read`/`per_index` as well as hyphenated forms.
Manufacturer adapter sequences are never inferred from an instrument name.

Successful startup adapter inference enables conservative partial matching for
that run. Explicit adapter configurations otherwise retain the separate
`--partial-pass-through` control.

## Matching safeguards and deliberate differences

A sequence-supported full adapter match uses an exact diverse 12–16-base anchor,
then masked Hamming verification of all known positions in chunks of at most 32.
It requires at least **20 known bases** by default, configurable with
`--min-full-adapter-known-bases`. Masked T positions add no evidence. A mismatch
budget never relaxes the exact anchor or treats an observed N as a known base.
Partial sequence matching is separately enabled with `--partial-pass-through`,
requires at least 20 informative **known** bases by default, and is exact at known
positions. Across multiple adapters it favors the most known evidence, then the
least clipping; full matching has precedence over partial matching.

Paired overlap considers relative alignments of biological R1 and reverse-complement
R2, including ordinary long inserts and read-through short inserts. Every position
must satisfy both absolute and percentage mismatch limits. A qualifying region
must contain at least three bases, six dinucleotides, and no base above 75%.
Multiple qualifying alignments preserve the pair. Overlap alone can establish
an adapter boundary when sequence matching offers no candidate. Where a full
candidate exists, it must agree with the inferred boundary; an unconfirmed full
match is never replaced by a weaker partial match elsewhere. Configured front
clipping is included in the physical boundary calculation. Different mate lengths
are allowed. Sequence/quality bounds always move together.

fastp's inspected overlap implementation permits a protected-prefix shortcut
that can accept excessive late mismatches and chooses the first qualifying
alignment. Plexless checks the whole overlap and rejects non-unique alignments.
This can retain contamination that fastp trims, especially in repeated sequence,
very short inserts, or low-complexity libraries. It is an intentional specificity
tradeoff. No inference can resolve DNA that is indistinguishable from an adapter
or a unique apparent overlap by sequence alone; library-specific validation remains
necessary. These regression tests do not replace validation on representative library data.

Plexless does not copy fastp's short four-base sequence matches, manufacturer
fallback catalog, protected-prefix mismatch shortcut, or probabilistic read
removal. Barcode correction remains subject to strict route-local geometry.
Merging selects the higher-quality observation at a disagreement and preserves
biological overhangs covered by the other mate after front clipping; header annotation
and report schemas differ. An optional one-gap result authorizes adapter bounds
only, with separate mate lengths; merge/correction remain ungapped. Indels in an
explicit full/partial adapter signature are not accepted.

## Workflow and pairing

Processing is original-read barcode assignment, explicit proximal/UMI clipping,
fixed/quality clipping and polyG, adapter/overlap analysis, polyX/maximum-length
clipping, optional correction/merge, filtering, retained QC, ordered compression.
No barcode assignment depends on corrected or filtered bases.

Without merging, if either mate fails filtering the primary pair is rejected.
Optional auxiliary files can retain the passing mate and failing reads; they do
not change the paired primary counts or retained-pair metric. With merging,
a successfully merged fragment is filtered as one consensus read. Unmerged
passing pairs remain paired, or are interleaved in merged output when
`--include-unmerged` is enabled. Rejected pairs are still rejected as pairs.
Existing resynchronization orphans retain their separate legacy routing; they
are not newly classified as QC-passing unpaired fragments.

Custom primary/merged paths require exactly one output sample. Multisample demux
uses per-sample default filenames and a merged stream per sample. Auxiliary
outputs intentionally combine samples. Stdout mirrors emitted retained records;
with merging it streams consensus and unmerged records in input order.
Splitting and stdout cannot be combined. Output path collisions are errors,
except the explicitly shared unpaired path. Inputs and existing outputs are
never overwritten. Native CRAM permits shared bounds/filter/QC operations, but
FASTQ-only UMI/header, encoding, merge, auxiliary, custom paths, split, and stream
features are rejected explicitly so original CRAM metadata remains meaningful.

## Bounded resources and costs

Normal adapter compilation happens once. Exact barcode lookup remains packed;
large correction neighborhoods switch to a packed-distance scan instead of an
exponential index. Adapter masks/chunks and overlap scratch are reused in each
worker batch. Queues, buffered output, compression work, and writer caches remain
bounded. A rare long read/signature can grow reusable scratch to that read's size.

Exact duplicate keys initially use bounded memory, then a transient SQLite BLOB
catalog with a bounded page cache. Temporary disk use scales with distinct
original fragments. The database is removed on close. Database/disk errors stop
the run and preserve the incomplete marker; reads are never silently classified
as duplicates on an I/O error. This optional feature adds ordered producer work
and can reduce throughput. The budget is an approximate catalog/cache budget,
not a process RSS ceiling; scratch, allocator metadata, queues, and SQLite have
additional overhead. Bloom reporting may overestimate duplicates as its bitset
fills; it never controls removal.

Overrepresentation uses 4,096 Space-Saving candidates per mate/stage and a bounded
stale-entry heap. It examines sampled substrings of lengths 10, 20, 40, 100, 150,
plus full reads up to 512 bases. Reports show at most 100 sequences and lower/upper
occurrence bounds including replacement error. It counts occurrences, not unique
reads, and differs from fastp's startup-candidate implementation. Counts are sampled,
not a complete frequency census. Post-filter sampling restarts per worker batch;
input sampling is global per mate. Selection and ordered aggregation are deterministic
for a fixed input and batching policy.

JSON/HTML detailed curves and overrepresentation add sequence scans only when
requested. Quality filters and retained summary QC share a base-level scan.
Read transformations allocate owned records only when needed by mutable/header or
auxiliary features. Exact dedup stores original sequences intentionally.

`--split-by-lines` is streaming and counts serialized FASTQ records. `--split`
adds a full counting pass over seekable input, then uses `ceil(input fragments/N)`
records per part per stream. Filtering/fanout can yield fewer than N parts; merged
interleaved and combined auxiliary streams can yield more. It is an input-based
size target, not an exact count of final files. Use a lines limit for strict
per-file record limits or streaming input. Splitting is performed in the ordered
buffered pipeline without decompressing and recompressing output.

## Examples

```bash
# Established Plexless behavior: explicit opt-in filtering, default paired overlap.
plexless demux --r1 R1.fq.gz --r2 R2.fq.gz \
  --r1-structure R1_10A33T --r2-structure R2_10B \
  --barcodes barcodes.tsv --samples samples.tsv --output demux --threads 8 \
  --technical-sequences technical.tsv --min-length 30 --max-n 5 --fastq-stats

# Already-demultiplexed T7: no proximal structure or barcode file.
R1_DISTAL=$(tr -d '\r\n' < r1_distal.txt)
plexless trim --reads t7.fq.gz --r1-distal-sequence "$R1_DISTAL" \
  --partial-pass-through --min-partial-adapter-overlap 20 \
  --cut-right --cut-window-size 4 --cut-mean-quality 20 \
  --min-length 30 --output t7-trim --json t7.json --html t7.html --threads 8

# Paired library: overlap inference, optional adapter discovery, profile filters.
plexless trim --r1 R1.fq.gz --r2 R2.fq.gz --fastp-defaults \
  --detect-adapter-for-pe --correction --output cleaned --threads 8

# UMI, exact duplicate removal, merging and explicit outputs.
plexless trim --r1 R1.fq.gz --r2 R2.fq.gz --umi --umi-loc per-read \
  --umi-len 8 --umi-skip 2 --dedup --duplication-memory-mib 256 \
  --merge --merged-out merged.fq.gz --out1 unmerged.R1.fq.gz \
  --out2 unmerged.R2.fq.gz --output metrics --fastq-stats --threads 8

# Streaming interleaved input, mirrored ordered stdout and split file mode.
cat interleaved.fq | plexless trim --stdin --interleaved-in --stdout \
  --output stream-metrics --threads 8 > cleaned.interleaved.fq
plexless trim --reads input.fq.gz --split-by-lines 4000000 \
  --output split --threads 8

# Standalone reports with sampled overrepresented sequences.
plexless trim --reads input.fq.gz --output qc --json qc.json --html qc.html \
  --overrepresentation-analysis --overrepresentation-sampling 20 --fastq-stats
```

## Validation and fastp comparisons

`tests/fixtures/fastp_v1_4_0.json` contains synthetic original inputs and sequence/
quality outputs generated by the official fastp 1.4.0 binary. Cargo integration
checks replay them without needing fastp installed. Regenerate/recheck with:

```bash
python3 scripts/compare_fastp.py --fastp /path/to/fastp \
  --plexless target/debug/plexless
```

The nine comparisons cover fixed clipping, front/tail/right windows, polyG,
polyX, UMI clipping, paired merging, and short-insert overlap clipping. Separate
tests cover full/partial adapters, technical orientation, large barcode mismatch
budgets, repeated/invalid/late-mismatch overlaps, correction quality thresholds,
exact dedup determinism/spill, reports, paired filtering/auxiliary streams,
interleaved stdin, and split gzip output. Header decorations are tested by
Plexless-specific tests rather than assumed identical to fastp.

For performance comparisons, use the same input files, quality encoding, adapter
knowledge, evidence and mismatch policies, filtering, duplicate mode, compression
level, output format, storage, and thread/CPU budget. First compare every retained
ID/sequence/quality and rejection count; report intentional differences. Measure
wall/CPU time, peak RSS, temporary disk, compressed output bytes, and reads/bases
processed on SE, ordinary PE, short-insert PE, long technical structures, high-N/
low-quality, many-sample fanout, and duplicate-rich workloads. Separate inference
startup and exact-dedup spill costs. There is no measured speedup claim here.
