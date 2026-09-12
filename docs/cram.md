# Unmapped CRAM support

Plexless uses rust-htslib/HTSlib as a format boundary around the same qualified
read-structure, decoder, hierarchical routing, trimming, and QC implementation
used for FASTQ. Record metadata never participates in barcode assignment.

## Supported matrix

| Input | Output | Status |
| --- | --- | --- |
| FASTQ/FASTQ.gz | FASTQ.gz | Existing behavior, unchanged |
| Unmapped CRAM | FASTQ.gz | Supported; metadata outside FASTQ is intentionally lost |
| Unmapped CRAM | Unmapped CRAM | Supported; metadata-preserving transformation contract below |
| FASTQ | CRAM | Not supported |
| Aligned CRAM | Any | Not supported |

CRAM input always requires one explicit `--read-mode single|paired` and one
explicit `--output-format fastq|cram`. Plexless does not infer either choice.

## Input validation

All records must be primary, unmapped records with no reference ID, position,
mapping quality, or CIGAR. Secondary and supplementary records are rejected.
SEQ must be present. Every numeric Phred score must be present and in `0..=93`;
the CRAM boundary converts it to Phred+33 for the existing core and converts it
back exactly for CRAM output.

Single mode rejects paired flags and mate fields. Paired mode requires the
paired and mate-unmapped flags, exactly one of READ1/READ2, and unmapped mate
fields. A QNAME group may contain one R1 and one R2 in either physical order,
or one primary mate, which becomes an existing Plexless R1/R2 orphan. Duplicate
primary R1 or R2 records are rejected.

Paired CRAM must carry authoritative `@HD SO:queryname` or `@HD GO:query`
metadata. Known incompatible or unspecified grouping is rejected before the
output directory is created. Plexless processes one adjacent QNAME group at a
time; it does not maintain an unbounded coordinate-sorted mate cache. The
header declaration is treated as the input contract and is not independently
re-sorted.

## Metadata preservation policy

CRAM output preserves the input header through HTSlib, including applicable
`@HD`, `@SQ`, `@RG`, existing `@PG`, `@CO`, and other valid records, then
appends a uniquely named Plexless `@PG` entry with program name, version, and
command line. Existing read groups are not reassigned. Output filenames express
the demultiplexing destination.

For assigned records, Plexless intentionally replaces SEQ and QUAL with the
trimmed biological suffix. For unassigned records and orphans, SEQ/QUAL and all
auxiliary tags are unchanged. The following table is the v1 policy when a
record or its mate is trimmed:

| Metadata | Action | Reason |
| --- | --- | --- |
| QNAME; paired/READ1/READ2/unmapped, duplicate, and QC-fail flags; RG, BC, QT, RX, QX, MI | Preserve | The demultiplexing transformation does not alter their defined value |
| Unknown/nonstandard auxiliary tags | Preserve | Plexless cannot infer private semantics |
| `OQ`, `BQ`, `E2`, `U2` | Rewrite | Trim the position-for-position sequence/quality string by the same prefix; a wrong type or pre-trim length fails the run |
| `R2`, `Q2` | Rewrite when the represented mate was trimmed | Replace with the transformed mate sequence or Phred+33 quality string |
| `MM`, `ML`, `MN` | Remove | Base-modification coordinates/probabilities cannot be safely repaired in v1 |
| `MD`, `NM`, `SA`, `MC`, `AS`, `XS`, `UQ`, `MQ`, `AM`, `SM`, `CM`, `NH`, `HI`, `IH`, `CC`, `CP`, `CG`, `H0`, `H1`, `H2`, `PQ`, `TS` | Remove | Current alignment/mapping-derived values are invalid or meaningless after an intentional sequence change |
| `OA`, `OC`, `OP`; `CS`, `CQ`, `FZ`; `FI`, `FS`, `TC`; other defined barcode/library/program tags | Preserve | These describe historical/original data, technology metadata, template structure, or provenance rather than a current per-base alignment |
| `PT` | Remove | Local-coordinate annotation can be invalidated by prefix trimming |

The preserve/rewrite/remove decision is aggregate-counted by tag in
`cram_metadata.tsv` and summarized once per run. It is never logged once per
record. Unknown tags are preserved by contract, but their private semantics
remain the producer's responsibility.

## Output and completion

Single-end output uses `<sample>.cram`. Paired mates share the same
`<sample>.cram` and preserve their physical relative order within the ordered
QNAME stream. Enabled unassigned output is `unassigned.cram`; there are no
separate R1/R2 CRAM files and no index.

CRAM streams cannot safely be closed and later reopened for append like
concatenated gzip. An HTSlib append/reopen experiment produced an unreadable
file with a missing terminal EOF block, so Plexless deliberately rejects that
design rather than treating CRAM like gzip. All active CRAM destinations remain
open for the run.

Before creating the output directory, Plexless reserves one descriptor for
every possible sample/unassigned CRAM destination plus 64 descriptors for
input, reports, and process headroom. On Unix it reads `RLIMIT_NOFILE`; when the
hard limit is sufficient but the soft limit is not, Plexless raises its own
soft limit to exactly the required value. An insufficient hard limit, or an
explicit `--max-open-files` smaller than the required writer count, fails before
processing and reports the required, soft, and hard values. Plexless never
lowers an existing limit. On platforms without Unix resource limits, a
conservative 256-writer preflight applies.

CRAM writer memory is approximately linear in the number of populated output
files. On the qualification host, actual Plexless runs with four records per
sample peaked at about 152 MiB for 384 outputs, 291 MiB for 750, and 575 MiB for
1,500. The 1,500-output case finalized in 0.113 seconds and passed independent
`samtools quickcheck`; workloads with substantially more outputs should budget
memory accordingly.

Normal production execution does not reopen or decode completed outputs.
Success requires every record write to return successfully, every writer to
return success from explicit `hts_close`, actual successful-write counts to
reconcile with assigned/unassigned and mate-aware QC counts, and all reports to
finish. Only then is `PLEXLESS_INCOMPLETE` removed. Decode, routing,
transformation, write, finalization, report, or reconciliation failure retains
the marker. Full output decoding, exact-record verification, and independent
format checks belong to the qualification tooling.

## Implementation and native requirements

The first backend is `rust-htslib 1.0.1`, which bundles HTSlib 1.19.1, with
default features disabled and only `bzip2` and `lzma` enabled. This keeps
remote/curl support out while retaining normal CRAM codec compatibility.
HTSlib is C software, so source builds require a C compiler, `make`,
`pkg-config`, Clang/libclang for bindings, and the zlib development library.
Plexless remains MIT OR Apache-2.0; rust-htslib is MIT-licensed and HTSlib uses
the MIT/Expat license.

`hts-sys` transitively requests static zlib. Plexless's checked-in Cargo
configuration sets `LIBZ_SYS_STATIC=0`, the supported `libz-sys` override, so
the existing flate2 path and HTSlib share the established system `libz` rather
than silently replacing FASTQ compression. This is both a compatibility and a
measured FASTQ-regression safeguard.

For CRAM input, `--threads` is the global staged-pipeline budget. Plexless calls
rust-htslib `Reader::set_threads(N)`; HTSlib defines `N` as extra background
decoder workers in addition to the calling reader. FASTQ output assigns one
quarter of budgets of four or more to CRAM decoding, capped at eight decoder
workers. CRAM output assigns one decoder worker because measurements show its
single ordered metadata-preserving writer is the bottleneck. One reader and one
ordered output writer are accounted first; remaining slots become Plexless
routing/output workers. At `--threads 1`, a dedicated serial path reads,
routes, and writes inline on the caller thread with no background decoder or
stage threads. A budget of two uses the staged pipeline with no HTSlib
background decoder and one reported thread of minimum stage overcommit. Budgets
of four and above fit exactly. FASTQ input thread planning is unchanged.

| `--threads` | Decode -> FASTQ | Plexless -> FASTQ | Decode -> CRAM | Plexless -> CRAM | Reader + writer | Overcommit |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0 | 0 | 0 | 0 | 1 (inline) | 0 |
| 2 | 0 | 1 | 0 | 1 | 2 | 1 |
| 4 | 1 | 1 | 1 | 1 | 2 | 0 |
| 8 | 2 | 4 | 1 | 5 | 2 | 0 |
| 12 | 3 | 7 | 1 | 9 | 2 | 0 |
| 16 | 4 | 10 | 1 | 13 | 2 | 0 |
| 24 | 6 | 16 | 1 | 21 | 2 | 0 |
| 32 | 8 | 22 | 1 | 29 | 2 | 0 |
| 64 | 8 | 54 | 1 | 61 | 2 | 0 |

Per-output CRAM writers remain single-threaded. rust-htslib exposes a shared
HTSlib thread-pool API, but attaching private pools to hundreds of writers would
violate the global budget, and writer-side threading was not needed to fix the
measured decoder starvation. Packed-base shortcuts remain deferred.

## Validation and known limitations

The deterministic validation generator can create equivalent FASTQ and
unmapped CRAM fixtures with `--with-cram`. `--include-cram` extends the existing
deep harness to compare classification counts, samples, trimming, SEQ, QUAL,
mate relationships, per-output ordering, logical hashes, and QC reports. CRAM
output is reopened through rust-htslib, decoded by a validation-only adapter,
and checked against the same truth. If available, `samtools quickcheck -u`
provides a separate reference-free CRAM integrity check.

Current limitations are intentional: no aligned input/output, no FASTQ to
CRAM conversion, no missing qualities, no coordinate-sorted mate collation,
no CRAM index, one ordered CRAM output writer stage, linear memory/descriptor
growth with output count, and removal rather than coordinate rewriting for
complex base-modification tags. A two-thread CRAM request retains the
one-thread staged-pipeline minimum overcommit described above.
