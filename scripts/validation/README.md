# FASTQ and unmapped-CRAM qualification

These scripts provide deterministic, self-contained correctness and
concurrency qualification for the existing FASTQ implementation. Generated
data and outputs default below `$HOME/plexless_validation_data`; they are not
written into the repository.

The FASTQ suite remains the reference oracle. Optional equivalent unmapped
CRAM fixtures enter the same routing core and are checked against the same
truth. CRAM fixture conversion is a repository-contained Rust example and is
not product FASTQ-to-CRAM functionality.

The suite treats `--threads 1 --output-mode direct` as the conservative
reference. Every other run is compared against that reference using complete
decompressed FASTQ SHA-256 values, exact records, and deterministic QC reports.
Compressed gzip bytes are deliberately not used as the equivalence contract.

## What is covered

The generated fixture families are:

| Fixture | Modes | Semantics |
| --- | --- | --- |
| `core` | SE and PE | exact assignment, one-mismatch correction, one observed `N`, unmatched, valid unrouted root, short structured reads, technical trimming, biological suffix preservation, repeated/cross-mate symbols, reverse-complement pieces, and R1/R2 orphans with resynchronization |
| `n2` | SE | two-mismatch correction and an observed `N` plus a nucleotide mismatch with `--max-mismatches 2` |
| `r1_only` | PE | asymmetric paired layout with only R1 structured and all of R2 biological |
| `r2_only` | PE | asymmetric paired layout with only reverse-complemented R2 structured and all of R1 biological |
| `negative` | SE and PE | truncated FASTQ, sequence/quality mismatch, corrupt gzip CRC, desynchronization beyond 1,024 records, correction-unsafe whitelist geometry, invalid barcode character, and malformed sample sheet |

Runtime ambiguity is not generated: valid Plexless configuration rejects the
unsafe sibling geometry required to create it. Existing Rust defense-in-depth
tests cover the decoder's ambiguous branch directly.

Truth is grouped by output filename and output index so the verifier remains
streaming for stress fixtures. Every truth row includes input ordinal, identity,
mate, semantic case, routing category, sample, full expected output header,
trimmed or untrimmed sequence, and quality.

With CRAM enabled, the same five positive fixture/mode combinations cover
single/paired interpretation, flag-based mate ordering, orphans, repeated and
cross-mate symbols, reverse-complement segments, asymmetric structures,
trimming, and all route categories. Focused Rust tests additionally cover both
physical mate orders, quality 0/93, aligned/secondary/supplementary rejection,
malformed flags and duplicate mates, incompatible headers, missing SEQ/QUAL,
corrupt CRAM, header preservation, and auxiliary-tag policy.

## Generate fixtures

Tiny fixtures contain 4,096 core records/events and are suitable for quick CI
or integration qualification:

```bash
python3 scripts/validation/make_plexless_validation_dataset.py \
  --profile tiny \
  --outdir "$HOME/plexless_validation_data/tiny"
```

Add equivalent reference-free, unmapped CRAM fixtures with `--with-cram`:

```bash
python3 scripts/validation/make_plexless_validation_dataset.py \
  --profile tiny --with-cram \
  --outdir "$HOME/plexless_validation_data/tiny-cram"
```

The converter preserves deterministic logical records, writes numeric CRAM
qualities, marks R1/R2 with SAM flags, and declares paired files queryname
sorted. External-source injection currently remains FASTQ-only because its
input ordering is not guaranteed by the generator.

The standard desktop fixture contains 1,000,000 SE records and 1,000,000 PE
events. Focused fixtures remain 4,096 records/pairs:

```bash
python3 scripts/validation/make_plexless_validation_dataset.py \
  --profile standard --with-cram \
  --outdir "$HOME/plexless_validation_data/standard-cram"
```

Stress generation is opt-in. It defaults to 10,000,000 core records/events and
streams directly to compressed files and truth spools:

```bash
python3 scripts/validation/make_plexless_validation_dataset.py \
  --profile stress \
  --outdir "$HOME/plexless_validation_data/stress"
```

Override SE and PE sizes independently:

```bash
python3 scripts/validation/make_plexless_validation_dataset.py \
  --profile stress --records 2000000 --pairs 12000000 \
  --outdir "$HOME/plexless_validation_data/custom-stress"
```

The same seed and options produce logically identical FASTQ and truth content.
The gzip writer also fixes header timestamps, so generated compressed bytes are
reproducible. Use `--seed` to select another deterministic dataset. Existing
directories are refused; pass `--replace` only when replacement is intended.

The original external-source injection workflow remains available:

```bash
python3 scripts/validation/make_plexless_validation_dataset.py \
  --r1 raw_R1.fastq.gz --r2 raw_R2.fastq.gz --limit 100000 \
  --outdir "$HOME/plexless_validation_data/injected"
```

## Run qualification

Build first:

```bash
cargo build --release
```

Fast CI qualification uses threads 1 and 4, both output modes, every focused
fixture, and all negative cases:

```bash
scripts/validation/run_plexless_deep_validation.sh \
  --preset ci --dataset "$HOME/plexless_validation_data/tiny"
```

For a dataset generated with `--with-cram`, add `--include-cram` to run both
CRAM output formats after the unchanged FASTQ qualification:

```bash
scripts/validation/run_plexless_deep_validation.sh \
  --preset ci --include-cram \
  --dataset "$HOME/plexless_validation_data/tiny-cram"
```

Normal deep validation uses 1, 2, 4, and 8 threads:

```bash
scripts/validation/run_plexless_deep_validation.sh \
  --preset normal --dataset "$HOME/plexless_validation_data/tiny"
```

The full i7-13700K desktop matrix uses 1, 2, 4, 8, 12, 16, and 24 threads for
both SE and PE core data and both direct/buffered strategies. Focused semantic
fixtures run at 1 and 24 threads:

```bash
scripts/validation/run_plexless_deep_validation.sh \
  --preset desktop --include-cram \
  --dataset "$HOME/plexless_validation_data/standard-cram"
```

`--threads` and `--output-modes` can override a preset. The first thread count
must be `1`, and `direct` must be present because that combination is the
reference.

### Repeatability / race qualification

`--repeat N` repeats both output strategies at the highest requested thread
count and compares every run with the one-thread reference. With CRAM enabled,
it also repeats CRAM-to-FASTQ and native CRAM output at the highest thread
count. Lower thread counts still run once:

```bash
scripts/validation/run_plexless_deep_validation.sh \
  --preset desktop --repeat 10 \
  --dataset "$HOME/plexless_validation_data/standard"
```

This targets nondeterministic ordering, missing/duplicate batches, wrong
routing, partial output, shutdown races, and intermittent deadlocks.

To repeat only CRAM paths while retaining the one-thread CRAM-to-FASTQ
reference, use the validation-only `--cram-only` selector:

```bash
scripts/validation/run_plexless_deep_validation.sh \
  --preset desktop --cram-only --threads '1 24' \
  --output-modes direct --repeat 10 \
  --dataset "$HOME/plexless_validation_data/standard-cram"
```

### Stress qualification

Use the stress preset only with an explicitly generated stress fixture:

```bash
scripts/validation/run_plexless_deep_validation.sh \
  --preset stress --repeat 3 \
  --dataset "$HOME/plexless_validation_data/stress"
```

The stress preset defaults to threads 1 and 24. Large stress qualification is
intentionally separate from CI and normal desktop work.

### CRAM writer fan-out qualification

The deterministic fan-out fixture opens real CRAM destinations without a
large biological dataset:

```bash
for samples in 384 750 1500; do
  scripts/validation/make_cram_writer_fanout_fixture.py \
    --samples "$samples" --records-per-sample 4 \
    --outdir "$HOME/plexless_validation_data/cram-fanout-$samples"
done
```

Run each with `--output-format cram`, `--read-mode single`, structure
`R1_8A2T`, and `--max-mismatches 0`. Measure with GNU
`time -f '%e\t%P\t%M'`. On Linux,
`prlimit --nofile=256:4096 -- COMMAND` can test automatic soft-limit raising in
a child process; a hard value below `samples + 64` tests fail-before-output
behavior without altering the invoking shell.

## Verification and artifacts

For every successful run, the harness requires:

- exact classification and optional short/orphan counts;
- exact output-file inventory;
- every expected record once, in promised per-output order;
- exact expected header, trimmed/untrimmed sequence, and quality;
- correct sample destination and paired-mate identity;
- exact `sample_metrics.tsv`, `fastq_stats.tsv`, and `barcode_stats.tsv`
  equivalence with the reference;
- Python gzip decompression and, when installed, independent `gzip -t`;
- matching decompressed logical FASTQ SHA-256 for all non-reference runs;
- absence of `PLEXLESS_INCOMPLETE` after success.

CRAM-to-FASTQ runs additionally require the exact same category, sample,
trimmed sequence, numeric-quality-equivalent FASTQ quality, mate relationship,
order, logical SHA-256, and reports as the CRAM one-thread reference. After the
Plexless process exits, CRAM output is decoded to validation FASTQ and compared
to truth; this reread is external qualification, not production completion
logic. The harness also checks `cram_metadata.tsv`. When samtools is installed,
`samtools quickcheck -u` is an optional independent integrity check. Its `-u`
mode is required because these raw CRAM fixtures intentionally have no reference
targets.

Negative operational runs must exit nonzero, must not print a success summary,
must retain `PLEXLESS_INCOMPLETE`, and must not write final sample metrics.
Startup configuration failures occur before output setup and therefore do not
create the marker. A missing mate is not a fatal case under current Plexless
semantics: it is an orphan, and the core PE truth verifies its unassigned output
and later resynchronization.

Each qualification gets a new timestamped directory under the dataset unless
`--results` is supplied. Existing result paths are refused. Important files:

- `run_metrics.tsv`: Plexless wall time, throughput, aggregate CPU, peak RSS
  (GNU `time` when available), output size, and separately measured external
  validation time;
- `*.logical-sha256.json`: per-file decompressed hash, bytes, records, and case
  counts;
- `*.stderr`: complete Plexless diagnostics and summary;
- `*.time.tsv`: raw wall-time and RSS measurement.

If `gzip` is absent, the independent check is clearly reported as skipped;
Python gzip parsing and logical hashing still run. See
[`AUDIT_DESIGN.md`](AUDIT_DESIGN.md) for the production manifest and hashing
cost analysis.

If samtools is absent, only the independent CRAM compatibility check is
skipped; the repository-contained rust-htslib decoder still fully decodes
every produced CRAM as an external qualification condition.

## Interpreting failures

- A truth mismatch names the output record, input ordinal, and semantic case;
  it indicates routing, trimming, quality, mate, duplication/loss, or ordering
  divergence.
- A logical SHA-256 mismatch with truth checks passing would indicate a verifier
  defect; both are intentionally required.
- A QC report mismatch means route counts or input-level statistics changed
  despite identical logical FASTQ.
- A retained marker after a nominally successful run is a completion failure.
- A missing marker after an operational failure means partial output could look
  complete and is a release blocker.
- A timeout or hang is not converted into success. Use an external job timeout
  for unattended stress runs and preserve the result directory for diagnosis.
