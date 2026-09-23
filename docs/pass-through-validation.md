# Pass-through validation for v0.6.0

The feature checkpoint is `2afbf07`; the baseline is release commit
`9897ff68caa985f2b1892e0962bcfe975b032953` (v0.5.0). All tests and measurements
use portable release settings, without `target-cpu=native`.

## Implementation

- `src/pass_through.rs`: startup signature compilation, exact rolling scanner,
  chunked verification, and shared biological coordinates.
- `src/routing.rs`: immutable per-sample/per-mate signature ownership, separate
  from routing decisions.
- `src/demux.rs`, `src/parallel.rs`: invoke detection only after assignment and
  carry both boundaries through serial, direct, buffered, and CRAM paths.
- `src/cram.rs`: apply identical boundaries to sequence, quality, and per-base
  metadata; update mate tags from independently trimmed mate sequences.
- `src/encoder.rs`'s existing `EncodedBarcode` is reused unchanged. The existing
  Hamming primitive in `src/decoder.rs` is shared with verification.
- `src/cli.rs` and `src/lib.rs`: default-on behavior, module registration, `--no-pass-through`, independent distal
  mismatch budget, and generic optional R1/R2 overrides.
- Tests, README, architecture/CRAM documentation, and the Python validation
  harness cover the new behavior. Existing Rust option literals gain the new
  default fields. `Cargo.toml`/`Cargo.lock` advance to 0.6.0.

The exact signature representation is:

```rust
struct PassThroughSignature {
    anchor: u64,
    anchor_mask: u64,
    anchor_len: usize,
    anchor_offset: usize,
    verification: Box<[VerificationChunk]>,
    verification_len: usize,
    mismatch_limit: u8,
}
struct VerificationChunk {
    expected: EncodedBarcode,
    known_mask: u64,
}
```

Compilation indexes canonical assigned barcode sequences using the same logical
extraction plans as demultiplexing. Default inference reverses physical segment
order and inverts each piece's orientation. PE uses the opposite mate; SE assumes
symmetric ends. Repeated cross-mate pieces and `(rc)` therefore retain their
existing meaning. Overrides are written in observed distal orientation and use
the existing segment grammar and barcode database.

Anchors are 12–16 contiguous known bases, longest first, then greatest
dinucleotide diversity. At least three bases and six dinucleotides must occur;
no base may exceed 75%. Unavailable ends produce a startup diagnostic. Technical
`T` positions retain their length and reversed position but are expected unknowns.

Only exact anchor hits invoke full verification. All known positions are checked
using masked packed XOR/Hamming comparison in 32-base chunks; total length is
not limited to one word. Default distal tolerance is one substitution, separately
configurable from barcode correction. Observed unknowns reset the anchor window
and reject candidates at known positions. No observed wildcard or indel matching
is used. Failed candidates continue scanning. Bounds are computed before the
loop, and the result is `anchor_start - anchor_offset`, not `anchor_start`.

Each mate gets independent `[bio_start, bio_end]` coordinates against its original
read. The same coordinates slice sequence and quality. Classification, pair
identity/order, output selection, compression, queues, and cancellation remain
unchanged. Input-level QC keeps its existing meaning; output base counts reflect
the cleaned reads. Worker counters report R1/SE and R2 detections, scanned
negatives, and distal bases removed. Disabled/unavailable ends are not negatives.

## Correctness qualification

- `cargo fmt --check`: passed.
- `cargo clippy --all-targets --all-features -- -D warnings`: passed.
- `cargo test --all-features`: **140 passed**, none failed or ignored.
- `git diff --check`: passed.
- `cargo build --release` and release example builds: passed.
- Production `--help`, `demux --help`, and `--version`: passed (`plexless 0.6.0`).

New Rust coverage includes legal boundary extremes, empty/short inserts, exact
biological output and quality coordinates, distal-start versus anchor-start,
false anchor followed by genuine distal, partial targets, threshold failures,
N/malformed/lowercase behavior, technical masks, >32-base verification, weak
anchors, a scalar-reference comparison, sample specificity, corrected assignment,
unassigned preservation, default/disabled behavior, override validation, repeated
pieces, segment order, and independent PE lengths. CRAM tests reopen native
output and verify per-base and mate metadata as well as sequence/quality.
An additional PE regression moves an observed N through each known distal
position, creating both broken anchors and exact anchors with failed verification;
later genuine distal structure must still be found and the other mate must
remain independent, in serial/parallel and direct/buffered modes.

The existing Python generator and complete desktop validation matrix passed on
the deterministic `tiny` fixture (4,096 core records/events and 4,096 per focused
fixture): **105 successful runs plus seven expected configuration failures**.
Threads: 1, 2, 4, 8, 12, 16, 24. Both FASTQ output modes, SE/PE, CRAM-to-FASTQ,
native CRAM, and three repetitions at the highest thread count were covered.
The harness checked classification, accounting, identities, order, sequence,
quality, deterministic reports, logical hashes, gzip integrity, and independent
samtools CRAM integrity. This is the full matrix on small fixtures, not the
separate opt-in 10-million-event stress preset.

The new Python truth matrix passed **1,240 exact output-record checks** across
FASTQ SE/PE, CRAM-to-FASTQ SE/PE, and reopened native CRAM SE/PE, enabled and
disabled, one/four threads, and both applicable output modes. Inserts of
50/75/100/150/200/250/300/400/500 bases use 350-base sequencing reach beyond the
proximal prefix. Covered distal structures yield the exact biological insert;
long inserts and negative controls remain unchanged. PE additionally exercises
neither, R1-only, R2-only, and both pass-through. Controls include genomic-like
random sequence, barcode-like substrings, exact anchor-like substrings, and
partial distal motifs. Exact equality detects residual distal contamination
and over-trimming; validation is independent of detector implementation.

Reproduce using the commands in [the validation guide](../scripts/validation/README.md).

## Performance

Intel Core i7-13700K, Rust 1.98.0, portable Linux GNU builds. Each run has
1,000,000 SE reads or 1,000,000 PE pairs, 377 bases per raw read, one assigned
sample, default auto/direct output, and gzip level 2. FASTQ input is uncompressed;
CRAM input uses the existing reference-free converter. Both input formats emit
FASTQ for timing; native CRAM output is separately correctness-qualified above.
A seeded 4,096-insert pool is reused across states. Positive data contain 50%
read-through per mate, with 100-base R1/SE and 135-base R2 inserts.

There are five timed repetitions after one warmup per state, rotating state
order. Wall time includes startup, decoding, assignment, scanning, serialization,
and compression; output validation/hashing is untimed. PE rates count individual
reads (2,000,000 per run). Percentage change is throughput versus the old baseline
on the same data, with negative numbers meaning regression. Positive-case speed
changes also include reduced output work. Default negative cases isolate the
permanent scanning cost; disabled cases expose unrelated pipeline differences.

The desktop was not isolated from ordinary system activity. Timing ranges are
reported because scheduling variation is visible, particularly with PE. Positive
throughput deltas do not prove that scanning is free or improves performance.
No claim of negligible overhead is made.

### Four threads

CPU affinity: 0,2,4,6,8,10,12,14 (performance cores). The measured no-read-through
SE throughput regressions were **19.2% FASTQ** and **12.1% CRAM** versus baseline.
Disabled negative-case changes were +1.0% and +2.6%, respectively. PE ranges
overlap substantially; see raw repeats before drawing broader conclusions.

| Input | Read-through | State | Median s | M reads/s | Throughput change | Min–max s |
| --- | --- | --- | ---: | ---: | ---: | --- |
| FASTQ SE | 0% | baseline | 0.996 | 1.004 | 0.0% | 0.931–1.082 |
| FASTQ SE | 0% | disabled | 0.986 | 1.014 | 1.0% | 0.944–1.076 |
| FASTQ SE | 0% | enabled | 1.232 | 0.812 | -19.2% | 1.195–1.286 |
| FASTQ SE | 50% | baseline | 0.888 | 1.126 | 0.0% | 0.842–0.910 |
| FASTQ SE | 50% | disabled | 0.924 | 1.082 | -3.9% | 0.851–0.948 |
| FASTQ SE | 50% | enabled | 1.042 | 0.960 | -14.8% | 0.913–1.063 |
| FASTQ PE | 0% | baseline | 2.926 | 0.684 | 0.0% | 2.761–3.475 |
| FASTQ PE | 0% | disabled | 2.906 | 0.688 | 0.7% | 2.275–3.047 |
| FASTQ PE | 0% | enabled | 2.836 | 0.705 | 3.2% | 2.761–2.984 |
| FASTQ PE | 50% | baseline | 2.232 | 0.896 | 0.0% | 1.926–2.494 |
| FASTQ PE | 50% | disabled | 2.472 | 0.809 | -9.7% | 2.176–3.262 |
| FASTQ PE | 50% | enabled | 2.453 | 0.815 | -9.0% | 2.387–2.843 |
| CRAM SE | 0% | baseline | 2.611 | 0.383 | 0.0% | 2.493–3.179 |
| CRAM SE | 0% | disabled | 2.544 | 0.393 | 2.6% | 2.464–2.565 |
| CRAM SE | 0% | enabled | 2.970 | 0.337 | -12.1% | 2.868–3.200 |
| CRAM SE | 50% | baseline | 2.034 | 0.492 | 0.0% | 1.985–2.126 |
| CRAM SE | 50% | disabled | 1.896 | 0.528 | 7.3% | 1.848–2.009 |
| CRAM SE | 50% | enabled | 2.136 | 0.468 | -4.8% | 2.084–2.190 |
| CRAM PE | 0% | baseline | 4.851 | 0.412 | 0.0% | 4.653–4.969 |
| CRAM PE | 0% | disabled | 4.696 | 0.426 | 3.3% | 4.611–4.718 |
| CRAM PE | 0% | enabled | 4.928 | 0.406 | -1.6% | 4.894–5.058 |
| CRAM PE | 50% | baseline | 3.785 | 0.528 | 0.0% | 3.700–4.248 |
| CRAM PE | 50% | disabled | 3.836 | 0.521 | -1.3% | 3.687–4.012 |
| CRAM PE | 50% | enabled | 3.924 | 0.510 | -3.5% | 3.816–4.030 |

Raw [summary and binary hashes](benchmarks/pass-through-v0.6.0-4t.json) and
[all timed repetitions](benchmarks/pass-through-v0.6.0-4t.csv) are retained.

### Twenty-four threads

CPU affinity: all 24 logical CPUs (0–23), matching the requested full CPU budget.
The same input files were reused from the four-thread run. The pass-through
truth matrix also passed at one/24 threads (1,240 exact output checks), including
real read-through in FASTQ SE/PE, CRAM-to-FASTQ SE/PE, and native CRAM SE/PE.
Every enabled benchmark run asserted exact positive and negative detection
counts: 500,000 detections per applicable mate in positive data and none in
negative data. Disabled runs asserted zero scanned mates.

| Input | Read-through | State | Median s | M reads/s | Throughput change | Min–max s |
| --- | --- | --- | ---: | ---: | ---: | --- |
| FASTQ SE | 0% | baseline | 1.901 | 0.526 | 0.0% | 0.922–3.031 |
| FASTQ SE | 0% | disabled | 1.635 | 0.612 | 16.2% | 1.317–1.881 |
| FASTQ SE | 0% | enabled | 1.756 | 0.570 | 8.3% | 1.208–2.139 |
| FASTQ SE | 50% | baseline | 1.317 | 0.759 | 0.0% | 0.913–2.643 |
| FASTQ SE | 50% | disabled | 2.158 | 0.463 | -39.0% | 0.761–2.486 |
| FASTQ SE | 50% | enabled | 1.509 | 0.663 | -12.7% | 1.086–3.842 |
| FASTQ PE | 0% | baseline | 1.229 | 1.627 | 0.0% | 1.044–1.292 |
| FASTQ PE | 0% | disabled | 1.289 | 1.551 | -4.6% | 1.207–1.334 |
| FASTQ PE | 0% | enabled | 1.247 | 1.603 | -1.5% | 1.188–1.282 |
| FASTQ PE | 50% | baseline | 1.337 | 1.496 | 0.0% | 1.197–1.439 |
| FASTQ PE | 50% | disabled | 1.241 | 1.612 | 7.8% | 1.031–1.582 |
| FASTQ PE | 50% | enabled | 1.307 | 1.530 | 2.3% | 1.263–1.872 |
| CRAM SE | 0% | baseline | 1.997 | 0.501 | 0.0% | 1.849–2.919 |
| CRAM SE | 0% | disabled | 2.209 | 0.453 | -9.6% | 1.805–2.579 |
| CRAM SE | 0% | enabled | 2.830 | 0.353 | -29.4% | 2.010–2.959 |
| CRAM SE | 50% | baseline | 2.077 | 0.482 | 0.0% | 1.964–2.537 |
| CRAM SE | 50% | disabled | 2.396 | 0.417 | -13.3% | 1.680–2.978 |
| CRAM SE | 50% | enabled | 2.198 | 0.455 | -5.5% | 1.942–2.856 |
| CRAM PE | 0% | baseline | 4.027 | 0.497 | 0.0% | 3.339–4.575 |
| CRAM PE | 0% | disabled | 3.551 | 0.563 | 13.4% | 3.323–4.035 |
| CRAM PE | 0% | enabled | 3.573 | 0.560 | 12.7% | 3.475–4.797 |
| CRAM PE | 50% | baseline | 3.761 | 0.532 | 0.0% | 3.334–7.330 |
| CRAM PE | 50% | disabled | 4.035 | 0.496 | -6.8% | 3.394–5.323 |
| CRAM PE | 50% | enabled | 4.022 | 0.497 | -6.5% | 3.926–4.097 |

Raw [summary and binary hashes](benchmarks/pass-through-v0.6.0-24t.json) and
[all timed repetitions](benchmarks/pass-through-v0.6.0-24t.csv) are retained.

The 24-thread no-read-through median changes versus baseline were +8.3% FASTQ SE,
-1.5% FASTQ PE, -29.4% CRAM SE, and +12.7% CRAM PE. Wide overlapping ranges and
substantial disabled-build variation mean these are observed medians, not reliable
isolated scanner-cost estimates. In particular, the apparent positive changes
must not be advertised as speedups. The remaining SE scan tax is material. The
feature is default-on as requested; `--no-pass-through` restores proximal-only
behavior and completely bypasses detection. These one-sample direct-output
measurements do not establish performance for high-multiplex buffered workloads.

## Limitations

SE symmetry and PE opposite-end geometry are explicit assumptions. Incorrect
geometry requires an override or disabling detection. Complete distal coverage
and an exact anchor are necessary; partial targets, anchor substitutions, indels,
and observed Ns at known positions deliberately favor false negatives. A true
biological occurrence indistinguishable from the expected complete structure
cannot be distinguished by sequence alone. Hardware `perf` events were unavailable
under the host's `perf_event_paranoid=4`; assembly inspection guided removal of
verification setup from the negative scan, reordering its common branch, and
tracking clean-window validity through iterator position instead of advancing
two additional counters per base.

## Release artifact

Version/tag: `0.6.0` / `v0.6.0`, on `main`. The archive follows the existing
Linux GNU release convention and contains the optimized `plexless` executable:

`plexless-v0.6.0-x86_64-unknown-linux-gnu.tar.gz`

Archive SHA-256 (also supplied as the adjacent `.sha256` release asset):

```text
4ed37fdc073317d2fabdaae286a827e223a004343e750114c66480baac815aa9
```

Binary SHA-256, identical to both benchmarked production binaries:

```text
2f60da982ebc89ad617c71da65ee215c65a5da297d9f238bb3d2aae27237a462
```

The packaged executable successfully runs `--help`, `demux --help`, and
`--version`. It is an x86-64 Linux GNU dynamically linked binary, consistent
with the preceding release, not a new musl/static build. No compiled executable
or generated FASTQ/CRAM fixture belongs in normal Git history.
