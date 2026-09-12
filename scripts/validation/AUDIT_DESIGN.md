# Production run-manifest and audit design

This is a design recommendation, not an implemented CLI or hot-path change.
The qualification work deliberately leaves production behavior unchanged until
the costs and artifact contract are agreed.

## Minimal always-on manifest

A successful run should eventually publish one small, versioned
`plexless_run_manifest.json` containing:

- schema version and Plexless version;
- normalized input filenames and SE/PE mode;
- read structures, mismatch limit, output strategy, compression level, thread
  budget, and relevant QC/output options;
- SHA-256 of the barcode catalog and sample sheet;
- total, assigned, unmatched, ambiguous, unrouted, short-read, orphan-R1, and
  orphan-R2 counts;
- every worklist sample's assigned fragment and mate counts;
- explicit `status: "complete"`.

The manifest should be written to a temporary file, flushed, and atomically
renamed only after writers, count reconciliation, and QC reports succeed. Until
then, `PLEXLESS_INCOMPLETE` remains the authoritative state. Publishing the
manifest and removing the marker must be ordered so an interrupted run can
never look complete.

The small configuration files can be read and SHA-256 hashed during startup.
The other proposed fields already exist in memory or CLI state. This design
adds no per-read locking and should be cheap enough for normal production, but
it still needs an implementation-specific benchmark before becoming default.

## Hash and digest cost classes

### A. Suitable for normal production

- Configuration SHA-256: small, bounded startup reads with useful provenance.
- Existing classification, per-sample, mate, batch, chunk, and unassigned
  reconciliation: already always-on and measured as part of the current
  implementation.
- Atomic successful-completion manifest: a tiny final serialization and write.

Input filenames are provenance, not content identity. Store normalized paths
without implying that they prove immutable input.

### B. Qualification or explicit `--audit`

- Input SHA-256: reading compressed input bytes adds a full I/O pass unless the
  digest is carefully integrated below the decompressor. A digest of compressed
  bytes identifies the artifact but not equivalent recompressions.
- Compressed output SHA-256: cheap to stream at the final writer relative to
  output I/O, but member boundaries and gzip headers make it unsuitable as the
  principal cross-mode equivalence value.
- Logical/decompressed output SHA-256: the best portable output-content
  identity, but an after-the-fact implementation rereads and decompresses all
  output. The current external verifier does exactly this so normal runs pay
  nothing.
- Deterministic assignment digest: valuable for detecting routing changes even
  when output writing is disabled, but it touches every fragment and needs a
  stable canonical schema.

## Assignment digest design

If implemented, each work item should contribute a canonical length-delimited
record containing at least:

```text
input ordinal | normalized read ID | mate state | route result |
sample ID when assigned | R1 trim length | R2 trim length
```

Workers should hash records inside each existing batch. The ordered aggregator
should combine `(batch_id, batch_digest)` in batch-ID order. This avoids one
global mutex-protected hasher in the per-read hot path and makes direct and
buffered output strategies comparable. The schema, byte order, string
normalization, and treatment of orphans must be versioned before the value is a
durable audit artifact.

SHA-256 is recommended for published configuration, input, output, and
assignment artifacts because it is widely supported and independently
auditable. The pure-Rust `sha2` crate would be a focused dependency if this is
implemented. A faster non-cryptographic hash could reduce CPU cost for an
internal race detector, but should not be labeled SHA-256 or used as the sole
portable provenance digest.

## Current qualification boundary

`verify_plexless_outputs.py` writes per-run logical hash manifests externally.
It hashes decompressed FASTQ, verifies every record against content-rich truth,
compares deterministic QC reports, and independently invokes `gzip -t` when
available. This intentionally belongs outside the production hot path.

Because this task adds no production integrity feature, normal-run throughput
overhead from these changes is exactly zero by construction: the compiled Rust
binary is unchanged. `run_metrics.tsv` measures the current production path;
qualification/verifier work runs after that timer and is not misreported as
Plexless throughput. A future manifest implementation must be benchmarked
against the same generated standard fixture before it is enabled by default.
