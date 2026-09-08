# plexless validation suite

Keep these files under `scripts/validation/`.

- `make_simplex_validation_dataset.py`: injects deterministic A/B/T barcode
  prefixes into paired FASTQ input and writes SE + PE validation fixtures.
- `verify_simplex_summary.py`: checks aggregate category counts.
- `verify_simplex_outputs.py`: checks every output read ID against truth.
- `run_simplex_deep_validation.sh`: runs SE/PE validation across thread counts.

Typical generation command:

```bash
python3 scripts/validation/make_simplex_validation_dataset.py \
  --r1 ~/simplex_validation_data/raw_R1.fastq.gz \
  --r2 ~/simplex_validation_data/raw_R2.fastq.gz \
  --outdir ~/simplex_validation_data/simplex_validation
```

Typical validation command:

```bash
DATASET="$HOME/simplex_validation_data/simplex_validation" \
THREAD_COUNTS="1 2 4 8" \
scripts/validation/run_simplex_deep_validation.sh
```
