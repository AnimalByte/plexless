# plexless validation suite

Keep these files under `scripts/validation/`.

- `make_plexless_validation_dataset.py`: injects deterministic A/B/T barcode
  prefixes into paired FASTQ input and writes SE + PE validation fixtures.
- `verify_plexless_summary.py`: checks aggregate category counts.
- `verify_plexless_outputs.py`: checks every output read ID against truth.
- `run_plexless_deep_validation.sh`: runs SE/PE validation across thread counts.

Typical generation command:

```bash
python3 scripts/validation/make_plexless_validation_dataset.py \
  --r1 ~/plexless_validation_data/raw_R1.fastq.gz \
  --r2 ~/plexless_validation_data/raw_R2.fastq.gz \
  --outdir ~/plexless_validation_data/plexless_validation
```

Typical validation command:

```bash
DATASET="$HOME/plexless_validation_data/plexless_validation" \
THREAD_COUNTS="1 2 4 8" \
scripts/validation/run_plexless_deep_validation.sh
```
