#!/usr/bin/env bash
set -euo pipefail

PLEXLESS_BIN="${PLEXLESS_BIN:-./target/release/plexless}"
DATASET="${DATASET:-$HOME/plexless_validation_data/plexless_validation}"
THREAD_COUNTS="${THREAD_COUNTS:-1 2 4 8}"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VERIFY_SUMMARY="$ROOT/verify_plexless_summary.py"
VERIFY_OUTPUTS="$ROOT/verify_plexless_outputs.py"

EXPECTED="$DATASET/expected_summary.tsv"
BARCODES="$DATASET/barcodes.tsv"
SAMPLES="$DATASET/samples.tsv"

if [[ ! -x "$PLEXLESS_BIN" ]]; then
    echo "plexless binary not executable: $PLEXLESS_BIN" >&2
    exit 1
fi

for mode in se pe; do
    echo
    echo "======================================"
    echo "DEEP VALIDATION MODE: $mode"
    echo "======================================"

    baseline=""

    for threads in $THREAD_COUNTS; do
        out="deep_validation_${mode}_${threads}t"
        log="deep_validation_${mode}_${threads}t.stderr"

        rm -rf "$out" "$log"

        echo
        echo "--- threads=$threads ---"

        if [[ "$mode" == "se" ]]; then
            "$PLEXLESS_BIN" \
                --threads "$threads" \
                demux \
                --reads "$DATASET/se.fastq.gz" \
                --structure R1_4A4B2T \
                --barcodes "$BARCODES" \
                --samples "$SAMPLES" \
                --write-unassigned \
                --output "$out" \
                2> >(tee "$log" >&2)
        else
            "$PLEXLESS_BIN" \
                --threads "$threads" \
                demux \
                --r1 "$DATASET/pe_R1.fastq.gz" \
                --r2 "$DATASET/pe_R2.fastq.gz" \
                --r1-structure R1_2A2B2T \
                --r2-structure R2_2A2B2T \
                --barcodes "$BARCODES" \
                --samples "$SAMPLES" \
                --write-unassigned \
                --output "$out" \
                2> >(tee "$log" >&2)
        fi

        if [[ -e "$out/PLEXLESS_INCOMPLETE" ]]; then
            echo "Successful run retained $out/PLEXLESS_INCOMPLETE" >&2
            exit 1
        fi

        python3 "$VERIFY_SUMMARY" "$EXPECTED" "$log"

        if [[ "$threads" == "1" ]]; then
            python3 "$VERIFY_OUTPUTS" \
                --dataset "$DATASET" \
                --output "$out" \
                --mode "$mode"
            baseline="$out"
        else
            if [[ -z "$baseline" ]]; then
                echo "The thread list must include 1 first so a baseline exists." >&2
                exit 1
            fi

            python3 "$VERIFY_OUTPUTS" \
                --dataset "$DATASET" \
                --output "$out" \
                --mode "$mode" \
                --baseline "$baseline"
        fi
    done
done

echo
echo "All deep SE/PE validation runs passed."
