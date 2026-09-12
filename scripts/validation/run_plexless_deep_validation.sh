#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PLEXLESS_BIN="${PLEXLESS_BIN:-./target/release/plexless}"
DATASET="${DATASET:-$HOME/plexless_validation_data/plexless_validation}"
THREAD_COUNTS="${THREAD_COUNTS:-1 2 4 8}"
OUTPUT_MODES="${OUTPUT_MODES:-direct buffered}"
REPEAT="${REPEAT:-1}"
PRESET="${PRESET:-normal}"
RESULTS_ROOT="${RESULTS_ROOT:-}"
RUN_NEGATIVE="${RUN_NEGATIVE:-1}"
INCLUDE_CRAM="${INCLUDE_CRAM:-0}"
CRAM_ONLY="${CRAM_ONLY:-0}"
VERIFY_SUMMARY="$ROOT/verify_plexless_summary.py"
VERIFY_OUTPUTS="$ROOT/verify_plexless_outputs.py"

usage() {
    echo "Usage: $0 [--preset ci|normal|desktop|stress] [--dataset DIR]" >&2
    echo "          [--threads '1 2 4 8'] [--output-modes 'direct buffered']" >&2
    echo "          [--repeat N] [--results DIR] [--skip-negative] [--include-cram|--cram-only]" >&2
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --preset) PRESET="$2"; shift 2 ;;
        --dataset) DATASET="$2"; shift 2 ;;
        --threads) THREAD_COUNTS="$2"; shift 2 ;;
        --output-modes) OUTPUT_MODES="$2"; shift 2 ;;
        --repeat) REPEAT="$2"; shift 2 ;;
        --results) RESULTS_ROOT="$2"; shift 2 ;;
        --skip-negative) RUN_NEGATIVE=0; shift ;;
        --include-cram) INCLUDE_CRAM=1; shift ;;
        --cram-only) INCLUDE_CRAM=1; CRAM_ONLY=1; RUN_NEGATIVE=0; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 2 ;;
    esac
done

case "$PRESET" in
    ci)
        [[ "$THREAD_COUNTS" != "1 2 4 8" ]] || THREAD_COUNTS="1 4"
        ;;
    normal) ;;
    desktop)
        [[ "$THREAD_COUNTS" != "1 2 4 8" ]] || THREAD_COUNTS="1 2 4 8 12 16 24"
        ;;
    stress)
        [[ "$THREAD_COUNTS" != "1 2 4 8" ]] || THREAD_COUNTS="1 24"
        ;;
    *) echo "Invalid preset: $PRESET" >&2; exit 2 ;;
esac

if [[ ! "$REPEAT" =~ ^[1-9][0-9]*$ ]]; then
    echo "--repeat must be a positive integer" >&2
    exit 2
fi
if [[ ! -x "$PLEXLESS_BIN" ]]; then
    echo "plexless binary not executable: $PLEXLESS_BIN" >&2
    exit 1
fi
if [[ ! -f "$DATASET/metadata.json" ]]; then
    echo "not a generated validation dataset: $DATASET" >&2
    exit 1
fi

read -r -a THREAD_ARRAY <<< "$THREAD_COUNTS"
read -r -a MODE_ARRAY <<< "$OUTPUT_MODES"
if [[ "${THREAD_ARRAY[0]:-}" != "1" ]]; then
    echo "thread list must start with 1 for the conservative reference" >&2
    exit 2
fi
for mode in "${MODE_ARRAY[@]}"; do
    if [[ "$mode" != "direct" && "$mode" != "buffered" ]]; then
        echo "invalid output mode: $mode" >&2
        exit 2
    fi
done
if [[ ! " ${MODE_ARRAY[*]} " =~ " direct " ]]; then
    echo "output modes must include direct for the 1-thread reference" >&2
    exit 2
fi

max_threads="${THREAD_ARRAY[${#THREAD_ARRAY[@]} - 1]}"
if [[ -z "$RESULTS_ROOT" ]]; then
    stamp="$(date -u +%Y%m%dT%H%M%SZ)"
    RESULTS_ROOT="$DATASET/qualification-$PRESET-$stamp"
fi
if [[ -e "$RESULTS_ROOT" ]]; then
    echo "results path already exists; choose a new --results path: $RESULTS_ROOT" >&2
    exit 1
fi
mkdir -p "$RESULTS_ROOT"
METRICS="$RESULTS_ROOT/run_metrics.tsv"
printf 'fixture\tmode\tthreads\toutput_mode\trepetition\tseconds\treads_or_events_per_second\taggregate_cpu\tpeak_rss_kib\toutput_bytes\texternal_validation_seconds\n' > "$METRICS"

fixture_args() {
    local fixture="$1"
    local mode="$2"
    local directory="$DATASET/$fixture"
    FIXTURE_CLI=(
        --barcodes "$directory/barcodes.tsv"
        --samples "$directory/samples.tsv"
        --write-unassigned
        --fastq-stats
    )
    case "$fixture:$mode" in
        core:se)
            FIXTURE_CLI+=(--reads "$directory/se.fastq.gz" --structure R1_6A6B2T --max-mismatches 1)
            ;;
        core:pe)
            FIXTURE_CLI+=(
                --r1 "$directory/pe_R1.fastq.gz"
                --r2 "$directory/pe_R2.fastq.gz"
                --r1-structure R1_3A3B2T
                --r2-structure 'R2_3A(rc)3B(rc)2T'
                --max-mismatches 1
            )
            ;;
        n2:se)
            FIXTURE_CLI+=(--reads "$directory/se.fastq.gz" --structure R1_6A2T --max-mismatches 2)
            ;;
        r1_only:pe)
            FIXTURE_CLI+=(
                --r1 "$directory/pe_R1.fastq.gz"
                --r2 "$directory/pe_R2.fastq.gz"
                --r1-structure R1_6A2T
                --max-mismatches 1
            )
            ;;
        r2_only:pe)
            FIXTURE_CLI+=(
                --r1 "$directory/pe_R1.fastq.gz"
                --r2 "$directory/pe_R2.fastq.gz"
                --r2-structure 'R2_6A(rc)2T'
                --max-mismatches 1
            )
            ;;
        *) echo "unsupported fixture/mode: $fixture:$mode" >&2; exit 2 ;;
    esac
}

cram_fixture_args() {
    local fixture="$1"
    local mode="$2"
    local output_format="$3"
    local directory="$DATASET/$fixture"
    local read_mode
    if [[ "$mode" == se ]]; then
        read_mode=single
    else
        read_mode=paired
    fi
    FIXTURE_CLI=(
        --cram "$directory/$mode.cram"
        --read-mode "$read_mode"
        --output-format "$output_format"
        --barcodes "$directory/barcodes.tsv"
        --samples "$directory/samples.tsv"
        --write-unassigned
        --fastq-stats
    )
    case "$fixture:$mode" in
        core:se) FIXTURE_CLI+=(--structure R1_6A6B2T --max-mismatches 1) ;;
        core:pe)
            FIXTURE_CLI+=(
                --r1-structure R1_3A3B2T
                --r2-structure 'R2_3A(rc)3B(rc)2T'
                --max-mismatches 1
            )
            ;;
        n2:se) FIXTURE_CLI+=(--structure R1_6A2T --max-mismatches 2) ;;
        r1_only:pe) FIXTURE_CLI+=(--r1-structure R1_6A2T --max-mismatches 1) ;;
        r2_only:pe) FIXTURE_CLI+=(--r2-structure 'R2_6A(rc)2T' --max-mismatches 1) ;;
        *) echo "unsupported CRAM fixture/mode: $fixture:$mode" >&2; exit 2 ;;
    esac
}

measure_run() {
    local time_file="$1"
    local log_file="$2"
    shift 2
    local status
    set +e
    if [[ -x /usr/bin/time ]]; then
        /usr/bin/time -f '%e\t%P\t%M' -o "$time_file" "$@" 2> "$log_file"
        status=$?
    else
        local started ended
        started="$(date +%s%N)"
        "$@" 2> "$log_file"
        status=$?
        ended="$(date +%s%N)"
        python3 - "$started" "$ended" > "$time_file" <<'PY'
import sys
print(f"{(int(sys.argv[2]) - int(sys.argv[1])) / 1_000_000_000:.6f}\tNA\tNA")
PY
    fi
    set -e
    cat "$log_file" >&2
    return "$status"
}

elapsed_seconds() {
    python3 - "$1" "$2" <<'PY'
import sys
print(f"{(int(sys.argv[2]) - int(sys.argv[1])) / 1_000_000_000:.6f}")
PY
}

run_positive() {
    local fixture="$1" mode="$2" threads="$3" output_mode="$4" repetition="$5" baseline="$6"
    local label="${fixture}_${mode}_${threads}t_${output_mode}_r${repetition}"
    local output="$RESULTS_ROOT/$label"
    local log="$RESULTS_ROOT/$label.stderr"
    local timing="$RESULTS_ROOT/$label.time.tsv"
    local hashes="$RESULTS_ROOT/$label.logical-sha256.json"
    fixture_args "$fixture" "$mode"
    echo
    echo "--- $label ---"
    measure_run "$timing" "$log" \
        "$PLEXLESS_BIN" --threads "$threads" demux \
        "${FIXTURE_CLI[@]}" --output-mode "$output_mode" --output "$output"
    if [[ -e "$output/PLEXLESS_INCOMPLETE" ]]; then
        echo "successful run retained PLEXLESS_INCOMPLETE: $output" >&2
        exit 1
    fi
    local validation_started validation_seconds
    validation_started="$(date +%s%N)"
    python3 "$VERIFY_SUMMARY" \
        "$DATASET/$fixture/expected_summary_${mode}.tsv" "$log"
    verify=(
        python3 "$VERIFY_OUTPUTS"
        --dataset "$DATASET/$fixture"
        --output "$output"
        --mode "$mode"
        --expect-fastq-stats
        --hash-manifest "$hashes"
    )
    if [[ -n "$baseline" ]]; then
        verify+=(--baseline "$baseline")
    fi
    "${verify[@]}"
    validation_seconds="$(elapsed_seconds "$validation_started" "$(date +%s%N)")"

    local seconds cpu peak total output_bytes rate
    IFS=$'\t' read -r seconds cpu peak < "$timing"
    total="$(awk -F '\t' '$1 == "total" {print $2}' "$DATASET/$fixture/expected_summary_${mode}.tsv")"
    output_bytes="$(find "$output" -maxdepth 1 -type f -printf '%s\n' | awk '{sum += $1} END {print sum + 0}')"
    rate="$(python3 - "$total" "$seconds" <<'PY'
import sys
total, seconds = int(sys.argv[1]), float(sys.argv[2])
print(f"{total / seconds:.2f}" if seconds else "inf")
PY
)"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$fixture" "$mode" "$threads" "$output_mode" "$repetition" \
        "$seconds" "$rate" "$cpu" "$peak" "$output_bytes" "$validation_seconds" >> "$METRICS"
    LAST_OUTPUT="$output"
}

run_matrix() {
    local fixture="$1" mode="$2"
    echo
    echo "======================================"
    echo "QUALIFICATION: $fixture $mode"
    echo "======================================"
    run_positive "$fixture" "$mode" 1 direct 1 ""
    local baseline="$LAST_OUTPUT"
    for threads in "${THREAD_ARRAY[@]}"; do
        for output_mode in "${MODE_ARRAY[@]}"; do
            if [[ "$threads" == 1 && "$output_mode" == direct ]]; then
                continue
            fi
            local repetitions=1
            if [[ "$threads" == "$max_threads" ]]; then
                repetitions="$REPEAT"
            fi
            for ((repetition = 1; repetition <= repetitions; repetition++)); do
                run_positive "$fixture" "$mode" "$threads" "$output_mode" "$repetition" "$baseline"
            done
        done
    done
}

run_focused() {
    local fixture="$1" mode="$2"
    echo
    echo "======================================"
    echo "FOCUSED SEMANTICS: $fixture $mode"
    echo "======================================"
    run_positive "$fixture" "$mode" 1 direct 1 ""
    local baseline="$LAST_OUTPUT"
    if [[ "$max_threads" != 1 ]]; then
        for output_mode in "${MODE_ARRAY[@]}"; do
            run_positive "$fixture" "$mode" "$max_threads" "$output_mode" 1 "$baseline"
        done
    elif [[ " ${MODE_ARRAY[*]} " =~ " buffered " ]]; then
        run_positive "$fixture" "$mode" 1 buffered 1 "$baseline"
    fi
}

run_cram_fastq_positive() {
    local fixture="$1" mode="$2" threads="$3" output_mode="$4" repetition="$5" baseline="$6"
    local label="${fixture}_${mode}_cram_to_fastq_${threads}t_${output_mode}_r${repetition}"
    local output="$RESULTS_ROOT/$label"
    local log="$RESULTS_ROOT/$label.stderr"
    local timing="$RESULTS_ROOT/$label.time.tsv"
    local hashes="$RESULTS_ROOT/$label.logical-sha256.json"
    cram_fixture_args "$fixture" "$mode" fastq
    echo
    echo "--- $label ---"
    measure_run "$timing" "$log" \
        "$PLEXLESS_BIN" --threads "$threads" demux \
        "${FIXTURE_CLI[@]}" --output-mode "$output_mode" --output "$output"
    if [[ -e "$output/PLEXLESS_INCOMPLETE" ]]; then
        echo "successful CRAM input run retained PLEXLESS_INCOMPLETE: $output" >&2
        exit 1
    fi
    local validation_started validation_seconds
    validation_started="$(date +%s%N)"
    python3 "$VERIFY_SUMMARY" \
        "$DATASET/$fixture/expected_summary_${mode}.tsv" "$log"
    verify=(
        python3 "$VERIFY_OUTPUTS"
        --dataset "$DATASET/$fixture"
        --output "$output"
        --mode "$mode"
        --cram-qnames
        --expect-fastq-stats
        --hash-manifest "$hashes"
    )
    if [[ -n "$baseline" ]]; then
        verify+=(--baseline "$baseline")
    fi
    "${verify[@]}"
    validation_seconds="$(elapsed_seconds "$validation_started" "$(date +%s%N)")"

    local seconds cpu peak total output_bytes rate
    IFS=$'\t' read -r seconds cpu peak < "$timing"
    total="$(awk -F '\t' '$1 == "total" {print $2}' "$DATASET/$fixture/expected_summary_${mode}.tsv")"
    output_bytes="$(find "$output" -maxdepth 1 -type f -printf '%s\n' | awk '{sum += $1} END {print sum + 0}')"
    rate="$(python3 - "$total" "$seconds" <<'PY'
import sys
total, seconds = int(sys.argv[1]), float(sys.argv[2])
print(f"{total / seconds:.2f}" if seconds else "inf")
PY
)"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$fixture" "cram_to_fastq_$mode" "$threads" "$output_mode" "$repetition" \
        "$seconds" "$rate" "$cpu" "$peak" "$output_bytes" "$validation_seconds" >> "$METRICS"
    LAST_OUTPUT="$output"
}

run_cram_fastq_matrix() {
    local fixture="$1" mode="$2"
    echo
    echo "======================================"
    echo "CRAM -> FASTQ QUALIFICATION: $fixture $mode"
    echo "======================================"
    run_cram_fastq_positive "$fixture" "$mode" 1 direct 1 ""
    local baseline="$LAST_OUTPUT"
    for threads in "${THREAD_ARRAY[@]}"; do
        for output_mode in "${MODE_ARRAY[@]}"; do
            if [[ "$threads" == 1 && "$output_mode" == direct ]]; then
                continue
            fi
            local repetitions=1
            if [[ "$threads" == "$max_threads" ]]; then
                repetitions="$REPEAT"
            fi
            for ((repetition = 1; repetition <= repetitions; repetition++)); do
                run_cram_fastq_positive "$fixture" "$mode" "$threads" "$output_mode" "$repetition" "$baseline"
            done
        done
    done
}

run_cram_fastq_focused() {
    local fixture="$1" mode="$2"
    run_cram_fastq_positive "$fixture" "$mode" 1 direct 1 ""
    local baseline="$LAST_OUTPUT"
    if [[ "$max_threads" != 1 ]]; then
        run_cram_fastq_positive "$fixture" "$mode" "$max_threads" direct 1 "$baseline"
    fi
}

run_cram_output_positive() {
    local fixture="$1" mode="$2" threads="$3" repetition="$4"
    local label="${fixture}_${mode}_cram_to_cram_${threads}t_r${repetition}"
    local output="$RESULTS_ROOT/$label"
    local decoded="$RESULTS_ROOT/$label-decoded-fastq"
    local log="$RESULTS_ROOT/$label.stderr"
    local timing="$RESULTS_ROOT/$label.time.tsv"
    local hashes="$RESULTS_ROOT/$label.logical-sha256.json"
    cram_fixture_args "$fixture" "$mode" cram
    echo
    echo "--- $label ---"
    measure_run "$timing" "$log" \
        "$PLEXLESS_BIN" --threads "$threads" demux "${FIXTURE_CLI[@]}" --output "$output"
    if [[ -e "$output/PLEXLESS_INCOMPLETE" ]]; then
        echo "successful CRAM output run retained PLEXLESS_INCOMPLETE: $output" >&2
        exit 1
    fi
    local validation_started validation_seconds
    validation_started="$(date +%s%N)"
    if [[ ! -f "$output/cram_metadata.tsv" ]] ||
        [[ "$(head -n 1 "$output/cram_metadata.tsv")" != $'Action\tTag\tCount' ]]; then
        echo "CRAM output run did not write a valid metadata action report: $output" >&2
        exit 1
    fi
    python3 "$VERIFY_SUMMARY" \
        "$DATASET/$fixture/expected_summary_${mode}.tsv" "$log"
    cargo run --quiet --release --example extract_cram_outputs_to_fastq -- \
        "$([[ "$mode" == se ]] && echo single || echo paired)" "$output" "$decoded"
    python3 "$VERIFY_OUTPUTS" \
        --dataset "$DATASET/$fixture" \
        --output "$output" \
        --fastq-directory "$decoded" \
        --mode "$mode" \
        --cram-qnames \
        --expect-fastq-stats \
        --hash-manifest "$hashes"
    if command -v samtools >/dev/null 2>&1; then
        mapfile -t cram_outputs < <(find "$output" -maxdepth 1 -type f -name '*.cram' -print | sort)
        samtools quickcheck -u -v "${cram_outputs[@]}"
        echo "PASS  independent samtools quickcheck      ${#cram_outputs[@]} files"
    else
        echo "SKIP  independent samtools quickcheck (samtools is not installed)"
    fi
    validation_seconds="$(elapsed_seconds "$validation_started" "$(date +%s%N)")"

    local seconds cpu peak total output_bytes rate
    IFS=$'\t' read -r seconds cpu peak < "$timing"
    total="$(awk -F '\t' '$1 == "total" {print $2}' "$DATASET/$fixture/expected_summary_${mode}.tsv")"
    output_bytes="$(find "$output" -maxdepth 1 -type f -printf '%s\n' | awk '{sum += $1} END {print sum + 0}')"
    rate="$(python3 - "$total" "$seconds" <<'PY'
import sys
total, seconds = int(sys.argv[1]), float(sys.argv[2])
print(f"{total / seconds:.2f}" if seconds else "inf")
PY
)"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$fixture" "cram_to_cram_$mode" "$threads" native "$repetition" \
        "$seconds" "$rate" "$cpu" "$peak" "$output_bytes" "$validation_seconds" >> "$METRICS"
}

run_cram_output_matrix() {
    local fixture="$1" mode="$2"
    echo
    echo "======================================"
    echo "CRAM -> CRAM QUALIFICATION: $fixture $mode"
    echo "======================================"
    for threads in "${THREAD_ARRAY[@]}"; do
        local repetitions=1
        if [[ "$threads" == "$max_threads" ]]; then
            repetitions="$REPEAT"
        fi
        for ((repetition = 1; repetition <= repetitions; repetition++)); do
            run_cram_output_positive "$fixture" "$mode" "$threads" "$repetition"
        done
    done
}

expect_failure() {
    local name="$1" marker_expected="$2" threads="$3"
    shift 3
    local output="$RESULTS_ROOT/negative_$name"
    local log="$RESULTS_ROOT/negative_$name.stderr"
    local timing="$RESULTS_ROOT/negative_$name.time.tsv"
    echo
    echo "--- expected failure: $name ---"
    if measure_run "$timing" "$log" "$PLEXLESS_BIN" --threads "$threads" demux "$@" --output "$output"; then
        echo "negative case unexpectedly succeeded: $name" >&2
        exit 1
    fi
    if grep -q '^Demultiplexing complete$' "$log"; then
        echo "negative case printed a success summary: $name" >&2
        exit 1
    fi
    if [[ "$marker_expected" == yes ]]; then
        if [[ ! -f "$output/PLEXLESS_INCOMPLETE" ]]; then
            echo "operational failure lacks PLEXLESS_INCOMPLETE: $name" >&2
            exit 1
        fi
        if [[ -e "$output/sample_metrics.tsv" ]]; then
            echo "operational failure wrote final sample metrics: $name" >&2
            exit 1
        fi
    elif [[ -e "$output/PLEXLESS_INCOMPLETE" ]]; then
        echo "startup validation unexpectedly created an incomplete marker: $name" >&2
        exit 1
    fi
    echo "PASS  expected failure: $name"
}

run_negative_suite() {
    local negative="$DATASET/negative"
    local n2="$DATASET/n2"
    local common=(--structure R1_6A2T --barcodes "$negative/barcodes.tsv" --samples "$negative/samples.tsv" --output-mode direct --write-unassigned)
    echo
    echo "======================================"
    echo "NEGATIVE / COMPLETION QUALIFICATION"
    echo "======================================"
    expect_failure truncated_fastq yes 1 --reads "$negative/truncated.fastq" "${common[@]}"
    expect_failure length_mismatch yes 1 --reads "$negative/length_mismatch.fastq" "${common[@]}"
    expect_failure corrupt_gzip yes 8 --reads "$negative/corrupt.fastq.gz" "${common[@]}"
    expect_failure desync_beyond_window yes 8 \
        --r1 "$negative/desync_R1.fastq.gz" \
        --r2 "$negative/desync_R2.fastq.gz" \
        --r1-structure R1_6A2T \
        --barcodes "$negative/barcodes.tsv" \
        --samples "$negative/samples.tsv" \
        --max-mismatches 1 --output-mode buffered --write-unassigned
    expect_failure unsafe_whitelist no 1 \
        --reads "$n2/se.fastq.gz" --structure R1_6A2T \
        --barcodes "$negative/unsafe_barcodes.tsv" \
        --samples "$negative/unsafe_samples.tsv" --max-mismatches 1
    expect_failure invalid_barcode_character no 1 \
        --reads "$n2/se.fastq.gz" --structure R1_6A2T \
        --barcodes "$negative/invalid_barcode.tsv" \
        --samples "$negative/samples.tsv" --max-mismatches 1
    expect_failure malformed_sample_sheet no 1 \
        --reads "$n2/se.fastq.gz" --structure R1_6A2T \
        --barcodes "$negative/barcodes.tsv" \
        --samples "$negative/malformed_samples.tsv" --max-mismatches 1
}

if [[ "$CRAM_ONLY" == 0 ]]; then
    run_matrix core se
    run_matrix core pe
    run_focused n2 se
    run_focused r1_only pe
    run_focused r2_only pe
fi

if [[ "$INCLUDE_CRAM" == 1 ]]; then
    for fixture_mode in core:se core:pe n2:se r1_only:pe r2_only:pe; do
        fixture="${fixture_mode%%:*}"
        mode="${fixture_mode##*:}"
        if [[ ! -f "$DATASET/$fixture/$mode.cram" ]]; then
            echo "missing CRAM fixture (regenerate with --with-cram): $DATASET/$fixture/$mode.cram" >&2
            exit 1
        fi
    done
    run_cram_fastq_matrix core se
    run_cram_fastq_matrix core pe
    run_cram_fastq_focused n2 se
    run_cram_fastq_focused r1_only pe
    run_cram_fastq_focused r2_only pe
    run_cram_output_matrix core se
    run_cram_output_matrix core pe
fi

if [[ "$RUN_NEGATIVE" == 1 ]]; then
    run_negative_suite
fi

echo
echo "All requested validation runs passed."
echo "Results: $RESULTS_ROOT"
echo "Metrics: $METRICS"
