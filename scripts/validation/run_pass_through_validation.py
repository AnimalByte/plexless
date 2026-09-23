#!/usr/bin/env python3
"""Deterministic pass-through truth validation and paired baseline timings.

Uses the repository CRAM fixture converter. Timings include process startup,
decoding, demultiplexing, trimming, and gzip writing; verification is untimed.
All large inputs/outputs stay in the explicitly supplied new work directory.
"""
import argparse
import csv
import gzip
import hashlib
import json
import platform
import random
import shutil
import statistics
import subprocess
import time
from pathlib import Path

A = "ACGTTGCACTGATCGAGGTACCTA"
B = "GCTACAGTTCGATACCGTGATCAG"
SIZES = (50, 75, 100, 150, 200, 250, 300, 400, 500)
REACH = 350


def rc(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def quality(length):
    return "".join(chr(33 + i % 60) for i in range(length))


def record(name, seq):
    return f"@{name}\n{seq}\n+\n{quality(len(seq))}\n"


def fixture(root, paired, positive, count, validation=False):
    rng = random.Random(71029)
    root.mkdir()
    (root / "barcodes.tsv").write_text(f"Set\tID\tSequence\nA\ta\t{A}\nB\tb\t{B}\n")
    (root / "samples.tsv").write_text("Sample\tA\tB\ns\ta\tb\n" if paired else "Sample\tA\ns\ta\n")
    known = [A, B] if paired else [A]
    expected = [[] for _ in known]
    handles = [(root / f"r{mate + 1}.fastq").open("w") for mate in range(len(known))]
    # Negative controls include known exact substrings with deliberately wrong
    # surrounding bases, partial motifs, and barcode-like sequence.
    for index in range(count):
        for mate, proximal in enumerate(known):
            opposite = known[1 - mate] if paired else proximal
            distal = "NNN" + rc(opposite)
            if validation:
                insert_len = SIZES[index % len(SIZES)]
                bio = "".join(rng.choices("ACGT", k=insert_len))
                if index >= len(SIZES):
                    # Each control retains an exact 16-base anchor-like motif
                    # but never the complete assigned distal structure.
                    bio = "G" * 30 + rc(opposite)[:16] + "C" * 30 + proximal[:8] + "T" * 30
                    suffix = bio
                    clean = suffix
                else:
                    suffix = (bio + distal + "A" * REACH)[:REACH]
                    clean = bio if insert_len + len(distal) <= REACH else suffix
            else:
                # Equal raw lengths in positive/negative benchmarks. 50% of
                # molecules reach distal structure; each PE mate is independent.
                insert_len = 100 + mate * 35 if positive and (index + mate) % 2 == 0 else 500
                bio = "".join(rng.choices("ACGT", k=insert_len))
                suffix = (bio + distal + "A" * REACH)[:REACH]
                clean = bio if insert_len + len(distal) <= REACH else suffix
            seq = proximal + "NNN" + suffix
            name = f"r{index:09d}" + (f"/{mate + 1}" if paired else "")
            handles[mate].write(record(name, seq))
            if validation:
                start = len(proximal) + 3
                expected[mate].append(f"@{name}\n{clean}\n+\n{quality(len(seq))[start:start + len(clean)]}\n")
    for handle in handles:
        handle.close()
    return ["".join(records) for records in expected]


def command(binary, root, output, paired, threads, disabled=False, cram=False):
    cmd = [str(binary), "--threads", str(threads), "demux", "--barcodes", str(root / "barcodes.tsv"),
           "--samples", str(root / "samples.tsv"), "--output", str(output), "--max-mismatches", "0"]
    if cram:
        cmd += ["--cram", str(root / "input.cram"), "--read-mode", "paired" if paired else "single", "--output-format", "fastq"]
    elif paired:
        cmd += ["--r1", str(root / "r1.fastq"), "--r2", str(root / "r2.fastq")]
    else:
        cmd += ["--reads", str(root / "r1.fastq")]
    if paired:
        cmd += ["--r1-structure", "R1_24A3T", "--r2-structure", "R2_24B3T"]
    else:
        cmd += ["--structure", "R1_24A3T"]
    if disabled:
        cmd += ["--no-pass-through"]
    return cmd


def digest_outputs(output):
    digest = hashlib.sha256()
    for path in sorted(output.glob("*.fastq.gz")):
        with gzip.open(path, "rb") as stream:
            while block := stream.read(1024 * 1024):
                digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=Path("target/release/plexless"))
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--records", type=int, default=200000)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--with-cram", action="store_true")
    args = parser.parse_args()
    if min(args.records, args.repeats, args.threads) < 1:
        parser.error("records, repeats, and threads must be positive")
    args.binary = args.binary.resolve()
    if args.baseline:
        args.baseline = args.baseline.resolve()
    root = args.work_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    converter = args.binary.parent / "examples/make_unmapped_cram_fixture"
    validation_cases = 0
    for paired in [False, True]:
        dataset = root / f"validation-{paired}"
        expected = fixture(dataset, paired, True, 18, validation=True)
        for threads in [1, args.threads]:
            for mode in ["direct", "buffered"]:
                output = root / "validation-output"
                cmd = command(args.binary, dataset, output, paired, threads) + ["--output-mode", mode]
                result = subprocess.run(cmd, capture_output=True, check=True)
                (root / f"validation-{paired}-{threads}-{mode}.stderr").write_bytes(result.stderr)
                for mate, truth in enumerate(expected):
                    filename = f"s_R{mate + 1}.fastq.gz" if paired else "s.fastq.gz"
                    with gzip.open(output / filename, "rt") as stream:
                        actual = stream.read()
                    if actual != truth:
                        raise AssertionError(f"Exact sequence/quality/name truth mismatch: {paired=} {threads=} {mode=} {mate=}")
                    validation_cases += len(SIZES) * 2
                shutil.rmtree(output)
    print(f"Exact synthetic truth checks passed: {validation_cases}", flush=True)
    if not args.baseline:
        return
    rows = []
    for paired in [False, True]:
        for positive in [False, True]:
            dataset = root / f"bench-{paired}-{positive}"
            fixture(dataset, paired, positive, args.records)
            if args.with_cram:
                cmd = [str(converter), "paired" if paired else "single", str(dataset / "r1.fastq")]
                if paired:
                    cmd.append(str(dataset / "r2.fastq"))
                subprocess.run(cmd + [str(dataset / "input.cram")], check=True, capture_output=True)
            for cram in [False, True] if args.with_cram else [False]:
                states = [("baseline", args.baseline, False), ("disabled", args.binary, True), ("enabled", args.binary, False)]
                hashes = {}
                # Untimed warmup per state, then rotate state order each round.
                for repetition in range(-1, args.repeats):
                    for offset in range(3):
                        state, binary, disabled = states[(offset + max(repetition, 0)) % 3]
                        output = root / "benchmark-output"
                        cmd = command(binary, dataset, output, paired, args.threads, disabled, cram)
                        start = time.perf_counter()
                        result = subprocess.run(cmd, capture_output=True, check=True)
                        seconds = time.perf_counter() - start
                        name = f"{'cram' if cram else 'fastq'}-{'pe' if paired else 'se'}-{positive}-{state}-{repetition}"
                        (root / f"{name}.stderr").write_bytes(result.stderr)
                        if repetition == 0:
                            hashes[state] = digest_outputs(output)
                        shutil.rmtree(output)
                        if repetition >= 0:
                            rows.append(dict(format="CRAM" if cram else "FASTQ", mode="PE" if paired else "SE", positive=positive,
                                             state=state, repetition=repetition, seconds=seconds,
                                             reads_per_second=args.records * (2 if paired else 1) / seconds))
                        print(f"{name}: {seconds:.4f}s", flush=True)
                assert hashes["baseline"] == hashes["disabled"], "Disabled path differs from baseline output"
                if not positive:
                    assert hashes["baseline"] == hashes["enabled"], "Negative control was incorrectly trimmed"
                else:
                    assert hashes["baseline"] != hashes["enabled"], "Positive control was not trimmed"
    with (root / "timings.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summaries = []
    for fmt in ("FASTQ", "CRAM") if args.with_cram else ("FASTQ",):
        for mode in ("SE", "PE"):
            for positive in (False, True):
                subset = [row for row in rows if (row["format"], row["mode"], row["positive"]) == (fmt, mode, positive)]
                baseline = statistics.median(row["seconds"] for row in subset if row["state"] == "baseline")
                for state in ("baseline", "disabled", "enabled"):
                    times = [row["seconds"] for row in subset if row["state"] == state]
                    median = statistics.median(times)
                    summaries.append(dict(format=fmt, mode=mode, positive=positive, state=state, median_seconds=median,
                                          reads_per_second=args.records * (2 if mode == "PE" else 1) / median,
                                          throughput_change_percent=(baseline / median - 1) * 100, min_seconds=min(times), max_seconds=max(times)))
    report = dict(platform=platform.platform(), threads=args.threads, records_or_pairs=args.records, repeats=args.repeats,
                  raw_read_length=REACH + 27, positive_fraction=0.5, validation_checks=validation_cases,
                  binary_sha256=hashlib.sha256(args.binary.read_bytes()).hexdigest(), baseline_sha256=hashlib.sha256(args.baseline.read_bytes()).hexdigest(), results=summaries)
    (root / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
