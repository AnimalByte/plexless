#!/usr/bin/env python3
"""
Deep validation of plexless output against truth.tsv.gz.

Checks:
  * every assigned read ID is in the expected sample
  * every non-assigned read ID is in unassigned output
  * per-sample order matches source/truth order
  * paired R1/R2 output IDs both match truth
  * no expected records are missing or duplicated
  * output file set is exactly what is expected
  * optionally, every decompressed FASTQ is byte-identical to a baseline
    run (normally --threads 1)
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
from pathlib import Path
from typing import Iterator, TextIO


SAMPLES = ("sample_1", "sample_2", "sample_3", "sample_4")


def core_read_id(header: str) -> str:
    token = header[1:] if header.startswith("@") else header
    token = token.split(maxsplit=1)[0]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def fastq_ids(path: Path) -> Iterator[str]:
    with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
        record_number = 0
        while True:
            header = handle.readline()
            if header == "":
                return

            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            record_number += 1

            if not seq or not plus or not qual:
                raise SystemExit(f"{path}: truncated FASTQ record {record_number}")

            header = header.rstrip("\r\n")
            seq = seq.rstrip("\r\n")
            plus = plus.rstrip("\r\n")
            qual = qual.rstrip("\r\n")

            if not header.startswith("@"):
                raise SystemExit(f"{path}: invalid header at record {record_number}")
            if not plus.startswith("+"):
                raise SystemExit(f"{path}: invalid plus line at record {record_number}")
            if len(seq) != len(qual):
                raise SystemExit(
                    f"{path}: sequence/quality length mismatch at record {record_number}"
                )

            yield core_read_id(header)


def load_truth(path: Path) -> tuple[dict[str, list[str]], list[str]]:
    assigned = {sample: [] for sample in SAMPLES}
    unassigned: list[str] = []

    with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
        header = handle.readline().rstrip("\r\n")
        expected_header = "Index\tReadID\tCategory\tExpectedSample\tA\tB"
        if header != expected_header:
            raise SystemExit(f"Unexpected truth.tsv.gz header: {header!r}")

        for line_number, line in enumerate(handle, start=2):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 6:
                raise SystemExit(
                    f"{path}: line {line_number} has {len(fields)} columns, expected 6"
                )

            _, read_id, category, expected_sample, _, _ = fields

            if category == "assigned":
                if expected_sample not in assigned:
                    raise SystemExit(
                        f"{path}: unknown expected sample {expected_sample!r}"
                    )
                assigned[expected_sample].append(read_id)
            elif category in ("unmatched", "ambiguous", "unrouted"):
                unassigned.append(read_id)
            else:
                raise SystemExit(f"{path}: unknown category {category!r}")

    return assigned, unassigned


def compare_ids(path: Path, expected: list[str], label: str) -> None:
    if not path.exists():
        raise SystemExit(f"Missing expected output file: {path}")

    observed_count = 0

    for observed_count, (observed, wanted) in enumerate(
        zip(fastq_ids(path), expected, strict=False),
        start=1,
    ):
        if observed != wanted:
            raise SystemExit(
                f"{label}: read-order/assignment mismatch at output record "
                f"{observed_count:,}: expected {wanted!r}, observed {observed!r}"
            )

    # zip() stops at the shortest input, so count the actual file independently
    actual_count = sum(1 for _ in fastq_ids(path))

    if actual_count != len(expected):
        raise SystemExit(
            f"{label}: expected {len(expected):,} records, observed {actual_count:,}"
        )

    print(f"PASS  {label:<28} {actual_count:>10,d} records")


def decompressed_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def expected_filenames(mode: str) -> set[str]:
    if mode == "se":
        return {f"{sample}.fastq.gz" for sample in SAMPLES} | {"unassigned.fastq.gz"}

    return (
        {f"{sample}_R1.fastq.gz" for sample in SAMPLES}
        | {f"{sample}_R2.fastq.gz" for sample in SAMPLES}
        | {"unassigned_R1.fastq.gz", "unassigned_R2.fastq.gz"}
    )


def check_file_set(output: Path, mode: str) -> set[str]:
    expected = expected_filenames(mode)
    observed = {path.name for path in output.glob("*.fastq.gz")}

    if observed != expected:
        missing = sorted(expected - observed)
        extra = sorted(observed - expected)
        parts = []
        if missing:
            parts.append(f"missing={missing}")
        if extra:
            parts.append(f"extra={extra}")
        raise SystemExit("Unexpected FASTQ output file set: " + ", ".join(parts))

    print(f"PASS  output file set             {len(observed):>10,d} files")
    return observed


def compare_to_baseline(output: Path, baseline: Path, filenames: set[str]) -> None:
    print("\nComparing decompressed FASTQ bytes with 1-thread baseline:")

    baseline_files = {path.name for path in baseline.glob("*.fastq.gz")}
    if baseline_files != filenames:
        raise SystemExit(
            "Baseline output file set differs from current output file set"
        )

    for filename in sorted(filenames):
        current_hash = decompressed_sha256(output / filename)
        baseline_hash = decompressed_sha256(baseline / filename)

        if current_hash != baseline_hash:
            raise SystemExit(
                f"Decompressed FASTQ differs from baseline: {filename}"
            )

        print(f"PASS  {filename}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--mode", choices=("se", "pe"), required=True)
    p.add_argument("--baseline", type=Path)
    args = p.parse_args()

    truth_path = args.dataset / "truth.tsv.gz"
    if not truth_path.exists():
        raise SystemExit(f"Missing truth file: {truth_path}")

    assigned, unassigned = load_truth(truth_path)

    print("\nDeep record-level validation:")
    filenames = check_file_set(args.output, args.mode)

    if args.mode == "se":
        for sample in SAMPLES:
            compare_ids(
                args.output / f"{sample}.fastq.gz",
                assigned[sample],
                sample,
            )

        compare_ids(
            args.output / "unassigned.fastq.gz",
            unassigned,
            "unassigned",
        )

    else:
        for sample in SAMPLES:
            expected = assigned[sample]
            compare_ids(
                args.output / f"{sample}_R1.fastq.gz",
                expected,
                f"{sample} R1",
            )
            compare_ids(
                args.output / f"{sample}_R2.fastq.gz",
                expected,
                f"{sample} R2",
            )

        compare_ids(
            args.output / "unassigned_R1.fastq.gz",
            unassigned,
            "unassigned R1",
        )
        compare_ids(
            args.output / "unassigned_R2.fastq.gz",
            unassigned,
            "unassigned R2",
        )

    assigned_total = sum(len(values) for values in assigned.values())
    print(
        f"\nPASS  truth routing totals: assigned={assigned_total:,}, "
        f"unassigned={len(unassigned):,}"
    )

    if args.baseline is not None:
        compare_to_baseline(args.output, args.baseline, filenames)

    print("\nAll record-level checks passed.")


if __name__ == "__main__":
    main()
