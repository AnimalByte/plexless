#!/usr/bin/env python3
"""Create a small deterministic CRAM fixture with many routed destinations.

This fixture is intentionally tiny in record count. It exists to qualify CRAM
writer descriptor and memory scaling without generating a large biological
dataset.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


BASES = "ACGT"


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be greater than zero")
    return parsed


def barcode(index: int, width: int = 8) -> str:
    digits = []
    for _ in range(width):
        digits.append(BASES[index % 4])
        index //= 4
    if index:
        raise ValueError("sample count exceeds the eight-base barcode space")
    return "".join(reversed(digits))


def biological_sequence(sample_index: int, repetition: int, length: int = 100) -> str:
    state = (sample_index + 1) * 0x9E3779B1 ^ (repetition + 1) * 0x85EBCA77
    sequence = []
    for _ in range(length):
        state = (1664525 * state + 1013904223) & 0xFFFFFFFF
        sequence.append(BASES[(state >> 30) & 3])
    return "".join(sequence)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=positive_int, required=True)
    parser.add_argument("--records-per-sample", type=positive_int, default=4)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--plexless-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.samples > 4**8:
        raise SystemExit("--samples exceeds the eight-base barcode space")
    if args.outdir.exists():
        raise SystemExit(f"output already exists: {args.outdir}")
    args.outdir.mkdir(parents=True)

    barcodes = args.outdir / "barcodes.tsv"
    samples = args.outdir / "samples.tsv"
    reads = args.outdir / "reads.fastq"
    cram = args.outdir / "reads.cram"

    with barcodes.open("wt", encoding="ascii", newline="") as handle:
        handle.write("Set\tID\tSequence\n")
        for index in range(args.samples):
            handle.write(f"A\tA{index + 1:05d}\t{barcode(index)}\n")

    with samples.open("wt", encoding="ascii", newline="") as handle:
        handle.write("Sample\tA\n")
        for index in range(args.samples):
            handle.write(f"sample_{index + 1:05d}\tA{index + 1:05d}\n")

    with reads.open("wt", encoding="ascii", newline="") as handle:
        ordinal = 0
        for index in range(args.samples):
            for repetition in range(args.records_per_sample):
                ordinal += 1
                sequence = barcode(index) + "TT" + biological_sequence(index, repetition)
                handle.write(
                    f"@fanout:{ordinal:09d}\n{sequence}\n+\n{'I' * len(sequence)}\n"
                )

    subprocess.run(
        [
            "cargo",
            "run",
            "--quiet",
            "--release",
            "--example",
            "make_unmapped_cram_fixture",
            "--",
            "single",
            str(reads),
            str(cram),
        ],
        cwd=args.plexless_root,
        check=True,
    )
    print(
        f"Created {args.samples} destinations and "
        f"{args.samples * args.records_per_sample} records in {args.outdir}"
    )


if __name__ == "__main__":
    main()
