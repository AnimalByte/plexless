#!/usr/bin/env python3
"""
Create a deterministic plexless validation dataset from paired FASTQ input.

Typical input is InSilicoSeq output:
    raw_R1.fastq.gz
    raw_R2.fastq.gz

Generated categories:
    assigned    70%
    unmatched   15%
    unrouted    15%
    ambiguous    0%

Ambiguous is intentionally 0 because plexless rejects correction-unsafe
whitelist geometry before demultiplexing.

Single-end structure:
    R1_4A4B2T

Paired-end structures:
    R1_2A2B2T
    R2_2A2B2T

The paired fixture splits logical A and B barcodes across R1 and R2 so the
test also verifies cross-mate barcode concatenation.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import IO, Iterator


A_BARCODES = [
    ("A01", "AAAA"),
    ("A02", "CCCC"),
    ("A03", "GGGG"),
    ("A04", "TTTT"),
    # Correction-safe but deliberately absent from every sample route.
    ("A05", "ACGT"),
]

B_BARCODES = [
    ("B01", "ACGT"),
    ("B02", "CATG"),
    ("B03", "GTAC"),
    ("B04", "TGCA"),
]

SAMPLES = [
    ("sample_1", "A01", "B01"),
    ("sample_2", "A02", "B02"),
    ("sample_3", "A03", "B03"),
    ("sample_4", "A04", "B04"),
]

# A valid unused root barcode produces the RoutingTree Unrouted terminal.
UNROUTED_COMBINATIONS = [
    ("A05", "B01"),
    ("A05", "B02"),
    ("A05", "B03"),
    ("A05", "B04"),
]

# Hamming distance >= 2 from every A whitelist barcode, so it cannot be
# rescued with --max-mismatches 1.
UNMATCHED_A = "ATAT"

SE_STRUCTURE = "R1_4A4B2T"
PE_R1_STRUCTURE = "R1_2A2B2T"
PE_R2_STRUCTURE = "R2_2A2B2T"

SE_TECH = "GG"
PE_R1_TECH = "GG"
PE_R2_TECH = "CC"

# Deterministic 20-record repeating schedule.
SCHEDULE = (
    ["assigned"] * 14
    + ["unmatched"] * 3
    + ["unrouted"] * 3
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--r1", type=Path, default=Path("raw_R1.fastq.gz"))
    p.add_argument("--r2", type=Path, default=Path("raw_R2.fastq.gz"))
    p.add_argument("--outdir", type=Path, default=Path("plexless_validation"))
    p.add_argument("--limit", type=int, default=None)
    p.add_argument(
        "--gzip-level",
        type=int,
        choices=range(10),
        default=2,
        metavar="0-9",
    )
    return p.parse_args()


def open_read(path: Path) -> IO[str]:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="ascii", newline="")
    return path.open("rt", encoding="ascii", newline="")


def open_write(path: Path, level: int) -> IO[str]:
    if path.suffix == ".gz":
        return gzip.open(
            path,
            "wt",
            encoding="ascii",
            newline="",
            compresslevel=level,
        )
    return path.open("wt", encoding="ascii", newline="")


def fastq_records(
    handle: IO[str],
    source: Path,
) -> Iterator[tuple[str, str, str, str]]:
    n = 0
    while True:
        h = handle.readline()
        if h == "":
            return

        s = handle.readline()
        p = handle.readline()
        q = handle.readline()
        n += 1

        if not s or not p or not q:
            raise ValueError(f"{source}: truncated FASTQ record {n}")

        h = h.rstrip("\r\n")
        s = s.rstrip("\r\n")
        p = p.rstrip("\r\n")
        q = q.rstrip("\r\n")

        if not h.startswith("@"):
            raise ValueError(f"{source}: record {n} header does not start with @")
        if not p.startswith("+"):
            raise ValueError(f"{source}: record {n} plus line does not start with +")
        if len(s) != len(q):
            raise ValueError(f"{source}: record {n} sequence/quality length mismatch")

        yield h, s, p, q


def core_id(header: str) -> str:
    token = header[1:] if header.startswith("@") else header
    token = token.split(maxsplit=1)[0]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def write_record(
    out: IO[str],
    h: str,
    s: str,
    p: str,
    q: str,
) -> None:
    out.write(f"{h}\n{s}\n{p}\n{q}\n")


def lookup(rows: list[tuple[str, str]], barcode_id: str) -> str:
    for current_id, seq in rows:
        if current_id == barcode_id:
            return seq
    raise KeyError(barcode_id)


def classify(index: int) -> str:
    return SCHEDULE[index % len(SCHEDULE)]


def write_metadata(
    outdir: Path,
    category_counts: dict[str, int],
    sample_counts: dict[str, int],
    total: int,
) -> None:
    with (outdir / "barcodes.tsv").open("w", encoding="ascii") as f:
        f.write("Set\tID\tSequence\n")
        for barcode_id, seq in A_BARCODES:
            f.write(f"A\t{barcode_id}\t{seq}\n")
        for barcode_id, seq in B_BARCODES:
            f.write(f"B\t{barcode_id}\t{seq}\n")

    with (outdir / "samples.tsv").open("w", encoding="ascii") as f:
        f.write("Sample\tA\tB\n")
        for sample, a_id, b_id in SAMPLES:
            f.write(f"{sample}\t{a_id}\t{b_id}\n")

    with (outdir / "expected_summary.tsv").open("w", encoding="ascii") as f:
        f.write("Category\tExpected\n")
        f.write(f"total\t{total}\n")
        for category in ("assigned", "unmatched", "ambiguous", "unrouted"):
            f.write(f"{category}\t{category_counts.get(category, 0)}\n")

    with (outdir / "expected_samples.tsv").open("w", encoding="ascii") as f:
        f.write("Sample\tExpectedAssigned\n")
        for sample, _, _ in SAMPLES:
            f.write(f"{sample}\t{sample_counts[sample]}\n")

    with (outdir / "metadata.txt").open("w", encoding="ascii") as f:
        f.write("plexless deterministic validation fixture\n")
        f.write(f"records_or_pairs\t{total}\n")
        f.write(f"single_end_structure\t{SE_STRUCTURE}\n")
        f.write(f"paired_r1_structure\t{PE_R1_STRUCTURE}\n")
        f.write(f"paired_r2_structure\t{PE_R2_STRUCTURE}\n")
        f.write("max_mismatches\t1\n")
        f.write("category_schedule\t14 assigned, 3 unmatched, 3 unrouted per 20\n")
        f.write("ambiguous_expected\t0\n")
        f.write(
            "ambiguous_note\tRuntime ambiguity is intentionally unreachable "
            "for a valid correction-safe whitelist under current plexless rules.\n"
        )


def main() -> None:
    args = parse_args()

    if args.limit is not None and args.limit <= 0:
        raise SystemExit("--limit must be > 0")
    if not args.r1.exists():
        raise SystemExit(f"Missing R1: {args.r1}")
    if not args.r2.exists():
        raise SystemExit(f"Missing R2: {args.r2}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    se_path = args.outdir / "se.fastq.gz"
    pe_r1_path = args.outdir / "pe_R1.fastq.gz"
    pe_r2_path = args.outdir / "pe_R2.fastq.gz"
    truth_path = args.outdir / "truth.tsv.gz"

    for path in (se_path, pe_r1_path, pe_r2_path, truth_path):
        if path.exists():
            path.unlink()

    category_counts = {
        "assigned": 0,
        "unmatched": 0,
        "ambiguous": 0,
        "unrouted": 0,
    }
    sample_counts = {sample: 0 for sample, _, _ in SAMPLES}

    processed = 0

    with (
        open_read(args.r1) as r1h,
        open_read(args.r2) as r2h,
        open_write(se_path, args.gzip_level) as se_out,
        open_write(pe_r1_path, args.gzip_level) as pe1_out,
        open_write(pe_r2_path, args.gzip_level) as pe2_out,
        gzip.open(
            truth_path,
            "wt",
            encoding="ascii",
            newline="",
            compresslevel=args.gzip_level,
        ) as truth,
    ):
        r1_iter = fastq_records(r1h, args.r1)
        r2_iter = fastq_records(r2h, args.r2)

        truth.write("Index\tReadID\tCategory\tExpectedSample\tA\tB\n")

        while args.limit is None or processed < args.limit:
            try:
                r1 = next(r1_iter)
            except StopIteration:
                r1 = None

            try:
                r2 = next(r2_iter)
            except StopIteration:
                r2 = None

            if r1 is None and r2 is None:
                break
            if r1 is None or r2 is None:
                raise ValueError("R1 and R2 contain different record counts")

            h1, s1, p1, q1 = r1
            h2, s2, p2, q2 = r2

            if core_id(h1) != core_id(h2):
                raise ValueError(
                    f"Pair {processed + 1} IDs differ: {h1!r} vs {h2!r}"
                )

            category = classify(processed)
            expected_sample = "."

            if category == "assigned":
                sample, a_id, b_id = SAMPLES[processed % len(SAMPLES)]
                a = lookup(A_BARCODES, a_id)
                b = lookup(B_BARCODES, b_id)
                expected_sample = sample
                sample_counts[sample] += 1

            elif category == "unmatched":
                _, _, b_id = SAMPLES[processed % len(SAMPLES)]
                a_id = "NO_MATCH"
                a = UNMATCHED_A
                b = lookup(B_BARCODES, b_id)

            elif category == "unrouted":
                a_id, b_id = UNROUTED_COMBINATIONS[
                    processed % len(UNROUTED_COMBINATIONS)
                ]
                a = lookup(A_BARCODES, a_id)
                b = lookup(B_BARCODES, b_id)

            else:
                raise AssertionError(category)

            category_counts[category] += 1

            se_prefix = a + b + SE_TECH
            write_record(
                se_out,
                h1,
                se_prefix + s1,
                p1,
                ("I" * len(se_prefix)) + q1,
            )

            pe1_prefix = a[:2] + b[:2] + PE_R1_TECH
            pe2_prefix = a[2:] + b[2:] + PE_R2_TECH

            write_record(
                pe1_out,
                h1,
                pe1_prefix + s1,
                p1,
                ("I" * len(pe1_prefix)) + q1,
            )
            write_record(
                pe2_out,
                h2,
                pe2_prefix + s2,
                p2,
                ("I" * len(pe2_prefix)) + q2,
            )

            truth.write(
                f"{processed + 1}\t{core_id(h1)}\t{category}\t"
                f"{expected_sample}\t{a_id}\t{b_id}\n"
            )

            processed += 1
            if processed % 100_000 == 0:
                print(f"Processed {processed:,} pairs...", flush=True)

    write_metadata(
        args.outdir,
        category_counts,
        sample_counts,
        processed,
    )

    print()
    print(f"Done: {processed:,} source pairs")
    print("Expected plexless summary:")
    for category in ("assigned", "unmatched", "ambiguous", "unrouted"):
        print(f"  {category:10s} {category_counts[category]:,}")
    print(f"  {'total':10s} {processed:,}")

    print()
    print("Expected assigned sample counts:")
    for sample, _, _ in SAMPLES:
        print(f"  {sample:10s} {sample_counts[sample]:,}")


if __name__ == "__main__":
    main()
