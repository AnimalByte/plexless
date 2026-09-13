#!/usr/bin/env python3
"""Generate deterministic FASTQ fixtures shared by Plexless, fqtk, and fastp.

The performance fixtures intentionally exercise only the common denominator:
one fixed-position R1 sample barcode, a two-base skipped technical segment,
paired biological templates, flat sample routing, and no unassigned reads.
"""

from __future__ import annotations

import argparse
import contextlib
import gzip
import hashlib
import json
import shutil
from pathlib import Path


BARCODE_BITS = 16
DATA_BITS = 11
DATA_POSITIONS = (3, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15)
PARITY_POSITIONS = (1, 2, 4, 8)
BASE_PAIRS = (
    ("A", "C"),
    ("G", "T"),
    ("C", "A"),
    ("T", "G"),
)
BASE_TRANSLATION = bytes(b"ACGT"[value & 3] for value in range(256))
QUALITY_TRANSLATION = bytes(33 + 30 + value % 11 for value in range(256))


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be greater than zero")
    return parsed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--pairs", type=positive_int, default=10_000_000)
    parser.add_argument(
        "--samples", type=positive_int, nargs="+", default=(96, 384, 1500)
    )
    parser.add_argument("--correctness-pairs", type=positive_int, default=4096)
    parser.add_argument("--correctness-samples", type=positive_int, default=8)
    parser.add_argument("--read-length", type=positive_int, default=150)
    parser.add_argument("--sequence-pool-size", type=positive_int, default=65_536)
    parser.add_argument("--seed", type=int, default=20_250_524)
    parser.add_argument(
        "--gzip-level", type=int, choices=range(1, 10), default=2, metavar="1-9"
    )
    parser.add_argument("--replace", action="store_true")
    return parser.parse_args()


def extended_hamming_barcode(value: int) -> str:
    """Return one word from a binary extended Hamming(16, 11, 4) code."""
    if not 0 <= value < 2**DATA_BITS:
        raise ValueError(f"barcode index {value} exceeds {2**DATA_BITS} code words")
    bits = [0] * (BARCODE_BITS + 1)
    for bit_index, position in enumerate(DATA_POSITIONS):
        bits[position] = (value >> bit_index) & 1
    for parity in PARITY_POSITIONS:
        bits[parity] = sum(
            bits[position]
            for position in range(1, 16)
            if position != parity and position & parity
        ) % 2
    bits[16] = sum(bits[1:16]) % 2
    return "".join(
        BASE_PAIRS[(position - 1) % len(BASE_PAIRS)][bits[position]]
        for position in range(1, 17)
    )


def hamming_distance(left: str, right: str) -> int:
    return sum(a != b for a, b in zip(left, right, strict=True))


def validate_barcodes(barcodes: list[str]) -> None:
    if len(set(barcodes)) != len(barcodes):
        raise RuntimeError("generated barcode code contains duplicates")
    minimum = BARCODE_BITS
    for index, left in enumerate(barcodes):
        for right in barcodes[index + 1 :]:
            minimum = min(minimum, hamming_distance(left, right))
            if minimum < 4:
                raise RuntimeError(
                    f"generated barcode code has distance {minimum}, expected at least 4"
                )


def deterministic_record_bytes(
    seed: int, mate: int, index: int, read_length: int
) -> tuple[bytes, bytes]:
    material = f"plexless-compare:{seed}:{mate}:{index}".encode("ascii")
    raw = hashlib.shake_256(material).digest(read_length * 2)
    sequence = raw[:read_length].translate(BASE_TRANSLATION)
    quality = raw[read_length:].translate(QUALITY_TRANSLATION)
    return sequence, quality


def write_sequence_pool(
    outdir: Path, seed: int, read_length: int, pool_size: int
) -> dict[str, str]:
    paths: dict[str, str] = {}
    handles = {}
    try:
        for mate in (1, 2):
            for field in ("sequence", "quality"):
                path = outdir / f"pool_R{mate}_{field}.bin"
                handles[(mate, field)] = path.open("wb")
                paths[f"r{mate}_{field}"] = path.name
        for index in range(pool_size):
            for mate in (1, 2):
                sequence, quality = deterministic_record_bytes(
                    seed, mate, index, read_length
                )
                handles[(mate, "sequence")].write(sequence)
                handles[(mate, "quality")].write(quality)
    finally:
        for handle in handles.values():
            handle.close()
    return paths


def load_sequence_pool(
    outdir: Path, paths: dict[str, str]
) -> dict[tuple[int, str], bytes]:
    return {
        (mate, field): (outdir / paths[f"r{mate}_{field}"]).read_bytes()
        for mate in (1, 2)
        for field in ("sequence", "quality")
    }


def write_routing_files(directory: Path, barcodes: list[str]) -> None:
    with (directory / "barcodes.tsv").open(
        "wt", encoding="ascii", newline=""
    ) as handle:
        handle.write("Set\tID\tSequence\n")
        for index, barcode in enumerate(barcodes, start=1):
            handle.write(f"A\tBC{index:04d}\t{barcode}\n")
    with (directory / "samples.tsv").open(
        "wt", encoding="ascii", newline=""
    ) as handle:
        handle.write("Sample\tA\n")
        for index in range(1, len(barcodes) + 1):
            handle.write(f"sample_{index:04d}\tBC{index:04d}\n")
    with (directory / "fqtk_samples.tsv").open(
        "wt", encoding="ascii", newline=""
    ) as handle:
        handle.write("sample_id\tbarcode\n")
        for index, barcode in enumerate(barcodes, start=1):
            handle.write(f"sample_{index:04d}\t{barcode}\n")


@contextlib.contextmanager
def open_deterministic_gzip(path: Path, level: int):
    with path.open("wb") as raw:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw, compresslevel=level, mtime=0
        ) as compressed:
            yield compressed


def pool_record(
    pool: dict[tuple[int, str], bytes],
    mate: int,
    pool_index: int,
    read_length: int,
) -> tuple[bytes, bytes]:
    start = pool_index * read_length
    end = start + read_length
    return pool[(mate, "sequence")][start:end], pool[(mate, "quality")][start:end]


def append_pair(
    r1_buffer: bytearray,
    r2_buffer: bytearray,
    ordinal: int,
    barcode: str,
    r1_sequence: bytes,
    r1_quality: bytes,
    r2_sequence: bytes,
    r2_quality: bytes,
) -> None:
    header = f"@plexcmp:{ordinal:012d}".encode("ascii")
    barcode_bytes = barcode.encode("ascii")
    r1_buffer.extend(header + b"/1\n")
    r1_buffer.extend(barcode_bytes + b"GG" + r1_sequence + b"\n+\n")
    r1_buffer.extend(b"I" * (BARCODE_BITS + 2) + r1_quality + b"\n")
    r2_buffer.extend(header + b"/2\n")
    r2_buffer.extend(r2_sequence + b"\n+\n" + r2_quality + b"\n")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(4 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def write_fastqs(
    directory: Path,
    pairs: int,
    barcodes: list[str],
    pool: dict[tuple[int, str], bytes],
    pool_size: int,
    read_length: int,
    gzip_level: int,
    correctness: bool,
) -> dict[str, object]:
    r1_path = directory / "R1.fastq.gz"
    r2_path = directory / "R2.fastq.gz"
    unmatched_barcode = extended_hamming_barcode(2**DATA_BITS - 1)
    if unmatched_barcode in barcodes:
        raise RuntimeError("reserved unmatched barcode is present in the sample set")
    case_counts = {"exact": 0, "mismatch": 0, "unmatched": 0}
    with open_deterministic_gzip(r1_path, gzip_level) as r1, open_deterministic_gzip(
        r2_path, gzip_level
    ) as r2:
        r1_buffer = bytearray()
        r2_buffer = bytearray()
        for zero_based in range(pairs):
            ordinal = zero_based + 1
            sample_index = zero_based % len(barcodes)
            barcode = barcodes[sample_index]
            case = "exact"
            if correctness:
                case = ("exact", "mismatch", "unmatched")[zero_based % 3]
                if case == "mismatch":
                    replacement = (
                        BASE_PAIRS[0][1]
                        if barcode[0] == BASE_PAIRS[0][0]
                        else BASE_PAIRS[0][0]
                    )
                    barcode = replacement + barcode[1:]
                elif case == "unmatched":
                    barcode = unmatched_barcode
            case_counts[case] += 1
            pool_index = zero_based % pool_size
            r1_sequence, r1_quality = pool_record(
                pool, 1, pool_index, read_length
            )
            r2_sequence, r2_quality = pool_record(
                pool, 2, pool_index, read_length
            )
            append_pair(
                r1_buffer,
                r2_buffer,
                ordinal,
                barcode,
                r1_sequence,
                r1_quality,
                r2_sequence,
                r2_quality,
            )
            if len(r1_buffer) >= 4 * 1024 * 1024:
                r1.write(r1_buffer)
                r2.write(r2_buffer)
                r1_buffer.clear()
                r2_buffer.clear()
            if ordinal % 1_000_000 == 0:
                print(
                    f"  {directory.name}: wrote {ordinal:,}/{pairs:,} pairs",
                    flush=True,
                )
        if r1_buffer:
            r1.write(r1_buffer)
            r2.write(r2_buffer)
    return {
        "r1": r1_path.name,
        "r2": r2_path.name,
        "r1_bytes": r1_path.stat().st_size,
        "r2_bytes": r2_path.stat().st_size,
        "r1_sha256": sha256_file(r1_path),
        "r2_sha256": sha256_file(r2_path),
        "case_counts": case_counts,
    }


def write_fixture(
    directory: Path,
    pairs: int,
    sample_count: int,
    all_barcodes: list[str],
    pool: dict[tuple[int, str], bytes],
    args: argparse.Namespace,
    correctness: bool,
) -> dict[str, object]:
    directory.mkdir(parents=True)
    barcodes = all_barcodes[:sample_count]
    write_routing_files(directory, barcodes)
    fastqs = write_fastqs(
        directory,
        pairs,
        barcodes,
        pool,
        args.sequence_pool_size,
        args.read_length,
        args.gzip_level,
        correctness,
    )
    manifest: dict[str, object] = {
        "schema": 1,
        "kind": "correctness" if correctness else "performance",
        "pairs": pairs,
        "samples": sample_count,
        "seed": args.seed,
        "read_length": args.read_length,
        "raw_r1_length": BARCODE_BITS + 2 + args.read_length,
        "raw_r2_length": args.read_length,
        "barcode_length": BARCODE_BITS,
        "barcode_code": "binary extended Hamming(16,11,4), position-wise mapped to A/C/G/T",
        "minimum_barcode_distance": 4,
        "technical_skip_length": 2,
        "plexless_structure": f"R1_{BARCODE_BITS}A2T",
        "fqtk_r1_structure": f"{BARCODE_BITS}B2S+T",
        "fqtk_r2_structure": "+T",
        "max_mismatches": 1,
        "fqtk_min_mismatch_delta": 2,
        "gzip_level": args.gzip_level,
        "sequence_pool_size": args.sequence_pool_size,
        **fastqs,
    }
    (directory / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="ascii"
    )
    return manifest


def main() -> None:
    args = parse_args()
    if max((*args.samples, args.correctness_samples)) > 1500:
        raise SystemExit("this benchmark intentionally supports at most 1500 samples")
    if args.sequence_pool_size > args.pairs:
        args.sequence_pool_size = args.pairs
    outdir = args.outdir.expanduser().resolve()
    if outdir.exists():
        if not args.replace:
            raise SystemExit(f"output already exists (use --replace): {outdir}")
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)

    all_barcodes = [extended_hamming_barcode(index) for index in range(1500)]
    validate_barcodes(all_barcodes)
    pool_paths = write_sequence_pool(
        outdir, args.seed, args.read_length, args.sequence_pool_size
    )
    pool = load_sequence_pool(outdir, pool_paths)
    fixtures = []
    performance = outdir / "performance"
    for sample_count in args.samples:
        print(
            f"Generating {args.pairs:,} performance pairs for {sample_count} samples...",
            flush=True,
        )
        directory = performance / f"samples_{sample_count:04d}"
        fixtures.append(
            str(directory.relative_to(outdir))
        )
        write_fixture(
            directory,
            args.pairs,
            sample_count,
            all_barcodes,
            pool,
            args,
            correctness=False,
        )

    correctness_dir = outdir / "correctness"
    print(
        f"Generating {args.correctness_pairs:,} correctness pairs for "
        f"{args.correctness_samples} samples...",
        flush=True,
    )
    write_fixture(
        correctness_dir,
        args.correctness_pairs,
        args.correctness_samples,
        all_barcodes,
        pool,
        args,
        correctness=True,
    )
    root_manifest = {
        "schema": 1,
        "fixtures": fixtures,
        "correctness_fixture": str(correctness_dir.relative_to(outdir)),
        "pool_files": pool_paths,
        "read_length": args.read_length,
        "sequence_pool_size": args.sequence_pool_size,
        "seed": args.seed,
    }
    (outdir / "manifest.json").write_text(
        json.dumps(root_manifest, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
    )
    print(f"Created deterministic comparison dataset at {outdir}")


if __name__ == "__main__":
    main()
