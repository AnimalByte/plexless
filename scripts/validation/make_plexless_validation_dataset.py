#!/usr/bin/env python3
"""Generate deterministic, self-contained plexless FASTQ qualification data.

Without --r1/--r2 this script synthesizes biological reads. Supplying both
options preserves the original behavior of injecting Plexless prefixes into
external paired FASTQ data. Generation is streaming: records and truth rows
are never accumulated in memory.
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import shutil
import subprocess
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import IO, Iterator


PROFILE_SIZES = {"tiny": 4_096, "standard": 1_000_000, "stress": 10_000_000}
SAMPLE_DISTRIBUTIONS = ("legacy", "balanced", "skewed", "extreme-skew")
DEFAULT_OUTDIR = Path.home() / "plexless_validation_data" / "plexless_validation"
TRUTH_HEADER = (
    "OutputFile\tOutputIndex\tOrdinal\tReadID\tMate\tCase\tCategory\t"
    "ExpectedSample\tExpectedHeader\tExpectedSequence\tExpectedQuality"
)

CORE_A = [
    ("A01", "AAAAAA"),
    ("A02", "CCCCCC"),
    ("A03", "GGGGGG"),
    ("A04", "TTTTTT"),
    ("A05", "ACGTAC"),  # valid root, deliberately unrouted
]
CORE_B = [
    ("B01", "ACGTAC"),
    ("B02", "CATGCA"),
    ("B03", "GTACGT"),
    ("B04", "TGCATG"),
]
CORE_SAMPLES = [
    ("sample_1", "A01", "B01"),
    ("sample_2", "A02", "B02"),
    ("sample_3", "A03", "B03"),
    ("sample_4", "A04", "B04"),
]
CORE_CASES = (
    "exact",
    "exact",
    "mismatch_1",
    "observed_n",
    "exact",
    "unmatched",
    "unrouted",
    "short",
    "exact",
    "mismatch_1",
    "observed_n",
    "unmatched",
    "orphan_r1",
    "exact",
    "orphan_r2",
    "exact",
)
CORE_SE_STRUCTURE = "R1_6A6B2T"
CORE_PE_R1_STRUCTURE = "R1_3A3B2T"
CORE_PE_R2_STRUCTURE = "R2_3A(rc)3B(rc)2T"


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be greater than zero")
    return parsed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", choices=PROFILE_SIZES, default="tiny")
    parser.add_argument("--records", type=positive_int)
    parser.add_argument("--pairs", type=positive_int)
    parser.add_argument(
        "--special-records",
        type=positive_int,
        default=4_096,
        help="records/pairs in each focused semantic fixture (default: 4096)",
    )
    parser.add_argument("--seed", type=int, default=20_250_524)
    parser.add_argument(
        "--sample-distribution",
        choices=SAMPLE_DISTRIBUTIONS,
        default="legacy",
        help=(
            "core assigned-read distribution; legacy preserves the original "
            "fixture, while the other choices support CRAM writer benchmarks"
        ),
    )
    parser.add_argument("--r1", type=Path)
    parser.add_argument("--r2", type=Path)
    parser.add_argument(
        "--limit",
        type=positive_int,
        help="legacy alias setting both --records and --pairs in injection mode",
    )
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--gzip-level", type=int, choices=range(10), default=2, metavar="0-9"
    )
    parser.add_argument(
        "--replace",
        action="store_true",
        help="replace an existing generated dataset directory",
    )
    parser.add_argument(
        "--with-cram",
        action="store_true",
        help="also create equivalent unmapped CRAM fixtures with the Rust helper",
    )
    return parser.parse_args()


@contextmanager
def deterministic_gzip_text(path: Path, level: int) -> Iterator[IO[str]]:
    """Write text gzip with mtime=0 and no path-dependent gzip header."""
    raw = path.open("wb")
    compressed = gzip.GzipFile(
        filename="", mode="wb", fileobj=raw, compresslevel=level, mtime=0
    )
    text = io.TextIOWrapper(compressed, encoding="ascii", newline="")
    try:
        yield text
    finally:
        text.close()


def open_read(path: Path) -> IO[str]:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="ascii", newline="")
    return path.open("rt", encoding="ascii", newline="")


def fastq_records(
    handle: IO[str], source: Path
) -> Iterator[tuple[str, str, str, str]]:
    record_number = 0
    while True:
        header = handle.readline()
        if header == "":
            return
        sequence = handle.readline()
        plus = handle.readline()
        quality = handle.readline()
        record_number += 1
        if not sequence or not plus or not quality:
            raise ValueError(f"{source}: truncated FASTQ record {record_number}")
        header = header.rstrip("\r\n")
        sequence = sequence.rstrip("\r\n")
        plus = plus.rstrip("\r\n")
        quality = quality.rstrip("\r\n")
        if not header.startswith("@"):
            raise ValueError(f"{source}: record {record_number} header lacks @")
        if not plus.startswith("+"):
            raise ValueError(f"{source}: record {record_number} plus line is invalid")
        if len(sequence) != len(quality):
            raise ValueError(
                f"{source}: record {record_number} sequence/quality length mismatch"
            )
        yield header, sequence, plus, quality


def core_id(header: str) -> str:
    token = header.removeprefix("@").split(maxsplit=1)[0]
    if token.endswith(("/1", "/2")):
        return token[:-2]
    return token


def write_fastq(handle: IO[str], record: tuple[str, str, str, str]) -> None:
    header, sequence, plus, quality = record
    handle.write(f"{header}\n{sequence}\n{plus}\n{quality}\n")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def lookup(rows: list[tuple[str, str]], barcode_id: str) -> str:
    return dict(rows)[barcode_id]


def mutate(sequence: str, positions: tuple[int, ...]) -> str:
    replacements = {"A": "C", "C": "G", "G": "T", "T": "A"}
    result = list(sequence)
    for position in positions:
        result[position] = replacements[result[position]]
    return "".join(result)


def synthetic_biology(index: int, mate: int, seed: int) -> tuple[str, str]:
    """Stable small PRNG with variable 88-151 bp inserts and realistic Q scores."""
    state = (
        (seed & 0xFFFFFFFFFFFFFFFF)
        ^ ((index + 1) * 0x9E3779B97F4A7C15)
        ^ (mate * 0xD1B54A32D192ED03)
    ) & 0xFFFFFFFFFFFFFFFF
    length = 88 + ((state >> 17) % 64)
    bases: list[str] = []
    qualities: list[str] = []
    alphabet = "ACGT"
    for cycle in range(length):
        state ^= state >> 12
        state ^= (state << 25) & 0xFFFFFFFFFFFFFFFF
        state ^= state >> 27
        state = (state * 0x2545F4914F6CDD1D) & 0xFFFFFFFFFFFFFFFF
        bases.append(alphabet[state & 3])
        phred = max(25, 40 - cycle // 24 - ((state >> 8) & 3))
        qualities.append(chr(33 + phred))
    return "".join(bases), "".join(qualities)


class TruthSpool:
    """Spool truth per output so final truth is stream-verifiable at any size."""

    def __init__(self, fixture_dir: Path, mode: str, level: int) -> None:
        self.fixture_dir = fixture_dir
        self.mode = mode
        self.level = level
        self.temp = tempfile.TemporaryDirectory(prefix=f"truth-{mode}-", dir=fixture_dir)
        self.handles: dict[str, IO[str]] = {}
        self.counts: dict[str, int] = {}

    def add(
        self,
        output_file: str,
        ordinal: int,
        read_id: str,
        mate: str,
        case: str,
        category: str,
        sample: str,
        header: str,
        sequence: str,
        quality: str,
    ) -> None:
        if any("\t" in value or "\n" in value for value in (header, sequence, quality)):
            raise ValueError("truth fields cannot contain tabs or newlines")
        if output_file not in self.handles:
            spool = Path(self.temp.name) / f"stream-{len(self.handles):04d}.tsv"
            self.handles[output_file] = spool.open("w", encoding="ascii", newline="")
            self.counts[output_file] = 0
        self.counts[output_file] += 1
        self.handles[output_file].write(
            f"{output_file}\t{self.counts[output_file]}\t{ordinal}\t{read_id}\t"
            f"{mate}\t{case}\t{category}\t{sample}\t{header}\t{sequence}\t"
            f"{quality}\n"
        )

    def finish(self) -> Path:
        spool_paths = {
            output: Path(handle.name) for output, handle in self.handles.items()
        }
        for handle in self.handles.values():
            handle.close()
        truth_path = self.fixture_dir / f"truth_{self.mode}.tsv.gz"
        with deterministic_gzip_text(truth_path, self.level) as output:
            output.write(TRUTH_HEADER + "\n")
            for filename in sorted(spool_paths):
                with spool_paths[filename].open("r", encoding="ascii", newline="") as source:
                    shutil.copyfileobj(source, output, length=1024 * 1024)
        self.temp.cleanup()
        return truth_path


def write_catalog(
    fixture_dir: Path,
    sets: dict[str, list[tuple[str, str]]],
    samples: list[tuple[str, ...]],
) -> None:
    with (fixture_dir / "barcodes.tsv").open("w", encoding="ascii") as handle:
        handle.write("Set\tID\tSequence\n")
        for symbol, rows in sets.items():
            for barcode_id, sequence in rows:
                handle.write(f"{symbol}\t{barcode_id}\t{sequence}\n")
    symbols = list(sets)
    with (fixture_dir / "samples.tsv").open("w", encoding="ascii") as handle:
        handle.write("Sample\t" + "\t".join(symbols) + "\n")
        for row in samples:
            handle.write("\t".join(row) + "\n")


def empty_counts() -> dict[str, int]:
    return {
        "assigned": 0,
        "unmatched": 0,
        "ambiguous": 0,
        "unrouted": 0,
        "short_reads": 0,
        "orphan_r1": 0,
        "orphan_r2": 0,
    }


def write_expected(
    fixture_dir: Path,
    mode: str,
    counts: dict[str, int],
    sample_counts: dict[str, int],
) -> None:
    total = sum(counts[key] for key in ("assigned", "unmatched", "ambiguous", "unrouted"))
    total += counts["orphan_r1"] + counts["orphan_r2"]
    with (fixture_dir / f"expected_summary_{mode}.tsv").open(
        "w", encoding="ascii"
    ) as handle:
        handle.write("Category\tExpected\n")
        handle.write(f"total\t{total}\n")
        for key in counts:
            handle.write(f"{key}\t{counts[key]}\n")
    with (fixture_dir / f"expected_samples_{mode}.tsv").open(
        "w", encoding="ascii"
    ) as handle:
        handle.write("Sample\tExpectedAssigned\n")
        for sample, count in sample_counts.items():
            handle.write(f"{sample}\t{count}\n")


def distribution_sample_index(ordinal: int, distribution: str) -> int:
    if distribution == "legacy":
        return ordinal % len(CORE_SAMPLES)
    if distribution == "balanced":
        return ordinal % len(CORE_SAMPLES)
    if distribution == "skewed":
        position = ordinal % 10
        return 0 if position < 6 else 1 + ((position - 6) % 3)
    if distribution == "extreme-skew":
        return 0 if ordinal % 10 < 9 else 1 + ((ordinal // 10) % 3)
    raise ValueError(f"unknown sample distribution: {distribution}")


def observed_core_barcodes(
    index: int, case: str, distribution: str, assigned_ordinal: int
) -> tuple[str, str, str]:
    selection_ordinal = index if distribution == "legacy" else assigned_ordinal
    sample, a_id, b_id = CORE_SAMPLES[
        distribution_sample_index(selection_ordinal, distribution)
    ]
    a = lookup(CORE_A, a_id)
    b = lookup(CORE_B, b_id)
    if case == "mismatch_1":
        a = mutate(a, (index % len(a),))
    elif case == "observed_n":
        position = index % len(b)
        b = b[:position] + "N" + b[position + 1 :]
    elif case == "unmatched":
        a = mutate(lookup(CORE_A, "A01"), (0, 1))
        sample = "."
    elif case == "unrouted":
        a = lookup(CORE_A, "A05")
        sample = "."
    return a, b, sample


def add_truth_record(
    truth: TruthSpool,
    output_file: str,
    ordinal: int,
    record: tuple[str, str, str, str],
    mate: str,
    case: str,
    category: str,
    sample: str,
) -> None:
    header, sequence, _, quality = record
    truth.add(
        output_file,
        ordinal,
        core_id(header),
        mate,
        case,
        category,
        sample,
        header,
        sequence,
        quality,
    )


def generate_core(
    fixture_dir: Path,
    records: int,
    pairs: int,
    seed: int,
    level: int,
    source_r1: Path | None,
    source_r2: Path | None,
    sample_distribution: str,
) -> dict[str, object]:
    fixture_dir.mkdir(parents=True)
    write_catalog(fixture_dir, {"A": CORE_A, "B": CORE_B}, CORE_SAMPLES)
    se_counts = empty_counts()
    pe_counts = empty_counts()
    se_samples = {sample: 0 for sample, _, _ in CORE_SAMPLES}
    pe_samples = dict(se_samples)
    se_truth = TruthSpool(fixture_dir, "se", level)
    pe_truth = TruthSpool(fixture_dir, "pe", level)

    source_handles: tuple[IO[str], IO[str]] | None = None
    if source_r1 is not None and source_r2 is not None:
        source_handles = (open_read(source_r1), open_read(source_r2))
        source_iterators = (
            fastq_records(source_handles[0], source_r1),
            fastq_records(source_handles[1], source_r2),
        )
    else:
        source_iterators = None

    with (
        deterministic_gzip_text(fixture_dir / "se.fastq.gz", level) as se_out,
        deterministic_gzip_text(fixture_dir / "pe_R1.fastq.gz", level) as r1_out,
        deterministic_gzip_text(fixture_dir / "pe_R2.fastq.gz", level) as r2_out,
    ):
        maximum = max(records, pairs)
        produced = 0
        se_assigned_ordinal = 0
        pe_assigned_ordinal = 0
        for index in range(maximum):
            if source_iterators is None:
                bio1, qual1 = synthetic_biology(index, 1, seed)
                bio2, qual2 = synthetic_biology(index, 2, seed)
                read_id = f"plexless:{seed}:{index + 1:012d}"
                raw1 = (f"@{read_id}/1 lane=1 synthetic", bio1, "+source", qual1)
                raw2 = (f"@{read_id}/2 lane=1 synthetic", bio2, "+source", qual2)
            else:
                try:
                    raw1 = next(source_iterators[0])
                    raw2 = next(source_iterators[1])
                except StopIteration:
                    break
                if core_id(raw1[0]) != core_id(raw2[0]):
                    raise ValueError(
                        f"source pair {index + 1} differs: {raw1[0]!r} vs {raw2[0]!r}"
                    )
            case = CORE_CASES[index % len(CORE_CASES)]

            if index < records:
                se_case = "exact" if case.startswith("orphan_") else case
                se_is_assigned = se_case in ("exact", "mismatch_1", "observed_n")
                a, b, sample = observed_core_barcodes(
                    index, se_case, sample_distribution, se_assigned_ordinal
                )
                if se_is_assigned:
                    se_assigned_ordinal += 1
                if se_case == "short":
                    se_record = (raw1[0], "ACG", raw1[2], "III")
                    category = "unmatched"
                    expected = se_record
                    se_counts["short_reads"] += 1
                else:
                    prefix = a + b + "GG"
                    se_record = (
                        raw1[0],
                        prefix + raw1[1],
                        raw1[2],
                        "I" * len(prefix) + raw1[3],
                    )
                    if se_is_assigned:
                        category = "assigned"
                        expected = (raw1[0], raw1[1], "+", raw1[3])
                        se_samples[sample] += 1
                    else:
                        category = se_case
                        expected = se_record
                se_counts[category] += 1
                write_fastq(se_out, se_record)
                output_name = (
                    f"{sample}.fastq.gz"
                    if category == "assigned"
                    else "unassigned.fastq.gz"
                )
                add_truth_record(
                    se_truth,
                    output_name,
                    index + 1,
                    expected,
                    "SE",
                    se_case,
                    category,
                    sample if category == "assigned" else ".",
                )

            if index < pairs:
                pe_is_assigned = case in ("exact", "mismatch_1", "observed_n")
                a, b, sample = observed_core_barcodes(
                    index, case, sample_distribution, pe_assigned_ordinal
                )
                if pe_is_assigned:
                    pe_assigned_ordinal += 1
                r1_prefix = a[:3] + b[:3] + "GG"
                r2_prefix = reverse_complement(a[3:]) + reverse_complement(b[3:]) + "CC"
                input_r1 = (
                    raw1[0],
                    r1_prefix + raw1[1],
                    raw1[2],
                    "I" * len(r1_prefix) + raw1[3],
                )
                input_r2 = (
                    raw2[0],
                    r2_prefix + raw2[1],
                    raw2[2],
                    "I" * len(r2_prefix) + raw2[3],
                )
                if case == "short":
                    input_r1 = (raw1[0], "ACG", raw1[2], "III")
                    category = "unmatched"
                    expected_r1, expected_r2 = input_r1, input_r2
                    pe_counts["short_reads"] += 1
                elif case in ("exact", "mismatch_1", "observed_n"):
                    category = "assigned"
                    expected_r1 = (raw1[0], raw1[1], "+", raw1[3])
                    expected_r2 = (raw2[0], raw2[1], "+", raw2[3])
                else:
                    category = case
                    expected_r1, expected_r2 = input_r1, input_r2

                if case == "orphan_r1":
                    write_fastq(r1_out, input_r1)
                    pe_counts["orphan_r1"] += 1
                    add_truth_record(
                        pe_truth,
                        "unassigned_R1.fastq.gz",
                        index + 1,
                        expected_r1,
                        "R1",
                        case,
                        "orphan_r1",
                        ".",
                    )
                elif case == "orphan_r2":
                    write_fastq(r2_out, input_r2)
                    pe_counts["orphan_r2"] += 1
                    add_truth_record(
                        pe_truth,
                        "unassigned_R2.fastq.gz",
                        index + 1,
                        expected_r2,
                        "R2",
                        case,
                        "orphan_r2",
                        ".",
                    )
                else:
                    write_fastq(r1_out, input_r1)
                    write_fastq(r2_out, input_r2)
                    pe_counts[category] += 1
                    if category == "assigned":
                        pe_samples[sample] += 1
                        out_r1 = f"{sample}_R1.fastq.gz"
                        out_r2 = f"{sample}_R2.fastq.gz"
                    else:
                        out_r1 = "unassigned_R1.fastq.gz"
                        out_r2 = "unassigned_R2.fastq.gz"
                        sample = "."
                    add_truth_record(
                        pe_truth,
                        out_r1,
                        index + 1,
                        expected_r1,
                        "R1",
                        case,
                        category,
                        sample,
                    )
                    add_truth_record(
                        pe_truth,
                        out_r2,
                        index + 1,
                        expected_r2,
                        "R2",
                        case,
                        category,
                        sample,
                    )
            produced = index + 1
            if produced % 100_000 == 0:
                print(f"core: generated {produced:,} logical records", flush=True)

    if source_handles is not None:
        for handle in source_handles:
            handle.close()
    if produced < maximum:
        raise ValueError(
            f"source FASTQs ended after {produced:,} pairs; requested {maximum:,}"
        )
    se_truth.finish()
    pe_truth.finish()
    write_expected(fixture_dir, "se", se_counts, se_samples)
    write_expected(fixture_dir, "pe", pe_counts, pe_samples)
    return {
        "name": "core",
        "se": {
            "records": records,
            "structure": CORE_SE_STRUCTURE,
            "max_mismatches": 1,
            "sample_distribution": sample_distribution,
        },
        "pe": {
            "events": pairs,
            "r1_structure": CORE_PE_R1_STRUCTURE,
            "r2_structure": CORE_PE_R2_STRUCTURE,
            "max_mismatches": 1,
            "sample_distribution": sample_distribution,
        },
    }


def generate_single_a_fixture(
    root: Path,
    name: str,
    count: int,
    seed: int,
    level: int,
    mode: str,
) -> dict[str, object]:
    fixture = root / name
    fixture.mkdir()
    barcodes = [("A01", "AAAAAA"), ("A02", "CCCCCC")]
    samples = [("sample_a", "A01"), ("sample_c", "A02")]
    write_catalog(fixture, {"A": barcodes}, samples)
    truth = TruthSpool(fixture, mode, level)
    counts = empty_counts()
    sample_counts = {"sample_a": 0, "sample_c": 0}

    if mode == "se":
        with deterministic_gzip_text(fixture / "se.fastq.gz", level) as output:
            for index in range(count):
                case = ("exact", "n_plus_mismatch", "mismatch_2", "unmatched")[
                    index % 4
                ]
                if case == "exact":
                    barcode = "AAAAAA" if index % 8 < 4 else "CCCCCC"
                elif case == "n_plus_mismatch":
                    barcode = "NAAAAT"
                elif case == "mismatch_2":
                    barcode = "CCAAAA"
                else:
                    barcode = "CCCAAA"
                bio, qual = synthetic_biology(index, 1, seed + 101)
                read_id = f"n2:{seed}:{index + 1:012d}"
                raw = (f"@{read_id} synthetic", barcode + "GG" + bio, "+", "I" * 8 + qual)
                write_fastq(output, raw)
                if case == "unmatched":
                    category, sample, expected = "unmatched", ".", raw
                    output_name = "unassigned.fastq.gz"
                else:
                    category = "assigned"
                    sample = "sample_c" if barcode == "CCCCCC" else "sample_a"
                    expected = (raw[0], bio, "+", qual)
                    output_name = f"{sample}.fastq.gz"
                    sample_counts[sample] += 1
                counts[category] += 1
                add_truth_record(
                    truth,
                    output_name,
                    index + 1,
                    expected,
                    "SE",
                    case,
                    category,
                    sample,
                )
        truth.finish()
        write_expected(fixture, mode, counts, sample_counts)
        return {
            "name": name,
            "se": {"records": count, "structure": "R1_6A2T", "max_mismatches": 2},
        }

    structure = "R1_6A2T" if name == "r1_only" else "R2_6A(rc)2T"
    with (
        deterministic_gzip_text(fixture / "pe_R1.fastq.gz", level) as r1_out,
        deterministic_gzip_text(fixture / "pe_R2.fastq.gz", level) as r2_out,
    ):
        for index in range(count):
            case = ("exact", "mismatch_1", "unmatched", "exact")[index % 4]
            expected_barcode = "AAAAAA" if index % 8 < 4 else "CCCCCC"
            observed = expected_barcode
            if case == "mismatch_1":
                observed = mutate(observed, (0,))
            elif case == "unmatched":
                observed = "GGGAAA"
            bio1, qual1 = synthetic_biology(index, 1, seed + 211)
            bio2, qual2 = synthetic_biology(index, 2, seed + 211)
            read_id = f"{name}:{seed}:{index + 1:012d}"
            base_r1 = (f"@{read_id}/1 synthetic", bio1, "+", qual1)
            base_r2 = (f"@{read_id}/2 synthetic", bio2, "+", qual2)
            if name == "r1_only":
                input_r1 = (
                    base_r1[0], observed + "GG" + bio1, "+", "I" * 8 + qual1
                )
                input_r2 = base_r2
            else:
                input_r1 = base_r1
                prefix = reverse_complement(observed) + "CC"
                input_r2 = (base_r2[0], prefix + bio2, "+", "I" * 8 + qual2)
            write_fastq(r1_out, input_r1)
            write_fastq(r2_out, input_r2)
            if case == "unmatched":
                category, sample = "unmatched", "."
                expected_r1, expected_r2 = input_r1, input_r2
                out_r1, out_r2 = "unassigned_R1.fastq.gz", "unassigned_R2.fastq.gz"
            else:
                category = "assigned"
                sample = "sample_a" if expected_barcode == "AAAAAA" else "sample_c"
                expected_r1, expected_r2 = base_r1, base_r2
                out_r1 = f"{sample}_R1.fastq.gz"
                out_r2 = f"{sample}_R2.fastq.gz"
                sample_counts[sample] += 1
            counts[category] += 1
            add_truth_record(
                truth, out_r1, index + 1, expected_r1, "R1", case, category, sample
            )
            add_truth_record(
                truth, out_r2, index + 1, expected_r2, "R2", case, category, sample
            )
    truth.finish()
    write_expected(fixture, mode, counts, sample_counts)
    return {
        "name": name,
        "pe": {
            "pairs": count,
            "r1_structure": structure if name == "r1_only" else None,
            "r2_structure": structure if name == "r2_only" else None,
            "max_mismatches": 1,
        },
    }


def generate_negative_fixtures(root: Path, level: int) -> None:
    negative = root / "negative"
    negative.mkdir()
    (negative / "barcodes.tsv").write_text(
        "Set\tID\tSequence\nA\tA01\tAAAAAA\n", encoding="ascii"
    )
    (negative / "samples.tsv").write_text(
        "Sample\tA\nsample_a\tA01\n", encoding="ascii"
    )
    (negative / "unsafe_barcodes.tsv").write_text(
        "Set\tID\tSequence\nA\tA01\tAAAAAA\nA\tA02\tAAAAAC\n",
        encoding="ascii",
    )
    (negative / "unsafe_samples.tsv").write_text(
        "Sample\tA\nsample_a\tA01\nsample_b\tA02\n", encoding="ascii"
    )
    (negative / "invalid_barcode.tsv").write_text(
        "Set\tID\tSequence\nA\tA01\tAAAAXA\n", encoding="ascii"
    )
    (negative / "malformed_samples.tsv").write_text(
        "Wrong\tA\nsample_a\tA01\n", encoding="ascii"
    )
    (negative / "truncated.fastq").write_text(
        "@truncated\nAAAAAAGGACGT\n+\n", encoding="ascii"
    )
    (negative / "length_mismatch.fastq").write_text(
        "@bad-length\nAAAAAAGGACGT\n+\nIIII\n", encoding="ascii"
    )
    valid_text = "@corrupt\nAAAAAAGGACGT\n+\nIIIIIIIIIIII\n"
    corrupt_path = negative / "corrupt.fastq.gz"
    with deterministic_gzip_text(corrupt_path, level) as handle:
        handle.write(valid_text)
    compressed = bytearray(corrupt_path.read_bytes())
    compressed[-8] ^= 0xFF
    corrupt_path.write_bytes(compressed)

    with (
        deterministic_gzip_text(negative / "desync_R1.fastq.gz", level) as r1,
        deterministic_gzip_text(negative / "desync_R2.fastq.gz", level) as r2,
    ):
        for index in range(1_030):
            seq = "AAAAAAGGACGT"
            qual = "I" * len(seq)
            write_fastq(r1, (f"@r1-only-{index}/1", seq, "+", qual))
            write_fastq(r2, (f"@r2-only-{index}/2", seq, "+", qual))


def write_manifest(
    outdir: Path,
    profile: str,
    records: int,
    pairs: int,
    special_records: int,
    seed: int,
    source_mode: bool,
    with_cram: bool,
    sample_distribution: str,
    fixtures: list[dict[str, object]],
) -> None:
    manifest = {
        "schema_version": 2,
        "profile": profile,
        "seed": seed,
        "source_mode": "injected" if source_mode else "synthetic",
        "sample_distribution": sample_distribution,
        "cram_fixtures": with_cram,
        "core_records": records,
        "core_pair_events": pairs,
        "special_records": special_records,
        "truth_order": "grouped by OutputFile, then OutputIndex",
        "fixtures": fixtures,
    }
    (outdir / "metadata.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="ascii"
    )


def generate_cram_fixtures(
    outdir: Path, fixtures: list[dict[str, object]]
) -> None:
    repo = Path(__file__).resolve().parents[2]
    for fixture in fixtures:
        fixture_dir = outdir / str(fixture["name"])
        if "se" in fixture:
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
                    str(fixture_dir / "se.fastq.gz"),
                    str(fixture_dir / "se.cram"),
                ],
                cwd=repo,
                check=True,
            )
            assert isinstance(fixture["se"], dict)
            fixture["se"]["cram"] = "se.cram"
        if "pe" in fixture:
            subprocess.run(
                [
                    "cargo",
                    "run",
                    "--quiet",
                    "--release",
                    "--example",
                    "make_unmapped_cram_fixture",
                    "--",
                    "paired",
                    str(fixture_dir / "pe_R1.fastq.gz"),
                    str(fixture_dir / "pe_R2.fastq.gz"),
                    str(fixture_dir / "pe.cram"),
                ],
                cwd=repo,
                check=True,
            )
            assert isinstance(fixture["pe"], dict)
            fixture["pe"]["cram"] = "pe.cram"


def main() -> None:
    args = parse_args()
    if (args.r1 is None) != (args.r2 is None):
        raise SystemExit("--r1 and --r2 must be supplied together")
    source_mode = args.r1 is not None
    if source_mode and args.with_cram:
        raise SystemExit("--with-cram currently requires self-contained synthetic generation")
    if source_mode:
        assert args.r1 is not None and args.r2 is not None
        if not args.r1.is_file() or not args.r2.is_file():
            raise SystemExit("--r1 and --r2 must name existing FASTQ files")
    if args.limit is not None and (args.records is not None or args.pairs is not None):
        raise SystemExit("--limit cannot be combined with --records or --pairs")

    preset = PROFILE_SIZES[args.profile]
    records = args.limit or args.records or preset
    pairs = args.limit or args.pairs or preset
    outdir = args.outdir.expanduser().resolve()
    if outdir.exists():
        if not args.replace:
            raise SystemExit(f"output already exists (use --replace): {outdir}")
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)

    try:
        fixtures = [
            generate_core(
                outdir / "core",
                records,
                pairs,
                args.seed,
                args.gzip_level,
                args.r1,
                args.r2,
                args.sample_distribution,
            )
        ]
        focused_count = min(args.special_records, max(records, pairs))
        fixtures.append(
            generate_single_a_fixture(
                outdir, "n2", focused_count, args.seed, args.gzip_level, "se"
            )
        )
        fixtures.append(
            generate_single_a_fixture(
                outdir, "r1_only", focused_count, args.seed, args.gzip_level, "pe"
            )
        )
        fixtures.append(
            generate_single_a_fixture(
                outdir, "r2_only", focused_count, args.seed, args.gzip_level, "pe"
            )
        )
        generate_negative_fixtures(outdir, args.gzip_level)
        if args.with_cram:
            generate_cram_fixtures(outdir, fixtures)
        write_manifest(
            outdir,
            args.profile,
            records,
            pairs,
            focused_count,
            args.seed,
            source_mode,
            args.with_cram,
            args.sample_distribution,
            fixtures,
        )
    except Exception:
        (outdir / "GENERATION_INCOMPLETE").write_text(
            "Dataset generation failed; remove or regenerate this directory.\n",
            encoding="ascii",
        )
        raise

    print(f"Generated deterministic {args.profile} fixture at {outdir}")
    print(f"  core SE records: {records:,}")
    print(f"  core PE events:  {pairs:,}")
    print(f"  focused records: {focused_count:,} per fixture")


if __name__ == "__main__":
    main()
