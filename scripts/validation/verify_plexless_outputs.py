#!/usr/bin/env python3
"""Verify logical FASTQ output exactly against deterministic fixture truth."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import subprocess
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Iterator, TextIO


TRUTH_HEADER = (
    "OutputFile\tOutputIndex\tOrdinal\tReadID\tMate\tCase\tCategory\t"
    "ExpectedSample\tExpectedHeader\tExpectedSequence\tExpectedQuality"
)
REPORTS = ("sample_metrics.tsv", "fastq_stats.tsv", "barcode_stats.tsv")


@dataclass(frozen=True)
class TruthRecord:
    output_file: str
    output_index: int
    ordinal: int
    read_id: str
    mate: str
    case: str
    category: str
    sample: str
    header: str
    sequence: str
    quality: str


@dataclass
class LogicalFile:
    records: int
    bytes: int
    sha256: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--fastq-directory",
        type=Path,
        help="directory containing validation-decoded FASTQ (defaults to --output)",
    )
    parser.add_argument("--mode", choices=("se", "pe"), required=True)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--hash-manifest", type=Path)
    parser.add_argument("--expect-fastq-stats", action="store_true")
    parser.add_argument("--skip-system-gzip", action="store_true")
    parser.add_argument(
        "--cram-qnames",
        action="store_true",
        help="expect CRAM-derived FASTQ headers normalized to @QNAME",
    )
    return parser.parse_args()


def core_read_id(header: str) -> str:
    token = header.removeprefix("@").split(maxsplit=1)[0]
    if token.endswith(("/1", "/2")):
        return token[:-2]
    return token


def truth_records(path: Path) -> Iterator[TruthRecord]:
    with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
        header = handle.readline().rstrip("\r\n")
        if header != TRUTH_HEADER:
            raise SystemExit(f"Unexpected truth header in {path}: {header!r}")
        previous_file = ""
        previous_index = 0
        for line_number, line in enumerate(handle, start=2):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 11:
                raise SystemExit(
                    f"{path}: truth line {line_number} has {len(fields)} fields, expected 11"
                )
            record = TruthRecord(
                output_file=fields[0],
                output_index=int(fields[1]),
                ordinal=int(fields[2]),
                read_id=fields[3],
                mate=fields[4],
                case=fields[5],
                category=fields[6],
                sample=fields[7],
                header=fields[8],
                sequence=fields[9],
                quality=fields[10],
            )
            if record.output_file < previous_file:
                raise SystemExit(f"{path}: truth output groups are not sorted")
            if record.output_file != previous_file:
                previous_file = record.output_file
                previous_index = 0
            previous_index += 1
            if record.output_index != previous_index:
                raise SystemExit(
                    f"{path}: {record.output_file} truth index is "
                    f"{record.output_index}, expected {previous_index}"
                )
            if len(record.sequence) != len(record.quality):
                raise SystemExit(
                    f"{path}: truth sequence/quality mismatch at line {line_number}"
                )
            if core_read_id(record.header) != record.read_id:
                raise SystemExit(f"{path}: truth read ID mismatch at line {line_number}")
            yield record


def strip_eol(line: bytes) -> bytes:
    return line.rstrip(b"\r\n")


def read_fastq_record(
    handle: BinaryIO,
    path: Path,
    record_number: int,
    digest: hashlib._Hash,
) -> tuple[str, str, str] | None:
    header = handle.readline()
    if header == b"":
        return None
    sequence = handle.readline()
    plus = handle.readline()
    quality = handle.readline()
    if not sequence or not plus or not quality:
        raise SystemExit(f"{path}: truncated FASTQ record {record_number}")
    for part in (header, sequence, plus, quality):
        digest.update(part)
    try:
        header_text = strip_eol(header).decode("ascii")
        sequence_text = strip_eol(sequence).decode("ascii")
        plus_text = strip_eol(plus).decode("ascii")
        quality_text = strip_eol(quality).decode("ascii")
    except UnicodeDecodeError as error:
        raise SystemExit(f"{path}: non-ASCII FASTQ record {record_number}: {error}")
    if not header_text.startswith("@"):
        raise SystemExit(f"{path}: record {record_number} header does not start with @")
    if plus_text != "+":
        raise SystemExit(f"{path}: record {record_number} plus line is not normalized '+'")
    if len(sequence_text) != len(quality_text):
        raise SystemExit(
            f"{path}: record {record_number} sequence/quality length mismatch"
        )
    return header_text, sequence_text, quality_text


def assert_record(
    path: Path,
    expected: TruthRecord,
    observed: tuple[str, str, str],
    cram_qnames: bool,
) -> None:
    header, sequence, quality = observed
    expected_header = f"@{expected.read_id}" if cram_qnames else expected.header
    if header != expected_header:
        raise SystemExit(
            f"{path}: header mismatch at record {expected.output_index}: "
            f"expected {expected_header!r}, observed {header!r}"
        )
    if core_read_id(header) != expected.read_id:
        raise SystemExit(
            f"{path}: read identity mismatch at record {expected.output_index}"
        )
    if sequence != expected.sequence:
        raise SystemExit(
            f"{path}: sequence mismatch at record {expected.output_index} "
            f"({expected.case}, input ordinal {expected.ordinal})"
        )
    if quality != expected.quality:
        raise SystemExit(
            f"{path}: quality mismatch at record {expected.output_index} "
            f"({expected.case}, input ordinal {expected.ordinal})"
        )


def finish_output(
    path: Path,
    handle: BinaryIO,
    digest: hashlib._Hash,
    records: int,
) -> LogicalFile:
    if handle.read(1):
        raise SystemExit(f"{path}: extra content follows expected FASTQ records")
    handle.close()
    # SHA-256 above is computed over the actual decompressed byte stream. Count
    # bytes independently here so the manifest is auditable without trusting
    # the truth fields.
    with gzip.open(path, "rb") as reader:
        byte_count = 0
        while chunk := reader.read(1024 * 1024):
            byte_count += len(chunk)
    return LogicalFile(records=records, bytes=byte_count, sha256=digest.hexdigest())


def verify_against_truth(
    truth_path: Path, output_dir: Path, cram_qnames: bool
) -> tuple[dict[str, LogicalFile], Counter[str], int]:
    logical: dict[str, LogicalFile] = {}
    cases: Counter[str] = Counter()
    expected_files: set[str] = set()
    current_name: str | None = None
    current_path: Path | None = None
    current_handle: BinaryIO | None = None
    current_digest: hashlib._Hash | None = None
    current_records = 0
    total_records = 0

    def close_current() -> None:
        nonlocal current_handle, current_digest, current_records
        if current_name is None or current_path is None:
            return
        assert current_handle is not None and current_digest is not None
        logical[current_name] = finish_output(
            current_path, current_handle, current_digest, current_records
        )
        print(f"PASS  {current_name:<34} {current_records:>12,d} records")
        current_handle = None
        current_digest = None
        current_records = 0

    try:
        for expected in truth_records(truth_path):
            if expected.output_file != current_name:
                close_current()
                current_name = expected.output_file
                current_path = output_dir / current_name
                expected_files.add(current_name)
                if not current_path.is_file():
                    raise SystemExit(f"Missing expected output file: {current_path}")
                current_handle = gzip.open(current_path, "rb")
                current_digest = hashlib.sha256()
            assert current_handle is not None and current_digest is not None
            current_records += 1
            total_records += 1
            cases[expected.case] += 1
            observed = read_fastq_record(
                current_handle, current_path, current_records, current_digest
            )
            if observed is None:
                raise SystemExit(
                    f"{current_path}: expected {expected.output_index} records, "
                    f"file ended after {current_records - 1}"
                )
            assert_record(current_path, expected, observed, cram_qnames)
        close_current()
    except (OSError, EOFError) as error:
        raise SystemExit(f"gzip/FASTQ read failed: {error}")

    observed_files = {path.name for path in output_dir.glob("*.fastq.gz")}
    if observed_files != expected_files:
        raise SystemExit(
            "Unexpected FASTQ output file set: "
            f"missing={sorted(expected_files - observed_files)}, "
            f"extra={sorted(observed_files - expected_files)}"
        )
    return logical, cases, total_records


def logical_sha256(path: Path) -> LogicalFile:
    digest = hashlib.sha256()
    byte_count = 0
    records = 0
    with gzip.open(path, "rb") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not sequence or not plus or not quality:
                raise SystemExit(f"{path}: truncated baseline FASTQ")
            for part in (header, sequence, plus, quality):
                digest.update(part)
                byte_count += len(part)
            records += 1
    return LogicalFile(records=records, bytes=byte_count, sha256=digest.hexdigest())


def compare_baseline(
    output: Path, baseline: Path, current: dict[str, LogicalFile], expect_stats: bool
) -> None:
    baseline_files = {path.name for path in baseline.glob("*.fastq.gz")}
    if baseline_files != set(current):
        raise SystemExit("Baseline FASTQ file set differs from current output")
    for filename in sorted(current):
        baseline_logical = logical_sha256(baseline / filename)
        if baseline_logical != current[filename]:
            raise SystemExit(f"Logical FASTQ differs from baseline: {filename}")
    report_names = REPORTS if expect_stats else REPORTS[:1]
    for report in report_names:
        current_path = output / report
        baseline_path = baseline / report
        if not current_path.is_file() or not baseline_path.is_file():
            raise SystemExit(f"Missing deterministic report for comparison: {report}")
        if current_path.read_bytes() != baseline_path.read_bytes():
            raise SystemExit(f"Deterministic report differs from baseline: {report}")
    print(f"PASS  baseline logical FASTQ SHA-256   {len(current):>12,d} files")
    print(f"PASS  baseline deterministic reports   {len(report_names):>12,d} files")


def verify_sample_metrics(dataset: Path, output: Path, mode: str) -> None:
    expected_path = dataset / f"expected_samples_{mode}.tsv"
    expected_lines = expected_path.read_text(encoding="ascii").splitlines()
    if not expected_lines or expected_lines[0] != "Sample\tExpectedAssigned":
        raise SystemExit(f"Unexpected expected-samples format: {expected_path}")
    expected = {
        sample: int(count)
        for sample, count in (line.split("\t") for line in expected_lines[1:])
    }
    metrics_path = output / "sample_metrics.tsv"
    lines = metrics_path.read_text(encoding="ascii").splitlines()
    if not lines or not lines[0].startswith("sample\tassigned_fragments\t"):
        raise SystemExit(f"Unexpected sample metrics format: {metrics_path}")
    observed = {}
    for line in lines[1:]:
        fields = line.split("\t")
        observed[fields[0]] = int(fields[1])
    if observed != expected:
        raise SystemExit(
            f"sample_metrics.tsv differs: expected={expected}, observed={observed}"
        )
    print(f"PASS  sample metrics reconciliation    {sum(observed.values()):>12,d} assigned")


def system_gzip_test(output: Path, files: dict[str, LogicalFile], skip: bool) -> None:
    if skip:
        print("SKIP  external gzip -t (requested)")
        return
    executable = shutil.which("gzip")
    if executable is None:
        print("SKIP  external gzip -t (gzip is not installed)")
        return
    for filename in sorted(files):
        result = subprocess.run(
            [executable, "-t", str(output / filename)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            raise SystemExit(
                f"gzip -t failed for {filename}: {result.stderr.strip()}"
            )
    print(f"PASS  independent gzip -t              {len(files):>12,d} files")


def write_manifest(
    path: Path,
    dataset: Path,
    output: Path,
    mode: str,
    files: dict[str, LogicalFile],
    cases: Counter[str],
    total_records: int,
) -> None:
    data = {
        "schema_version": 1,
        "dataset": str(dataset.resolve()),
        "output": str(output.resolve()),
        "mode": mode,
        "records_verified": total_records,
        "cases": dict(sorted(cases.items())),
        "logical_fastq": {
            name: {
                "records": value.records,
                "decompressed_bytes": value.bytes,
                "sha256": value.sha256,
            }
            for name, value in sorted(files.items())
        },
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="ascii")


def main() -> None:
    args = parse_args()
    truth_path = args.dataset / f"truth_{args.mode}.tsv.gz"
    if not truth_path.is_file():
        raise SystemExit(f"Missing truth file: {truth_path}")
    if not args.output.is_dir():
        raise SystemExit(f"Missing output directory: {args.output}")
    fastq_directory = args.fastq_directory or args.output
    if not fastq_directory.is_dir():
        raise SystemExit(f"Missing FASTQ directory: {fastq_directory}")
    if (args.output / "PLEXLESS_INCOMPLETE").exists():
        raise SystemExit("Output still contains PLEXLESS_INCOMPLETE")

    expected_reports = REPORTS if args.expect_fastq_stats else REPORTS[:1]
    for report in expected_reports:
        if not (args.output / report).is_file():
            raise SystemExit(f"Missing expected report: {args.output / report}")

    print("Deep logical FASTQ validation:")
    files, cases, total_records = verify_against_truth(
        truth_path, fastq_directory, args.cram_qnames
    )
    verify_sample_metrics(args.dataset, args.output, args.mode)
    system_gzip_test(fastq_directory, files, args.skip_system_gzip)
    if args.baseline is not None:
        compare_baseline(
            args.output, args.baseline, files, args.expect_fastq_stats
        )
    if args.hash_manifest is not None:
        write_manifest(
            args.hash_manifest,
            args.dataset,
            args.output,
            args.mode,
            files,
            cases,
            total_records,
        )
        print(f"PASS  logical hash manifest             {args.hash_manifest}")
    print(
        f"PASS  complete accounting/order/content {total_records:>12,d} records"
    )
    print("Cases: " + ", ".join(f"{key}={value:,}" for key, value in sorted(cases.items())))


if __name__ == "__main__":
    main()
