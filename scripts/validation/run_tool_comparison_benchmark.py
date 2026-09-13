#!/usr/bin/env python3
"""Benchmark Plexless vs fqtk, with fastp as a non-demux I/O reference."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import multiprocessing
import os
import platform
import random
import shutil
import statistics
import subprocess
import time
import tomllib
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


REPO = Path(__file__).resolve().parents[2]
DEFAULT_ENV = Path.home() / ".local/share/miniforge3/envs/plexless-compare/bin"
TIME_FORMAT = "%e\t%P\t%M"
VALIDATION_POOLS: dict[tuple[int, str], bytes] = {}
VALIDATION_FIXTURE: Fixture | None = None
METRIC_FIELDS = (
    "phase",
    "tool",
    "samples",
    "threads",
    "repetition",
    "pairs",
    "seconds",
    "pairs_per_second",
    "cpu_percent",
    "peak_rss_kib",
    "fastq_output_bytes",
    "total_output_bytes",
    "fastq_output_files",
    "total_output_files",
    "validation_seconds",
    "system_gzip_tested",
    "canonical_sha256",
)


@dataclass(frozen=True)
class Fixture:
    directory: Path
    pairs: int
    samples: int
    read_length: int
    raw_r1_length: int
    barcode_length: int
    skip_length: int
    pool_size: int
    plexless_structure: str
    fqtk_r1_structure: str
    fqtk_r2_structure: str
    max_mismatches: int
    fqtk_min_mismatch_delta: int

    @classmethod
    def load(cls, directory: Path) -> "Fixture":
        manifest = json.loads((directory / "manifest.json").read_text())
        return cls(
            directory=directory,
            pairs=int(manifest["pairs"]),
            samples=int(manifest["samples"]),
            read_length=int(manifest["read_length"]),
            raw_r1_length=int(manifest["raw_r1_length"]),
            barcode_length=int(manifest["barcode_length"]),
            skip_length=int(manifest["technical_skip_length"]),
            pool_size=int(manifest["sequence_pool_size"]),
            plexless_structure=str(manifest["plexless_structure"]),
            fqtk_r1_structure=str(manifest["fqtk_r1_structure"]),
            fqtk_r2_structure=str(manifest["fqtk_r2_structure"]),
            max_mismatches=int(manifest["max_mismatches"]),
            fqtk_min_mismatch_delta=int(manifest["fqtk_min_mismatch_delta"]),
        )


@dataclass(frozen=True)
class FastqRecord:
    header: bytes
    sequence: bytes
    quality: bytes


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be greater than zero")
    return parsed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument(
        "--plexless-bin", type=Path, default=REPO / "target/release/plexless"
    )
    parser.add_argument("--fqtk-bin", type=Path, default=DEFAULT_ENV / "fqtk")
    parser.add_argument("--fastp-bin", type=Path, default=DEFAULT_ENV / "fastp")
    parser.add_argument(
        "--sample-counts", nargs="+", type=positive_int, default=(96, 384, 1500)
    )
    parser.add_argument(
        "--shared-threads", nargs="+", type=positive_int, default=(5, 8, 12, 16, 24)
    )
    parser.add_argument(
        "--plexless-extra-threads", nargs="+", type=positive_int, default=(1, 2, 4)
    )
    parser.add_argument("--repeats", type=positive_int, default=5)
    parser.add_argument("--compression-level", type=int, choices=range(1, 10), default=2)
    parser.add_argument("--fastp-input-samples", type=positive_int, default=384)
    parser.add_argument(
        "--validation-processes",
        type=positive_int,
        default=min(16, os.cpu_count() or 1),
    )
    parser.add_argument(
        "--phases",
        nargs="+",
        choices=("correctness", "comparison", "plexless-scaling", "fastp"),
        default=("correctness", "comparison", "plexless-scaling", "fastp"),
    )
    parser.add_argument("--skip-system-gzip", action="store_true")
    parser.add_argument("--keep-outputs", action="store_true")
    return parser.parse_args()


def run_text(command: list[str], cwd: Path | None = None) -> str:
    return subprocess.run(
        command,
        cwd=cwd,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    ).stdout.strip()


def executable_version(path: Path, arguments: list[str]) -> str:
    return run_text([str(path), *arguments]).splitlines()[0]


def collect_environment(args: argparse.Namespace) -> dict[str, object]:
    cargo = tomllib.loads((REPO / "Cargo.toml").read_text())
    return {
        "captured_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "platform": platform.platform(),
        "kernel": platform.release(),
        "cpu": run_text(["lscpu"]),
        "memory": run_text(["free", "-h"]),
        "filesystem": run_text(
            ["findmnt", "-T", str(args.dataset), "-no", "SOURCE,FSTYPE,OPTIONS"]
        ),
        "plexless_version": cargo["package"]["version"],
        "plexless_git": run_text(["git", "describe", "--tags", "--always", "--dirty"], REPO),
        "plexless_ldd": run_text(["ldd", str(args.plexless_bin)]),
        "fqtk_version": executable_version(args.fqtk_bin, ["--version"]),
        "fqtk_ldd": run_text(["ldd", str(args.fqtk_bin)]),
        "fastp_version": executable_version(args.fastp_bin, ["--version"]),
        "fastp_ldd": run_text(["ldd", str(args.fastp_bin)]),
        "compression_level": args.compression_level,
        "cache_policy": (
            "warm: both compressed inputs are read in full immediately before every "
            "measured process; warming and validation are excluded from runtime"
        ),
        "fqtk_output": "BGZF through pooled-writer/bgzf/libdeflater",
        "plexless_output": "multi-member gzip through flate2 and zlib-rs",
        "fastp_output": "gzip through fastp with dynamic libdeflate/ISA-L available",
    }


def load_pools(dataset: Path) -> dict[tuple[int, str], bytes]:
    root = json.loads((dataset / "manifest.json").read_text())
    paths = root["pool_files"]
    return {
        (mate, field): (dataset / paths[f"r{mate}_{field}"]).read_bytes()
        for mate in (1, 2)
        for field in ("sequence", "quality")
    }


def pool_record(
    pools: dict[tuple[int, str], bytes], fixture: Fixture, mate: int, ordinal: int
) -> tuple[bytes, bytes]:
    pool_index = (ordinal - 1) % fixture.pool_size
    start = pool_index * fixture.read_length
    end = start + fixture.read_length
    return pools[(mate, "sequence")][start:end], pools[(mate, "quality")][start:end]


def fastq_records(path: Path) -> Iterator[FastqRecord]:
    with gzip.open(path, "rb") as handle:
        record_number = 0
        while True:
            header = handle.readline()
            if not header:
                return
            record_number += 1
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not sequence or not plus or not quality:
                raise RuntimeError(f"{path}: truncated FASTQ record {record_number}")
            header = header.rstrip(b"\r\n")
            sequence = sequence.rstrip(b"\r\n")
            plus = plus.rstrip(b"\r\n")
            quality = quality.rstrip(b"\r\n")
            if not header.startswith(b"@"):
                raise RuntimeError(f"{path}: malformed header at record {record_number}")
            if plus != b"+":
                raise RuntimeError(f"{path}: malformed plus line at record {record_number}")
            if len(sequence) != len(quality):
                raise RuntimeError(
                    f"{path}: sequence/quality mismatch at record {record_number}"
                )
            yield FastqRecord(header, sequence, quality)


def canonical_ordinal(header: bytes, barcode_length: int) -> int:
    token = header.removeprefix(b"@").split(maxsplit=1)[0]
    suffix = token.rsplit(b":", 1)
    if (
        len(suffix) == 2
        and len(suffix[1]) == barcode_length
        and set(suffix[1]).issubset(set(b"ACGTNRYSWKMBDHV"))
    ):
        token = suffix[0]
    if token.endswith((b"/1", b"/2")):
        token = token[:-2]
    if len(token) != 20 or not token.startswith(b"plexcmp:") or not token[8:].isdigit():
        raise RuntimeError(f"cannot canonicalize FASTQ header {header!r}")
    return int(token[8:])


def expected_data_files(tool: str, output: Path, fixture: Fixture) -> list[Path]:
    paths = []
    for sample in range(1, fixture.samples + 1):
        for mate in (1, 2):
            if tool == "plexless":
                name = f"sample_{sample:04d}_R{mate}.fastq.gz"
            elif tool == "fqtk":
                name = f"sample_{sample:04d}.R{mate}.fq.gz"
            else:
                raise ValueError(tool)
            paths.append(output / name)
    return paths


def update_canonical_digest(
    digest: hashlib._Hash,
    ordinal: int,
    mate: int,
    sequence: bytes,
    quality: bytes,
) -> None:
    digest.update(f"{ordinal}\tR{mate}\t".encode("ascii"))
    digest.update(sequence)
    digest.update(b"\t")
    digest.update(quality)
    digest.update(b"\n")


def configure_validation_workers(
    pools: dict[tuple[int, str], bytes], fixture: Fixture
) -> None:
    global VALIDATION_FIXTURE, VALIDATION_POOLS
    VALIDATION_POOLS = pools
    VALIDATION_FIXTURE = fixture


def validate_demux_file_task(
    task: tuple[Path, int, int]
) -> tuple[int, int, int, str]:
    path, sample, mate = task
    fixture = VALIDATION_FIXTURE
    if fixture is None:
        raise RuntimeError("validation worker was not configured")
    digest = hashlib.sha256()
    count = 0
    for record in fastq_records(path):
        ordinal = sample + count * fixture.samples
        if ordinal > fixture.pairs:
            raise RuntimeError(f"{path}: contains extra record {count + 1}")
        observed_ordinal = canonical_ordinal(record.header, fixture.barcode_length)
        if observed_ordinal != ordinal:
            raise RuntimeError(
                f"{path}: record {count + 1} ID is {observed_ordinal}, expected {ordinal}"
            )
        sequence, quality = pool_record(VALIDATION_POOLS, fixture, mate, ordinal)
        if record.sequence != sequence:
            raise RuntimeError(f"{path}: sequence mismatch for ordinal {ordinal}")
        if record.quality != quality:
            raise RuntimeError(f"{path}: quality mismatch for ordinal {ordinal}")
        update_canonical_digest(digest, ordinal, mate, sequence, quality)
        count += 1
    expected_count = (fixture.pairs - sample) // fixture.samples + 1
    if count != expected_count:
        raise RuntimeError(f"{path}: {count} records, expected {expected_count}")
    return sample, mate, count, digest.hexdigest()


def validate_fastp_file_task(task: tuple[Path, int]) -> tuple[int, int, str]:
    path, mate = task
    fixture = VALIDATION_FIXTURE
    if fixture is None:
        raise RuntimeError("validation worker was not configured")
    digest = hashlib.sha256()
    count = 0
    for count, record in enumerate(fastq_records(path), start=1):
        if canonical_ordinal(record.header, fixture.barcode_length) != count:
            raise RuntimeError(f"{path}: identity mismatch at record {count}")
        sequence, quality = pool_record(VALIDATION_POOLS, fixture, mate, count)
        if record.sequence != sequence or record.quality != quality:
            raise RuntimeError(f"{path}: content mismatch at record {count}")
        update_canonical_digest(digest, count, mate, sequence, quality)
    if count != fixture.pairs:
        raise RuntimeError(f"{path}: {count} records, expected {fixture.pairs}")
    return mate, count, digest.hexdigest()


def parallel_map(function, tasks: list, processes: int):
    if processes == 1 or len(tasks) == 1:
        return [function(task) for task in tasks]
    context = multiprocessing.get_context("fork")
    with context.Pool(processes=min(processes, len(tasks))) as pool:
        return pool.map(function, tasks)


def validate_demux_output(
    tool: str,
    output: Path,
    fixture: Fixture,
    pools: dict[tuple[int, str], bytes],
    processes: int,
) -> tuple[str, list[Path]]:
    expected = expected_data_files(tool, output, fixture)
    data_paths = list(expected)
    if tool == "fqtk":
        unmatched = [output / "unmatched.R1.fq.gz", output / "unmatched.R2.fq.gz"]
        for path in unmatched:
            if next(fastq_records(path), None) is not None:
                raise RuntimeError(f"{path}: expected an empty unmatched output")
        data_paths.extend(unmatched)
    expected_names = {path.name for path in data_paths}
    suffix = "*.fastq.gz" if tool == "plexless" else "*.fq.gz"
    actual_names = {path.name for path in output.glob(suffix)}
    if actual_names != expected_names:
        raise RuntimeError(
            f"{tool}: FASTQ inventory mismatch: missing={sorted(expected_names-actual_names)[:5]}, "
            f"extra={sorted(actual_names-expected_names)[:5]}"
        )
    tasks = []
    for sample in range(1, fixture.samples + 1):
        for mate in (1, 2):
            path = expected[(sample - 1) * 2 + mate - 1]
            tasks.append((path, sample, mate))
    configure_validation_workers(pools, fixture)
    results = parallel_map(validate_demux_file_task, tasks, processes)
    digest = hashlib.sha256()
    total = 0
    for sample, mate, count, file_digest in results:
        digest.update(f"{sample}\tR{mate}\t".encode("ascii"))
        digest.update(bytes.fromhex(file_digest))
        total += count
    if total != fixture.pairs * 2:
        raise RuntimeError(f"{tool}: {total} mate records, expected {fixture.pairs * 2}")
    return digest.hexdigest(), data_paths


def validate_fastp_output(
    output: Path,
    fixture: Fixture,
    pools: dict[tuple[int, str], bytes],
    processes: int,
) -> tuple[str, list[Path]]:
    paths = [output / "R1.fastq.gz", output / "R2.fastq.gz"]
    configure_validation_workers(pools, fixture)
    results = parallel_map(
        validate_fastp_file_task,
        [(path, mate) for mate, path in enumerate(paths, start=1)],
        processes,
    )
    digest = hashlib.sha256()
    for mate, _, file_digest in results:
        digest.update(f"R{mate}\t".encode("ascii"))
        digest.update(bytes.fromhex(file_digest))
    return digest.hexdigest(), paths


def system_gzip_test(paths: list[Path]) -> None:
    if shutil.which("gzip") is None:
        raise RuntimeError("system gzip is unavailable")
    for start in range(0, len(paths), 128):
        subprocess.run(
            ["gzip", "-t", "--", *map(str, paths[start : start + 128])], check=True
        )


def warm_inputs(fixture: Fixture) -> None:
    for name in ("R1.fastq.gz", "R2.fastq.gz"):
        with (fixture.directory / name).open("rb") as handle:
            while handle.read(8 * 1024 * 1024):
                pass


def plexless_command(
    args: argparse.Namespace,
    fixture: Fixture,
    output: Path,
    threads: int,
    write_unassigned: bool = False,
) -> list[str]:
    command = [
        str(args.plexless_bin),
        "--threads",
        str(threads),
        "demux",
        "--r1",
        str(fixture.directory / "R1.fastq.gz"),
        "--r2",
        str(fixture.directory / "R2.fastq.gz"),
        "--r1-structure",
        fixture.plexless_structure,
        "--barcodes",
        str(fixture.directory / "barcodes.tsv"),
        "--samples",
        str(fixture.directory / "samples.tsv"),
        "--output",
        str(output),
        "--compression-level",
        str(args.compression_level),
        "--output-mode",
        "auto",
        "--max-mismatches",
        str(fixture.max_mismatches),
    ]
    if write_unassigned:
        command.append("--write-unassigned")
    return command


def fqtk_command(
    args: argparse.Namespace,
    fixture: Fixture,
    output: Path,
    threads: int,
) -> list[str]:
    return [
        str(args.fqtk_bin),
        "demux",
        "--inputs",
        str(fixture.directory / "R1.fastq.gz"),
        str(fixture.directory / "R2.fastq.gz"),
        "--read-structures",
        fixture.fqtk_r1_structure,
        fixture.fqtk_r2_structure,
        "--sample-metadata",
        str(fixture.directory / "fqtk_samples.tsv"),
        "--output",
        str(output),
        "--max-mismatches",
        str(fixture.max_mismatches),
        "--min-mismatch-delta",
        str(fixture.fqtk_min_mismatch_delta),
        "--threads",
        str(threads),
        "--compression-level",
        str(args.compression_level),
    ]


def fastp_command(
    args: argparse.Namespace,
    fixture: Fixture,
    output: Path,
    threads: int,
) -> list[str]:
    return [
        str(args.fastp_bin),
        "--in1",
        str(fixture.directory / "R1.fastq.gz"),
        "--in2",
        str(fixture.directory / "R2.fastq.gz"),
        "--out1",
        str(output / "R1.fastq.gz"),
        "--out2",
        str(output / "R2.fastq.gz"),
        "--json",
        str(output / "fastp.json"),
        "--html",
        str(output / "fastp.html"),
        "--disable_adapter_trimming",
        "--disable_quality_filtering",
        "--disable_length_filtering",
        "--disable_trim_poly_g",
        "--dont_eval_duplication",
        "--trim_front1",
        str(fixture.barcode_length + fixture.skip_length),
        "--trim_front2",
        "0",
        "--thread",
        str(threads),
        "--compression",
        str(args.compression_level),
    ]


def output_stats(output: Path, data_paths: list[Path]) -> tuple[int, int, int, int]:
    all_files = [path for path in output.rglob("*") if path.is_file()]
    return (
        sum(path.stat().st_size for path in data_paths),
        sum(path.stat().st_size for path in all_files),
        len(data_paths),
        len(all_files),
    )


def parse_time_file(path: Path) -> tuple[float, float, int]:
    fields = path.read_text().strip().split("\t")
    if len(fields) != 3:
        raise RuntimeError(f"unexpected GNU time output in {path}: {fields}")
    return float(fields[0]), float(fields[1].rstrip("%")), int(fields[2])


def completed_keys(metrics_path: Path) -> set[tuple[str, str, int, int, int]]:
    if not metrics_path.is_file():
        return set()
    with metrics_path.open(newline="") as handle:
        return {
            (
                row["phase"],
                row["tool"],
                int(row["samples"]),
                int(row["threads"]),
                int(row["repetition"]),
            )
            for row in csv.DictReader(handle, delimiter="\t")
        }


def append_metric(path: Path, metric: dict[str, object]) -> None:
    exists = path.is_file()
    with path.open("a", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=METRIC_FIELDS, delimiter="\t")
        if not exists:
            writer.writeheader()
        writer.writerow(metric)
        handle.flush()
        os.fsync(handle.fileno())


def append_command(path: Path, payload: dict[str, object]) -> None:
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")


def execute_job(
    args: argparse.Namespace,
    phase: str,
    tool: str,
    fixture: Fixture,
    threads: int,
    repetition: int,
    pools: dict[tuple[int, str], bytes],
    results: Path,
) -> dict[str, object]:
    label = f"{phase}_{tool}_{fixture.samples}s_{threads}t_r{repetition}"
    output = results / "work" / label
    if output.exists():
        shutil.rmtree(output)
    output.mkdir(parents=True)
    logs = results / "logs"
    logs.mkdir(exist_ok=True)
    time_path = logs / f"{label}.time.tsv"
    stdout_path = logs / f"{label}.stdout"
    stderr_path = logs / f"{label}.stderr"

    if tool == "plexless":
        command = plexless_command(args, fixture, output, threads)
    elif tool == "fqtk":
        command = fqtk_command(args, fixture, output, threads)
    elif tool == "fastp":
        command = fastp_command(args, fixture, output, threads)
    else:
        raise ValueError(tool)

    warm_inputs(fixture)
    timed = ["/usr/bin/time", "-f", TIME_FORMAT, "-o", str(time_path), *command]
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        result = subprocess.run(timed, stdout=stdout, stderr=stderr)
    if result.returncode != 0:
        raise RuntimeError(
            f"{label} failed with exit {result.returncode}; output retained at {output}"
        )
    seconds, cpu_percent, peak_rss = parse_time_file(time_path)

    validation_start = time.perf_counter()
    if tool in ("plexless", "fqtk"):
        canonical_digest, data_paths = validate_demux_output(
            tool, output, fixture, pools, args.validation_processes
        )
    else:
        canonical_digest, data_paths = validate_fastp_output(
            output, fixture, pools, args.validation_processes
        )
    gzip_tested = repetition == 1 and not args.skip_system_gzip
    if gzip_tested:
        system_gzip_test(data_paths)
    validation_seconds = time.perf_counter() - validation_start
    fastq_bytes, total_bytes, fastq_files, total_files = output_stats(
        output, data_paths
    )
    payload = {
        "phase": phase,
        "tool": tool,
        "samples": fixture.samples,
        "threads": threads,
        "repetition": repetition,
        "pairs": fixture.pairs,
        "seconds": f"{seconds:.6f}",
        "pairs_per_second": f"{fixture.pairs / seconds:.3f}",
        "cpu_percent": f"{cpu_percent:.1f}",
        "peak_rss_kib": peak_rss,
        "fastq_output_bytes": fastq_bytes,
        "total_output_bytes": total_bytes,
        "fastq_output_files": fastq_files,
        "total_output_files": total_files,
        "validation_seconds": f"{validation_seconds:.6f}",
        "system_gzip_tested": str(gzip_tested).lower(),
        "canonical_sha256": canonical_digest,
    }
    append_command(
        results / "commands.jsonl",
        {"label": label, "command": command, "metric": payload},
    )
    if not args.keep_outputs:
        shutil.rmtree(output)
    return payload


def correctness_case(zero_based: int) -> str:
    return ("exact", "mismatch", "unmatched")[zero_based % 3]


def expected_correctness_ordinals(
    fixture: Fixture, sample: int | None
) -> list[int]:
    ordinals = []
    for zero_based in range(fixture.pairs):
        case = correctness_case(zero_based)
        expected_sample = zero_based % fixture.samples + 1
        if (sample is None and case == "unmatched") or (
            sample == expected_sample and case != "unmatched"
        ):
            ordinals.append(zero_based + 1)
    return ordinals


def validate_correctness_file(
    path: Path,
    tool: str,
    fixture: Fixture,
    pools: dict[tuple[int, str], bytes],
    mate: int,
    ordinals: list[int],
    unmatched: bool,
    input_records: dict[int, list[FastqRecord]],
) -> tuple[list[bytes], list[int]]:
    raw_headers = []
    lengths = []
    records = list(fastq_records(path))
    if len(records) != len(ordinals):
        raise RuntimeError(f"{path}: {len(records)} records, expected {len(ordinals)}")
    for record, ordinal in zip(records, ordinals, strict=True):
        if canonical_ordinal(record.header, fixture.barcode_length) != ordinal:
            raise RuntimeError(f"{path}: identity mismatch for ordinal {ordinal}")
        sequence, quality = pool_record(pools, fixture, mate, ordinal)
        if unmatched and tool == "plexless":
            source = input_records[mate][ordinal - 1]
            sequence, quality = source.sequence, source.quality
        if record.sequence != sequence or record.quality != quality:
            raise RuntimeError(f"{path}: content mismatch for ordinal {ordinal}")
        raw_headers.append(record.header)
        lengths.append(len(record.sequence))
    return raw_headers, lengths


def run_correctness(
    args: argparse.Namespace,
    fixture: Fixture,
    pools: dict[tuple[int, str], bytes],
    results: Path,
) -> dict[str, object]:
    root = results / "correctness"
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)
    outputs = {"plexless": root / "plexless", "fqtk": root / "fqtk"}
    commands = {
        "plexless": plexless_command(args, fixture, outputs["plexless"], 5, True),
        "fqtk": fqtk_command(args, fixture, outputs["fqtk"], 5),
    }
    for tool, command in commands.items():
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    input_records = {
        mate: list(fastq_records(fixture.directory / f"R{mate}.fastq.gz"))
        for mate in (1, 2)
    }

    gzip_paths = []
    fqtk_r1_header = b""
    fqtk_r2_header = b""
    unmatched_lengths: dict[str, dict[str, list[int]]] = defaultdict(dict)
    for tool, output in outputs.items():
        for sample in range(1, fixture.samples + 1):
            ordinals = expected_correctness_ordinals(fixture, sample)
            for mate in (1, 2):
                name = (
                    f"sample_{sample:04d}_R{mate}.fastq.gz"
                    if tool == "plexless"
                    else f"sample_{sample:04d}.R{mate}.fq.gz"
                )
                path = output / name
                headers, _ = validate_correctness_file(
                    path,
                    tool,
                    fixture,
                    pools,
                    mate,
                    ordinals,
                    False,
                    input_records,
                )
                gzip_paths.append(path)
                if tool == "fqtk" and headers:
                    if mate == 1 and not fqtk_r1_header:
                        fqtk_r1_header = headers[0]
                    elif mate == 2 and not fqtk_r2_header:
                        fqtk_r2_header = headers[0]
        unmatched_ordinals = expected_correctness_ordinals(fixture, None)
        for mate in (1, 2):
            name = (
                f"unassigned_R{mate}.fastq.gz"
                if tool == "plexless"
                else f"unmatched.R{mate}.fq.gz"
            )
            path = output / name
            headers, lengths = validate_correctness_file(
                path,
                tool,
                fixture,
                pools,
                mate,
                unmatched_ordinals,
                True,
                input_records,
            )
            gzip_paths.append(path)
            unmatched_lengths[tool][f"R{mate}"] = sorted(set(lengths))
            if tool == "fqtk" and headers and not fqtk_r1_header:
                fqtk_r1_header = headers[0]
    if not args.skip_system_gzip:
        system_gzip_test(gzip_paths)
    outcome = {
        "pairs": fixture.pairs,
        "samples": fixture.samples,
        "case_counts": {
            case: sum(correctness_case(index) == case for index in range(fixture.pairs))
            for case in ("exact", "mismatch", "unmatched")
        },
        "assigned_equivalence": "pass",
        "canonicalization": (
            "take the first whitespace-delimited header token, remove one trailing "
            ":<16-IUPAC-base> fqtk annotation when present, remove a terminal /1 or /2, "
            "require plexcmp:<12 decimal digits>; infer mate from the output-file context"
        ),
        "fqtk_example_r1_header": fqtk_r1_header.decode("ascii"),
        "fqtk_example_r2_header": fqtk_r2_header.decode("ascii"),
        "unmatched_lengths": unmatched_lengths,
        "semantic_difference": (
            "fqtk trims B/S segments from unmatched R1; Plexless preserves the complete "
            "unassigned R1 including the barcode and skipped technical bases"
        ),
    }
    (results / "correctness.json").write_text(
        json.dumps(outcome, indent=2, sort_keys=True) + "\n"
    )
    if not args.keep_outputs:
        shutil.rmtree(root)
    return outcome


def read_metrics(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_summary(metrics_path: Path, results: Path) -> None:
    rows = read_metrics(metrics_path)
    groups: dict[tuple[str, str, int, int], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[
            (row["phase"], row["tool"], int(row["samples"]), int(row["threads"]))
        ].append(row)
    summary_fields = (
        "phase",
        "tool",
        "samples",
        "threads",
        "repeats",
        "median_seconds",
        "min_seconds",
        "max_seconds",
        "median_pairs_per_second",
        "median_cpu_percent",
        "median_peak_rss_kib",
        "median_fastq_output_bytes",
        "fastq_output_files",
    )
    summary_rows = []
    for key, members in sorted(groups.items()):
        phase, tool, samples, threads = key
        seconds = [float(row["seconds"]) for row in members]
        summary_rows.append(
            {
                "phase": phase,
                "tool": tool,
                "samples": samples,
                "threads": threads,
                "repeats": len(members),
                "median_seconds": f"{statistics.median(seconds):.6f}",
                "min_seconds": f"{min(seconds):.6f}",
                "max_seconds": f"{max(seconds):.6f}",
                "median_pairs_per_second": f"{statistics.median(float(row['pairs_per_second']) for row in members):.3f}",
                "median_cpu_percent": f"{statistics.median(float(row['cpu_percent']) for row in members):.1f}",
                "median_peak_rss_kib": f"{statistics.median(int(row['peak_rss_kib']) for row in members):.0f}",
                "median_fastq_output_bytes": f"{statistics.median(int(row['fastq_output_bytes']) for row in members):.0f}",
                "fastq_output_files": members[0]["fastq_output_files"],
            }
        )
    with (results / "summary.tsv").open("w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)


def fixture_path(dataset: Path, sample_count: int) -> Path:
    return dataset / "performance" / f"samples_{sample_count:04d}"


def main() -> None:
    args = parse_args()
    args.dataset = args.dataset.expanduser().resolve()
    args.results = args.results.expanduser().resolve()
    for executable in (args.plexless_bin, args.fqtk_bin, args.fastp_bin):
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise SystemExit(f"executable not found: {executable}")
    if not (args.dataset / "manifest.json").is_file():
        raise SystemExit(f"not a comparison dataset: {args.dataset}")
    required_sample_counts = set(args.sample_counts)
    if "fastp" in args.phases:
        required_sample_counts.add(args.fastp_input_samples)
    fixtures = {
        samples: Fixture.load(fixture_path(args.dataset, samples))
        for samples in required_sample_counts
    }
    if any(fixture.pairs != next(iter(fixtures.values())).pairs for fixture in fixtures.values()):
        raise SystemExit("all performance fixtures must contain the same pair count")
    args.results.mkdir(parents=True, exist_ok=True)
    (args.results / "work").mkdir(exist_ok=True)
    environment_path = args.results / "environment.json"
    if not environment_path.exists():
        environment_path.write_text(
            json.dumps(collect_environment(args), indent=2, sort_keys=True) + "\n"
        )
    pools = load_pools(args.dataset)

    if "correctness" in args.phases and not (args.results / "correctness.json").is_file():
        correctness = Fixture.load(args.dataset / "correctness")
        print("Running common-denominator correctness probe...", flush=True)
        outcome = run_correctness(args, correctness, pools, args.results)
        print(f"Correctness: {outcome['assigned_equivalence']}", flush=True)

    jobs: list[tuple[str, str, Fixture, int, int]] = []
    if "comparison" in args.phases:
        for samples in args.sample_counts:
            for threads in args.shared_threads:
                for repetition in range(1, args.repeats + 1):
                    tools = ["plexless", "fqtk"]
                    random.Random(samples * 100_000 + threads * 100 + repetition).shuffle(tools)
                    jobs.extend(
                        ("comparison", tool, fixtures[samples], threads, repetition)
                        for tool in tools
                    )
    if "plexless-scaling" in args.phases:
        for samples in args.sample_counts:
            for threads in args.plexless_extra_threads:
                for repetition in range(1, args.repeats + 1):
                    jobs.append(
                        (
                            "plexless-scaling",
                            "plexless",
                            fixtures[samples],
                            threads,
                            repetition,
                        )
                    )
    if "fastp" in args.phases:
        fixture = fixtures[args.fastp_input_samples]
        for threads in args.shared_threads:
            for repetition in range(1, args.repeats + 1):
                jobs.append(("fastp-reference", "fastp", fixture, threads, repetition))

    metrics_path = args.results / "metrics.tsv"
    done = completed_keys(metrics_path)
    remaining = [
        job
        for job in jobs
        if (job[0], job[1], job[2].samples, job[3], job[4]) not in done
    ]
    print(f"Benchmark jobs remaining: {len(remaining)} of {len(jobs)}", flush=True)
    for index, (phase, tool, fixture, threads, repetition) in enumerate(
        remaining, start=1
    ):
        print(
            f"[{index}/{len(remaining)}] {tool} {fixture.samples} samples "
            f"{threads} threads repetition {repetition}",
            flush=True,
        )
        metric = execute_job(
            args,
            phase,
            tool,
            fixture,
            threads,
            repetition,
            pools,
            args.results,
        )
        append_metric(metrics_path, metric)
        print(
            f"  {metric['seconds']} s, {metric['pairs_per_second']} pairs/s, "
            f"validation {metric['validation_seconds']} s",
            flush=True,
        )
    if metrics_path.is_file():
        write_summary(metrics_path, args.results)
    print(f"Benchmark artifacts: {args.results}")


if __name__ == "__main__":
    main()
