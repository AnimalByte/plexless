#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

KEYS = ("total", "assigned", "unmatched", "ambiguous", "unrouted")


def parse_expected(path: Path) -> dict[str, int]:
    expected: dict[str, int] = {}
    lines = path.read_text().splitlines()

    if not lines or lines[0] != "Category\tExpected":
        raise SystemExit(f"Unexpected expected-summary format: {path}")

    for line in lines[1:]:
        key, value = line.split("\t")
        expected[key] = int(value)

    return expected


def parse_summary(path: Path) -> dict[str, int]:
    text = path.read_text()
    observed: dict[str, int] = {}

    for key in KEYS:
        match = re.search(
            rf"^\s*{re.escape(key)}:\s+(\d+)\s*$",
            text,
            flags=re.MULTILINE,
        )
        if not match:
            raise SystemExit(f"Could not find '{key}' in simplex summary: {path}")
        observed[key] = int(match.group(1))

    return observed


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("expected", type=Path)
    p.add_argument("summary", type=Path)
    args = p.parse_args()

    expected = parse_expected(args.expected)
    observed = parse_summary(args.summary)

    failed = False
    print("Category      Expected      Observed      Result")
    print("------------------------------------------------")

    for key in KEYS:
        e = expected[key]
        o = observed[key]
        ok = e == o
        failed |= not ok
        print(f"{key:10s} {e:12,d} {o:13,d}      {'PASS' if ok else 'FAIL'}")

    if failed:
        raise SystemExit(1)

    print("\nSummary matches expected truth.")


if __name__ == "__main__":
    main()
