#!/usr/bin/env python3
"""Compute mean Phred+33 base qualities from FASTQ inputs."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
import sys


def open_fastq(path: Path):
    """Yield a text handle for FASTQ or FASTQ.GZ files."""
    mode = "rt"
    if path.suffixes and path.suffixes[-1].lower() == ".gz":
        return gzip.open(path, mode=mode, newline="")
    return open(path, mode=mode, newline="")


def phred33_sum(quality_line: str) -> int:
    """Return the Phred+33 score for one quality line."""
    return sum((ord(ch) - 33) for ch in quality_line)


def file_quality(path: Path) -> tuple[int, int, int]:
    """Return (total quality, total bases, total reads) for a FASTQ file."""
    total_score = 0
    total_bases = 0
    reads = 0

    with open_fastq(path) as handle:
        for idx, line in enumerate(handle, start=1):
            if idx % 4 != 0:
                continue
            quality = line.rstrip()
            if not quality:
                continue
            total_score += phred33_sum(quality)
            total_bases += len(quality)
            reads += 1

    return total_score, total_bases, reads


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Measure the mean Phred+33 base quality for one or more FASTQ files "
            "(supports gzipped FASTQ)."
        )
    )
    parser.add_argument(
        "fastq",
        type=Path,
        nargs="+",
        help="FASTQ/FASTQ.gz file(s) to summarize",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    total_quality = 0
    total_bases = 0
    total_reads = 0

    for fastq_path in args.fastq:
        if not fastq_path.exists():
            print(f"error: {fastq_path} does not exist", file=sys.stderr)
            return 1

        quality_sum, base_count, reads = file_quality(fastq_path)
        if base_count == 0:
            print(f"warning: {fastq_path} contains no quality lines", file=sys.stderr)
            continue

        mean_quality = quality_sum / base_count
        print(
            f"{fastq_path}: {reads} reads, {base_count} bases, "
            f"mean quality {mean_quality:.2f}"
        )

        total_quality += quality_sum
        total_bases += base_count
        total_reads += reads

    if total_bases:
        overall_mean = total_quality / total_bases
        print(
            f"overall: {total_reads} reads, {total_bases} bases, "
            f"mean quality {overall_mean:.2f}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
