#!/usr/bin/env python3
"""Find regions of uneven read depth, and each replicon's depth relative to the chromosome.

Reads mosdepth's fixed-window depths (``--by``). A window below ``--low`` or above
``--high`` times its contig's median is flagged, and adjacent flagged windows of the same
kind are merged into one region:

- a **low** region suggests a misassembly or an under-represented stretch;
- a **high** region (about 2x) suggests a collapsed repeat, such as rRNA operons.

Circular contigs are mapped as they are: a read spanning the origin is split into a
primary and a supplementary alignment, and mosdepth counts both, so depth at the ends of a
circular contig needs no correction.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import statistics
import sys
from pathlib import Path


def read_windows(path: Path) -> dict[str, list[tuple[int, int, float]]]:
    """mosdepth ``regions.bed.gz`` as {contig: [(start, end, depth), ...]}."""
    windows: dict[str, list[tuple[int, int, float]]] = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 4:
                continue
            windows.setdefault(fields[0], []).append(
                (int(fields[1]), int(fields[2]), float(fields[3]))
            )
    return windows


def read_contigs(path: Path | None) -> dict[str, str]:
    """Replicon type of each contig, from classify_replicons.py's table."""
    if path is None or not path.exists():
        return {}
    with path.open() as handle:
        return {
            row["contig"]: row["replicon_type"]
            for row in csv.DictReader(handle, delimiter="\t")
        }


def flagged_regions(
    windows: list[tuple[int, int, float]], median: float, low: float, high: float
) -> list[dict[str, object]]:
    """Merge runs of adjacent windows that are all low, or all high."""
    regions: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    depths: list[float] = []

    def close() -> None:
        if current is not None:
            current["mean_depth"] = round(statistics.fmean(depths), 2)
            current["depth_ratio"] = round(current["mean_depth"] / median, 3)
            regions.append(current)

    for start, end, depth in windows:
        kind = None
        if depth < low * median:
            kind = "low"
        elif depth > high * median:
            kind = "high"

        if current is not None and kind == current["kind"] and start == current["end"]:
            current["end"] = end
            depths.append(depth)
            continue
        close()
        current, depths = (
            ({"kind": kind, "start": start, "end": end}, [depth])
            if kind
            else (None, [])
        )
    close()
    return regions


def summarise(
    windows: dict[str, list[tuple[int, int, float]]],
    replicon_types: dict[str, str],
    low: float,
    high: float,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Per-contig depth rows and flagged-region rows."""
    medians = {
        contig: statistics.median(depth for _, _, depth in rows)
        for contig, rows in windows.items()
        if rows
    }
    chromosomes = [
        contig
        for contig, kind in replicon_types.items()
        if kind == "chromosome" and contig in medians
    ]
    # The longest chromosome is the reference when a genome has several.
    reference = max(chromosomes, key=lambda c: windows[c][-1][1], default=None)
    reference_depth = medians[reference] if reference else None

    depth_rows: list[dict[str, object]] = []
    region_rows: list[dict[str, object]] = []
    for contig, rows in windows.items():
        median = medians[contig]
        ratio = round(median / reference_depth, 3) if reference_depth else ""
        depth_rows.append(
            {
                "contig": contig,
                "replicon_type": replicon_types.get(contig, ""),
                "length": rows[-1][1],
                "median_depth": median,
                "depth_ratio": ratio,
            }
        )
        if median > 0:
            for region in flagged_regions(rows, median, low, high):
                region_rows.append({"contig": contig} | region)
    return depth_rows, region_rows


def write_tsv(path: Path, sample: str, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["sample", *fields], delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows({"sample": sample} | row for row in rows)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--windows", type=Path, required=True, help="regions.bed.gz")
    parser.add_argument("--contigs", type=Path, help="<id>.contigs.tsv")
    parser.add_argument("--low", type=float, default=0.5)
    parser.add_argument("--high", type=float, default=2.0)
    parser.add_argument("--regions-output", type=Path, required=True)
    parser.add_argument("--depth-output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    depth_rows, region_rows = summarise(
        read_windows(args.windows), read_contigs(args.contigs), args.low, args.high
    )
    write_tsv(
        args.depth_output,
        args.sample,
        ["contig", "replicon_type", "length", "median_depth", "depth_ratio"],
        depth_rows,
    )
    write_tsv(
        args.regions_output,
        args.sample,
        ["contig", "kind", "start", "end", "mean_depth", "depth_ratio"],
        region_rows,
    )
    if region_rows:
        print(f"{len(region_rows)} uneven-depth region(s)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
