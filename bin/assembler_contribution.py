#!/usr/bin/env python3
"""Which assemblers' contigs made it into each QC-pass Autocycler cluster.

Reads Autocycler v0.7.0's own records in ``autocycler_out/``:

- ``clustering/clustering.tsv``: every input contig, its file and its QC-pass cluster;
- ``clustering/qc_pass/cluster_NNN/2_trimmed.gfa``: the contigs ``trim`` kept, as ``P``
  lines tagged ``FN:Z:<file>`` and ``HD:Z:<header>``.

Input files are named ``<assembler>_<subset>.fasta`` by NORMALISE_HEADERS. An assembler
that routinely reaches few subsets of a cluster, or is routinely trimmed out, is a
candidate for removal from the defaults.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


def assembler_and_subset(filename: str) -> tuple[str, str]:
    stem = Path(filename).name.split(".")[0]
    assembler, _, subset = stem.rpartition("_")
    return (assembler, subset) if assembler else (stem, "")


def trimmed_contigs(gfa: Path) -> set[tuple[str, str]]:
    """(file, contig name) pairs that survived ``autocycler trim``."""
    kept = set()
    if not gfa.exists():
        return kept
    for line in gfa.read_text().splitlines():
        if not line.startswith("P\t"):
            continue
        tags = dict(
            (field[:2], field[5:])
            for field in line.split("\t")[4:]
            if re.match(r"^[A-Z]{2}:Z:", field)
        )
        if "FN" in tags and "HD" in tags:
            kept.add((Path(tags["FN"]).name, tags["HD"].split()[0]))
    return kept


def contribution(autocycler_dir: Path) -> list[dict[str, object]]:
    clustering = autocycler_dir / "clustering" / "clustering.tsv"
    if not clustering.exists():
        return []
    with clustering.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    subsets: dict[str, set[str]] = defaultdict(set)
    members: dict[int, list[tuple[str, str, str, str]]] = defaultdict(list)
    for row in rows:
        filename = Path(row["file_name"]).name
        assembler, subset = assembler_and_subset(filename)
        subsets[assembler].add(subset)
        if row["passing_clusters"].isdigit():
            members[int(row["passing_clusters"])].append(
                (assembler, subset, filename, row["contig_name"])
            )

    results = []
    for cluster, contigs in sorted(members.items()):
        gfa = (
            autocycler_dir
            / "clustering"
            / "qc_pass"
            / f"cluster_{cluster:03d}"
            / "2_trimmed.gfa"
        )
        kept = trimmed_contigs(gfa)
        for assembler in sorted(subsets):
            mine = [c for c in contigs if c[0] == assembler]
            results.append(
                {
                    "assembler": assembler,
                    "cluster": cluster,
                    "subsets_total": len(subsets[assembler]),
                    "subsets_in_cluster": len({c[1] for c in mine}),
                    "contigs_in_cluster": len(mine),
                    # Without a trimmed graph (trim failed) nothing can be counted.
                    "contigs_trimmed_out": sum((c[2], c[3]) not in kept for c in mine)
                    if gfa.exists()
                    else "",
                }
            )
    return results


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--autocycler-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    fields = [
        "assembler",
        "cluster",
        "subsets_total",
        "subsets_in_cluster",
        "contigs_in_cluster",
        "contigs_trimmed_out",
    ]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["sample", *fields], delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(
            {"sample": args.sample} | row for row in contribution(args.autocycler_dir)
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
