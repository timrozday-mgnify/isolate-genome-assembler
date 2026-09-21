#!/usr/bin/env python3
"""Gene-level integrity: the IDEEL protein length test and the rRNA collapse check.

**IDEEL.** Each Bakta protein's length is divided by the length of its best DIAMOND hit.
An indel error in the assembly shifts the reading frame and truncates the protein, so a
high fraction of proteins below ``--min-ratio`` of their hit points to indel errors.

**rRNA collapse.** Depth over each annotated rRNA locus is compared with its contig's
median depth. Identical rRNA operons that the assembler merged into one show up as
roughly n-times depth over the single copy left.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import statistics
import sys
from pathlib import Path


def read_hits(path: Path | None) -> dict[str, tuple[str, float]]:
    """Best hit per protein from DIAMOND outfmt 6 ``qseqid sseqid qlen slen ...``.

    DIAMOND writes each query's hits best first, so the first line per query wins.
    """
    hits: dict[str, tuple[str, float]] = {}
    if path is None or not path.exists():
        return hits
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < 4 or fields[0] in hits:
            continue
        query_length, subject_length = float(fields[2]), float(fields[3])
        if subject_length > 0:
            hits[fields[0]] = (fields[1], query_length / subject_length)
    return hits


def count_proteins(path: Path | None) -> int:
    if path is None or not path.exists():
        return 0
    return sum(line.startswith(">") for line in path.read_text().splitlines())


def rrna_loci(path: Path | None) -> list[tuple[str, int, int, str]]:
    """(contig, start, end, product) for each rRNA feature in a Bakta GFF3."""
    loci = []
    if path is None or not path.exists():
        return loci
    for line in path.read_text().splitlines():
        if line.startswith("##FASTA"):
            break
        fields = line.split("\t")
        if line.startswith("#") or len(fields) < 9 or fields[2] != "rRNA":
            continue
        attributes = dict(
            item.split("=", 1) for item in fields[8].split(";") if "=" in item
        )
        product = attributes.get("product") or attributes.get("Name", "rRNA")
        loci.append((fields[0], int(fields[3]) - 1, int(fields[4]), product))
    return loci


def read_windows(path: Path | None) -> dict[str, list[tuple[int, int, float]]]:
    windows: dict[str, list[tuple[int, int, float]]] = {}
    if path is None or not path.exists():
        return windows
    with gzip.open(path, "rt") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 4:
                windows.setdefault(fields[0], []).append(
                    (int(fields[1]), int(fields[2]), float(fields[3]))
                )
    return windows


def rrna_depths(
    loci: list[tuple[str, int, int, str]],
    windows: dict[str, list[tuple[int, int, float]]],
) -> list[dict[str, object]]:
    medians = {
        contig: statistics.median(d for _, _, d in rows)
        for contig, rows in windows.items()
        if rows
    }
    rows = []
    for contig, start, end, product in loci:
        overlapping = [
            d for s, e, d in windows.get(contig, []) if s < end and e > start
        ]
        if not overlapping or not medians.get(contig):
            continue
        depth = statistics.fmean(overlapping)
        rows.append(
            {
                "contig": contig,
                "start": start,
                "end": end,
                "product": product,
                "depth": round(depth, 2),
                "depth_ratio": round(depth / medians[contig], 3),
            }
        )
    return rows


def write_tsv(path: Path, sample: str, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["sample", *fields], delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows({"sample": sample} | row for row in rows)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--proteins", type=Path, help="Bakta .faa")
    parser.add_argument("--hits", type=Path, help="DIAMOND outfmt 6")
    parser.add_argument("--gff", type=Path, help="Bakta .gff3")
    parser.add_argument("--windows", type=Path, help="mosdepth regions.bed.gz")
    parser.add_argument("--min-ratio", type=float, default=0.9)
    parser.add_argument("--output", type=Path, required=True, help="per protein")
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--rrna-output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    hits = read_hits(args.hits)
    write_tsv(
        args.output,
        args.sample,
        [
            {"protein": protein, "hit": hit, "ratio": round(ratio, 4)}
            for protein, (hit, ratio) in hits.items()
        ],
        ["protein", "hit", "ratio"],
    )

    truncated = sum(ratio < args.min_ratio for _, ratio in hits.values())
    write_tsv(
        args.summary,
        args.sample,
        [
            {
                "proteins": count_proteins(args.proteins),
                "proteins_with_hit": len(hits),
                "truncated": truncated,
                "truncated_fraction": round(truncated / len(hits), 4) if hits else "",
            }
        ],
        ["proteins", "proteins_with_hit", "truncated", "truncated_fraction"],
    )

    write_tsv(
        args.rrna_output,
        args.sample,
        rrna_depths(rrna_loci(args.gff), read_windows(args.windows)),
        ["contig", "start", "end", "product", "depth", "depth_ratio"],
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
