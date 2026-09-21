#!/usr/bin/env python3
"""Turn a rotated assembly into the finished, named replicons for one sample.

Nothing is dropped silently. Every contig is classified, named and written to
``<id>.contigs.tsv`` with the flags it earned; ``--drop-flagged`` is what actually moves a
flagged contig out of the final FASTA and into ``<id>.removed_contigs.fasta``.

Flags are measurements, not verdicts: ``qc_gates.py`` decides what a flagged contig means
for the sample's status.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

CHROMOSOME = "chromosome"
PLASMID = "plasmid"
UNPLACED = "unplaced"


def load_fasta(path: Path | None) -> list[tuple[str, str]]:
    """Read a FASTA into (header, sequence) pairs, joining multi-line sequences."""
    if path is None or not path.exists():
        return []
    records: list[tuple[str, str]] = []
    header = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header = line[1:]
            chunks = []
        elif header is not None:
            chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return [record for record in records if record[1]]


def tag(header: str, name: str) -> str | None:
    """Read a ``key=value`` tag out of a contig header."""
    match = re.search(rf"\b{name}=([^\s]+)", header)
    return match.group(1) if match else None


def self_overlapping(paf: Path | None) -> set[str]:
    """Contigs whose first 10 kb still aligns to their own last 10 kb.

    Autocycler's `trim` removes that duplication; the Flye fallback path is the one that
    can still carry it, and a circular contig that does is over-long by the overlap.
    """
    overlapping: set[str] = set()
    if paf is None or not paf.exists():
        return overlapping
    for line in paf.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        query, target = fields[0], fields[5]
        if query.endswith("_start") and target.endswith("_end"):
            contig = query[: -len("_start")]
            if contig == target[: -len("_end")]:
                overlapping.add(contig)
    return overlapping


class Contig:
    """One assembled replicon and everything known about it."""

    def __init__(self, header: str, sequence: str) -> None:
        self.original_name = header.split()[0]
        self.sequence = sequence
        self.length = len(sequence)
        self.circular = (tag(header, "circular") or "").lower() == "true"
        depth = tag(header, "depth")
        self.depth = float(depth) if depth else None
        self.flags: list[str] = []
        self.replicon_type = UNPLACED
        self.name = self.original_name

    def header(self) -> str:
        parts = [self.name, f"length={self.length}"]
        if self.depth is not None:
            parts.append(f"depth={self.depth}")
        parts.append(f"circular={'true' if self.circular else 'false'}")
        return " ".join(parts)


def classify(contigs: list[Contig], chromosome_min_len: int) -> None:
    """Call each contig a chromosome, a plasmid, or neither.

    The longest contig is the chromosome even when it is shorter than the threshold: a
    sample whose assembly fell apart still has a chromosome, and the length flag is what
    says the assembly is poor.
    """
    if not contigs:
        return
    longest = max(contigs, key=lambda contig: contig.length)
    for contig in contigs:
        if contig is longest or contig.length >= chromosome_min_len:
            contig.replicon_type = CHROMOSOME
        elif contig.circular:
            contig.replicon_type = PLASMID
        else:
            contig.replicon_type = UNPLACED


def flag(
    contigs: list[Contig],
    min_contig_len: int,
    min_depth_ratio: float,
    unrotated: set[str],
    overlapping: set[str],
) -> None:
    """Attach every flag a contig earns. Depth is relative to the chromosome."""
    chromosome_depths = [
        c.depth for c in contigs if c.replicon_type == CHROMOSOME and c.depth is not None
    ]
    chromosome_depth = max(chromosome_depths) if chromosome_depths else None

    for contig in contigs:
        if not contig.circular:
            contig.flags.append("linear")
        if contig.length < min_contig_len:
            contig.flags.append("short")
        if (
            chromosome_depth
            and contig.depth is not None
            and contig.replicon_type != CHROMOSOME
            and contig.depth / chromosome_depth < min_depth_ratio
        ):
            contig.flags.append("low_depth")
        if contig.original_name in unrotated:
            contig.flags.append("unrotated")
        if contig.original_name in overlapping:
            contig.flags.append("end_overlap")


def rename(contigs: list[Contig], sample: str) -> None:
    """Name the replicons ``<id>_chromosome`` and ``<id>_plasmid_1..n``, longest first."""
    counters = {CHROMOSOME: 0, PLASMID: 0, UNPLACED: 0}
    totals = {kind: sum(1 for c in contigs if c.replicon_type == kind) for kind in counters}

    for contig in sorted(contigs, key=lambda c: c.length, reverse=True):
        kind = contig.replicon_type
        counters[kind] += 1
        if kind == CHROMOSOME and totals[CHROMOSOME] == 1:
            contig.name = f"{sample}_chromosome"
        else:
            contig.name = f"{sample}_{kind}_{counters[kind]}"


def write_outputs(contigs: list[Contig], args: argparse.Namespace) -> None:
    """Write the final FASTA, the removed contigs, and the per-contig table."""
    kept = [c for c in contigs if not (args.drop_flagged and c.flags)]
    removed = [c for c in contigs if args.drop_flagged and c.flags]

    order = sorted(
        kept,
        key=lambda c: ({CHROMOSOME: 0, PLASMID: 1, UNPLACED: 2}[c.replicon_type], -c.length),
    )
    args.output.write_text("".join(f">{c.header()}\n{c.sequence}\n" for c in order))
    args.removed.write_text("".join(f">{c.header()}\n{c.sequence}\n" for c in removed))

    with args.table.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["sample", "contig", "original_contig", "replicon_type", "length", "depth",
             "circular", "flags", "in_final_assembly"]
        )
        for contig in sorted(contigs, key=lambda c: c.length, reverse=True):
            writer.writerow(
                [
                    args.sample,
                    contig.name,
                    contig.original_name,
                    contig.replicon_type,
                    contig.length,
                    "" if contig.depth is None else contig.depth,
                    "true" if contig.circular else "false",
                    ",".join(contig.flags),
                    "false" if contig in removed else "true",
                ]
            )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--rotated", type=Path, required=True, help="dnaapler's reoriented FASTA")
    parser.add_argument("--unrotated", type=Path, help="dnaapler's failed_to_reorient FASTA")
    parser.add_argument("--end-overlaps", type=Path, help="PAF of contig ends against each other")
    parser.add_argument("--min-contig-len", type=int, default=1000)
    parser.add_argument("--chromosome-min-len", type=int, default=1_000_000)
    parser.add_argument("--min-depth-ratio", type=float, default=0.1)
    parser.add_argument("--drop-flagged", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--removed", type=Path, required=True)
    parser.add_argument("--table", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    unrotated_records = load_fasta(args.unrotated)
    contigs = [
        Contig(header, sequence)
        for header, sequence in load_fasta(args.rotated) + unrotated_records
    ]
    if not contigs:
        print(f"no contigs found for {args.sample}", file=sys.stderr)

    unrotated = {header.split()[0] for header, _ in unrotated_records}
    classify(contigs, args.chromosome_min_len)
    flag(contigs, args.min_contig_len, args.min_depth_ratio, unrotated, self_overlapping(args.end_overlaps))
    rename(contigs, args.sample)
    write_outputs(contigs, args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
