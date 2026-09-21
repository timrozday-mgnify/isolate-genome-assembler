#!/usr/bin/env python3
"""Check the final assembly against Plassembler's own view of the sample's plasmids.

Plassembler is run a second time on the full read set, independently of the consensus.
Each plasmid it reports is matched against the finished contigs by skani containment: a
plasmid that is not there is **missing**, which is a warn the report shows so a person can
decide whether to add it back.

HiFi library prep removes short fragments, so plasmids under ~10 kb are under-represented
in the reads to begin with. That caveat belongs in the report next to this table.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def read_tsv(path: Path | None) -> list[dict[str, str]]:
    """Read a TSV into a list of rows, tolerating a missing or empty file."""
    if path is None or not path.exists() or path.stat().st_size == 0:
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value: object, default: float = 0.0) -> float:
    """Coerce a field to float, falling back to ``default`` for blanks and junk."""
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return default


def first_word(value: object) -> str:
    """The first whitespace-separated token, or empty for a blank or missing field."""
    return str(value or "").split()[0] if str(value or "").strip() else ""


def best_matches(
    rows: list[dict[str, str]], min_identity: float, min_coverage: float
) -> dict[str, dict[str, str]]:
    """The best final contig for each Plassembler plasmid, keyed by the plasmid's name.

    skani's query is the Plassembler plasmid and the reference is the final assembly, so
    `Align_fraction_query` is how much of the plasmid the assembly actually contains.
    """
    matches: dict[str, dict[str, str]] = {}
    for row in rows:
        plasmid = first_word(row.get("Query_name"))
        if not plasmid:
            continue
        identity = number(row.get("ANI"))
        coverage = number(row.get("Align_fraction_query"))
        if identity < min_identity or coverage < min_coverage:
            continue
        best = matches.get(plasmid)
        if best is None or coverage > number(best.get("Align_fraction_query")):
            matches[plasmid] = row
    return matches


def audit(args: argparse.Namespace) -> list[dict[str, object]]:
    """One row per plasmid Plassembler reported, recovered or missing."""
    matches = best_matches(read_tsv(args.skani), args.min_identity, args.min_coverage)

    rows: list[dict[str, object]] = []
    for plasmid in read_tsv(args.plassembler_summary):
        # Plassembler names its plasmids by number in a `contig` column.
        name = plasmid.get("contig") or plasmid.get("Contig") or ""
        match = matches.get(name)
        rows.append(
            {
                "sample": args.sample,
                "plassembler_contig": name,
                "length": plasmid.get("length") or plasmid.get("Length") or "",
                "copy_number": plasmid.get("copy_number_short")
                or plasmid.get("copy_number_long")
                or "",
                "plsdb_hit": plasmid.get("PLSDB_hit") or plasmid.get("plsdb_hit") or "",
                "status": "recovered" if match else "missing",
                "matched_contig": first_word((match or {}).get("Ref_name")),
                "identity": (match or {}).get("ANI", ""),
                "coverage": (match or {}).get("Align_fraction_query", ""),
            }
        )
    return rows


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument(
        "--plassembler-summary", type=Path, help="plassembler_summary.tsv"
    )
    parser.add_argument(
        "--skani", type=Path, help="skani dist of the plasmids vs the assembly"
    )
    parser.add_argument(
        "--min-identity", type=float, default=95.0, help="Minimum skani ANI"
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=90.0,
        help="Minimum percent of the plasmid the assembly must contain",
    )
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    rows = audit(args)
    fields = [
        "sample",
        "plassembler_contig",
        "length",
        "copy_number",
        "plsdb_hit",
        "status",
        "matched_contig",
        "identity",
        "coverage",
    ]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    missing = sum(1 for row in rows if row["status"] == "missing")
    if missing:
        print(
            f"{missing} Plassembler plasmid(s) missing from the assembly",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
