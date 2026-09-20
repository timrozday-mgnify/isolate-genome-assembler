#!/usr/bin/env python3
"""Combine the sylph GTDB profile, the human query and the human read fraction.

Writes one JSON of measurements per sample plus the `--remove_human` decision, which is
the only threshold this script acts on: everything else is gated by ``qc_gates.py``.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

RANK_PREFIXES = ("d__", "p__", "c__", "o__", "f__", "g__", "s__")


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


def lineage_of(clade: str) -> dict[str, str]:
    """Split a sylph-tax clade string into its GTDB ranks."""
    return {part[:3]: part for part in clade.split("|") if part[:3] in RANK_PREFIXES}


def species_table(sylph_tax_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    """Species-level rows of a sylph-tax profile, most abundant first."""
    species = [
        {
            "clade": row["clade_name"],
            "species": lineage_of(row["clade_name"]).get("s__", ""),
            "relative_abundance": number(row.get("relative_abundance")),
            "sequence_abundance": number(row.get("sequence_abundance")),
        }
        for row in sylph_tax_rows
        if row.get("clade_name", "").split("|")[-1].startswith("s__")
    ]
    return sorted(species, key=lambda row: row["sequence_abundance"], reverse=True)


def taxon_matches(expected: str | None, clade: str) -> bool | None:
    """Whether the samplesheet's expected taxon appears in the dominant lineage."""
    if not expected:
        return None
    return expected.strip() in lineage_of(clade).values()


def summarise(args: argparse.Namespace) -> dict[str, object]:
    """Build the measurement bundle for one sample."""
    species = species_table(read_tsv(args.sylph_tax))
    total_abundance = sum(row["sequence_abundance"] for row in species)
    dominant = species[0] if species else None
    secondary = species[1] if len(species) > 1 else None

    human_rows = read_tsv(args.human_query)
    human_fraction = number(
        next(iter(read_tsv(args.human_fraction)), {}).get("human_fraction"),
        default=0.0,
    )

    # sylph reports abundances as percentages that sum to the classified fraction.
    return {
        "sample": args.sample,
        "species": species,
        "dominant_species": dominant["species"] if dominant else "",
        "dominant_abundance": (dominant["sequence_abundance"] / 100)
        if dominant
        else 0.0,
        "secondary_species": secondary["species"] if secondary else "",
        "secondary_abundance": (secondary["sequence_abundance"] / 100)
        if secondary
        else 0.0,
        "unknown_fraction": max(0.0, 1.0 - total_abundance / 100),
        "sylph_profile_rows": len(read_tsv(args.sylph_profile)),
        "human_query_hits": len(human_rows),
        "human_query_max_ani": max(
            (number(row.get("Adjusted_ANI")) for row in human_rows), default=0.0
        ),
        "human_fraction": human_fraction,
        "expected_taxon": args.expected_taxon or "",
        "expected_taxon_matches": taxon_matches(
            args.expected_taxon, dominant["clade"] if dominant else ""
        ),
    }


def human_removal_decision(
    summary: dict[str, object], mode: str, max_human_fraction: float
) -> bool:
    """`auto` removes human reads only when human is actually detected."""
    if mode == "true":
        return True
    if mode == "false":
        return False
    return (
        bool(summary["human_query_hits"])
        or summary["human_fraction"] >= max_human_fraction
    )  # type: ignore[operator]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--sylph-profile", type=Path)
    parser.add_argument("--sylph-tax", type=Path)
    parser.add_argument("--human-query", type=Path)
    parser.add_argument("--human-fraction", type=Path)
    parser.add_argument("--expected-taxon")
    parser.add_argument(
        "--remove-human", choices=["auto", "true", "false"], default="auto"
    )
    parser.add_argument("--max-human-fraction", type=float, default=0.001)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--remove-human-decision", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    summary = summarise(args)
    remove_human = human_removal_decision(
        summary, args.remove_human, args.max_human_fraction
    )
    summary["remove_human_mode"] = args.remove_human
    summary["remove_human_applied"] = remove_human
    args.output.write_text(json.dumps(summary, indent=2) + "\n")
    args.remove_human_decision.write_text("true\n" if remove_human else "false\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
