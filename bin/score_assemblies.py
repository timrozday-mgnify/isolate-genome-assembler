#!/usr/bin/env python3
"""Rank one sample's candidate assemblies against each other, reference-free.

Every signal comes from a check the pipeline already runs, so scoring costs no new tool.
There is no weighted score: hard filters remove an assembly that cannot be right, and the
survivors are compared rung by rung, first difference winning. A weighted sum would hide
why one assembly beat another; this way the whole comparison is in the table it writes.

The order of the rungs is the argument: structural correctness, then structural
completeness, then base accuracy. An assembly with one misjoin and QV 60 is worse than a
slightly noisier one that is correctly joined, because polishing fixes the second and
nothing downstream fixes the first.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

# Plassembler assembles plasmids only, so as a whole-genome candidate it is all filter and
# no information. It runs on the full read set in FINISHING, where its answer is used.
EXCLUDED = {"plassembler": "plasmid-only assembler, audited in FINISHING instead"}

FIELDS = [
    "sample",
    "assembler",
    "contigs",
    "circular_contigs",
    "n50",
    "total_length",
    "size_ratio",
    "merqury_qv",
    "merqury_completeness",
    "unmapped_read_fraction",
    "clipping_confirmed",
    "filtered",
    "rank",
    "selected",
]


def number(value: object) -> float | None:
    """A float, or None when the field is missing or not a number."""
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


def count(value: object) -> int | None:
    """A whole number, for the columns that are counts rather than measurements."""
    number_value = number(value)
    return None if number_value is None else int(number_value)


def rows(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    """A headed TSV as dictionaries, or nothing when it is missing or empty."""
    if not path.exists() or not path.read_text().strip():
        return []
    return list(csv.DictReader(path.read_text().splitlines(), delimiter=delimiter))


def by_assembler(directory: Path, sample: str) -> dict[str, Path]:
    """Files named ``<sample>.<assembler>.<whatever>``, keyed by assembler."""
    found: dict[str, Path] = {}
    if not directory.is_dir():
        return found
    for path in sorted(directory.rglob("*")):
        if path.is_file() and path.name.startswith(f"{sample}."):
            assembler = path.name[len(sample) + 1 :].split(".")[0]
            found.setdefault(assembler, path)
    return found


def headless(path: Path | None, column: int) -> float | None:
    """One column of the first line of a headerless table, as a number."""
    if path is None or not path.exists():
        return None
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) > column:
            return number(fields[column])
    return None


def measure(
    assembler: str,
    sample: str,
    stats: Path | None,
    circularity: Path | None,
    qv: Path | None,
    completeness: Path | None,
    mapping: Path | None,
    clipping: Path | None,
    genome_size: float | None,
) -> dict[str, object]:
    """Everything the comparison reads about one candidate, as one table row."""
    stat = (rows(stats) or [{}])[0]
    total_length = number(stat.get("sum_len"))
    circular = [row for row in rows(circularity) if row.get("circular") == "true"]
    return {
        "sample": sample,
        "assembler": assembler,
        "contigs": count(stat.get("num_seqs")),
        "circular_contigs": len(circular),
        "n50": count(stat.get("N50")),
        "total_length": count(total_length),
        "size_ratio": (
            round(total_length / genome_size, 4)
            if total_length and genome_size
            else None
        ),
        # Merqury's <prefix>.qv is assembly, asm-only k-mers, total k-mers, QV, error;
        # its completeness.stats is assembly, set, found, total, percent.
        "merqury_qv": headless(qv, 3),
        "merqury_completeness": headless(completeness, 4),
        "unmapped_read_fraction": number(
            (rows(mapping) or [{}])[0].get("unmapped_read_fraction")
        ),
        "clipping_confirmed": sum(
            row.get("verdict") == "confirmed" for row in rows(clipping)
        ),
        "filtered": "",
        "rank": "",
        "selected": "",
    }


def filter_reason(row: dict[str, object], size_tolerance: float, min_qv: float) -> str:
    """Why this candidate cannot be selected, or an empty string when it can."""
    if row["assembler"] in EXCLUDED:
        return EXCLUDED[row["assembler"]]
    if not row["total_length"]:
        return "empty assembly"
    ratio = row["size_ratio"]
    if ratio is not None and abs(ratio - 1) > size_tolerance:
        return f"size_ratio {ratio:.3f} outside 1 +/- {size_tolerance}"
    qv = row["merqury_qv"]
    if qv is not None and qv < min_qv:
        return f"merqury_qv {qv:.1f} below {min_qv}"
    return ""


def order(row: dict[str, object]) -> tuple:
    """The ordered comparison, as a sort key: first difference wins.

    Missing measurements sort last within their rung rather than winning by absence.
    """
    return (
        row["clipping_confirmed"] if row["clipping_confirmed"] is not None else 1e9,
        -(row["circular_contigs"] or 0),
        row["contigs"] if row["contigs"] is not None else 1e9,
        -(row["merqury_completeness"] or 0),
        -(row["merqury_qv"] or 0),
        row["assembler"],
    )


def rank(
    candidates: list[dict[str, object]], size_tolerance: float, min_qv: float
) -> list[dict[str, object]]:
    """Filter, rank the survivors, and mark the winner.

    A filtered candidate is ranked last rather than dropped: it is still evidence, and a
    sample where everything is filtered must still deliver something.
    """
    for row in candidates:
        row["filtered"] = filter_reason(row, size_tolerance, min_qv)
    ranked = sorted(candidates, key=lambda row: (bool(row["filtered"]), order(row)))
    for position, row in enumerate(ranked, start=1):
        row["rank"] = position
        row["selected"] = "true" if position == 1 else "false"
    return ranked


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True, help="Sample id")
    for name, help_text in [
        ("stats", "SEQKIT_STATS tables"),
        ("circularity", "CIRCULARISE tables"),
        ("qv", "Merqury <prefix>.qv files"),
        ("completeness", "Merqury completeness.stats files"),
        ("mapping", "MAP_READS tables"),
        ("clipping", "CLIPPING_PILEUPS tables"),
    ]:
        parser.add_argument(
            f"--{name}",
            type=Path,
            default=Path(name),
            help=f"Directory of per-candidate {help_text}",
        )
    parser.add_argument(
        "--genome-size", type=float, help="Genome size the size ratio is taken against"
    )
    parser.add_argument("--size-tolerance", type=float, default=0.1)
    parser.add_argument("--min-qv", type=float, default=40.0)
    parser.add_argument("--output", type=Path, required=True, help="TSV to write")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    found = {
        name: by_assembler(getattr(args, name), args.sample)
        for name in [
            "stats",
            "circularity",
            "qv",
            "completeness",
            "mapping",
            "clipping",
        ]
    }
    assemblers = sorted({a for files in found.values() for a in files})
    candidates = [
        measure(
            assembler,
            args.sample,
            *[
                found[name].get(assembler)
                for name in [
                    "stats",
                    "circularity",
                    "qv",
                    "completeness",
                    "mapping",
                    "clipping",
                ]
            ],
            genome_size=args.genome_size,
        )
        for assembler in assemblers
    ]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rank(candidates, args.size_tolerance, args.min_qv))
    if not candidates:
        print(f"no candidate assemblies found for {args.sample}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
