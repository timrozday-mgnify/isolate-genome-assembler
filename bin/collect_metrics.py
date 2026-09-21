#!/usr/bin/env python3
"""Gather every sample's measurements into one tidy ``run_summary/`` bundle.

The report reads nothing else, so it renders the same from a real run as from the test
fixture. The input is a manifest of ``sample<TAB>kind<TAB>path`` lines. A kind is the
name of the table the file joins: most files are TSVs that are simply concatenated, with
a ``sample`` column added when they lack one. The few that are not TSVs have a parser
below, and ``image-*`` kinds are copied into ``images/``.

Only measurements are written; the pass/warn/fail calls come from ``<id>.qc.json``.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import sys
from collections import defaultdict
from pathlib import Path

STATUSES = ["pass", "warn", "fail"]

PAF_COLUMNS = [
    "query",
    "query_length",
    "query_start",
    "query_end",
    "strand",
    "target",
    "target_length",
    "target_start",
    "target_end",
    "matches",
    "block_length",
    "mapq",
]


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def read_tsv(path: Path) -> list[dict[str, str]]:
    """A headed TSV; a leading ``#`` on the header (Flye's assembly_info) is dropped."""
    if not path.exists() or path.stat().st_size == 0:
        return []
    with open_text(path) as handle:
        lines = handle.read().splitlines()
    if lines and lines[0].startswith("#"):
        lines[0] = lines[0][1:]
    return list(csv.DictReader(lines, delimiter="\t"))


def headless(path: Path, columns: list[str]) -> list[dict[str, str]]:
    """A TSV without a header, as rows of the given columns."""
    if not path.exists():
        return []
    with open_text(path) as handle:
        return [
            dict(zip(columns, line.rstrip("\n").split("\t")))
            for line in handle
            if line.strip()
        ]


def key_values(path: Path, separator: str) -> list[dict[str, str]]:
    """``key<sep>value`` lines (Inspector, Bakta) as one row."""
    row = {}
    for line in path.read_text().splitlines():
        key, found, value = line.partition(separator)
        if found and key.strip() and value.strip():
            row[key.strip()] = value.strip()
    return [row] if row else []


def qc_rows(path: Path) -> list[dict[str, object]]:
    """``<id>.qc.json`` in long format, one row per check."""
    result = json.loads(path.read_text())
    return [
        {
            "check": name,
            "value": check.get("value"),
            "warn": check.get("threshold", {}).get("warn"),
            "fail": check.get("threshold", {}).get("fail"),
            "status": check.get("status"),
            "message": check.get("message"),
        }
        for name, check in result.get("checks", {}).items()
    ]


def contamination_rows(path: Path) -> list[dict[str, object]]:
    """The scalar fields of the contamination summary; species go to their own table."""
    summary = json.loads(path.read_text())
    return [{k: v for k, v in summary.items() if not isinstance(v, (list, dict))}]


def species_rows(path: Path) -> list[dict[str, object]]:
    return json.loads(path.read_text()).get("species", [])


# kind -> (output table, parser). Anything not listed is a headed TSV.
PARSERS = {
    "qc": ("qc_gates", qc_rows),
    "contamination": ("contamination", contamination_rows),
    "inspector": ("inspector", lambda p: key_values(p, "\t")),
    "bakta": ("bakta", lambda p: key_values(p, ":")),
    "merqury": (
        "merqury",
        lambda p: headless(p, ["assembly", "asm_only_kmers", "kmers", "qv", "error"]),
    ),
    "merqury-completeness": (
        "merqury_completeness",
        lambda p: headless(p, ["assembly", "set", "found", "total", "completeness"]),
    ),
    "coverage": (
        "coverage",
        lambda p: headless(p, ["contig", "start", "end", "depth"]),
    ),
    "dotplot": ("dotplot", lambda p: headless(p, PAF_COLUMNS)),
}
EXTRA_TABLES = {"contamination": ("contamination_species", species_rows)}


def first(rows: list[dict], column: str) -> object:
    return rows[0].get(column, "") if rows else ""


def sample_rows(tables: dict[str, list[dict]], samples: list[str]) -> list[dict]:
    """One status-board row per sample, from the tables already collected."""

    def of(table: str, sample: str) -> list[dict]:
        return [row for row in tables.get(table, []) if row["sample"] == sample]

    board = []
    for sample in samples:
        contigs = [
            row
            for row in of("contigs", sample)
            if row.get("in_final_assembly") != "false"
        ]
        chromosomes = [
            row for row in contigs if row.get("replicon_type") == "chromosome"
        ]
        plasmids = [row for row in contigs if row.get("replicon_type") == "plasmid"]
        gtdbtk = [
            row.get("classification", "")
            for row in of("gtdbtk", sample)
            if row.get("classification", "")
            and not row["classification"].startswith("Unclassified")
        ]
        species = (gtdbtk[0].split(";")[-1] if gtdbtk else "") or first(
            of("contamination", sample), "dominant_species"
        )
        # The same rule as qc_gates.py: the worst measured check; not_measured never counts.
        worst = [
            row["status"] for row in of("qc_gates", sample) if row["status"] in STATUSES
        ]
        board.append(
            {
                "sample": sample,
                "status": max(worst, key=STATUSES.index, default=""),
                "depth": first(of("read_qc", sample), "depth"),
                "species": species,
                "chromosomes": len(chromosomes),
                "plasmids": len(plasmids),
                "all_circular": bool(contigs)
                and all(row.get("circular") == "true" for row in contigs),
                "total_length": sum(
                    int(float(row.get("length") or 0)) for row in contigs
                ),
                "merqury_qv": first(of("merqury", sample), "qv"),
                "checkm2_completeness": first(of("checkm2", sample), "Completeness"),
                "checkm2_contamination": first(of("checkm2", sample), "Contamination"),
                "assembly_source": first(
                    of("assembly_source", sample), "assembly_source"
                ),
            }
        )
    return board


def write_tsv(path: Path, rows: list[dict]) -> None:
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        # Lower-case booleans, as every other table in the pipeline writes them.
        writer.writerows(
            {k: str(v).lower() if isinstance(v, bool) else v for k, v in row.items()}
            for row in rows
        )


def collect(manifest: Path, outdir: Path) -> dict[str, list[dict]]:
    tables: dict[str, list[dict]] = defaultdict(list)
    images = []
    samples = []
    (outdir / "images").mkdir(parents=True, exist_ok=True)

    for entry in headless(manifest, ["sample", "kind", "path"]):
        sample, kind, path = entry["sample"], entry["kind"], Path(entry["path"])
        if sample not in samples:
            samples.append(sample)
        if kind.startswith("image-"):
            target = (
                Path("images") / f"{sample}.{kind.removeprefix('image-')}{path.suffix}"
            )
            shutil.copyfile(path, outdir / target)
            images.append(
                {
                    "sample": sample,
                    "kind": kind.removeprefix("image-"),
                    "path": str(target),
                }
            )
            continue
        parsers = [PARSERS.get(kind, (kind.replace("-", "_"), read_tsv))]
        parsers += [EXTRA_TABLES[kind]] if kind in EXTRA_TABLES else []
        for table, parser in parsers:
            tables[table] += [{"sample": sample} | row for row in parser(path)]

    tables["images"] = images
    tables["samples"] = sample_rows(tables, sorted(samples))
    for table, rows in tables.items():
        write_tsv(outdir / f"{table}.tsv", rows)
    return tables


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest", type=Path, required=True, help="sample, kind, path TSV"
    )
    parser.add_argument(
        "--software-versions", type=Path, help="process, tool, version TSV"
    )
    parser.add_argument("--run-info", type=Path, help="JSON: pipeline, run and params")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)

    tables = collect(args.manifest, args.outdir)
    if args.software_versions:
        write_tsv(
            args.outdir / "software_versions.tsv",
            headless(args.software_versions, ["process", "tool", "version"]),
        )
    if args.run_info:
        shutil.copyfile(args.run_info, args.outdir / "run_info.json")
    print(f"collected {len(tables['samples'])} sample(s)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
