#!/usr/bin/env python3
"""Apply the pass/warn/fail thresholds to every measurement taken for one sample.

Each check becomes one entry of ``<id>.qc.json``: ``{value, threshold, status, message}``.
The sample's overall status is its worst entry. A measurement that is missing (a check
that did not run, or a tool that wrote nothing) is ``not_measured`` and never changes the
overall status, so a partial run still gets an honest verdict for what it did measure.

The thresholds are data, not code: they arrive as JSON, converted by the workflow from
``assets/qc_thresholds.yml`` or ``--qc_thresholds``.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

STATUS_ORDER = {"pass": 0, "warn": 1, "fail": 2}


def read_tsv(path: Path | None) -> list[dict[str, str]]:
    """Read a TSV into a list of rows, tolerating a missing or empty file."""
    if path is None or not path.exists() or path.stat().st_size == 0:
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_json(path: Path | None) -> dict:
    if path is None or not path.exists() or path.stat().st_size == 0:
        return {}
    return json.loads(path.read_text())


def number(value: object) -> float | None:
    """Coerce a field to float; None for blanks, junk and NaN."""
    try:
        result = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return None if result != result else result


def first(path: Path | None, column: str) -> float | None:
    """A numeric column of a one-row TSV."""
    rows = read_tsv(path)
    return number(rows[0].get(column)) if rows else None


def count_rows(path: Path | None, **where: str) -> int | None:
    """Rows of a TSV matching every ``column=value``; None when there is no file at all."""
    if path is None or not path.exists():
        return None
    return sum(
        all(row.get(key) == value for key, value in where.items())
        for row in read_tsv(path)
    )


# --- Readers, one per input. Each returns None when the measurement is unavailable. ---


def consensus_unresolved(path: Path | None) -> float | None:
    """Did the consensus resolve? Not the same as whether it was the one delivered.

    Under `--assembly_selection always_score` a resolved consensus can still lose to a
    better-scoring single assembly, which is not a failure of the consensus, so the
    verdict Autocycler wrote is read first and the source only stands in for it.
    """
    rows = read_tsv(path)
    if not rows:
        return None
    resolved = (rows[0].get("fully_resolved") or "").lower()
    if resolved in ("true", "false"):
        return 0.0 if resolved == "true" else 1.0
    return 0.0 if rows[0].get("assembly_source") == "autocycler" else 1.0


def chromosome_not_circular(path: Path | None) -> float | None:
    rows = [row for row in read_tsv(path) if row.get("replicon_type") == "chromosome"]
    if not rows:
        return None
    return float(sum(row.get("circular") != "true" for row in rows))


def plasmid_depth_ratio(path: Path | None) -> float | None:
    ratios = [
        number(row.get("depth_ratio"))
        for row in read_tsv(path)
        if row.get("replicon_type") == "plasmid"
    ]
    ratios = [ratio for ratio in ratios if ratio is not None]
    return min(ratios) if ratios else None


def circular_flye_contigs(path: Path | None) -> float | None:
    """Circular contigs in Flye's assembly_info.txt of the unmapped reads."""
    if path is None or not path.exists():
        return None
    count = 0
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if line.startswith("#") or len(fields) < 4:
            continue
        count += fields[3].strip() == "Y"
    return float(count)


def inspector_structural_errors(path: Path | None) -> float | None:
    """`Structural error` from Inspector's summary_statistics."""
    if path is None or not path.exists():
        return None
    match = re.search(r"^Structural error\t(\d+)", path.read_text(), re.MULTILINE)
    return float(match.group(1)) if match else None


def merqury_qv(path: Path | None) -> float | None:
    """Merqury's ``<prefix>.qv``: assembly, asm-only k-mers, total k-mers, QV, error."""
    if path is None or not path.exists():
        return None
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) >= 4:
            return number(fields[3])
    return None


def lineage_ranks(lineage: str) -> set[str]:
    """``d__Bacteria;...;s__Escherichia coli`` (or ``|``-separated) as a set of ranks."""
    return {
        part.strip()
        for part in re.split(r"[;|]", lineage)
        if re.match(r"^[dpcofgs]__.+", part.strip())
    }


def gtdbtk_classification(paths: list[Path]) -> str:
    """The classification of the one genome GTDB-Tk was given, bacterial or archaeal."""
    for path in paths:
        for row in read_tsv(path):
            classification = row.get("classification", "")
            if classification and not classification.startswith("Unclassified"):
                return classification
    return ""


def taxon_disagreement(
    contamination: dict, gtdbtk: str, expected: str | None
) -> float | None:
    """1 when the available taxon calls disagree, 0 when they agree.

    sylph and GTDB-Tk are compared at the species rank. ``expected_taxon`` may name any
    rank, so it must appear in each lineage. Fewer than two sources is not measurable.
    """
    species = contamination.get("species") or []
    sylph = lineage_ranks(species[0]["clade"] if species else "")
    classified = [ranks for ranks in (sylph, lineage_ranks(gtdbtk)) if ranks]
    if len(classified) + bool(expected) < 2:
        return None

    agree = True
    if len(classified) == 2:
        agree = {r for r in classified[0] if r.startswith("s__")} == {
            r for r in classified[1] if r.startswith("s__")
        }
    if expected:
        agree = agree and all(expected.strip() in ranks for ranks in classified)
    return 0.0 if agree else 1.0


def measure(args: argparse.Namespace) -> dict[str, float | None]:
    """Every measurement the thresholds can refer to, keyed by check name."""
    contamination = read_json(args.contamination)

    def contamination_value(key: str) -> float | None:
        return number(contamination.get(key)) if contamination else None

    # sylph places nothing at all for a species GTDB does not hold, so an empty profile
    # says the isolate is novel, not that it is contaminated.
    placed = bool(contamination.get("species"))

    return {
        "read_depth": first(args.read_qc, "depth"),
        "genome_size_disagreement": first(args.read_qc, "genome_size_disagreement"),
        "dominant_species_abundance": contamination_value("dominant_abundance")
        if placed
        else None,
        "secondary_species_abundance": contamination_value("secondary_abundance")
        if placed
        else None,
        "human_read_fraction": contamination_value("human_fraction"),
        "taxon_disagreement": taxon_disagreement(
            contamination, gtdbtk_classification(args.gtdbtk), args.expected_taxon
        ),
        "consensus_unresolved": consensus_unresolved(args.assembly_source),
        "chromosome_not_circular": chromosome_not_circular(args.contigs),
        "plasmids_missing": count_rows(args.plasmid_audit, status="missing"),
        "possible_missing_replicons": circular_flye_contigs(args.unmapped_assembly),
        "unmapped_read_fraction": first(args.mapping, "unmapped_read_fraction"),
        "coverage_regions": count_rows(args.coverage_regions),
        # Only the pile-ups a permissive re-alignment could not extend through.
        "clipping_pileups": count_rows(args.clipping, verdict="confirmed"),
        "inspector_structural_errors": inspector_structural_errors(args.inspector),
        "variants_high_af": first(args.variants, "high_af"),
        "merqury_qv": merqury_qv(args.merqury),
        "ideel_truncated_fraction": first(args.ideel, "truncated_fraction"),
        "checkm2_completeness": first(args.checkm2, "Completeness"),
        "checkm2_contamination": first(args.checkm2, "Contamination"),
        "busco_duplicated": first(args.busco, "Duplicated"),
        "plasmid_depth_ratio": plasmid_depth_ratio(args.contig_depth),
    }


def gate(name: str, value: float | None, rule: dict) -> dict[str, object]:
    """Compare one value with its rule. Fail is checked before warn."""
    direction = rule.get("direction", "above")
    if direction not in ("above", "below"):
        raise ValueError(f"{name}: direction must be 'above' or 'below'")
    threshold = {key: rule[key] for key in ("warn", "fail") if key in rule}
    entry: dict[str, object] = {"value": value, "threshold": threshold}

    if value is None:
        return entry | {"status": "not_measured", "message": f"{name} not measured"}

    def trips(limit: float) -> bool:
        return value < limit if direction == "below" else value >= limit

    for status in ("fail", "warn"):
        if status in threshold and trips(threshold[status]):
            comparison = "<" if direction == "below" else ">="
            message = f"{name} {value:g} {comparison} {threshold[status]:g}"
            return entry | {"status": status, "message": message}
    return entry | {"status": "pass", "message": ""}


def evaluate(measurements: dict, thresholds: dict) -> dict[str, object]:
    unknown = sorted(set(thresholds) - set(measurements))
    if unknown:
        raise SystemExit(f"thresholds name unknown check(s): {', '.join(unknown)}")

    checks = {
        name: gate(name, measurements[name], thresholds[name]) for name in thresholds
    }
    statuses = [str(c["status"]) for c in checks.values()]
    overall = max(
        (s for s in statuses if s in STATUS_ORDER),
        key=lambda status: STATUS_ORDER[status],
        default="pass",
    )
    return {"status": overall, "checks": checks}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--thresholds", type=Path, required=True, help="JSON")
    parser.add_argument("--expected-taxon")
    parser.add_argument("--read-qc", type=Path)
    parser.add_argument("--contamination", type=Path)
    parser.add_argument("--assembly-source", type=Path)
    parser.add_argument("--contigs", type=Path)
    parser.add_argument("--plasmid-audit", type=Path)
    parser.add_argument("--unmapped-assembly", type=Path, help="Flye assembly_info.txt")
    parser.add_argument("--mapping", type=Path)
    parser.add_argument("--coverage-regions", type=Path)
    parser.add_argument("--contig-depth", type=Path)
    parser.add_argument("--clipping", type=Path)
    parser.add_argument("--inspector", type=Path, help="summary_statistics")
    parser.add_argument("--variants", type=Path, help="variant summary TSV")
    parser.add_argument("--merqury", type=Path, help="<prefix>.qv")
    parser.add_argument("--ideel", type=Path, help="IDEEL summary TSV")
    parser.add_argument("--checkm2", type=Path)
    parser.add_argument("--busco", type=Path, help="BUSCO batch summary")
    parser.add_argument("--gtdbtk", type=Path, action="append", default=[])
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    result = evaluate(measure(args), json.loads(args.thresholds.read_text()))
    args.output.write_text(
        json.dumps({"sample": args.sample} | result, indent=2) + "\n"
    )
    print(f"{args.sample}: {result['status']}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
