#!/usr/bin/env python3
"""Collapse a sample's read QC measurements into one tidy `read_qc.tsv` row.

Only measurements are written. Thresholds live in ``assets/qc_thresholds.yml`` and are
applied by ``qc_gates.py``, so this script never decides pass/warn/fail.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

SIZE_SUFFIXES = {"k": 1_000, "m": 1_000_000, "g": 1_000_000_000}

COLUMNS = [
    "sample",
    "reads_in",
    "reads_kept",
    "reads_below_min_qv",
    "bases",
    "read_n50",
    "read_length_mean",
    "read_length_min",
    "read_length_max",
    "q20_fraction",
    "q30_fraction",
    "gc_mean",
    "gc_sd",
    "duplicate_id_fraction",
    "duplicate_sequence_fraction",
    "adapter_fraction",
    "genome_size_declared",
    "genome_size_autocycler",
    "genome_size_genomescope",
    "genome_size_used",
    "genome_size_source",
    "genome_size_disagreement",
    "heterozygosity",
    "depth",
]


def parse_size(value: str | None) -> float | None:
    """Parse a samplesheet genome size such as ``5.2m`` into bases."""
    if not value:
        return None
    match = re.fullmatch(r"(\d+(?:\.\d+)?)\s*([kmg]?)", value.strip(), re.IGNORECASE)
    if not match:
        raise ValueError(f"cannot parse genome size {value!r}")
    return float(match.group(1)) * SIZE_SUFFIXES.get(match.group(2).lower(), 1)


def read_tsv_row(path: Path | None) -> dict[str, str]:
    """Return the first data row of a single-row TSV, or an empty mapping."""
    if path is None or not path.exists():
        return {}
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return rows[0] if rows else {}


def read_seqkit_stats(path: Path | None) -> dict[str, str]:
    """Normalise `seqkit stats --tabular --all` headers into snake_case keys."""
    row = read_tsv_row(path)
    return {
        key.strip().lower().replace("(%)", "").strip().replace(" ", "_"): value
        for key, value in row.items()
        if key
    }


def gc_moments(path: Path | None) -> tuple[float | None, float | None]:
    """Mean and standard deviation of per-read GC% from the histogram."""
    if path is None or not path.exists():
        return None, None
    total = 0
    weighted = 0.0
    weighted_squares = 0.0
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gc = float(row["gc_percent_bin"])
            count = int(row["read_count"])
            total += count
            weighted += gc * count
            weighted_squares += gc * gc * count
    if total == 0:
        return None, None
    mean = weighted / total
    variance = max(weighted_squares / total - mean * mean, 0.0)
    return mean, variance**0.5


def parse_genomescope(path: Path | None) -> tuple[float | None, float | None]:
    """Genome haploid length and heterozygosity from a GenomeScope2 summary."""
    if path is None or not path.exists():
        return None, None
    text = path.read_text()
    length = re.search(r"Genome Haploid Length\s+([\d,]+)\s*bp", text)
    heterozygosity = re.search(r"Heterozygosity\s+([\d.]+)%", text)
    return (
        float(length.group(1).replace(",", "")) if length else None,
        float(heterozygosity.group(1)) if heterozygosity else None,
    )


def first_number(path: Path | None) -> float | None:
    """First integer in a file, which is how Autocycler reports a genome size."""
    if path is None or not path.exists():
        return None
    match = re.search(r"\d+", path.read_text())
    return float(match.group()) if match else None


def ratio(numerator: object, denominator: object) -> float | None:
    """Safe division that tolerates missing or non-numeric inputs."""
    try:
        top, bottom = float(numerator), float(denominator)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return top / bottom if bottom else None


def summarise(args: argparse.Namespace) -> dict[str, object]:
    """Build the single output row from the per-tool measurement files."""
    normalisation = read_tsv_row(args.normalisation)
    stats = read_seqkit_stats(args.seqkit_stats)
    duplicates = read_tsv_row(args.duplicates)
    adapters = read_tsv_row(args.adapters)
    gc_mean, gc_sd = gc_moments(args.gc_hist)
    genomescope_size, heterozygosity = parse_genomescope(args.genomescope_summary)
    autocycler_size = first_number(args.genome_size_autocycler)
    declared_size = parse_size(args.declared_genome_size)

    used, source = next(
        (
            (size, name)
            for size, name in (
                (declared_size, "samplesheet"),
                (autocycler_size, "autocycler"),
                (genomescope_size, "genomescope2"),
            )
            if size
        ),
        (None, "none"),
    )

    disagreement = None
    if autocycler_size and genomescope_size:
        disagreement = abs(autocycler_size - genomescope_size) / max(
            autocycler_size, genomescope_size
        )

    bases = normalisation.get("bases_kept") or stats.get("sum_len")

    return {
        "sample": args.sample,
        "reads_in": normalisation.get("reads_in"),
        "reads_kept": normalisation.get("reads_kept"),
        "reads_below_min_qv": normalisation.get("reads_below_min_qv"),
        "bases": bases,
        "read_n50": stats.get("n50"),
        "read_length_mean": stats.get("avg_len"),
        "read_length_min": stats.get("min_len"),
        "read_length_max": stats.get("max_len"),
        "q20_fraction": ratio(stats.get("q20"), 100),
        "q30_fraction": ratio(stats.get("q30"), 100),
        "gc_mean": gc_mean,
        "gc_sd": gc_sd,
        "duplicate_id_fraction": ratio(
            duplicates.get("duplicate_ids"), duplicates.get("reads")
        ),
        "duplicate_sequence_fraction": ratio(
            duplicates.get("duplicate_sequences"), duplicates.get("reads")
        ),
        "adapter_fraction": adapters.get("adapter_fraction"),
        "genome_size_declared": declared_size,
        "genome_size_autocycler": autocycler_size,
        "genome_size_genomescope": genomescope_size,
        "genome_size_used": used,
        "genome_size_source": source,
        "genome_size_disagreement": disagreement,
        "heterozygosity": heterozygosity,
        "depth": ratio(bases, used),
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--normalisation", type=Path)
    parser.add_argument("--seqkit-stats", type=Path)
    parser.add_argument("--gc-hist", type=Path)
    parser.add_argument("--duplicates", type=Path)
    parser.add_argument("--adapters", type=Path)
    parser.add_argument("--genome-size-autocycler", type=Path)
    parser.add_argument("--genomescope-summary", type=Path)
    parser.add_argument("--declared-genome-size")
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    row = summarise(args)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerow(
            {key: "" if value is None else value for key, value in row.items()}
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
