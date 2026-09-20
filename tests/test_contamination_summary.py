"""Tests for bin/contamination_summary.py."""

import importlib.util
import json
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[1] / "bin" / "contamination_summary.py"
_spec = importlib.util.spec_from_file_location("contamination_summary", BIN)
contamination_summary = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(contamination_summary)

ECOLI = "d__Bacteria|p__Pseudomonadota|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Escherichia|s__Escherichia coli"
KLEBSIELLA = "d__Bacteria|p__Pseudomonadota|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Klebsiella|s__Klebsiella pneumoniae"


def write_tax(tmp_path: Path, rows: list[tuple[str, float, float]]) -> Path:
    path = tmp_path / "tax.tsv"
    lines = ["clade_name\trelative_abundance\tsequence_abundance"]
    lines += [f"{clade}\t{relative}\t{sequence}" for clade, relative, sequence in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


def summarise(tmp_path: Path, **overrides) -> dict:
    defaults = {
        "sample": "isolate01",
        "sylph_profile": None,
        "sylph_tax": None,
        "human_query": None,
        "human_fraction": None,
        "expected_taxon": None,
    }
    namespace = type("Args", (), {**defaults, **overrides})
    return contamination_summary.summarise(namespace)


def test_species_are_ranked_and_the_unknown_fraction_is_the_remainder(
    tmp_path: Path,
) -> None:
    tax = write_tax(tmp_path, [(KLEBSIELLA, 4.0, 3.0), (ECOLI, 92.0, 95.0)])
    summary = summarise(tmp_path, sylph_tax=tax)
    assert summary["dominant_species"] == "s__Escherichia coli"
    assert summary["dominant_abundance"] == pytest.approx(0.95)
    assert summary["secondary_species"] == "s__Klebsiella pneumoniae"
    assert summary["secondary_abundance"] == pytest.approx(0.03)
    assert summary["unknown_fraction"] == pytest.approx(0.02)


def test_non_species_rows_are_ignored(tmp_path: Path) -> None:
    genus_only = ECOLI.rsplit("|", 1)[0]
    tax = write_tax(tmp_path, [(genus_only, 100.0, 100.0)])
    summary = summarise(tmp_path, sylph_tax=tax)
    assert summary["species"] == []
    assert summary["unknown_fraction"] == pytest.approx(1.0)


@pytest.mark.parametrize(
    ("expected", "matches"),
    [
        ("g__Escherichia", True),
        ("s__Escherichia coli", True),
        ("g__Klebsiella", False),
        (None, None),
    ],
)
def test_expected_taxon_is_matched_at_any_rank(
    tmp_path: Path, expected, matches
) -> None:
    tax = write_tax(tmp_path, [(ECOLI, 100.0, 100.0)])
    assert (
        summarise(tmp_path, sylph_tax=tax, expected_taxon=expected)[
            "expected_taxon_matches"
        ]
        is matches
    )


def test_missing_inputs_give_an_empty_summary(tmp_path: Path) -> None:
    summary = summarise(tmp_path, sylph_tax=tmp_path / "absent.tsv")
    assert summary["dominant_species"] == ""
    assert summary["human_fraction"] == 0.0


@pytest.mark.parametrize(
    ("mode", "hits", "fraction", "expected"),
    [
        ("auto", 0, 0.0, False),
        ("auto", 1, 0.0, True),
        ("auto", 0, 0.001, True),
        ("auto", 0, 0.0009, False),
        ("true", 0, 0.0, True),
        ("false", 5, 0.5, False),
    ],
)
def test_human_removal_decision(mode, hits, fraction, expected) -> None:
    summary = {"human_query_hits": hits, "human_fraction": fraction}
    assert (
        contamination_summary.human_removal_decision(summary, mode, 0.001) is expected
    )


def test_main_writes_the_json_and_the_decision_file(tmp_path: Path) -> None:
    tax = write_tax(tmp_path, [(ECOLI, 100.0, 100.0)])
    human_fraction = tmp_path / "human_fraction.tsv"
    human_fraction.write_text("reads\thuman_reads\thuman_fraction\n1000\t5\t0.005\n")
    output = tmp_path / "summary.json"
    decision = tmp_path / "decision.txt"

    assert (
        contamination_summary.main(
            [
                "--sample",
                "isolate01",
                "--sylph-tax",
                str(tax),
                "--human-fraction",
                str(human_fraction),
                "--expected-taxon",
                "g__Escherichia",
                "--output",
                str(output),
                "--remove-human-decision",
                str(decision),
            ]
        )
        == 0
    )

    summary = json.loads(output.read_text())
    assert summary["expected_taxon_matches"] is True
    assert summary["remove_human_applied"] is True
    assert decision.read_text().strip() == "true"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
