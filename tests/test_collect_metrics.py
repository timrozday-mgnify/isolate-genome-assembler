"""Tests for bin/collect_metrics.py, through the fixture generator that drives it."""

import csv
import filecmp
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    "generate_run_summary", ROOT / "tests" / "data" / "generate_run_summary.py"
)
generator = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(generator)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_checked_in_fixture_is_current(tmp_path: Path) -> None:
    """The report's fixture is what the current scripts produce; regenerate if not."""
    generator.generate(tmp_path / "run_summary")
    comparison = filecmp.dircmp(tmp_path / "run_summary", generator.FIXTURE)
    assert not comparison.left_only and not comparison.right_only
    _, mismatch, errors = filecmp.cmpfiles(
        tmp_path / "run_summary",
        generator.FIXTURE,
        comparison.common_files,
        shallow=False,
    )
    assert not mismatch and not errors, (
        "run: python3 tests/data/generate_run_summary.py"
    )


def test_status_board_takes_the_worst_measured_check(tmp_path: Path) -> None:
    generator.generate(tmp_path)
    board = {row["sample"]: row for row in rows(tmp_path / "samples.tsv")}
    assert {s: row["status"] for s, row in board.items()} == {
        "iso_pass": "pass",
        "iso_warn": "warn",
        "iso_fail": "fail",
    }
    assert board["iso_pass"]["plasmids"] == "2"
    assert board["iso_pass"]["total_length"] == "264000"
    assert board["iso_pass"]["all_circular"] == "true"
    assert board["iso_fail"]["assembly_source"] == "fallback_flye"


def test_non_tsv_inputs_are_parsed(tmp_path: Path) -> None:
    generator.generate(tmp_path)
    assert rows(tmp_path / "merqury.tsv")[0]["qv"] == "62.1"
    assert rows(tmp_path / "inspector.tsv")[0]["Structural error"] == "0"
    assert rows(tmp_path / "bakta.tsv")[0]["pseudogenes"] == "4"
    assert rows(tmp_path / "coverage.tsv")[0].keys() >= {"contig", "start", "depth"}
    species = rows(tmp_path / "contamination_species.tsv")
    assert [r["species"] for r in species if r["sample"] == "iso_fail"] == [
        "s__Escherichia coli",
        "s__Klebsiella pneumoniae",
    ]
    # Headerless inputs gain a sample column; headed ones keep their own.
    assert rows(tmp_path / "checkm2.tsv")[0]["sample"] == "iso_pass"


def test_images_are_copied_and_listed(tmp_path: Path) -> None:
    generator.generate(tmp_path)
    images = rows(tmp_path / "images.tsv")
    assert {row["sample"] for row in images} == set(generator.SAMPLES)
    assert all((tmp_path / row["path"]).read_bytes() == generator.PNG for row in images)
