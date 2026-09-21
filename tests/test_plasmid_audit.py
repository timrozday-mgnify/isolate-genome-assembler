"""Tests for bin/plasmid_audit.py."""

import csv
import importlib.util
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "bin" / "plasmid_audit.py"
_spec = importlib.util.spec_from_file_location("plasmid_audit", BIN)
plasmid_audit = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(plasmid_audit)

SUMMARY = (
    "contig\tlength\tcopy_number_long\tPLSDB_hit\n"
    "1\t8000\t2.4\tNZ_CP012345\n"
    "2\t4000\t1.1\tNone\n"
)

SKANI_HEADER = "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"


def skani_row(query: str, reference: str, ani: float, coverage: float) -> str:
    return (
        f"final.fasta\tplasmids.fasta\t{ani}\t50.0\t{coverage}\t{reference}\t{query}\n"
    )


def run(tmp_path: Path, summary: str | None, skani: str) -> list[dict[str, str]]:
    output = tmp_path / "audit.tsv"
    argv = ["--sample", "isolate01", "--output", str(output)]

    if summary is not None:
        summary_path = tmp_path / "plassembler_summary.tsv"
        summary_path.write_text(summary)
        argv += ["--plassembler-summary", str(summary_path)]

    skani_path = tmp_path / "skani.tsv"
    skani_path.write_text(skani)
    argv += ["--skani", str(skani_path)]

    plasmid_audit.main(argv)
    with output.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_a_plasmid_found_in_the_assembly_is_recovered(tmp_path: Path) -> None:
    skani = SKANI_HEADER + skani_row("1", "isolate01_plasmid_1", 99.9, 99.0)
    rows = run(tmp_path, SUMMARY, skani)

    assert rows[0]["plassembler_contig"] == "1"
    assert rows[0]["status"] == "recovered"
    assert rows[0]["matched_contig"] == "isolate01_plasmid_1"
    assert rows[0]["plsdb_hit"] == "NZ_CP012345"
    assert rows[0]["copy_number"] == "2.4"


def test_a_plasmid_with_no_skani_hit_is_missing(tmp_path: Path) -> None:
    skani = SKANI_HEADER + skani_row("1", "isolate01_plasmid_1", 99.9, 99.0)
    rows = run(tmp_path, SUMMARY, skani)

    assert rows[1]["plassembler_contig"] == "2"
    assert rows[1]["status"] == "missing"
    assert rows[1]["matched_contig"] == ""


def test_a_hit_below_the_coverage_threshold_does_not_count_as_recovered(
    tmp_path: Path,
) -> None:
    # A plasmid the assembly only partly contains is exactly the case the audit exists to
    # catch, so a high-identity, low-coverage hit must not pass.
    skani = SKANI_HEADER + skani_row("1", "isolate01_chromosome", 99.9, 40.0)
    rows = run(tmp_path, SUMMARY, skani)
    assert rows[0]["status"] == "missing"


def test_a_hit_below_the_identity_threshold_does_not_count(tmp_path: Path) -> None:
    skani = SKANI_HEADER + skani_row("1", "isolate01_plasmid_1", 80.0, 99.0)
    rows = run(tmp_path, SUMMARY, skani)
    assert rows[0]["status"] == "missing"


def test_the_best_hit_wins_when_a_plasmid_matches_several_contigs(
    tmp_path: Path,
) -> None:
    skani = (
        SKANI_HEADER
        + skani_row("1", "isolate01_chromosome", 99.0, 91.0)
        + skani_row("1", "isolate01_plasmid_1", 99.9, 99.5)
    )
    rows = run(tmp_path, SUMMARY, skani)
    assert rows[0]["matched_contig"] == "isolate01_plasmid_1"


def test_a_sample_with_no_plasmids_gives_an_empty_table_not_an_error(
    tmp_path: Path,
) -> None:
    assert run(tmp_path, "contig\tlength\n", SKANI_HEADER) == []


def test_a_missing_plassembler_summary_is_tolerated(tmp_path: Path) -> None:
    assert run(tmp_path, None, SKANI_HEADER) == []
