"""Tests for bin/classify_replicons.py."""

import csv
import importlib.util
from pathlib import Path
from typing import NamedTuple

BIN = Path(__file__).resolve().parents[1] / "bin" / "classify_replicons.py"
_spec = importlib.util.spec_from_file_location("classify_replicons", BIN)
classify_replicons = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(classify_replicons)


class Result(NamedTuple):
    """Everything classify_replicons.py wrote for one sample."""

    final: str
    removed: str
    rows: list[dict[str, str]]


def run(tmp_path: Path, rotated: str, **extra: object) -> Result:
    """Run the script over one rotated FASTA and read back everything it wrote."""
    rotated_path = tmp_path / "rotated.fasta"
    rotated_path.write_text(rotated)

    output = tmp_path / "final.fasta"
    removed = tmp_path / "removed.fasta"
    table = tmp_path / "contigs.tsv"
    argv = [
        "--sample", "isolate01",
        "--rotated", str(rotated_path),
        "--output", str(output),
        "--removed", str(removed),
        "--table", str(table),
    ]
    for name, value in extra.items():
        flag = f"--{name.replace('_', '-')}"
        argv += [flag] if value is True else [flag, str(value)]

    classify_replicons.main(argv)
    with table.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return Result(output.read_text(), removed.read_text(), rows)


CHROMOSOME_AND_PLASMID = (
    ">contig_1 length=2000000 depth=40 circular=true\n" + "A" * 2_000_000 + "\n"
    ">contig_2 length=5000 depth=120 circular=true\n" + "C" * 5_000 + "\n"
)


def test_the_longest_contig_is_the_chromosome_and_circular_ones_are_plasmids(tmp_path: Path) -> None:
    result = run(tmp_path, CHROMOSOME_AND_PLASMID)
    names = [(row["contig"], row["replicon_type"]) for row in result.rows]
    assert names == [
        ("isolate01_chromosome", "chromosome"),
        ("isolate01_plasmid_1", "plasmid"),
    ]


def test_plasmids_are_numbered_by_length_longest_first(tmp_path: Path) -> None:
    fasta = (
        ">c1 length=2000000 circular=true\n" + "A" * 2_000_000 + "\n"
        ">small length=3000 circular=true\n" + "C" * 3_000 + "\n"
        ">big length=9000 circular=true\n" + "G" * 9_000 + "\n"
    )
    rows = run(tmp_path, fasta).rows
    by_original = {row["original_contig"]: row["contig"] for row in rows}
    assert by_original["big"] == "isolate01_plasmid_1"
    assert by_original["small"] == "isolate01_plasmid_2"


def test_several_long_contigs_are_all_chromosomes_and_numbered(tmp_path: Path) -> None:
    fasta = (
        ">c1 length=2000000 circular=true\n" + "A" * 2_000_000 + "\n"
        ">c2 length=1500000 circular=true\n" + "C" * 1_500_000 + "\n"
    )
    rows = run(tmp_path, fasta).rows
    assert [row["contig"] for row in rows] == [
        "isolate01_chromosome_1",
        "isolate01_chromosome_2",
    ]


def test_short_linear_and_low_depth_contigs_are_flagged_but_kept(tmp_path: Path) -> None:
    fasta = (
        ">c1 length=2000000 depth=40 circular=true\n" + "A" * 2_000_000 + "\n"
        ">c2 length=500 depth=1 circular=false\n" + "C" * 500 + "\n"
    )
    result = run(tmp_path, fasta)
    flagged = result.rows[1]
    assert set(flagged["flags"].split(",")) == {"linear", "short", "low_depth"}
    assert flagged["in_final_assembly"] == "true"
    assert result.final.count(">") == 2
    assert result.removed == ""


def test_drop_flagged_moves_them_out_of_the_final_assembly(tmp_path: Path) -> None:
    fasta = (
        ">c1 length=2000000 depth=40 circular=true\n" + "A" * 2_000_000 + "\n"
        ">c2 length=500 depth=1 circular=false\n" + "C" * 500 + "\n"
    )
    result = run(tmp_path, fasta, drop_flagged=True)
    assert result.final.count(">") == 1
    assert result.removed.count(">") == 1
    assert result.rows[1]["in_final_assembly"] == "false"


def test_the_chromosome_is_never_flagged_for_its_own_depth(tmp_path: Path) -> None:
    fasta = ">c1 length=2000000 depth=40 circular=true\n" + "A" * 2_000_000 + "\n"
    assert run(tmp_path, fasta).rows[0]["flags"] == ""


def test_a_contig_that_failed_rotation_is_flagged(tmp_path: Path) -> None:
    unrotated = tmp_path / "unrotated.fasta"
    unrotated.write_text(">mystery length=4000 circular=true\n" + "T" * 4_000 + "\n")
    result = run(tmp_path, CHROMOSOME_AND_PLASMID, unrotated=unrotated)

    row = next(row for row in result.rows if row["original_contig"] == "mystery")
    assert "unrotated" in row["flags"]
    # It is still part of the assembly: an unrotated contig is a naming problem, not a
    # reason to lose sequence.
    assert row["in_final_assembly"] == "true"
    assert result.final.count(">") == 3


def test_a_contig_whose_ends_still_overlap_is_flagged(tmp_path: Path) -> None:
    paf = tmp_path / "ends.paf"
    paf.write_text(
        "contig_2_start\t10000\t0\t9000\t+\tcontig_2_end\t10000\t1000\t10000\t8900\t9000\t60\n"
    )
    rows = run(tmp_path, CHROMOSOME_AND_PLASMID, end_overlaps=paf).rows
    assert "end_overlap" in rows[1]["flags"]
    assert rows[0]["flags"] == ""


def test_an_alignment_between_two_different_contigs_is_not_a_self_overlap(tmp_path: Path) -> None:
    paf = tmp_path / "ends.paf"
    paf.write_text(
        "contig_1_start\t10000\t0\t9000\t+\tcontig_2_end\t10000\t1000\t10000\t8900\t9000\t60\n"
    )
    rows = run(tmp_path, CHROMOSOME_AND_PLASMID, end_overlaps=paf).rows
    assert all(row["flags"] == "" for row in rows)


def test_an_empty_assembly_writes_headers_and_no_contigs(tmp_path: Path) -> None:
    result = run(tmp_path, "")
    assert result.rows == []
    assert result.final == ""
