"""Tests for bin/circularise.py."""

import csv
import importlib.util
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "bin" / "circularise.py"
_spec = importlib.util.spec_from_file_location("circularise", BIN)
circularise = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(circularise)

WINDOW = 1000


def paf_line(
    contig: str,
    qstart: int,
    qend: int,
    tstart: int,
    tend: int,
    strand: str = "+",
    identity: float = 1.0,
    window: int = WINDOW,
) -> str:
    """One CONTIG_ENDS record: the contig's start window aligned to its end window."""
    block = qend - qstart
    return "\t".join(
        str(field)
        for field in [
            f"{contig}_start", window, qstart, qend, strand,
            f"{contig}_end", window, tstart, tend,
            int(block * identity), block, 60,
        ]
    )


def run(tmp_path: Path, sequence: str, *records: str) -> tuple[str, dict[str, str]]:
    """Run the script over one contig and read back the FASTA and its table row."""
    assembly = tmp_path / "in.fasta"
    assembly.write_text(f">c1 length={len(sequence)} depth=30.0 circular=true\n{sequence}\n")
    paf = tmp_path / "ends.paf"
    paf.write_text("".join(f"{record}\n" for record in records))
    output, table = tmp_path / "out.fasta", tmp_path / "circularity.tsv"
    assert circularise.main([
        "--assembly", str(assembly), "--overlaps", str(paf),
        "--output", str(output), "--table", str(table),
    ]) == 0
    rows = list(csv.DictReader(table.read_text().splitlines(), delimiter="\t"))
    return output.read_text(), rows[0]


def sequence(length: int, overlap: int = 0) -> str:
    """A contig of `length` bases whose last `overlap` bases repeat its first."""
    body = "".join("ACGT"[i % 4] for i in range(length))
    return body[: length - overlap] + body[:overlap] if overlap else body


def test_a_clean_wrap_is_cut_once_and_called_circular(tmp_path: Path) -> None:
    seq = sequence(5000, overlap=300)
    fasta, row = run(tmp_path, seq, paf_line("c1", 0, 300, WINDOW - 300, WINDOW))
    assert fasta.splitlines()[1] == seq[:-300]
    assert row["overlap"] == "300"
    assert (row["length_before"], row["length_after"]) == ("5000", "4700")
    assert row["circular"] == "true"
    assert row["flag"] == ""
    assert "length=4700" in fasta.splitlines()[0]
    assert "circular=true" in fasta.splitlines()[0]


def test_identity_below_the_threshold_is_left_alone(tmp_path: Path) -> None:
    seq = sequence(5000)
    fasta, row = run(
        tmp_path, seq, paf_line("c1", 0, 300, WINDOW - 300, WINDOW, identity=0.80)
    )
    assert fasta.splitlines()[1] == seq
    assert row["circular"] == "false"
    assert row["flag"] == ""


def test_a_match_starting_inside_the_contig_is_an_internal_repeat(tmp_path: Path) -> None:
    seq = sequence(5000)
    fasta, row = run(tmp_path, seq, paf_line("c1", 400, 700, WINDOW - 300, WINDOW))
    assert fasta.splitlines()[1] == seq
    assert row["flag"] == "internal_repeat"
    assert row["circular"] == "false"


def test_a_match_stopping_short_of_the_end_is_an_internal_repeat(tmp_path: Path) -> None:
    seq = sequence(5000)
    fasta, row = run(tmp_path, seq, paf_line("c1", 0, 300, 200, 500))
    assert fasta.splitlines()[1] == seq
    assert row["flag"] == "internal_repeat"


def test_a_reverse_hit_is_an_inverted_repeat_not_a_wrap(tmp_path: Path) -> None:
    seq = sequence(5000)
    fasta, row = run(
        tmp_path, seq, paf_line("c1", 0, 300, WINDOW - 300, WINDOW, strand="-")
    )
    assert fasta.splitlines()[1] == seq
    assert row["flag"] == "internal_repeat"


def test_a_window_filling_match_is_flagged_and_not_trimmed(tmp_path: Path) -> None:
    seq = sequence(5000)
    fasta, row = run(tmp_path, seq, paf_line("c1", 0, WINDOW, 0, WINDOW))
    assert fasta.splitlines()[1] == seq
    assert row["flag"] == "overlap_exceeds_window"
    assert row["overlap"] == f"{WINDOW}+"
    assert row["circular"] == "false"


def test_a_short_contig_with_halved_windows_gives_no_spurious_wrap(tmp_path: Path) -> None:
    # CONTIG_ENDS halves the window on a contig shorter than 2 * window, so the two
    # windows are disjoint and a short contig no longer self-aligns for free.
    seq = sequence(900)
    fasta, row = run(tmp_path, seq, "")  # a blank PAF line, as an empty PAF leaves
    assert fasta.splitlines()[1] == seq
    assert row["circular"] == "false"
    assert row["flag"] == ""


def test_a_contig_with_no_alignment_is_reported_unchanged(tmp_path: Path) -> None:
    seq = sequence(5000)
    fasta, row = run(tmp_path, seq)
    assert fasta.splitlines()[1] == seq
    assert row["identity"] == ""
    assert row["circular"] == "false"
