"""Tests for bin/score_assemblies.py."""

import csv
import importlib.util
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "bin" / "score_assemblies.py"
_spec = importlib.util.spec_from_file_location("score_assemblies", BIN)
score_assemblies = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(score_assemblies)

SAMPLE = "isolate01"
GENOME_SIZE = 5_000_000


def candidate(
    tmp_path: Path,
    assembler: str,
    contigs: int = 1,
    circular: int = 1,
    total_length: int = GENOME_SIZE,
    qv: float = 60.0,
    completeness: float = 99.0,
    clipping: int = 0,
) -> None:
    """Write one candidate's check outputs, as the per-candidate processes name them."""

    def write(directory: str, suffix: str, text: str) -> None:
        folder = tmp_path / directory
        folder.mkdir(exist_ok=True)
        (folder / f"{SAMPLE}.{assembler}{suffix}").write_text(text)

    write(
        "stats",
        ".tsv",
        "file\tformat\ttype\tnum_seqs\tsum_len\tN50\n"
        f"x\tFASTA\tDNA\t{contigs}\t{total_length}\t{total_length}\n",
    )
    write(
        "circularity",
        ".circularity.tsv",
        "contig\tlength_before\toverlap\tlength_after\tidentity\tcircular\tflag\n"
        + "".join(
            f"c{i}\t1\t0\t1\t\t{'true' if i <= circular else 'false'}\t\n"
            for i in range(1, contigs + 1)
        ),
    )
    write("qv", ".qv", f"asm\t1\t2\t{qv}\t0.0001\n")
    write("completeness", ".completeness.stats", f"asm\tall\t1\t2\t{completeness}\n")
    write(
        "mapping",
        ".mapping.tsv",
        "sample\treads\tunmapped_reads\tunmapped_read_fraction\n"
        f"{SAMPLE}\t100\t1\t0.01\n",
    )
    write(
        "clipping",
        ".clipping.tsv",
        "sample\tcontig\tverdict\n"
        + "".join(f"{SAMPLE}\tc1\tconfirmed\n" for _ in range(clipping)),
    )


def run(tmp_path: Path, **options: object) -> list[dict[str, str]]:
    """Score whatever candidates have been written, and read back the table."""
    output = tmp_path / "scores.tsv"
    argv = [
        "--sample",
        SAMPLE,
        "--genome-size",
        str(GENOME_SIZE),
        "--output",
        str(output),
    ]
    for name in ["stats", "circularity", "qv", "completeness", "mapping", "clipping"]:
        argv += [f"--{name}", str(tmp_path / name)]
    for name, value in options.items():
        argv += [f"--{name.replace('_', '-')}", str(value)]
    assert score_assemblies.main(argv) == 0
    return list(csv.DictReader(output.read_text().splitlines(), delimiter="\t"))


def winner(rows: list[dict[str, str]]) -> str:
    return next(row["assembler"] for row in rows if row["selected"] == "true")


def test_a_too_small_assembly_is_filtered(tmp_path: Path) -> None:
    candidate(tmp_path, "flye")
    candidate(tmp_path, "raven", total_length=int(GENOME_SIZE * 0.8))
    rows = {row["assembler"]: row for row in run(tmp_path)}
    assert "size_ratio" in rows["raven"]["filtered"]
    assert rows["flye"]["filtered"] == ""
    assert winner(rows.values()) == "flye"


def test_a_low_qv_assembly_is_filtered(tmp_path: Path) -> None:
    candidate(tmp_path, "flye")
    candidate(tmp_path, "raven", qv=20.0)
    rows = {row["assembler"]: row for row in run(tmp_path)}
    assert "merqury_qv" in rows["raven"]["filtered"]


def test_plassembler_is_never_a_candidate(tmp_path: Path) -> None:
    candidate(tmp_path, "plassembler")
    candidate(tmp_path, "flye")
    rows = {row["assembler"]: row for row in run(tmp_path)}
    assert "plasmid-only" in rows["plassembler"]["filtered"]
    assert winner(rows.values()) == "flye"


def test_a_filtered_candidate_still_wins_when_nothing_else_survives(
    tmp_path: Path,
) -> None:
    candidate(tmp_path, "flye", qv=20.0)
    assert winner(run(tmp_path)) == "flye"


def test_confirmed_clipping_decides_before_anything_else(tmp_path: Path) -> None:
    # hifiasm is better on every later rung, and still loses on the misjoin.
    candidate(tmp_path, "hifiasm", clipping=2, qv=70.0, completeness=99.9)
    candidate(tmp_path, "flye", clipping=0, qv=50.0, completeness=98.0)
    assert winner(run(tmp_path)) == "flye"


def test_circularity_decides_when_clipping_ties(tmp_path: Path) -> None:
    candidate(tmp_path, "flye", contigs=2, circular=0, qv=70.0)
    candidate(tmp_path, "raven", contigs=2, circular=1, qv=50.0)
    assert winner(run(tmp_path)) == "raven"


def test_contig_count_decides_when_circularity_ties(tmp_path: Path) -> None:
    candidate(tmp_path, "flye", contigs=5, circular=1, qv=70.0)
    candidate(tmp_path, "raven", contigs=2, circular=1, qv=50.0)
    assert winner(run(tmp_path)) == "raven"


def test_completeness_decides_when_structure_ties(tmp_path: Path) -> None:
    candidate(tmp_path, "flye", completeness=97.0, qv=70.0)
    candidate(tmp_path, "raven", completeness=99.0, qv=50.0)
    assert winner(run(tmp_path)) == "raven"


def test_qv_decides_when_completeness_ties(tmp_path: Path) -> None:
    candidate(tmp_path, "flye", qv=50.0)
    candidate(tmp_path, "raven", qv=55.0)
    assert winner(run(tmp_path)) == "raven"


def test_identical_candidates_tie_break_alphabetically(tmp_path: Path) -> None:
    candidate(tmp_path, "raven")
    candidate(tmp_path, "flye")
    rows = run(tmp_path)
    assert [row["assembler"] for row in rows] == ["flye", "raven"]
    assert winner(rows) == "flye"


def test_the_ranking_is_written_out_in_full(tmp_path: Path) -> None:
    candidate(tmp_path, "flye", clipping=1)
    candidate(tmp_path, "raven")
    rows = run(tmp_path)
    assert [row["rank"] for row in rows] == ["1", "2"]
    assert rows[0]["assembler"] == "raven"
    assert rows[0]["size_ratio"] == "1.0"
    assert rows[0]["contigs"] == "1"
    assert rows[0]["unmapped_read_fraction"] == "0.01"
    assert rows[1]["clipping_confirmed"] == "1"
