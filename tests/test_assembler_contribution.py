"""Tests for bin/assembler_contribution.py, on an Autocycler v0.7.0-shaped directory."""

import csv
import importlib.util
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "bin" / "assembler_contribution.py"
_spec = importlib.util.spec_from_file_location("assembler_contribution", BIN)
assembler_contribution = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(assembler_contribution)

CLUSTERING = (
    "node_name\tpassing_clusters\tall_clusters\tsequence_id\tfile_name\tcontig_name\t"
    "length\ttrusted\tcluster_weight\tconsensus_weight\n"
    "1__flye_01.fasta__flye_01_1__5000000_bp\t1\t1\t1\tflye_01.fasta\tflye_01_1\t5000000\tfalse\t1\t1\n"
    "2__flye_02.fasta__flye_02_1__5000000_bp\t1\t1\t2\tflye_02.fasta\tflye_02_1\t5000000\tfalse\t1\t1\n"
    "3__raven_01.fasta__raven_01_1__5000000_bp\t1\t1\t3\traven_01.fasta\traven_01_1\t5000000\tfalse\t1\t1\n"
    "4__raven_02.fasta__raven_02_1__9000_bp\tnone\t2\t4\traven_02.fasta\traven_02_1\t9000\tfalse\t1\t1\n"
)

# trim kept flye_01, flye_02 and dropped raven_01.
TRIMMED = (
    "H\tVN:Z:1.0\n"
    "P\t1\t1+\t*\tLN:i:5000000\tFN:Z:flye_01.fasta\tHD:Z:flye_01_1 circular=true\tCL:i:1\n"
    "P\t2\t1+\t*\tLN:i:5000000\tFN:Z:flye_02.fasta\tHD:Z:flye_02_1\tCL:i:1\n"
)


def run(tmp_path: Path, trimmed: str | None = TRIMMED) -> list[dict[str, str]]:
    autocycler = tmp_path / "autocycler_out"
    cluster = autocycler / "clustering" / "qc_pass" / "cluster_001"
    cluster.mkdir(parents=True)
    (autocycler / "clustering" / "clustering.tsv").write_text(CLUSTERING)
    if trimmed is not None:
        (cluster / "2_trimmed.gfa").write_text(trimmed)
    output = tmp_path / "contribution.tsv"
    assembler_contribution.main(
        ["--sample", "s1", "--autocycler-dir", str(autocycler), "--output", str(output)]
    )
    with output.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_each_assembler_gets_a_row_per_pass_cluster(tmp_path: Path) -> None:
    rows = {row["assembler"]: row for row in run(tmp_path)}
    assert set(rows) == {"flye", "raven"}
    assert rows["flye"]["subsets_in_cluster"] == "2"
    assert rows["flye"]["subsets_total"] == "2"
    # raven_02 only reached a failed cluster, so one of raven's two subsets is here.
    assert rows["raven"]["subsets_in_cluster"] == "1"
    assert rows["raven"]["subsets_total"] == "2"


def test_contigs_dropped_by_trim_are_counted(tmp_path: Path) -> None:
    rows = {row["assembler"]: row for row in run(tmp_path)}
    assert rows["flye"]["contigs_trimmed_out"] == "0"
    assert rows["raven"]["contigs_trimmed_out"] == "1"


def test_without_a_trimmed_graph_trimming_is_left_blank(tmp_path: Path) -> None:
    rows = run(tmp_path, trimmed=None)
    assert all(row["contigs_trimmed_out"] == "" for row in rows)


def test_a_stub_directory_gives_an_empty_table(tmp_path: Path) -> None:
    output = tmp_path / "contribution.tsv"
    assembler_contribution.main(
        ["--sample", "s1", "--autocycler-dir", str(tmp_path), "--output", str(output)]
    )
    assert output.read_text().startswith("sample\tassembler")
    assert len(output.read_text().splitlines()) == 1
