"""Choosing between the Autocycler consensus and the Flye fallback."""

import importlib.util
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[1] / "bin" / "select_assembly.py"
_spec = importlib.util.spec_from_file_location("select_assembly", BIN)
sa = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sa)


def write(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


def run(tmp_path: Path, **kwargs: Path) -> tuple[str, str]:
    output = tmp_path / "chosen.fasta"
    summary = tmp_path / "chosen.tsv"
    argv = ["--sample", "isolate01", "--output", str(output), "--summary", str(summary)]
    for name, value in kwargs.items():
        argv += [f"--{name}", str(value)]
    sa.main(argv)
    return output.read_text(), summary.read_text()


def test_a_fully_resolved_consensus_is_used(tmp_path: Path) -> None:
    consensus = write(tmp_path / "consensus.fasta", ">cluster_1\nACGT\n")
    metrics = write(tmp_path / "consensus.yaml", "consensus_assembly_fully_resolved: true\n")
    fallback = write(tmp_path / "flye.fasta", ">flye_full_1\nAAAA\n")

    assembly, summary = run(tmp_path, consensus=consensus, metrics=metrics, fallback=fallback)
    assert assembly == ">cluster_1\nACGT\n"
    assert summary.splitlines()[1].split("\t")[:3] == ["isolate01", "autocycler", "true"]


def test_an_unresolved_consensus_falls_back_to_flye(tmp_path: Path) -> None:
    consensus = write(tmp_path / "consensus.fasta", ">cluster_1\nACGT\n")
    metrics = write(tmp_path / "consensus.yaml", "consensus_assembly_fully_resolved: false\n")
    fallback = write(tmp_path / "flye.fasta", ">flye_full_1\nAAAA\n")

    assembly, summary = run(tmp_path, consensus=consensus, metrics=metrics, fallback=fallback)
    assert assembly == ">flye_full_1\nAAAA\n"
    row = summary.splitlines()[1].split("\t")
    assert row[1:3] == ["fallback_flye", "false"]
    assert "not fully resolved" in row[3]


def test_a_crashed_autocycler_leaves_no_verdict_and_falls_back(tmp_path: Path) -> None:
    metrics = tmp_path / "missing.yaml"
    fallback = write(tmp_path / "flye.fasta", ">flye_full_1\nAAAA\n")

    _, summary = run(tmp_path, metrics=metrics, fallback=fallback)
    row = summary.splitlines()[1].split("\t")
    assert row[1:3] == ["fallback_flye", ""]
    assert "no consensus metrics" in row[3]


def test_an_empty_consensus_is_not_mistaken_for_an_assembly(tmp_path: Path) -> None:
    consensus = write(tmp_path / "consensus.fasta", "")
    metrics = write(tmp_path / "consensus.yaml", "consensus_assembly_fully_resolved: true\n")
    fallback = write(tmp_path / "flye.fasta", ">flye_full_1\nAAAA\n")

    assembly, summary = run(tmp_path, consensus=consensus, metrics=metrics, fallback=fallback)
    assert assembly == ">flye_full_1\nAAAA\n"
    assert "empty" in summary.splitlines()[1]


def test_no_usable_assembly_at_all_is_an_error(tmp_path: Path) -> None:
    metrics = write(tmp_path / "consensus.yaml", "consensus_assembly_fully_resolved: false\n")
    fallback = write(tmp_path / "flye.fasta", "")

    with pytest.raises(SystemExit):
        run(tmp_path, metrics=metrics, fallback=fallback)
