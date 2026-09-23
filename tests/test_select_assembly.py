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
    metrics = write(
        tmp_path / "consensus.yaml", "consensus_assembly_fully_resolved: true\n"
    )
    fallback = write(tmp_path / "flye.fasta", ">flye_full_1\nAAAA\n")

    assembly, summary = run(
        tmp_path, consensus=consensus, metrics=metrics, fallback=fallback
    )
    assert assembly == ">cluster_1\nACGT\n"
    assert summary.splitlines()[1].split("\t")[:3] == [
        "isolate01",
        "autocycler",
        "true",
    ]


def test_an_unresolved_consensus_falls_back_to_flye(tmp_path: Path) -> None:
    consensus = write(tmp_path / "consensus.fasta", ">cluster_1\nACGT\n")
    metrics = write(
        tmp_path / "consensus.yaml", "consensus_assembly_fully_resolved: false\n"
    )
    fallback = write(tmp_path / "flye.fasta", ">flye_full_1\nAAAA\n")

    assembly, summary = run(
        tmp_path, consensus=consensus, metrics=metrics, fallback=fallback
    )
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
    metrics = write(
        tmp_path / "consensus.yaml", "consensus_assembly_fully_resolved: true\n"
    )
    fallback = write(tmp_path / "flye.fasta", ">flye_full_1\nAAAA\n")

    assembly, summary = run(
        tmp_path, consensus=consensus, metrics=metrics, fallback=fallback
    )
    assert assembly == ">flye_full_1\nAAAA\n"
    assert "empty" in summary.splitlines()[1]


def test_no_usable_assembly_at_all_is_an_error(tmp_path: Path) -> None:
    metrics = write(
        tmp_path / "consensus.yaml", "consensus_assembly_fully_resolved: false\n"
    )
    fallback = write(tmp_path / "flye.fasta", "")

    with pytest.raises(SystemExit):
        run(tmp_path, metrics=metrics, fallback=fallback)


def scores(tmp_path: Path, *ranked: str) -> Path:
    """A score table listing the candidates best first."""
    header = "sample\tassembler\tcontigs\tcircular_contigs\tmerqury_qv\tclipping_confirmed\trank\tselected\n"
    rows = "".join(
        f"isolate01\t{assembler}\t1\t1\t58.2\t0\t{rank}\t{str(rank == 1).lower()}\n"
        for rank, assembler in enumerate(ranked, start=1)
    )
    return write(tmp_path / "full_assemblies.tsv", header + rows)


def candidates(tmp_path: Path, **assemblies: str) -> list[str]:
    """`--candidates` arguments for the given assembler -> sequence pairs."""
    arguments = []
    for assembler, seq in assemblies.items():
        path = write(tmp_path / f"{assembler}.fasta", f">{assembler}_1\n{seq}\n")
        arguments.append(f"{assembler}={path}")
    return arguments


def run_with(tmp_path: Path, *extra: str, **kwargs: Path) -> tuple[str, str]:
    output = tmp_path / "chosen.fasta"
    summary = tmp_path / "chosen.tsv"
    argv = ["--sample", "isolate01", "--output", str(output), "--summary", str(summary)]
    for name, value in kwargs.items():
        argv += [f"--{name}", str(value)]
    argv += list(extra)
    sa.main(argv)
    return output.read_text(), summary.read_text()


def test_the_top_ranked_candidate_wins_under_always_score(tmp_path: Path) -> None:
    metrics = write(
        tmp_path / "consensus.yaml", "consensus_assembly_fully_resolved: true\n"
    )
    consensus = write(tmp_path / "consensus.fasta", ">cluster_1\nACGT\n")
    assembly, summary = run_with(
        tmp_path,
        "--selection", "always_score",
        "--candidates", *candidates(tmp_path, hifiasm="GGGG", autocycler="ACGT"),
        consensus=consensus,
        metrics=metrics,
        scores=scores(tmp_path, "hifiasm", "autocycler"),
    )
    row = summary.splitlines()[1].split("\t")
    assert assembly == ">hifiasm_1\nGGGG\n"
    assert row[1] == "fallback_hifiasm"
    assert "not like for like" in row[3]


def test_a_resolved_consensus_that_scores_best_keeps_its_label(tmp_path: Path) -> None:
    metrics = write(
        tmp_path / "consensus.yaml", "consensus_assembly_fully_resolved: true\n"
    )
    consensus = write(tmp_path / "consensus.fasta", ">cluster_1\nACGT\n")
    _, summary = run_with(
        tmp_path,
        "--selection", "always_score",
        "--candidates", *candidates(tmp_path, autocycler="ACGT", hifiasm="GGGG"),
        consensus=consensus,
        metrics=metrics,
        scores=scores(tmp_path, "autocycler", "hifiasm"),
    )
    assert summary.splitlines()[1].split("\t")[1] == "autocycler"


def test_an_unresolved_consensus_takes_the_best_full_read_assembly(tmp_path: Path) -> None:
    metrics = write(
        tmp_path / "consensus.yaml", "consensus_assembly_fully_resolved: false\n"
    )
    consensus = write(tmp_path / "consensus.fasta", ">cluster_1\nACGT\n")
    fallback = write(tmp_path / "flye.fasta", ">flye_full_1\nAAAA\n")
    assembly, summary = run_with(
        tmp_path,
        "--selection", "score",
        "--candidates", *candidates(tmp_path, autocycler="ACGT", raven="TTTT"),
        consensus=consensus,
        metrics=metrics,
        fallback=fallback,
        scores=scores(tmp_path, "autocycler", "raven"),
    )
    row = summary.splitlines()[1].split("\t")
    # The consensus is top-ranked and still loses: it did not resolve.
    assert assembly == ">raven_1\nTTTT\n"
    assert row[1] == "fallback_raven"
    assert "best full-read assembly: raven" in row[3]


def test_selection_flye_ignores_the_scores_entirely(tmp_path: Path) -> None:
    metrics = write(
        tmp_path / "consensus.yaml", "consensus_assembly_fully_resolved: false\n"
    )
    fallback = write(tmp_path / "flye.fasta", ">flye_full_1\nAAAA\n")
    assembly, summary = run_with(
        tmp_path,
        "--selection", "flye",
        "--candidates", *candidates(tmp_path, hifiasm="GGGG"),
        metrics=metrics,
        fallback=fallback,
        scores=scores(tmp_path, "hifiasm"),
    )
    assert assembly == ">flye_full_1\nAAAA\n"
    assert summary.splitlines()[1].split("\t")[1] == "fallback_flye"


def test_flye_is_still_the_last_resort_when_nothing_was_scored(tmp_path: Path) -> None:
    metrics = write(
        tmp_path / "consensus.yaml", "consensus_assembly_fully_resolved: false\n"
    )
    fallback = write(tmp_path / "flye.fasta", ">flye_full_1\nAAAA\n")
    assembly, summary = run_with(
        tmp_path, "--selection", "always_score", metrics=metrics, fallback=fallback
    )
    assert assembly == ">flye_full_1\nAAAA\n"
    assert summary.splitlines()[1].split("\t")[1] == "fallback_flye"


def test_a_scored_candidate_whose_file_is_empty_is_skipped(tmp_path: Path) -> None:
    metrics = write(
        tmp_path / "consensus.yaml", "consensus_assembly_fully_resolved: false\n"
    )
    empty = write(tmp_path / "hifiasm.fasta", "")
    assembly, summary = run_with(
        tmp_path,
        "--selection", "always_score",
        "--candidates", f"hifiasm={empty}", *candidates(tmp_path, raven="TTTT"),
        metrics=metrics,
        scores=scores(tmp_path, "hifiasm", "raven"),
    )
    assert assembly == ">raven_1\nTTTT\n"
    assert summary.splitlines()[1].split("\t")[1] == "fallback_raven"


def test_an_unresolved_consensus_does_not_win_on_score(tmp_path: Path) -> None:
    # Ranking is reference-free and cannot see an unresolved bridge, so Autocycler's own
    # verdict still keeps the consensus out.
    metrics = write(
        tmp_path / "consensus.yaml", "consensus_assembly_fully_resolved: false\n"
    )
    consensus = write(tmp_path / "consensus.fasta", ">cluster_1\nACGT\n")
    assembly, summary = run_with(
        tmp_path,
        "--selection", "always_score",
        "--candidates", *candidates(tmp_path, autocycler="ACGT", raven="TTTT"),
        consensus=consensus,
        metrics=metrics,
        scores=scores(tmp_path, "autocycler", "raven"),
    )
    assert assembly == ">raven_1\nTTTT\n"
    assert summary.splitlines()[1].split("\t")[1] == "fallback_raven"
