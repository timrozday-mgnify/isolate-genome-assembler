"""Tests for bin/ideel_ratios.py: the IDEEL length test and the rRNA depth check."""

import csv
import gzip
import importlib.util
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "bin" / "ideel_ratios.py"
_spec = importlib.util.spec_from_file_location("ideel_ratios", BIN)
ideel_ratios = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ideel_ratios)


def read(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(tmp_path: Path, **inputs: str) -> dict[str, list[dict[str, str]]]:
    argv = ["--sample", "s1"]
    for name, text in inputs.items():
        path = tmp_path / name
        if name == "windows":
            with gzip.open(path, "wt") as handle:
                handle.write(text)
        else:
            path.write_text(text)
        argv += [f"--{name}", str(path)]
    outputs = {k: tmp_path / f"{k}.tsv" for k in ("output", "summary", "rrna-output")}
    for key, path in outputs.items():
        argv += [f"--{key}", str(path)]
    ideel_ratios.main(argv)
    return {key: read(path) for key, path in outputs.items()}


PROTEINS = ">p1\nM\n>p2\nM\n>p3\nM\n>p4\nM\n"
# qseqid sseqid qlen slen: p2 is truncated to half its hit; p1's second hit is ignored.
HITS = (
    "p1\tsp|A\t300\t300\np1\tsp|B\t300\t100\np2\tsp|C\t150\t300\np3\tsp|D\t295\t300\n"
)


def test_ratios_use_each_proteins_best_hit_only(tmp_path: Path) -> None:
    result = run(tmp_path, proteins=PROTEINS, hits=HITS)
    ratios = {row["protein"]: float(row["ratio"]) for row in result["output"]}
    assert ratios == {"p1": 1.0, "p2": 0.5, "p3": 0.9833}


def test_the_truncated_fraction_is_over_proteins_with_a_hit(tmp_path: Path) -> None:
    summary = run(tmp_path, proteins=PROTEINS, hits=HITS)["summary"][0]
    assert summary["proteins"] == "4"
    assert summary["proteins_with_hit"] == "3"
    assert summary["truncated"] == "1"
    assert summary["truncated_fraction"] == "0.3333"


def test_no_hits_leaves_the_fraction_unmeasured(tmp_path: Path) -> None:
    summary = run(tmp_path, proteins=PROTEINS, hits="")["summary"][0]
    assert summary["truncated_fraction"] == ""


GFF = (
    "##gff-version 3\n"
    "c1\tBakta\trRNA\t4001\t5500\t.\t+\t.\tID=r1;product=16S ribosomal RNA\n"
    "c1\tBakta\tCDS\t1\t900\t.\t+\t0\tID=g1;product=thing\n"
    "##FASTA\n>c1\nACGT\n"
)


def test_a_collapsed_rrna_locus_shows_as_raised_depth(tmp_path: Path) -> None:
    depths = [50] * 4 + [300, 300] + [50] * 4
    windows = "".join(
        f"c1\t{i * 1000}\t{(i + 1) * 1000}\t{d}\n" for i, d in enumerate(depths)
    )
    rrna = run(tmp_path, gff=GFF, windows=windows)["rrna-output"]
    assert len(rrna) == 1
    assert rrna[0]["product"] == "16S ribosomal RNA"
    assert rrna[0]["depth_ratio"] == "6.0"
