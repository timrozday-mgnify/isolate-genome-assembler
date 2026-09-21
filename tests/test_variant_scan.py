"""Tests for bin/variant_scan.py."""

import csv
import importlib.util
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "bin" / "variant_scan.py"
_spec = importlib.util.spec_from_file_location("variant_scan", BIN)
variant_scan = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(variant_scan)

# Position 11-20 (1-based) is a 10 bp run of A.
ASSEMBLY = ">c1\nCGTCGTCGTC" + "A" * 10 + "CGTCGTCGTC\n"
HEADER = (
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
)


def record(pos: int, ref: str, alt: str, ad: str) -> str:
    return f"c1\t{pos}\t.\t{ref}\t{alt}\t0\t.\t.\tPL:AD\t0,0,0:{ad}\n"


def run(tmp_path: Path, records: str) -> tuple[list[dict], dict]:
    vcf = tmp_path / "s1.vcf"
    vcf.write_text(HEADER + records)
    assembly = tmp_path / "s1.fasta"
    assembly.write_text(ASSEMBLY)
    output, summary = tmp_path / "variants.tsv", tmp_path / "summary.tsv"
    variant_scan.main(
        [
            "--sample", "s1",
            "--vcf", str(vcf),
            "--assembly", str(assembly),
            "--output", str(output),
            "--summary", str(summary),
        ]
    )  # fmt: skip
    with output.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    with summary.open() as handle:
        return rows, next(csv.DictReader(handle, delimiter="\t"))


def test_a_majority_alt_allele_is_a_likely_error(tmp_path: Path) -> None:
    rows, summary = run(tmp_path, record(3, "T", "G,<*>", "2,28,0"))
    assert rows[0]["af_bin"] == "high"
    assert rows[0]["af"] == "0.933"
    assert summary["high_af"] == "1"


def test_a_minority_allele_between_20_and_50_percent_is_mixed(tmp_path: Path) -> None:
    rows, summary = run(tmp_path, record(3, "T", "G,<*>", "20,10,0"))
    assert rows[0]["af_bin"] == "mixed"
    assert summary["mixed_af"] == "1" and summary["high_af"] == "0"


def test_background_noise_below_20_percent_is_dropped(tmp_path: Path) -> None:
    rows, _ = run(tmp_path, record(3, "T", "G,<*>", "27,3,0"))
    assert rows == []


def test_exactly_half_counts_as_high(tmp_path: Path) -> None:
    rows, _ = run(tmp_path, record(3, "T", "G,<*>", "15,15,0"))
    assert rows[0]["af_bin"] == "high"


def test_shallow_sites_are_skipped(tmp_path: Path) -> None:
    rows, _ = run(tmp_path, record(3, "T", "G,<*>", "1,8,0"))
    assert rows == []


def test_a_placeholder_only_site_is_not_a_variant(tmp_path: Path) -> None:
    rows, _ = run(tmp_path, record(3, "T", "<*>", "30,0"))
    assert rows == []


def test_an_indel_in_a_long_homopolymer_is_tallied_separately(tmp_path: Path) -> None:
    # Anchor base at position 10 (C), deleting one A from the run that follows.
    rows, summary = run(tmp_path, record(10, "CA", "C,<*>", "5,25,0"))
    assert rows[0]["type"] == "indel"
    assert rows[0]["homopolymer_length"] == "10"
    assert rows[0]["homopolymer_indel"] == "true"
    assert summary["high_af_homopolymer_indels"] == "1"


def test_an_snv_is_never_a_homopolymer_indel(tmp_path: Path) -> None:
    rows, _ = run(tmp_path, record(15, "A", "C,<*>", "5,25,0"))
    assert rows[0]["homopolymer_indel"] == "false"


def test_the_strongest_alt_allele_is_the_one_reported(tmp_path: Path) -> None:
    rows, _ = run(tmp_path, record(3, "T", "G,C,<*>", "10,4,16,0"))
    assert rows[0]["alt"] == "C"
