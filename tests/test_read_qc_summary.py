"""Tests for bin/read_qc_summary.py."""

import csv
import importlib.util
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[1] / "bin" / "read_qc_summary.py"
_spec = importlib.util.spec_from_file_location("read_qc_summary", BIN)
read_qc_summary = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(read_qc_summary)


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("5.2m", 5_200_000),
        ("500k", 500_000),
        ("2G", 2_000_000_000),
        ("4600000", 4_600_000),
        (None, None),
    ],
)
def test_parse_size(value, expected) -> None:
    assert read_qc_summary.parse_size(value) == expected


def test_parse_size_rejects_junk() -> None:
    with pytest.raises(ValueError, match="cannot parse genome size"):
        read_qc_summary.parse_size("about five megabases")


def test_gc_moments_weighs_by_read_count(tmp_path: Path) -> None:
    histogram = tmp_path / "gc.tsv"
    histogram.write_text("gc_percent_bin\tread_count\n40\t1\n60\t3\n")
    mean, sd = read_qc_summary.gc_moments(histogram)
    assert mean == pytest.approx(55.0)
    assert sd == pytest.approx(8.6602540, rel=1e-6)


def test_gc_moments_tolerates_a_missing_file(tmp_path: Path) -> None:
    assert read_qc_summary.gc_moments(tmp_path / "absent.tsv") == (None, None)


def test_parse_genomescope(tmp_path: Path) -> None:
    summary = tmp_path / "summary.txt"
    summary.write_text(
        "GenomeScope version 2.0\n"
        "property                      min          max\n"
        "Heterozygosity                0.0123%      0.0456%\n"
        "Genome Haploid Length         4,600,123 bp 4,700,000 bp\n"
    )
    assert read_qc_summary.parse_genomescope(summary) == (4_600_123.0, 0.0123)


def write_measurements(tmp_path: Path) -> dict[str, Path]:
    """The minimum set of upstream files the summary reads."""
    files = {
        "normalisation": "reads_in\treads_kept\treads_below_min_qv\tbases_kept\tmin_read_qv\n1100\t1000\t100\t500000\t20\n",
        "seqkit_stats": "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\tN50\tQ20(%)\tQ30(%)\n"
        "r.fastq.gz\tFASTQ\tDNA\t1000\t500000\t100\t500\t2000\t600\t99.5\t95.0\n",
        "gc_hist": "gc_percent_bin\tread_count\n50\t1000\n",
        "duplicates": "reads\tduplicate_ids\tduplicate_sequences\n1000\t10\t5\n",
        "adapters": "reads\tadapter_reads\tadapter_fraction\n1000\t2\t0.002\n",
    }
    paths = {}
    for name, content in files.items():
        path = tmp_path / f"{name}.tsv"
        path.write_text(content)
        paths[name] = path
    paths["genome_size_autocycler"] = tmp_path / "gs.txt"
    paths["genome_size_autocycler"].write_text("100000\n")
    paths["genomescope_summary"] = tmp_path / "summary.txt"
    paths["genomescope_summary"].write_text(
        "Genome Haploid Length 125,000 bp\nHeterozygosity 0.01%\n"
    )
    return paths


def run(tmp_path: Path, extra: list[str]) -> dict[str, str]:
    paths = write_measurements(tmp_path)
    output = tmp_path / "read_qc.tsv"
    argv = [
        "--sample",
        "isolate01",
        "--normalisation",
        str(paths["normalisation"]),
        "--seqkit-stats",
        str(paths["seqkit_stats"]),
        "--gc-hist",
        str(paths["gc_hist"]),
        "--duplicates",
        str(paths["duplicates"]),
        "--adapters",
        str(paths["adapters"]),
        "--genome-size-autocycler",
        str(paths["genome_size_autocycler"]),
        "--genomescope-summary",
        str(paths["genomescope_summary"]),
        "--output",
        str(output),
        *extra,
    ]
    assert read_qc_summary.main(argv) == 0
    with output.open() as handle:
        return next(iter(csv.DictReader(handle, delimiter="\t")))


def test_summary_prefers_the_samplesheet_genome_size(tmp_path: Path) -> None:
    row = run(tmp_path, ["--declared-genome-size", "200k"])
    assert row["genome_size_source"] == "samplesheet"
    assert float(row["genome_size_used"]) == 200_000
    assert float(row["depth"]) == pytest.approx(2.5)


def test_summary_falls_back_to_autocycler_and_flags_disagreement(
    tmp_path: Path,
) -> None:
    row = run(tmp_path, [])
    assert row["genome_size_source"] == "autocycler"
    assert float(row["genome_size_used"]) == 100_000
    # 100 kb vs 125 kb is a 20% disagreement, the point at which the gate warns.
    assert float(row["genome_size_disagreement"]) == pytest.approx(0.2)
    assert float(row["depth"]) == pytest.approx(5.0)
    assert float(row["duplicate_id_fraction"]) == pytest.approx(0.01)
    assert float(row["q20_fraction"]) == pytest.approx(0.995)
    assert float(row["read_n50"]) == 600


def test_summary_writes_every_column_even_with_no_inputs(tmp_path: Path) -> None:
    output = tmp_path / "read_qc.tsv"
    assert read_qc_summary.main(["--sample", "empty", "--output", str(output)]) == 0
    with output.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert list(rows[0]) == read_qc_summary.COLUMNS
    assert rows[0]["genome_size_source"] == "none"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
