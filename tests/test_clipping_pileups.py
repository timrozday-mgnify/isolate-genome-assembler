"""Tests for bin/clipping_pileups.py, on a BAM built in the test."""

import csv
import importlib.util
from pathlib import Path

import pytest

pysam = pytest.importorskip("pysam")

BIN = Path(__file__).resolve().parents[1] / "bin" / "clipping_pileups.py"
_spec = importlib.util.spec_from_file_location("clipping_pileups", BIN)
clipping_pileups = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(clipping_pileups)

CONTIG_LENGTH = 20_000


def make_bam(tmp_path: Path, alignments: list[tuple[int, str]]) -> Path:
    """A sorted, indexed BAM of (start, cigar) alignments on one contig."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "c1", "LN": CONTIG_LENGTH}],
    }
    unsorted = tmp_path / "unsorted.bam"
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as bam:
        for i, (start, cigar) in enumerate(alignments):
            read = pysam.AlignedSegment(bam.header)
            read.query_name = f"read_{i}"
            read.reference_id = 0
            read.reference_start = start
            read.cigarstring = cigar
            read.mapping_quality = 60
            length = read.infer_query_length()
            read.query_sequence = "A" * length
            read.query_qualities = pysam.qualitystring_to_array("I" * length)
            bam.write(read)
    sorted_bam = tmp_path / "sorted.bam"
    pysam.sort("-o", str(sorted_bam), str(unsorted))
    pysam.index(str(sorted_bam))
    return sorted_bam


def run(tmp_path: Path, alignments: list[tuple[int, str]]) -> list[dict[str, str]]:
    output = tmp_path / "clipping.tsv"
    bam = make_bam(tmp_path, alignments)
    clipping_pileups.main(
        ["--sample", "s1", "--bam", str(bam), "--output", str(output)]
    )
    with output.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def spanning(n: int) -> list[tuple[int, str]]:
    """n unclipped reads spanning position 10,000."""
    return [(8_000, "4000M")] * n


def test_a_pileup_of_clipped_reads_is_reported(tmp_path: Path) -> None:
    clipped = [(10_000, "1000S3000M")] * 6
    rows = run(tmp_path, clipped + spanning(10))
    assert len(rows) == 1
    assert rows[0]["start"] == "10000"
    assert rows[0]["clipped_reads"] == "6"


def test_right_hand_clips_count_at_the_alignment_end(tmp_path: Path) -> None:
    clipped = [(7_000, "3000M1000S")] * 6
    rows = run(tmp_path, clipped + spanning(10))
    assert [row["start"] for row in rows] == ["10000"]


def test_too_few_clipped_reads_is_not_a_pileup(tmp_path: Path) -> None:
    assert run(tmp_path, [(10_000, "1000S3000M")] * 4 + spanning(4)) == []


def test_a_small_fraction_of_deep_coverage_is_not_a_pileup(tmp_path: Path) -> None:
    # 5 clipped reads among 40 is 12.5%, under the 20% default.
    assert run(tmp_path, [(10_000, "1000S3000M")] * 5 + spanning(35)) == []


def test_short_clips_are_ignored(tmp_path: Path) -> None:
    assert run(tmp_path, [(10_000, "400S3000M")] * 10) == []


def test_clips_at_a_contig_end_are_ignored(tmp_path: Path) -> None:
    # Origin-spanning reads on a circular contig are clipped exactly here.
    assert run(tmp_path, [(0, "2000S3000M")] * 10) == []
