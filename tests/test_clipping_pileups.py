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


CONTIGS = ("c1", "c2")


def make_bam(tmp_path: Path, alignments: list[tuple]) -> Path:
    """A sorted, indexed BAM of (start, cigar[, contig[, name]]) alignments.

    Giving a name makes the record a supplementary alignment of that read, which is how
    the other half of a read cut at a misjoin appears.
    """
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c, "LN": CONTIG_LENGTH} for c in CONTIGS],
    }
    unsorted = tmp_path / "unsorted.bam"
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as bam:
        for i, spec in enumerate(alignments):
            start, cigar = spec[0], spec[1]
            contig = spec[2] if len(spec) > 2 else CONTIGS[0]
            read = pysam.AlignedSegment(bam.header)
            read.query_name = spec[3] if len(spec) > 3 else f"read_{i}"
            read.is_supplementary = len(spec) > 3
            read.reference_id = CONTIGS.index(contig)
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


def run(
    tmp_path: Path,
    alignments: list[tuple[int, str]],
    extend: list[tuple[int, str]] | None = None,
) -> list[dict[str, str]]:
    output = tmp_path / "clipping.tsv"
    argv = ["--sample", "s1", "--bam", str(make_bam(tmp_path, alignments))]
    if extend is not None:
        directory = tmp_path / "extend"
        directory.mkdir()
        argv += ["--extend-bam", str(make_bam(directory, extend))]
    clipping_pileups.main([*argv, "--output", str(output)])
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


def test_a_pileup_the_permissive_alignment_extends_through_is_resolved(
    tmp_path: Path,
) -> None:
    # The same reads, mapped permissively, carry on through 10,000 instead of clipping:
    # a local difference the aligner gave up on, not a misjoin.
    clipped = [(10_000, "1000S3000M")] * 6
    rows = run(tmp_path, clipped + spanning(10), extend=spanning(16))
    assert [row["verdict"] for row in rows] == ["resolved"]
    assert rows[0]["clipped_reads"] == "6"
    assert rows[0]["extend_clipped_reads"] == "0"


def test_a_pileup_that_survives_the_permissive_alignment_is_confirmed(
    tmp_path: Path,
) -> None:
    clipped = [(10_000, "1000S3000M")] * 6
    rows = run(tmp_path, clipped + spanning(10), extend=clipped + spanning(10))
    assert [row["verdict"] for row in rows] == ["confirmed"]


def test_without_a_permissive_alignment_nothing_is_cleared(tmp_path: Path) -> None:
    rows = run(tmp_path, [(10_000, "1000S3000M")] * 6 + spanning(10))
    assert [row["verdict"] for row in rows] == ["confirmed"]
    assert rows[0]["extend_clipped_reads"] == ""


def tails(contig: str, start: int, names: list[str]) -> list[tuple]:
    """A supplementary alignment at contig:start for each named read."""
    return [(start, "1000M", contig, name) for name in names]


def test_tails_that_land_together_name_the_join_the_assembly_missed(
    tmp_path: Path,
) -> None:
    clipped = [(10_000, "1000S3000M")] * 6
    names = [f"read_{i}" for i in range(6)]
    rows = run(tmp_path, clipped + tails("c2", 5_000, names) + spanning(10))
    assert len(rows) == 1
    assert rows[0]["tail_reads"] == "6"
    assert rows[0]["tail_target"] == "c2:0"
    assert rows[0]["tail_agree"] == "6"


def test_tails_that_land_nowhere_are_sequence_outside_the_assembly(
    tmp_path: Path,
) -> None:
    rows = run(tmp_path, [(10_000, "1000S3000M")] * 6 + spanning(10))
    assert rows[0]["tail_reads"] == "0"
    assert rows[0]["tail_target"] == ""
    assert rows[0]["tail_agree"] == "0"


def test_an_alignment_beside_the_clip_is_not_a_tail(tmp_path: Path) -> None:
    # The aligner splitting one alignment around a large indel puts the other piece a
    # few kb away, not elsewhere in the assembly.
    clipped = [(10_000, "1000S3000M")] * 6
    names = [f"read_{i}" for i in range(6)]
    rows = run(tmp_path, clipped + tails("c1", 13_000, names) + spanning(10))
    assert rows[0]["tail_reads"] == "0"


def test_tails_are_found_wherever_they_sit_in_the_file(tmp_path: Path) -> None:
    # The pile-up is on the second contig and the tails on the first, so a tail scan
    # that carried on from wherever the clip scan left the file would miss them.
    clipped = [(10_000, "1000S3000M", "c2")] * 6
    names = [f"read_{i}" for i in range(6)]
    rows = run(
        tmp_path,
        clipped + tails("c1", 5_000, names) + [(8_000, "4000M", "c2")] * 10,
    )
    assert len(rows) == 1
    assert rows[0]["tail_reads"] == "6"
    assert rows[0]["tail_target"] == "c1:0"
