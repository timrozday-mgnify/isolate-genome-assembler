"""Tests for bin/coverage_regions.py."""

import csv
import gzip
import importlib.util
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "bin" / "coverage_regions.py"
_spec = importlib.util.spec_from_file_location("coverage_regions", BIN)
coverage_regions = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(coverage_regions)

CONTIGS = (
    "sample\tcontig\treplicon_type\n"
    "s1\ts1_chromosome\tchromosome\n"
    "s1\ts1_plasmid_1\tplasmid\n"
)


def windows(contig: str, depths: list, size: int = 1000) -> str:
    return "".join(
        f"{contig}\t{i * size}\t{(i + 1) * size}\t{depth}\n"
        for i, depth in enumerate(depths)
    )


def run(tmp_path: Path, bed: str) -> tuple[list[dict], list[dict]]:
    bed_path = tmp_path / "s1.regions.bed.gz"
    with gzip.open(bed_path, "wt") as handle:
        handle.write(bed)
    contigs = tmp_path / "contigs.tsv"
    contigs.write_text(CONTIGS)
    regions, depth = tmp_path / "regions.tsv", tmp_path / "depth.tsv"
    coverage_regions.main(
        [
            "--sample", "s1",
            "--windows", str(bed_path),
            "--contigs", str(contigs),
            "--regions-output", str(regions),
            "--depth-output", str(depth),
        ]
    )  # fmt: skip

    def read(path: Path) -> list[dict]:
        with path.open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    return read(regions), read(depth)


def test_even_depth_has_no_regions(tmp_path: Path) -> None:
    regions, _ = run(tmp_path, windows("s1_chromosome", [50] * 10))
    assert regions == []


def test_adjacent_low_windows_merge_into_one_region(tmp_path: Path) -> None:
    depths = [50] * 4 + [10, 12] + [50] * 4
    regions, _ = run(tmp_path, windows("s1_chromosome", depths))
    assert len(regions) == 1
    assert regions[0]["kind"] == "low"
    assert (regions[0]["start"], regions[0]["end"]) == ("4000", "6000")
    assert regions[0]["mean_depth"] == "11.0"


def test_a_collapsed_repeat_is_a_high_region(tmp_path: Path) -> None:
    depths = [50] * 4 + [110] + [50] * 4
    regions, _ = run(tmp_path, windows("s1_chromosome", depths))
    assert [r["kind"] for r in regions] == ["high"]


def test_low_and_high_next_to_each_other_stay_separate(tmp_path: Path) -> None:
    depths = [50] * 4 + [10, 120] + [50] * 4
    regions, _ = run(tmp_path, windows("s1_chromosome", depths))
    assert [r["kind"] for r in regions] == ["low", "high"]


def test_the_threshold_itself_is_not_flagged(tmp_path: Path) -> None:
    # Exactly 0.5x and 2x the median are the edges of normal, not outside it.
    depths = [50] * 4 + [25, 100] + [50] * 4
    regions, _ = run(tmp_path, windows("s1_chromosome", depths))
    assert regions == []


def test_plasmid_depth_is_relative_to_the_chromosome(tmp_path: Path) -> None:
    bed = windows("s1_chromosome", [50] * 10) + windows("s1_plasmid_1", [150] * 3)
    regions, depth = run(tmp_path, bed)
    by_contig = {row["contig"]: row for row in depth}
    assert by_contig["s1_chromosome"]["depth_ratio"] == "1.0"
    assert by_contig["s1_plasmid_1"]["depth_ratio"] == "3.0"
    assert by_contig["s1_plasmid_1"]["replicon_type"] == "plasmid"
    # A plasmid is judged against its own median, so high copy number is not a region.
    assert regions == []


def test_an_empty_mosdepth_stub_gives_empty_tables(tmp_path: Path) -> None:
    regions, depth = run(tmp_path, "\n")
    assert regions == [] and depth == []
