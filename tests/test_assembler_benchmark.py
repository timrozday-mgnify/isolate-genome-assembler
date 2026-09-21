"""dev/assembler_benchmark.py: truth matching, edit distance and per-arm cost."""

import importlib.util
import random
from pathlib import Path

import pytest

pytest.importorskip("mappy")
pytest.importorskip("edlib")
pytest.importorskip("yaml")

SCRIPT = Path(__file__).resolve().parents[1] / "dev" / "assembler_benchmark.py"
_spec = importlib.util.spec_from_file_location("assembler_benchmark", SCRIPT)
ab = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ab)


def random_seq(length: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def test_rotated_reverse_contig_is_scored_and_missing_plasmid_is_not(
    tmp_path: Path,
) -> None:
    rng = random.Random(1)
    chromosome, plasmid = random_seq(20_000, rng), random_seq(3_000, rng)
    truth_path = tmp_path / "truth.fasta"
    truth_path.write_text(
        f">chr1 Example chromosome\n{chromosome}\n>p1 Example plasmid p1\n{plasmid}\n"
    )

    # Three substitutions, then rotated and reverse-complemented as an assembler might.
    mutated = list(chromosome)
    for position in (100, 5_000, 19_000):
        mutated[position] = "A" if mutated[position] != "A" else "C"
    contig = "".join(mutated)
    contig = ab.revcomp(contig[7_000:] + contig[:7_000])
    junk = random_seq(5_000, rng)

    rows, extra_n, extra_len = ab.score_assembly(
        [("1 length=20000 circular=true", contig), ("2 length=5000", junk)],
        ab.read_fasta(truth_path),
        truth_path,
    )

    chromosome_row, plasmid_row = rows
    assert chromosome_row["recovered"] and chromosome_row["circular"]
    assert chromosome_row["edit_distance"] == 3
    assert chromosome_row["replicon_type"] == "chromosome"
    assert plasmid_row["replicon_type"] == "plasmid"
    assert not plasmid_row["recovered"]
    assert (extra_n, extra_len) == (1, 5_000)


def test_arm_cost_counts_only_the_arms_assemblers_and_subsets() -> None:
    def task(process: str, tag: str, realtime: str, cpu: str) -> dict:
        return {
            "process": f"WF:ASSEMBLY:{process}",
            "tag": tag,
            "realtime": realtime,
            "%cpu": cpu,
        }

    trace = [
        task("FLYE", "s_1_01", "1h", "400%"),
        task("MINIASM_OVERLAP", "s_1_01", "30m", "100%"),
        task("MINIASM", "s_1_01", "1h 30m", "100%"),
        task("CANU", "s_1_01", "10h", "800%"),  # not in the arm
        task("FLYE", "s_1_05", "5h", "100%"),  # subset outside the arm
        task("FLYE", "s_2_01", "5h", "100%"),  # another sample
        task("FLYE_FULL", "s_1_full", "5h", "100%"),  # shared by every arm
    ]
    consensus = [{"tag": "s_1.fast", "realtime": "6m", "%cpu": "200%"}]
    arm = {"arm": "fast", "assemblers": ["flye", "miniasm"], "subsets": ["01", "02"]}

    cpu, wall = ab.arm_cost(trace, consensus, "s_1", arm)

    assert cpu == pytest.approx(4 + 0.5 + 1.5 + 0.2)
    assert wall == pytest.approx(2 + 0.1)  # the miniasm chain, then the consensus


def test_nextflow_durations_and_sizes_parse() -> None:
    assert ab.hours("1d 2h 30m") == pytest.approx(26.5)
    assert ab.hours("120ms") == pytest.approx(0.12 / 3600)
    assert ab.hours("-") == 0
    assert ab.gigabytes("512 MB") == pytest.approx(0.512)
