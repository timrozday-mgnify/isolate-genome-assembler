"""Tests for bin/qc_gates.py: every threshold edge, and how missing data is treated."""

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("qc_gates", ROOT / "bin" / "qc_gates.py")
qc_gates = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(qc_gates)


def default_thresholds() -> dict:
    """assets/qc_thresholds.yml, parsed the way the workflow's snakeyaml would."""
    yaml = pytest.importorskip("yaml")
    return yaml.safe_load((ROOT / "assets" / "qc_thresholds.yml").read_text())


def write(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


def run(tmp_path: Path, thresholds: dict, *extra: str) -> dict:
    rules = write(tmp_path / "thresholds.json", json.dumps(thresholds))
    output = tmp_path / "qc.json"
    qc_gates.main(
        ["--sample", "s1", "--thresholds", str(rules), "--output", str(output), *extra]
    )
    return json.loads(output.read_text())


# --- gate(): the comparison itself ---


@pytest.mark.parametrize(
    ("value", "status"), [(50, "pass"), (49.9, "warn"), (25, "warn"), (24.9, "fail")]
)
def test_a_below_rule_trips_strictly_under_its_threshold(value, status) -> None:
    rule = {"direction": "below", "warn": 50, "fail": 25}
    assert qc_gates.gate("read_depth", value, rule)["status"] == status


@pytest.mark.parametrize(
    ("value", "status"),
    [(0.0099, "pass"), (0.01, "warn"), (0.0499, "warn"), (0.05, "fail")],
)
def test_an_above_rule_trips_at_its_threshold(value, status) -> None:
    rule = {"direction": "above", "warn": 0.01, "fail": 0.05}
    assert qc_gates.gate("secondary", value, rule)["status"] == status


def test_a_rule_with_only_a_fail_threshold_never_warns() -> None:
    rule = {"direction": "below", "fail": 0.9}
    assert qc_gates.gate("dominant", 0.91, rule)["status"] == "pass"
    assert qc_gates.gate("dominant", 0.89, rule)["status"] == "fail"


def test_a_missing_value_is_not_measured_rather_than_passing() -> None:
    entry = qc_gates.gate("merqury_qv", None, {"direction": "below", "warn": 50})
    assert entry["status"] == "not_measured"


def test_a_tripped_gate_says_what_tripped_it() -> None:
    entry = qc_gates.gate("merqury_qv", 45.0, {"direction": "below", "warn": 50})
    assert entry["message"] == "merqury_qv 45 < 50"
    assert entry["threshold"] == {"warn": 50}


def test_an_unknown_direction_is_an_error() -> None:
    with pytest.raises(ValueError):
        qc_gates.gate("x", 1.0, {"direction": "sideways", "warn": 1})


# --- evaluate(): the overall status ---


def test_overall_status_is_the_worst_entry(tmp_path: Path) -> None:
    rules = {
        "read_depth": {"direction": "below", "warn": 50, "fail": 25},
        "merqury_qv": {"direction": "below", "warn": 50, "fail": 40},
    }
    read_qc = write(tmp_path / "read_qc.tsv", "sample\tdepth\ns1\t40\n")
    result = run(tmp_path, rules, "--read-qc", str(read_qc))
    assert result["checks"]["read_depth"]["status"] == "warn"
    assert result["checks"]["merqury_qv"]["status"] == "not_measured"
    assert result["status"] == "warn"


def test_a_sample_with_nothing_measured_passes(tmp_path: Path) -> None:
    result = run(tmp_path, default_thresholds())
    assert result["status"] == "pass"
    assert all(c["status"] == "not_measured" for c in result["checks"].values())


def test_a_threshold_for_an_unknown_check_is_rejected(tmp_path: Path) -> None:
    with pytest.raises(SystemExit):
        run(tmp_path, {"not_a_check": {"direction": "above", "warn": 1}})


def test_the_default_thresholds_cover_every_check() -> None:
    args = qc_gates.parse_args(["--sample", "s", "--thresholds", "t", "--output", "o"])
    assert set(default_thresholds()) == set(qc_gates.measure(args))


# --- measure(): each reader ---


def test_the_fallback_assembly_counts_as_unresolved(tmp_path: Path) -> None:
    source = write(
        tmp_path / "source.tsv", "sample\tassembly_source\ns1\tfallback_hifiasm\n"
    )
    assert qc_gates.consensus_unresolved(source) == 1.0
    write(source, "sample\tassembly_source\ns1\tautocycler\n")
    assert qc_gates.consensus_unresolved(source) == 0.0


def test_a_resolved_consensus_that_lost_on_score_is_not_unresolved(
    tmp_path: Path,
) -> None:
    # always_score can deliver a better-scoring single assembly over a consensus that
    # resolved perfectly well; the gate measures the consensus, not the choice.
    source = write(
        tmp_path / "source.tsv",
        "sample\tassembly_source\tfully_resolved\ns1\tfallback_hifiasm\ttrue\n",
    )
    assert qc_gates.consensus_unresolved(source) == 0.0
    write(
        source,
        "sample\tassembly_source\tfully_resolved\ns1\tfallback_hifiasm\tfalse\n",
    )
    assert qc_gates.consensus_unresolved(source) == 1.0


def test_a_linear_chromosome_is_counted(tmp_path: Path) -> None:
    contigs = write(
        tmp_path / "contigs.tsv",
        "contig\treplicon_type\tcircular\n"
        "s1_chromosome\tchromosome\tfalse\n"
        "s1_plasmid_1\tplasmid\ttrue\n",
    )
    assert qc_gates.chromosome_not_circular(contigs) == 1.0


def test_the_lowest_plasmid_depth_ratio_is_used(tmp_path: Path) -> None:
    depth = write(
        tmp_path / "depth.tsv",
        "contig\treplicon_type\tdepth_ratio\n"
        "c\tchromosome\t1.0\np1\tplasmid\t3.2\np2\tplasmid\t0.4\n",
    )
    assert qc_gates.plasmid_depth_ratio(depth) == 0.4


def test_no_plasmids_means_the_depth_ratio_is_not_measured(tmp_path: Path) -> None:
    depth = write(tmp_path / "depth.tsv", "contig\treplicon_type\tdepth_ratio\n")
    assert qc_gates.plasmid_depth_ratio(depth) is None


def test_an_empty_region_table_is_zero_not_missing(tmp_path: Path) -> None:
    regions = write(tmp_path / "regions.tsv", "sample\tcontig\tkind\n")
    assert qc_gates.count_rows(regions) == 0
    assert qc_gates.count_rows(tmp_path / "absent.tsv") is None


def test_circular_contigs_are_counted_from_flye_assembly_info(tmp_path: Path) -> None:
    info = write(
        tmp_path / "assembly_info.txt",
        "#seq_name\tlength\tcov.\tcirc.\ncontig_1\t5000\t30\tY\ncontig_2\t900\t4\tN\n",
    )
    assert qc_gates.circular_flye_contigs(info) == 1.0


def test_inspector_structural_errors_are_read_from_its_summary(tmp_path: Path) -> None:
    summary = write(
        tmp_path / "summary_statistics",
        "Statics of contigs:\nStructural error\t3\nExpansion\t1\n\nQV\t52.1\n",
    )
    assert qc_gates.inspector_structural_errors(summary) == 3.0


def test_merqury_qv_is_the_fourth_column(tmp_path: Path) -> None:
    qv = write(tmp_path / "s1.qv", "s1\t12\t5000000\t57.3\t1.8e-06\n")
    assert qc_gates.merqury_qv(qv) == 57.3


def test_empty_stub_outputs_are_not_measured(tmp_path: Path) -> None:
    for reader in (qc_gates.merqury_qv, qc_gates.inspector_structural_errors):
        assert reader(write(tmp_path / "empty", "")) is None
    assert qc_gates.first(write(tmp_path / "empty.tsv", ""), "Completeness") is None


# --- taxon agreement ---

ECOLI = "d__Bacteria;p__Pseudomonadota;g__Escherichia;s__Escherichia coli"
SYLPH_ECOLI = {"species": [{"clade": ECOLI.replace(";", "|")}]}


def test_sylph_and_gtdbtk_agreeing_on_species_is_agreement() -> None:
    assert qc_gates.taxon_disagreement(SYLPH_ECOLI, ECOLI, None) == 0.0


def test_different_species_is_disagreement() -> None:
    other = ECOLI.replace("s__Escherichia coli", "s__Escherichia fergusonii")
    assert qc_gates.taxon_disagreement(SYLPH_ECOLI, other, None) == 1.0


def test_expected_taxon_can_name_a_genus() -> None:
    assert qc_gates.taxon_disagreement(SYLPH_ECOLI, ECOLI, "g__Escherichia") == 0.0
    assert qc_gates.taxon_disagreement(SYLPH_ECOLI, ECOLI, "g__Salmonella") == 1.0


def test_one_source_alone_cannot_disagree() -> None:
    assert qc_gates.taxon_disagreement(SYLPH_ECOLI, "", None) is None
    assert qc_gates.taxon_disagreement({}, "", "g__Escherichia") is None


def test_gtdbtk_summaries_skip_the_empty_domain(tmp_path: Path) -> None:
    archaea = write(tmp_path / "ar53.tsv", "")
    bacteria = write(
        tmp_path / "bac120.tsv", f"user_genome\tclassification\ns1\t{ECOLI}\n"
    )
    assert qc_gates.gtdbtk_classification([archaea, bacteria]) == ECOLI


def test_a_novel_species_leaves_the_species_gates_unmeasured(tmp_path: Path) -> None:
    rules = {
        "dominant_species_abundance": {"direction": "below", "fail": 0.9},
        "secondary_species_abundance": {"direction": "above", "warn": 0.01},
    }
    empty = {"species": [], "dominant_abundance": 0.0, "secondary_abundance": 0.0}
    contamination = write(tmp_path / "c.json", json.dumps(empty))
    result = run(tmp_path, rules, "--contamination", str(contamination))
    assert result["checks"]["dominant_species_abundance"]["status"] == "not_measured"
    assert result["status"] == "pass"


def test_a_placed_species_is_still_gated(tmp_path: Path) -> None:
    rules = {"dominant_species_abundance": {"direction": "below", "fail": 0.9}}
    profile = {"species": [{"clade": ECOLI}], "dominant_abundance": 0.8}
    contamination = write(tmp_path / "c.json", json.dumps(profile))
    result = run(tmp_path, rules, "--contamination", str(contamination))
    assert result["checks"]["dominant_species_abundance"]["status"] == "fail"
