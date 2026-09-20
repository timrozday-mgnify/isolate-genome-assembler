"""The per-assembler rules ported from Autocycler v0.7.0's src/helper.rs."""

import gzip
import importlib.util
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "bin" / "normalise_headers.py"
_spec = importlib.util.spec_from_file_location("normalise_headers", BIN)
nh = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(nh)


def headers(path: Path) -> list[str]:
    return [line[1:] for line in path.read_text().splitlines() if line.startswith(">")]


def normalise(tmp_path: Path, assembler: str, subset: str = "01") -> Path:
    output = tmp_path / "out.fasta"
    nh.main(
        [
            "--assembler",
            assembler,
            "--subset",
            subset,
            "--input-dir",
            str(tmp_path),
            "--output",
            str(output),
        ]
    )
    return output


def test_flye_takes_circularity_and_depth_from_assembly_info(tmp_path: Path) -> None:
    (tmp_path / "assembly.fasta").write_text(">contig_1\nACGTACGT\n>contig_2\nAAAA\n")
    (tmp_path / "assembly_info.txt").write_text(
        "#seq_name\tlength\tcov.\tcirc.\ncontig_1\t8\t31\tY\ncontig_2\t4\t9\tN\n"
    )

    assert headers(normalise(tmp_path, "flye")) == [
        "flye_01_1 length=8 depth=31 circular=true",
        "flye_01_2 length=4 depth=9",
    ]


def test_canu_drops_repeats_and_trims_circular_contigs(tmp_path: Path) -> None:
    (tmp_path / "canu.contigs.fasta").write_text(
        ">tig00000001 len=8 suggestCircular=yes trim=2-6\nAACCGGTT\n"
        ">tig00000002 len=4 suggestRepeat=yes\nACGT\n"
        ">tig00000003 len=4 suggestBubble=yes\nACGT\n"
    )
    (tmp_path / "canu.contigs.layout.tigInfo").write_text(
        "#tigID\tlength\tcoverage\n1\t8\t42\n"
    )

    output = normalise(tmp_path, "canu")
    assert headers(output) == ["canu_01_1 length=4 depth=42 circular=true"]
    assert output.read_text().splitlines()[1] == "CCGG"


def test_gfa_segments_become_contigs_with_circular_and_depth_tags(
    tmp_path: Path,
) -> None:
    (tmp_path / "hifiasm.bp.p_ctg.gfa").write_text(
        "H\tVN:Z:1.0\n"
        "S\tptg000001c\tACGTACGT\tLN:i:8\tdp:f:30\n"
        "S\tptg000002l\tAAAA\tLN:i:4\n"
        "L\tptg000001c\t+\tptg000001c\t+\t0M\n"
    )

    assert headers(normalise(tmp_path, "hifiasm")) == [
        "hifiasm_01_1 length=8 depth=30 circular=true",
        "hifiasm_01_2 length=4",
    ]


def test_metamdbg_gzipped_contigs_are_read(tmp_path: Path) -> None:
    with gzip.open(tmp_path / "contigs.fasta.gz", "wt") as handle:
        handle.write(">ctg1\nACGT\nACGT\n")

    assert headers(normalise(tmp_path, "metamdbg")) == ["metamdbg_01_1 length=8"]


def test_plassembler_circular_plasmids_get_double_cluster_weight(
    tmp_path: Path,
) -> None:
    (tmp_path / "plassembler_plasmids.fasta").write_text(
        ">1 len=4 circular=true depth=120\nACGT\n>2 len=4\nAAAA\n"
    )

    assert headers(normalise(tmp_path, "plassembler")) == [
        "plassembler_01_1 length=4 depth=120 circular=true Autocycler_cluster_weight=2",
        "plassembler_01_2 length=4",
    ]


def test_a_failed_assembler_leaves_an_empty_file_not_an_error(tmp_path: Path) -> None:
    assert normalise(tmp_path, "raven").read_text() == ""
