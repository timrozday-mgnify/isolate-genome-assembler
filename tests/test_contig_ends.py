"""Tests for the end-window awk program in modules/local/contig_ends.

The program is run as shipped, extracted from the module, so the escaping the module
needs and the logic the test checks cannot drift apart.
"""

import re
import subprocess
from pathlib import Path

MODULE = (
    Path(__file__).resolve().parents[1]
    / "modules"
    / "local"
    / "contig_ends"
    / "main.nf"
)


def awk_program() -> str:
    """The awk program from the module, with Nextflow's escaping undone."""
    body = re.search(
        r"awk -v window=\S+ '(.*?)' \$\{assembly\}", MODULE.read_text(), re.S
    )
    assert body, "the awk program moved; update this test"
    return body.group(1).replace("\\$", "$").replace("\\\\n", "\\n")


def ends(tmp_path: Path, window: int, **contigs: str) -> dict[str, str]:
    """Run the program over one FASTA and read back the window records it wrote."""
    program = tmp_path / "ends.awk"
    program.write_text(awk_program())
    fasta = tmp_path / "in.fasta"
    fasta.write_text("".join(f">{name}\n{seq}\n" for name, seq in contigs.items()))
    out = subprocess.run(
        ["awk", "-v", f"window={window}", "-f", str(program), str(fasta)],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.splitlines()
    return dict(zip([line[1:] for line in out[::2]], out[1::2]))


def test_long_contig_uses_the_full_window(tmp_path: Path) -> None:
    seq = "ACGT" * 10
    records = ends(tmp_path, 10, c1=seq)
    assert records["c1_start"] == seq[:10]
    assert records["c1_end"] == seq[-10:]


def test_short_contig_windows_do_not_overlap(tmp_path: Path) -> None:
    # Shorter than 2 * window: the windows used to share sequence and self-align for a
    # trivial reason, flagging every such contig as end-overlapping.
    seq = "ACGTACGTACGTACG"  # 15 bases, window 10
    records = ends(tmp_path, 10, c1=seq)
    assert records["c1_start"] + records["c1_end"] == seq[:7] + seq[-7:]
    assert len(records["c1_start"]) == len(records["c1_end"]) == 7


def test_odd_length_contig_leaves_the_middle_base_out(tmp_path: Path) -> None:
    records = ends(tmp_path, 10, c1="ACGTA")
    assert records["c1_start"] == "AC"
    assert records["c1_end"] == "TA"


def test_single_base_contig_is_skipped(tmp_path: Path) -> None:
    assert ends(tmp_path, 10, c1="A", c2="ACGT") == {"c2_start": "AC", "c2_end": "GT"}
