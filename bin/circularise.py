#!/usr/bin/env python3
"""Trim the duplicated sequence a circular contig carries at its ends.

`autocycler trim` does this for the consensus, but it is cluster- and graph-shaped, so it
cannot be pointed at a single FASTA. This does the same one job from the self-alignment
`CONTIG_ENDS` already produces: a `<contig>_start` record that aligns to that contig's own
`<contig>_end` record, running from the very start of the contig to the very end of it, is
the wrap of a circular contig, and the sequence it covers is written twice.

Circularity here is derived from the sequence, not from the assembler's tag: Raven,
metaMDBG and LJA report no circularity at all, so a tag-based call would rank them last
whatever they assembled.

Nothing is cut on a guess. A match that is not anchored at both ends is an internal repeat
and is left alone; a match that fills the alignment window is not trimmed either, because
the true overlap is longer than the window and so unknowable from this PAF.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

# How far a match may start from the contig's first base, or stop short of its last, and
# still count as anchored there. Minimap2's ends wander by a few bases; 50 is well inside
# any real overlap and far outside a rounding difference.
SLOP = 50

FIELDS = [
    "contig",
    "length_before",
    "overlap",
    "length_after",
    "identity",
    "circular",
    "flag",
]


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """Read a FASTA into (header, sequence) pairs, joining multi-line sequences."""
    records: list[tuple[str, str]] = []
    header = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header = line[1:]
            chunks = []
        elif header is not None:
            chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


class Wrap:
    """A `<contig>_start` to `<contig>_end` self-alignment, judged."""

    def __init__(self, fields: list[str]) -> None:
        self.contig = fields[0][: -len("_start")]
        self.window = int(fields[1])
        self.qstart, self.qend = int(fields[2]), int(fields[3])
        self.strand = fields[4]
        self.tend = int(fields[8])
        block = int(fields[10])
        self.identity = int(fields[9]) / block if block else 0.0
        self.overlap = self.qend - self.qstart

    def anchored(self) -> bool:
        """Does the match run from the contig's first base to its last?

        A reverse hit is an inverted repeat rather than a wrap, and a match that starts
        inside the contig or stops short of its end is an internal repeat near the ends.
        """
        return (
            self.strand == "+"
            and self.qstart <= SLOP
            and self.tend >= self.window - SLOP
        )

    def fills_window(self) -> bool:
        """Is the overlap at least as long as the window that was aligned?

        Then the real overlap extends past what CONTIG_ENDS cut out, so its length cannot
        be read off this alignment and trimming the window would silently under-trim.
        """
        return self.overlap >= self.window - SLOP


def wraps(paf: Path, min_identity: float) -> dict[str, Wrap]:
    """The best-scoring candidate wrap for each contig, by identity then overlap."""
    best: dict[str, Wrap] = {}
    if not paf.exists():
        return best
    for line in paf.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        if not (fields[0].endswith("_start") and fields[5].endswith("_end")):
            continue
        if fields[0][: -len("_start")] != fields[5][: -len("_end")]:
            continue
        wrap = Wrap(fields)
        if wrap.identity < min_identity:
            continue
        current = best.get(wrap.contig)
        if current is None or (wrap.identity, wrap.overlap) > (
            current.identity,
            current.overlap,
        ):
            best[wrap.contig] = wrap
    return best


def retag(header: str, length: int, circular: bool) -> str:
    """Rewrite the length and circularity tags a normalised header carries."""
    parts = [
        part
        for part in header.split()
        if not part.startswith(("length=", "circular="))
    ]
    parts.insert(1, f"length={length}")
    parts.append(f"circular={'true' if circular else 'false'}")
    return " ".join(parts)


def circularise(
    records: list[tuple[str, str]], found: dict[str, Wrap]
) -> tuple[list[tuple[str, str]], list[dict[str, object]]]:
    """Trim every unambiguous wrap, and record what was done to every contig."""
    trimmed: list[tuple[str, str]] = []
    rows: list[dict[str, object]] = []
    for header, sequence in records:
        name = header.split()[0]
        wrap = found.get(name)
        row = {
            "contig": name,
            "length_before": len(sequence),
            "overlap": 0,
            "length_after": len(sequence),
            "identity": "",
            "circular": "false",
            "flag": "",
        }
        if wrap is not None:
            row["identity"] = f"{wrap.identity:.4f}"
            if not wrap.anchored():
                row["flag"] = "internal_repeat"
            elif wrap.fills_window():
                # Naming the cause, because the fix is a larger --contig_end_window.
                row["overlap"] = f"{wrap.window}+"
                row["flag"] = "overlap_exceeds_window"
            else:
                # Trim the tail, not the head: the start coordinate stays put, which is
                # what dnaapler's rotation downstream assumes.
                sequence = sequence[: len(sequence) - wrap.overlap]
                row["overlap"] = wrap.overlap
                row["length_after"] = len(sequence)
                row["circular"] = "true"
        trimmed.append((retag(header, len(sequence), row["circular"] == "true"), sequence))
        rows.append(row)
    return trimmed, rows


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly", type=Path, required=True, help="FASTA to trim")
    parser.add_argument(
        "--overlaps", type=Path, required=True, help="PAF written by CONTIG_ENDS"
    )
    parser.add_argument("--output", type=Path, required=True, help="Trimmed FASTA")
    parser.add_argument(
        "--table", type=Path, required=True, help="Per-contig circularity TSV"
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.95,
        help="Identity a start-end match must reach before anything is cut",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    records = read_fasta(args.assembly)
    trimmed, rows = circularise(records, wraps(args.overlaps, args.min_identity))
    with args.output.open("w") as handle:
        for header, sequence in trimmed:
            handle.write(f">{header}\n{sequence}\n")
    with args.table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
