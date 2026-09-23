#!/usr/bin/env python3
"""Find positions where many reads are clipped at once: candidate misjoins.

A read that disagrees with the assembly past some point is soft- or hard-clipped there.
One such read is noise; a pile-up of them at the same place, making up a real fraction of
the depth, means the assembly probably joined two things that do not belong together.

Clips within ``END_MARGIN`` of a contig end are ignored. A read spanning the origin of a
circular contig is split there into two clipped alignments, and a read running off the
end of a linear contig is clipped by definition; neither says anything about a misjoin.

Two further measurements say how bad each pile-up is.

*How hard the aligner tried.* Clipping alone does not mean the reads stop matching: an
aligner also gives up at a local difference it could have crossed, such as a long indel
or a diverged patch, because Z-drop truncates the alignment there. So each pile-up is
re-measured in a second, permissive alignment of the same reads (``--extend-bam``,
mapped with Z-drop and the end bonus raised) that extends through what it can. A pile-up
that survives that is a real misjoin; one that disappears was a local difference that
merely ended the alignment. Without ``--extend-bam`` every pile-up is reported as
``confirmed``, so a missing second alignment never quietly clears a sample.

*Where the clipped tails go.* A read cut at a misjoin usually still aligns somewhere, as
a supplementary alignment more than ``TAIL_MARGIN`` away. If the tails agree on one
destination, that is the join the assembly should have made, and ``tail_target`` names
it. If the tails land nowhere, the sequence past the clip is not in the assembly at all:
a missing replicon, contamination, or adapter rather than a misplaced join.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam

END_MARGIN = 1000
# ponytail: clips are grouped into fixed bins, so a pile-up straddling a bin edge is
# split in two; cluster by distance if that ever hides a real misjoin.
BIN = 100
# How far away another alignment of the same read has to be to count as a tail rather
# than the aligner splitting one alignment around a large indel. Tails are also grouped
# at this size, so reads landing near each other agree on one target.
TAIL_MARGIN = 10_000
CLIP_OPS = (4, 5)  # soft clip, hard clip


def clip_positions(alignment: pysam.AlignedSegment, min_clip: int) -> list[int]:
    """Reference positions at which this alignment is clipped by at least min_clip."""
    cigar = alignment.cigartuples or []
    positions = []
    if cigar and cigar[0][0] in CLIP_OPS and cigar[0][1] >= min_clip:
        positions.append(alignment.reference_start)
    if cigar and cigar[-1][0] in CLIP_OPS and cigar[-1][1] >= min_clip:
        positions.append(alignment.reference_end)
    return positions


def clip_reads(
    bam: pysam.AlignmentFile, min_clip: int
) -> dict[tuple[str, int], list[str]]:
    """Names of the reads clipped in each (contig, bin), ignoring clips at a contig end.

    One entry per clip, so a read clipped at both ends of the same bin appears twice.
    """
    lengths = dict(zip(bam.references, bam.lengths))
    clips: dict[tuple[str, int], list[str]] = defaultdict(list)
    bam.reset()  # a whole-file pass has to start at the first record, not wherever we are
    for alignment in bam.fetch(until_eof=True):
        if alignment.is_unmapped or alignment.is_secondary:
            continue
        contig = str(alignment.reference_name)
        for position in clip_positions(alignment, min_clip):
            if END_MARGIN <= position <= lengths[contig] - END_MARGIN:
                clips[(contig, position // BIN)].append(str(alignment.query_name))
    return clips


def bin_depth(bam: pysam.AlignmentFile, contig: str, bin_index: int) -> int:
    length = dict(zip(bam.references, bam.lengths))[contig]
    start = bin_index * BIN
    return bam.count(contig, start, min(start + BIN, length))


def find_pileups(
    bam: pysam.AlignmentFile,
    clips: dict[tuple[str, int], list[str]],
    min_reads: int,
    min_fraction: float,
) -> list[dict[str, object]]:
    pileups = []
    for (contig, bin_index), names in sorted(clips.items()):
        clipped = len(names)
        if clipped < min_reads:
            continue
        depth = bin_depth(bam, contig, bin_index)
        fraction = clipped / max(depth, clipped, 1)
        if fraction >= min_fraction:
            pileups.append(
                {
                    "contig": contig,
                    "start": bin_index * BIN,
                    "end": bin_index * BIN + BIN,
                    "clipped_reads": clipped,
                    "depth": depth,
                    "clipped_fraction": round(fraction, 3),
                }
            )
    return pileups


def add_tails(
    bam: pysam.AlignmentFile,
    pileups: list[dict[str, object]],
    clips: dict[tuple[str, int], list[str]],
) -> None:
    """Record where else the clipped reads align, in place.

    ``tail_reads`` counts the clipped reads with an alignment elsewhere, ``tail_target``
    is the commonest such place and ``tail_agree`` how many of them landed there.
    """
    wanted = {
        name for row in pileups for name in clips[(row["contig"], row["start"] // BIN)]
    }
    placements: dict[str, list[tuple[str, int]]] = defaultdict(list)
    if wanted:
        bam.reset()
        for alignment in bam.fetch(until_eof=True):
            if alignment.is_unmapped or alignment.is_secondary:
                continue
            name = str(alignment.query_name)
            if name in wanted:
                placements[name].append(
                    (str(alignment.reference_name), alignment.reference_start)
                )

    for row in pileups:
        targets: Counter[tuple[str, int]] = Counter()
        landed = 0
        for name in set(clips[(row["contig"], row["start"] // BIN)]):
            elsewhere = {
                (contig, position // TAIL_MARGIN)
                for contig, position in placements[name]
                if contig != row["contig"] or abs(position - row["start"]) > TAIL_MARGIN
            }
            landed += bool(elsewhere)
            targets.update(elsewhere)
        row["tail_reads"] = landed
        (contig, block), agree = targets.most_common(1)[0] if targets else (("", 0), 0)
        row["tail_target"] = f"{contig}:{block * TAIL_MARGIN}" if contig else ""
        row["tail_agree"] = agree


def confirm(
    pileups: list[dict[str, object]],
    bam: pysam.AlignmentFile,
    min_clip: int,
    min_reads: int,
    min_fraction: float,
) -> None:
    """Re-measure each pile-up in the permissive alignment, in place.

    A pile-up still over both thresholds once the aligner has been pushed to extend
    through the position is ``confirmed``; one that is not is ``resolved``, meaning the
    reads do carry on matching the assembly and the clip marked a local difference
    rather than a misjoin.
    """
    clips = clip_reads(bam, min_clip)
    for row in pileups:
        key = (row["contig"], row["start"] // BIN)
        clipped = len(clips[key])
        depth = bin_depth(bam, *key)
        fraction = clipped / max(depth, clipped, 1)
        row["extend_clipped_reads"] = clipped
        row["extend_depth"] = depth
        row["extend_clipped_fraction"] = round(fraction, 3)
        row["verdict"] = (
            "confirmed"
            if clipped >= min_reads and fraction >= min_fraction
            else "resolved"
        )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--bam", type=Path, required=True, help="sorted, indexed")
    parser.add_argument(
        "--extend-bam",
        type=Path,
        help="the same reads mapped permissively; sorted, indexed",
    )
    parser.add_argument("--min-clip", type=int, default=500)
    parser.add_argument("--min-reads", type=int, default=5)
    parser.add_argument("--min-fraction", type=float, default=0.2)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    with pysam.AlignmentFile(str(args.bam)) as bam:
        clips = clip_reads(bam, args.min_clip)
        pileups = find_pileups(bam, clips, args.min_reads, args.min_fraction)
        add_tails(bam, pileups, clips)

    # No second alignment: nothing has been ruled out, so nothing is cleared.
    for row in pileups:
        row |= {
            "extend_clipped_reads": "",
            "extend_depth": "",
            "extend_clipped_fraction": "",
            "verdict": "confirmed",
        }
    if args.extend_bam:
        with pysam.AlignmentFile(str(args.extend_bam)) as extend:
            confirm(pileups, extend, args.min_clip, args.min_reads, args.min_fraction)

    fields = [
        "contig",
        "start",
        "end",
        "clipped_reads",
        "depth",
        "clipped_fraction",
        "tail_reads",
        "tail_target",
        "tail_agree",
        "extend_clipped_reads",
        "extend_depth",
        "extend_clipped_fraction",
        "verdict",
    ]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["sample", *fields], delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows({"sample": args.sample} | row for row in pileups)
    confirmed = sum(row["verdict"] == "confirmed" for row in pileups)
    if pileups:
        print(
            f"{len(pileups)} clipping pile-up(s), {confirmed} confirmed",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
