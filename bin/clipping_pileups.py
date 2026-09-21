#!/usr/bin/env python3
"""Find positions where many reads are clipped at once: candidate misjoins.

A read that disagrees with the assembly past some point is soft- or hard-clipped there.
One such read is noise; a pile-up of them at the same place, making up a real fraction of
the depth, means the assembly probably joined two things that do not belong together.

Clips within ``END_MARGIN`` of a contig end are ignored. A read spanning the origin of a
circular contig is split there into two clipped alignments, and a read running off the
end of a linear contig is clipped by definition; neither says anything about a misjoin.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

import pysam

END_MARGIN = 1000
# ponytail: clips are grouped into fixed bins, so a pile-up straddling a bin edge is
# split in two; cluster by distance if that ever hides a real misjoin.
BIN = 100
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


def find_pileups(
    bam: pysam.AlignmentFile, min_clip: int, min_reads: int, min_fraction: float
) -> list[dict[str, object]]:
    lengths = dict(zip(bam.references, bam.lengths))
    clips: Counter[tuple[str, int]] = Counter()
    for alignment in bam.fetch(until_eof=True):
        if alignment.is_unmapped or alignment.is_secondary:
            continue
        contig = str(alignment.reference_name)
        for position in clip_positions(alignment, min_clip):
            if END_MARGIN <= position <= lengths[contig] - END_MARGIN:
                clips[(contig, position // BIN)] += 1

    pileups = []
    for (contig, bin_index), clipped in sorted(clips.items()):
        if clipped < min_reads:
            continue
        start = bin_index * BIN
        depth = bam.count(contig, start, min(start + BIN, lengths[contig]))
        fraction = clipped / max(depth, clipped)
        if fraction >= min_fraction:
            pileups.append(
                {
                    "contig": contig,
                    "start": start,
                    "end": start + BIN,
                    "clipped_reads": clipped,
                    "depth": depth,
                    "clipped_fraction": round(fraction, 3),
                }
            )
    return pileups


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--bam", type=Path, required=True, help="sorted, indexed")
    parser.add_argument("--min-clip", type=int, default=500)
    parser.add_argument("--min-reads", type=int, default=5)
    parser.add_argument("--min-fraction", type=float, default=0.2)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    with pysam.AlignmentFile(str(args.bam)) as bam:
        pileups = find_pileups(bam, args.min_clip, args.min_reads, args.min_fraction)

    fields = ["contig", "start", "end", "clipped_reads", "depth", "clipped_fraction"]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["sample", *fields], delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows({"sample": args.sample} | row for row in pileups)
    if pileups:
        print(f"{len(pileups)} clipping pile-up(s)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
