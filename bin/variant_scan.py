#!/usr/bin/env python3
"""Bin the sites where reads disagree with the assembly by allele frequency.

Input is a ``bcftools mpileup -a FORMAT/AD`` VCF, not a genotype call: a haploid caller
reports the majority allele and so hides exactly the minority alleles this check is for.
The alternative-allele frequency comes straight from the read counts:

- **AF >= 0.5**: most reads disagree with the assembly, so it is a likely assembly error.
  The target is zero.
- **0.2 <= AF < 0.5**: a mixed strain, or reads from a collapsed repeat.

Indels in homopolymers of at least ``--homopolymer`` bases are tallied separately: they
are HiFi's known systematic error, so they say more about the reads than the assembly.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path

HIGH_AF = 0.5
MIXED_AF = 0.2


def load_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    name = None
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                sequences[name] = []
            elif name is not None:
                sequences[name].append(line.strip())
    return {name: "".join(parts).upper() for name, parts in sequences.items()}


def homopolymer_length(sequence: str, index: int) -> int:
    """Length of the run of identical bases that covers 0-based ``index``."""
    if not 0 <= index < len(sequence):
        return 0
    base = sequence[index]
    start = index
    while start > 0 and sequence[start - 1] == base:
        start -= 1
    end = index
    while end + 1 < len(sequence) and sequence[end + 1] == base:
        end += 1
    return end - start + 1


def read_sites(path: Path):
    """Yield (contig, pos, ref, alts, allele depths) for each VCF record with AD."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            keys = fields[8].split(":")
            values = dict(zip(keys, fields[9].split(":")))
            if "AD" not in values:
                continue
            depths = [int(v) if v.isdigit() else 0 for v in values["AD"].split(",")]
            yield fields[0], int(fields[1]), fields[3], fields[4].split(","), depths


def scan(
    vcf: Path, sequences: dict[str, str], min_depth: int, homopolymer: int
) -> list[dict[str, object]]:
    variants = []
    for contig, pos, ref, alts, depths in read_sites(vcf):
        total = sum(depths)
        # <*> is mpileup's placeholder for "any other allele", never a real one.
        real = [(d, alt) for d, alt in zip(depths[1:], alts) if alt not in ("<*>", ".")]
        if total < min_depth or not real:
            continue
        alt_depth, alt = max(real)
        af = alt_depth / total
        if af < MIXED_AF:
            continue

        indel = len(ref) != len(alt)
        # An indel's VCF position is the base before it, so the run starts one base on.
        run = homopolymer_length(sequences.get(contig, ""), pos) if indel else 0
        variants.append(
            {
                "contig": contig,
                "position": pos,
                "ref": ref,
                "alt": alt,
                "depth": total,
                "alt_depth": alt_depth,
                "af": round(af, 3),
                "af_bin": "high" if af >= HIGH_AF else "mixed",
                "type": "indel" if indel else "snv",
                "homopolymer_length": run,
                "homopolymer_indel": "true"
                if indel and run >= homopolymer
                else "false",
            }
        )
    return variants


def summary(variants: list[dict[str, object]]) -> dict[str, int]:
    def count(bin_: str, homopolymer: bool | None = None) -> int:
        return sum(
            v["af_bin"] == bin_
            and (
                homopolymer is None or (v["homopolymer_indel"] == "true") == homopolymer
            )
            for v in variants
        )

    return {
        "high_af": count("high"),
        "high_af_homopolymer_indels": count("high", True),
        "mixed_af": count("mixed"),
        "mixed_af_homopolymer_indels": count("mixed", True),
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--assembly", type=Path, required=True)
    parser.add_argument("--min-depth", type=int, default=10)
    parser.add_argument("--homopolymer", type=int, default=8)
    parser.add_argument("--output", type=Path, required=True, help="one row per site")
    parser.add_argument("--summary", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    variants = scan(
        args.vcf, load_fasta(args.assembly), args.min_depth, args.homopolymer
    )
    fields = [
        "contig",
        "position",
        "ref",
        "alt",
        "depth",
        "alt_depth",
        "af",
        "af_bin",
        "type",
        "homopolymer_length",
        "homopolymer_indel",
    ]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["sample", *fields], delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows({"sample": args.sample} | row for row in variants)

    counts = summary(variants)
    with args.summary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["sample", *counts], delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow({"sample": args.sample} | counts)
    return 0


if __name__ == "__main__":
    sys.exit(main())
