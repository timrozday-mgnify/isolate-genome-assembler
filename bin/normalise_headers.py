#!/usr/bin/env python3
"""Turn one assembler's native output into an Autocycler input assembly.

The per-assembler rules are ported from Autocycler v0.7.0 ``src/helper.rs`` so the
consensus sees exactly what ``autocycler helper`` would have given it: circularity and
depth tags where the assembler reports them, Canu's repeat/bubble contigs dropped and its
circular contigs trimmed, GFA segments converted to FASTA.

On top of the helper's behaviour the contigs are renamed ``<assembler>_<subset>_<n>`` so a
consensus contig can be traced back to the assembly it came from.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path

# Plassembler's circular contigs get double weight in clustering, so a small plasmid only
# it found still passes Autocycler's cluster QC.
PLASSEMBLER_CLUSTER_WEIGHT = 2

CANU_TRIM = re.compile(r"trim=(\d+)-(\d+)")


class Contig:
    """One assembled sequence plus the tags Autocycler reads from its header."""

    def __init__(
        self,
        sequence: str,
        circular: bool = False,
        depth: str | None = None,
        cluster_weight: int | None = None,
    ) -> None:
        self.sequence = sequence
        self.circular = circular
        self.depth = depth
        self.cluster_weight = cluster_weight

    def header(self, name: str) -> str:
        parts = [name, f"length={len(self.sequence)}"]
        if self.depth is not None:
            parts.append(f"depth={self.depth}")
        if self.circular:
            parts.append("circular=true")
        if self.cluster_weight is not None:
            parts.append(f"Autocycler_cluster_weight={self.cluster_weight}")
        return " ".join(parts)


def open_text(path: Path):
    """Open a plain or gzipped text file."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def load_fasta(path: Path) -> list[tuple[str, str]]:
    """Read a FASTA into (header, sequence) pairs, joining multi-line sequences."""
    records: list[tuple[str, str]] = []
    header = None
    chunks: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:]
                chunks = []
            elif header is not None:
                chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return [record for record in records if record[1]]


def gfa_contigs(path: Path) -> list[Contig]:
    """Convert a GFA's segment lines to contigs, as helper.rs's gfa_to_fasta does."""
    contigs = []
    with open_text(path) as handle:
        for line in handle:
            if not line.startswith("S"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or not fields[2]:
                continue
            name, sequence = fields[1], fields[2]
            depth = None
            for field in fields[3:]:
                for prefix in ("dp:f:", "rd:i:"):
                    if field.startswith(prefix):
                        depth = field[len(prefix) :]
                        break
                if depth is not None:
                    break
            # Raven and miniasm mark a circular segment with a trailing 'c'.
            contigs.append(Contig(sequence, circular=name.endswith("c"), depth=depth))
    return contigs


def plain_contigs(path: Path) -> list[Contig]:
    """Read contigs from a FASTA, keeping any circular/depth tags already in the header."""
    contigs = []
    for header, sequence in load_fasta(path):
        lowered = header.lower()
        depth = None
        match = re.search(r"depth=([0-9.]+)", lowered)
        if match:
            depth = match.group(1)
        contigs.append(Contig(sequence, circular="circular=true" in lowered, depth=depth))
    return contigs


def flye_contigs(fasta: Path, assembly_info: Path | None) -> list[Contig]:
    """Flye reports circularity and depth in assembly_info.txt, not in the FASTA header."""
    info: dict[str, tuple[bool, str]] = {}
    if assembly_info is not None and assembly_info.exists():
        for line in assembly_info.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            columns = line.split("\t")
            if len(columns) < 4:
                continue
            info[columns[0]] = (columns[3] == "Y", columns[2])

    contigs = []
    for header, sequence in load_fasta(fasta):
        name = header.split()[0]
        circular, depth = info.get(name, (False, None))
        contigs.append(Contig(sequence, circular=circular, depth=depth))
    return contigs


def canu_depths(tig_info: Path | None) -> dict[str, str]:
    """Map Canu contig names to their depth from *.contigs.layout.tigInfo."""
    depths: dict[str, str] = {}
    if tig_info is None or not tig_info.exists():
        return depths
    for line in tig_info.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        columns = line.split("\t")
        if len(columns) < 3:
            continue
        try:
            tig_id = int(columns[0])
        except ValueError:
            continue
        depths[f"tig{tig_id:08d}"] = columns[2]
    return depths


def canu_contigs(fasta: Path, tig_info: Path | None) -> list[Contig]:
    """Drop Canu's repeat and bubble contigs and trim the overlap off circular ones."""
    depths = canu_depths(tig_info)
    contigs = []
    for header, sequence in load_fasta(fasta):
        if "suggestRepeat=yes" in header or "suggestBubble=yes" in header:
            continue
        circular = "suggestCircular=yes" in header
        if circular:
            match = CANU_TRIM.search(header)
            if match:
                start, end = int(match.group(1)), int(match.group(2))
                if start < end <= len(sequence):
                    sequence = sequence[start:end]
        name = header.split()[0]
        contigs.append(Contig(sequence, circular=circular, depth=depths.get(name)))
    return contigs


def plassembler_contigs(fasta: Path) -> list[Contig]:
    """Plassembler's circular plasmids carry extra weight in Autocycler's clustering."""
    contigs = plain_contigs(fasta)
    for contig in contigs:
        if contig.circular:
            contig.cluster_weight = PLASSEMBLER_CLUSTER_WEIGHT
    return contigs


def first_match(directory: Path, *patterns: str) -> Path | None:
    """Return the first file matching any of the glob patterns, in the order given."""
    for pattern in patterns:
        matches = sorted(directory.glob(pattern))
        if matches:
            return matches[0]
    return None


def collect(assembler: str, directory: Path) -> list[Contig]:
    """Find the assembler's output in `directory` and read it with that tool's rules."""
    if assembler == "flye":
        fasta = first_match(directory, "assembly.fasta")
        return flye_contigs(fasta, first_match(directory, "assembly_info.txt")) if fasta else []

    if assembler == "canu":
        fasta = first_match(directory, "*.contigs.fasta")
        return canu_contigs(fasta, first_match(directory, "*.tigInfo")) if fasta else []

    if assembler == "plassembler":
        fasta = first_match(directory, "plassembler_plasmids.fasta")
        return plassembler_contigs(fasta) if fasta else []

    if assembler in ("hifiasm", "miniasm"):
        gfa = first_match(directory, "*.gfa")
        return gfa_contigs(gfa) if gfa else []

    # Raven writes its contigs to stdout and metaMDBG to a gzipped FASTA; neither tags
    # circularity there, which matches what `autocycler helper` passes on.
    fasta = first_match(directory, "*.fasta", "*.fasta.gz", "*.fa", "*.fa.gz")
    return plain_contigs(fasta) if fasta else []


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembler", required=True, help="Assembler that produced the input")
    parser.add_argument("--subset", required=True, help="Subsampled read set id, e.g. 01")
    parser.add_argument("--input-dir", type=Path, default=Path(), help="Directory of native output")
    parser.add_argument("--output", type=Path, required=True, help="Normalised FASTA to write")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    contigs = collect(args.assembler, args.input_dir)
    with args.output.open("w") as handle:
        for index, contig in enumerate(contigs, start=1):
            name = f"{args.assembler}_{args.subset}_{index}"
            handle.write(f">{contig.header(name)}\n{contig.sequence}\n")
    if not contigs:
        print(f"no contigs found for {args.assembler} subset {args.subset}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
