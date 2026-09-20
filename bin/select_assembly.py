#!/usr/bin/env python3
"""Choose between the Autocycler consensus and the Flye fallback for one sample.

The consensus wins whenever ``autocycler combine`` reported
``consensus_assembly_fully_resolved: true`` and left a non-empty FASTA. Anything else —
clustering refused to run, a cluster failed to resolve, the process fell over — falls back
to the Flye assembly of the full read set, recorded as ``assembly_source=fallback_flye``.

No threshold lives here: the fallback is a warn, and that call is ``qc_gates.py``'s.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

RESOLVED = re.compile(r"^\s*consensus_assembly_fully_resolved:\s*(\S+)", re.MULTILINE)


def has_sequence(path: Path | None) -> bool:
    """True when the FASTA exists and holds at least one non-empty sequence."""
    if path is None or not path.exists():
        return False
    return any(line.strip() and not line.startswith(">") for line in path.read_text().splitlines())


def fully_resolved(metrics: Path | None) -> bool | None:
    """Read Autocycler's own verdict; None when it never got as far as writing one."""
    if metrics is None or not metrics.exists():
        return None
    match = RESOLVED.search(metrics.read_text())
    if match is None:
        return None
    return match.group(1).lower() == "true"


def select(consensus: Path | None, metrics: Path | None, fallback: Path | None) -> tuple[Path, str, str]:
    """Return the assembly to use, its source label, and why it was chosen."""
    resolved = fully_resolved(metrics)
    if resolved and consensus is not None and has_sequence(consensus):
        return consensus, "autocycler", "consensus assembly fully resolved"

    if resolved is None:
        reason = "autocycler produced no consensus metrics"
    elif not resolved:
        reason = "consensus assembly not fully resolved"
    else:
        reason = "consensus assembly is empty"

    if fallback is not None and has_sequence(fallback):
        return fallback, "fallback_flye", reason
    raise SystemExit(f"no assembly available: {reason}, and the Flye fallback is empty too")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--consensus", type=Path, help="autocycler consensus_assembly.fasta")
    parser.add_argument("--metrics", type=Path, help="autocycler consensus_assembly.yaml")
    parser.add_argument("--fallback", type=Path, help="Flye assembly of the full read set")
    parser.add_argument("--output", type=Path, required=True, help="Selected assembly FASTA")
    parser.add_argument("--summary", type=Path, required=True, help="One-row TSV of the choice")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    chosen, source, reason = select(args.consensus, args.metrics, args.fallback)
    args.output.write_text(chosen.read_text())
    resolved = fully_resolved(args.metrics)
    args.summary.write_text(
        "sample\tassembly_source\tfully_resolved\treason\n"
        f"{args.sample}\t{source}\t{'' if resolved is None else str(resolved).lower()}\t{reason}\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
