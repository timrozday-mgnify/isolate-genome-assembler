#!/usr/bin/env python3
"""Choose which candidate assembly one sample is finished from.

The candidates are the Autocycler consensus and every assembler's full-read assembly,
ranked against each other by ``score_assemblies.py``. ``--selection`` says how much weight
the consensus gets by right:

``always_score``
    Take the top-ranked candidate, consensus included. The consensus has to earn its place
    on every sample rather than by assumption.
``score``
    The consensus when ``autocycler combine`` reported
    ``consensus_assembly_fully_resolved: true`` and left a non-empty FASTA; otherwise the
    best-scoring full-read assembly.
``flye``
    The consensus on the same condition, otherwise Flye, whatever the scores say. Kept so
    a run made before scoring existed can be reproduced exactly.

Anything chosen over the consensus is recorded as ``assembly_source=fallback_<assembler>``.
Flye's full-read assembly is the last resort when there are no scores at all.

No threshold lives here: whether a fallback is a warn is ``qc_gates.py``'s call.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

RESOLVED = re.compile(r"^\s*consensus_assembly_fully_resolved:\s*(\S+)", re.MULTILINE)

# What the Autocycler consensus is called in the score table, where it is one candidate
# among the assemblers.
CONSENSUS = "autocycler"


def has_sequence(path: Path | None) -> bool:
    """True when the FASTA exists and holds at least one non-empty sequence."""
    if path is None or not path.exists():
        return False
    return any(
        line.strip() and not line.startswith(">")
        for line in path.read_text().splitlines()
    )


def fully_resolved(metrics: Path | None) -> bool | None:
    """Read Autocycler's own verdict; None when it never got as far as writing one."""
    if metrics is None or not metrics.exists():
        return None
    match = RESOLVED.search(metrics.read_text())
    if match is None:
        return None
    return match.group(1).lower() == "true"


def ranked(scores: Path | None) -> list[dict[str, str]]:
    """The scored candidates, best first; empty when nothing was scored."""
    if scores is None or not scores.exists() or not scores.read_text().strip():
        return []
    rows = list(csv.DictReader(scores.read_text().splitlines(), delimiter="\t"))
    return sorted(rows, key=lambda row: int(row.get("rank") or 0))


def evidence(row: dict[str, str]) -> str:
    """The part of a candidate's score a person would want to see in the reason."""
    parts = [
        f"{row.get('circular_contigs') or '0'} circular of {row.get('contigs') or '?'} contigs",
        f"QV {row.get('merqury_qv') or '?'}",
        f"{row.get('clipping_confirmed') or '0'} confirmed clipping pile-ups",
    ]
    return ", ".join(parts)


def best(
    candidates: dict[str, Path], scores: Path | None, allow_consensus: bool
) -> tuple[str, dict[str, str]] | None:
    """The top-ranked candidate that was actually delivered as a file."""
    for row in ranked(scores):
        assembler = row.get("assembler", "")
        if assembler == CONSENSUS and not allow_consensus:
            continue
        path = candidates.get(assembler)
        if path is not None and has_sequence(path):
            return assembler, row
    return None


def select(
    consensus: Path | None,
    metrics: Path | None,
    fallback: Path | None,
    candidates: dict[str, Path] | None = None,
    scores: Path | None = None,
    selection: str = "always_score",
) -> tuple[Path, str, str]:
    """Return the assembly to use, its source label, and why it was chosen."""
    candidates = candidates or {}
    resolved = fully_resolved(metrics)
    consensus_usable = bool(resolved) and consensus is not None and has_sequence(consensus)

    if selection == "always_score":
        # The consensus competes like any other candidate -- but only when Autocycler
        # resolved it. An unresolved consensus has unresolved bridges in it, which no
        # reference-free score sees, so it is not delivered on the strength of one.
        winner = best(candidates, scores, allow_consensus=consensus_usable)
        if winner is not None:
            assembler, row = winner
            if assembler == CONSENSUS:
                return candidates[assembler], "autocycler", (
                    f"highest-scoring assembly ({evidence(row)}); "
                    "the consensus also resolved"
                    if resolved
                    else f"highest-scoring assembly ({evidence(row)})"
                )
            return candidates[assembler], f"fallback_{assembler}", (
                f"{assembler} scored above the consensus ({evidence(row)}); "
                "the consensus has had autocycler resolve applied and a single assembly "
                "cannot, so the comparison is not like for like"
            )

    if consensus_usable:
        return consensus, "autocycler", "consensus assembly fully resolved"

    if resolved is None:
        reason = "autocycler produced no consensus metrics"
    elif not resolved:
        reason = "consensus assembly not fully resolved"
    else:
        reason = "consensus assembly is empty"

    if selection in ("score", "always_score"):
        winner = best(candidates, scores, allow_consensus=False)
        if winner is not None:
            assembler, row = winner
            return candidates[assembler], f"fallback_{assembler}", (
                f"{reason}; best full-read assembly: {assembler} ({evidence(row)})"
            )

    if fallback is not None and has_sequence(fallback):
        return fallback, "fallback_flye", reason
    raise SystemExit(
        f"no assembly available: {reason}, and the Flye fallback is empty too"
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument(
        "--consensus", type=Path, help="autocycler consensus_assembly.fasta"
    )
    parser.add_argument(
        "--metrics", type=Path, help="autocycler consensus_assembly.yaml"
    )
    parser.add_argument(
        "--fallback", type=Path, help="Flye assembly of the full read set"
    )
    parser.add_argument(
        "--candidates",
        nargs="*",
        default=[],
        metavar="ASSEMBLER=FASTA",
        help="Scored candidate assemblies, the consensus among them",
    )
    parser.add_argument("--scores", type=Path, help="score_assemblies.py table")
    parser.add_argument(
        "--selection",
        choices=["always_score", "score", "flye"],
        default="always_score",
        help="How much weight the consensus gets by right",
    )
    parser.add_argument(
        "--output", type=Path, required=True, help="Selected assembly FASTA"
    )
    parser.add_argument(
        "--summary", type=Path, required=True, help="One-row TSV of the choice"
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    candidates = dict(
        entry.split("=", 1) for entry in args.candidates if "=" in entry
    )
    chosen, source, reason = select(
        args.consensus,
        args.metrics,
        args.fallback,
        {name: Path(path) for name, path in candidates.items()},
        args.scores,
        args.selection,
    )
    args.output.write_text(chosen.read_text())
    resolved = fully_resolved(args.metrics)
    args.summary.write_text(
        "sample\tassembly_source\tfully_resolved\treason\n"
        f"{args.sample}\t{source}\t{'' if resolved is None else str(resolved).lower()}\t{reason}\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
