#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = ["edlib>=1.3", "mappy>=2.24", "pyyaml>=6"]
# ///
"""Score the Phase 6 benchmark arms against the truth genomes.

Reads what dev/benchmark.nf and the pipeline wrote, and writes
dev/assembler_benchmark.csv (one row per sample x arm x truth replicon) and
dev/assembler_benchmark.md (the summary tables the defaults are chosen from).

For each sample and arm the delivered assembly is scored: the arm's Autocycler
consensus when it is fully resolved, otherwise the full-read Flye fallback, as the
pipeline would. Each truth replicon is matched to the contig that aligns to it best;
that contig is oriented and rotated to the truth's start and its global edit distance
taken with edlib. A replicon counts as recovered when that distance is within 1% of
its length. Cost comes from the pipeline trace (the arm's assembler tasks) plus the
arm's own consensus task.

    uv run dev/assembler_benchmark.py --benchmark benchmark --results results
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path

import edlib
import mappy
import yaml

DEV = Path(__file__).resolve().parent
# Plasmids under this length are the ones HiFi library prep loses and the plan tracks.
SMALL_PLASMID = 10_000
# A replicon further than this fraction of its length from the truth is not recovered.
MAX_DIVERGENCE = 0.01
# Pipeline processes whose tasks belong to one assembler on one subset.
ASSEMBLER_PROCESSES = {
    "FLYE": "flye",
    "HIFIASM": "hifiasm",
    "RAVEN": "raven",
    "CANU": "canu",
    "METAMDBG": "metamdbg",
    "PLASSEMBLER": "plassembler",
    "MYLOASM": "myloasm",
    "LJA": "lja",
    "MINIASM_OVERLAP": "miniasm",
    "MINIASM": "miniasm",
    "MINIPOLISH": "miniasm",
}
DURATION_UNITS = {"d": 86400, "h": 3600, "m": 60, "s": 1, "ms": 0.001}
MEMORY_UNITS = {"B": 1, "KB": 1e3, "MB": 1e6, "GB": 1e9, "TB": 1e12}


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """(header, sequence) pairs from a plain or gzipped FASTA; missing file -> []."""
    if not path.exists():
        return []
    opener = gzip.open if path.suffix == ".gz" else open
    records: list[tuple[str, list[str]]] = []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                records.append((line[1:], []))
            elif records:
                records[-1][1].append(line)
    return [(h, "".join(s).upper()) for h, s in records if s]


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def hours(duration: str) -> float:
    """Nextflow's '1h 2m 3.4s' / '120ms' as hours; '-' or '' as 0."""
    total = 0.0
    for value, unit in re.findall(r"([\d.]+)(ms|d|h|m|s)", duration or ""):
        total += float(value) * DURATION_UNITS[unit]
    return total / 3600


def gigabytes(memory: str) -> float:
    match = re.match(r"([\d.]+)\s*([KMGT]?B)", memory or "")
    return float(match.group(1)) * MEMORY_UNITS[match.group(2)] / 1e9 if match else 0.0


def cpu_fraction(percent: str) -> float:
    try:
        return float(percent.rstrip("%")) / 100
    except ValueError:
        return 0.0


def read_trace(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as handle:
        return [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["status"] in ("COMPLETED", "FAILED", "CACHED")
        ]


def short_process(row: dict[str, str]) -> str:
    return row["process"].rsplit(":", 1)[-1]


def arm_cost(
    trace: list[dict[str, str]],
    consensus_trace: list[dict[str, str]],
    sample: str,
    arm: dict,
) -> tuple[float, float]:
    """CPU-hours and critical-path wall hours of one arm on one sample.

    The assemblers run in parallel, so wall time is the slowest assembler chain (miniasm
    is three tasks in a row) plus the consensus. Failed attempts cost CPU too, so they
    count. The shared subsample and Flye fallback tasks are left out: every arm has them.
    """
    cpu = 0.0
    chains: dict[tuple[str, str], float] = defaultdict(float)
    for row in trace:
        assembler = ASSEMBLER_PROCESSES.get(short_process(row))
        tag_sample, _, subset = row["tag"].rpartition("_")
        if assembler not in arm["assemblers"] or tag_sample != sample:
            continue
        if assembler is None or subset not in arm["subsets"]:
            continue
        cpu += hours(row["realtime"]) * cpu_fraction(row["%cpu"])
        chains[(assembler, subset)] += hours(row["realtime"])

    consensus_wall = 0.0
    for row in consensus_trace:
        if row["tag"] == f"{sample}.{arm['arm']}":
            cpu += hours(row["realtime"]) * cpu_fraction(row["%cpu"])
            consensus_wall = max(consensus_wall, hours(row["realtime"]))
    return cpu, max(chains.values(), default=0.0) + consensus_wall


def rotate_to(contig: str, truth: str) -> str:
    """Rotate a contig so it starts where the truth does, if the truth's start is found."""
    probe = truth[: min(1000, len(truth))]
    hit = edlib.align(
        probe, contig + contig, mode="HW", task="locations", k=len(probe) // 10
    )
    if hit["editDistance"] < 0:
        return contig
    start = hit["locations"][0][0] % len(contig)
    return contig[start:] + contig[:start]


def score_assembly(
    contigs: list[tuple[str, str]], truth: list[tuple[str, str]], truth_path: Path
) -> tuple[list[dict], int, int]:
    """Per-replicon scores, plus the count and length of contigs matching no replicon."""
    aligner = mappy.Aligner(str(truth_path), preset="asm5")
    replicon_names = [header.split()[0] for header, _ in truth]
    # contig index -> (replicon, aligned bases, strand of the longest hit)
    best: dict[int, tuple[str, int, int]] = {}
    for index, (_, sequence) in enumerate(contigs):
        matched: dict[str, int] = defaultdict(int)
        longest: dict[str, tuple[int, int]] = {}
        for hit in aligner.map(sequence):
            matched[hit.ctg] += hit.mlen
            if hit.mlen > longest.get(hit.ctg, (0, 0))[0]:
                longest[hit.ctg] = (hit.mlen, hit.strand)
        if matched:
            replicon = max(matched, key=lambda r: matched[r])
            best[index] = (replicon, matched[replicon], longest[replicon][1])

    rows = []
    used = set()
    for name, (header, truth_seq) in zip(replicon_names, truth):
        candidates = [i for i, (r, _, _) in best.items() if r == name]
        row = {
            "replicon": name,
            "replicon_type": "plasmid" if "plasmid" in header.lower() else "chromosome",
            "replicon_length": len(truth_seq),
            "contig_length": "",
            "circular": False,
            "recovered": False,
            "edit_distance": "",
            "qv": "",
        }
        if candidates:
            index = max(candidates, key=lambda i: best[i][1])
            used.add(index)
            header_c, contig = contigs[index]
            if best[index][2] < 0:
                contig = revcomp(contig)
            distance = edlib.align(
                rotate_to(contig, truth_seq),
                truth_seq,
                mode="NW",
                task="distance",
                k=int(MAX_DIVERGENCE * len(truth_seq)),
            )["editDistance"]
            row["contig_length"] = len(contig)
            row["circular"] = "circular=true" in header_c.lower()
            if distance >= 0:
                row["recovered"] = True
                row["edit_distance"] = distance
                # Zero errors is reported as one, so the QV is a lower bound.
                row["qv"] = round(
                    -10 * math.log10(max(distance, 1) / len(truth_seq)), 1
                )
        rows.append(row)

    extra = [c for i, c in enumerate(contigs) if i not in used]
    return rows, len(extra), sum(len(s) for _, s in extra)


def fully_resolved(metrics: Path) -> bool:
    if not metrics.exists():
        return False
    data = yaml.safe_load(metrics.read_text()) or {}
    return bool(data.get("consensus_assembly_fully_resolved"))


def load_arms(path: Path) -> list[dict]:
    with path.open() as handle:
        return [
            {
                "arm": row["arm"],
                "assemblers": row["assemblers"].split(","),
                "subsets": [f"{i:02d}" for i in range(1, int(row["subsets"]) + 1)],
            }
            for row in csv.DictReader(handle, delimiter="\t")
        ]


def depth_label(sample: str) -> str:
    match = re.search(r"_(\d+)x$", sample)
    return f"{match.group(1)}x" if match else "real"


def score(args: argparse.Namespace) -> list[dict]:
    samples = yaml.safe_load((args.benchmark / "samples.yml").read_text())
    arms = load_arms(args.arms)
    trace = read_trace(args.results / "pipeline_info" / "trace.txt")
    consensus_trace = read_trace(
        args.benchmark / "pipeline_info" / "trace_consensus.txt"
    )

    rows = []
    for sample in samples:
        sample_id = sample["id"]
        truth_path = args.benchmark / sample["reference"]
        truth = read_fasta(truth_path)
        fallback = (
            args.results / "assemblies" / sample_id / "flye_full" / "flye_full.fasta"
        )
        deliveries: list[tuple[str, Path, tuple[dict, bool] | None]] = [
            (
                "pipeline",
                args.results
                / "assemblies"
                / sample_id
                / "final"
                / f"{sample_id}.fasta",
                None,
            )
        ]
        for arm in arms:
            prefix = (
                args.benchmark / "consensus" / sample_id / f"{sample_id}.{arm['arm']}"
            )
            resolved = fully_resolved(Path(f"{prefix}.consensus.yaml"))
            assembly = Path(f"{prefix}.consensus.fasta") if resolved else fallback
            deliveries.append((arm["arm"], assembly, (arm, resolved)))

        for arm_name, assembly, arm_info in deliveries:
            replicons, extra_n, extra_len = score_assembly(
                read_fasta(assembly), truth, truth_path
            )
            if arm_info:
                arm, resolved = arm_info
                cpu, wall = arm_cost(trace, consensus_trace, sample_id, arm)
                source = "consensus" if resolved else "fallback_flye"
            else:
                resolved, cpu, wall, source = "", "", "", "pipeline_final"
            for replicon in replicons:
                rows.append(
                    {
                        "sample": sample_id,
                        "depth": depth_label(sample_id),
                        "arm": arm_name,
                        "fully_resolved": resolved,
                        "assembly_source": source,
                        **replicon,
                        "extra_contigs": extra_n,
                        "extra_length": extra_len,
                        "cpu_hours": round(cpu, 3) if cpu != "" else "",
                        "wall_hours": round(wall, 3) if wall != "" else "",
                    }
                )
    return rows


def percent(numerator: int, denominator: int) -> str:
    return (
        f"{100 * numerator / denominator:.0f}% ({numerator}/{denominator})"
        if denominator
        else "-"
    )


def markdown_table(header: list[str], rows: list[list]) -> list[str]:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return lines + [""]


def arm_summary(rows: list[dict], arm: str) -> list:
    arm_rows = [r for r in rows if r["arm"] == arm]
    per_sample = {r["sample"]: r for r in arm_rows}
    small = [
        r
        for r in arm_rows
        if r["replicon_type"] == "plasmid" and r["replicon_length"] < SMALL_PLASMID
    ]
    recovered = [r for r in arm_rows if r["recovered"]]
    return [
        arm,
        percent(
            sum(r["fully_resolved"] is True for r in per_sample.values()),
            len(per_sample),
        ),
        percent(len(recovered), len(arm_rows)),
        percent(sum(r["recovered"] for r in small), len(small)),
        percent(sum(r["circular"] for r in recovered), len(recovered)),
        sum(r["edit_distance"] for r in recovered),
        sum(r["extra_contigs"] for r in per_sample.values()),
        f"{statistics.mean(r['cpu_hours'] for r in per_sample.values()):.1f}",
        f"{statistics.mean(r['wall_hours'] for r in per_sample.values()):.1f}",
    ]


def process_resources(trace: list[dict[str, str]]) -> list[list]:
    """Per-process peaks from the pipeline trace, for setting per-process resources."""
    by_process: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in trace:
        by_process[short_process(row)].append(row)
    table = []
    for process, tasks in sorted(by_process.items()):
        realtimes = [hours(t["realtime"]) for t in tasks]
        table.append(
            [
                process,
                len(tasks),
                sum(t["status"] == "FAILED" for t in tasks),
                f"{max(gigabytes(t['peak_rss']) for t in tasks):.1f}",
                f"{statistics.median(realtimes):.2f}",
                f"{max(realtimes):.2f}",
                f"{max(cpu_fraction(t['%cpu']) for t in tasks):.1f}",
            ]
        )
    return table


def write_markdown(
    rows: list[dict], arms: list[dict], results: Path, path: Path
) -> None:
    arm_names = [a["arm"] for a in arms]
    lines = [
        "# Assembler benchmark results",
        "",
        "Generated by `dev/assembler_benchmark.py`; do not edit by hand. The choices made",
        "from these tables are recorded under Phase 6 in `docs/implementation_plan.md`.",
        "",
        "Recovered: the best-matching contig is within 1% edit distance of the truth",
        f"replicon. Small plasmids are under {SMALL_PLASMID // 1000} kb. Edit distance is summed",
        "over recovered replicons. CPU and wall hours are per sample; wall is the slowest",
        "assembler chain plus the consensus, as on an unloaded cluster.",
        "",
        "## Arms, all samples",
        "",
    ]
    header = [
        "arm",
        "fully resolved",
        "replicons recovered",
        "small plasmids recovered",
        "circular (of recovered)",
        "edit distance",
        "extra contigs",
        "CPU h",
        "wall h",
    ]
    lines += markdown_table(header, [arm_summary(rows, arm) for arm in arm_names])

    for depth in sorted(
        {r["depth"] for r in rows}, key=lambda d: (d == "real", len(d), d)
    ):
        depth_rows = [r for r in rows if r["depth"] == depth]
        lines += [f"## Arms, {depth} samples", ""]
        lines += markdown_table(
            header, [arm_summary(depth_rows, arm) for arm in arm_names]
        )

    lines += [
        "## Pipeline run (9 assemblers x 6 subsets) against truth",
        "",
        "What the benchmark's own pipeline run delivered, beside its gate calls: a gate",
        "that passes a sample with unrecovered replicons or many errors is too loose.",
        "",
    ]
    status = {}
    gates: dict[str, list[str]] = defaultdict(list)
    summary_dir = results / "report" / "run_summary"
    if (summary_dir / "samples.tsv").exists():
        with (summary_dir / "samples.tsv").open() as handle:
            status = {r["sample"]: r for r in csv.DictReader(handle, delimiter="\t")}
    if (summary_dir / "qc_gates.tsv").exists():
        with (summary_dir / "qc_gates.tsv").open() as handle:
            for r in csv.DictReader(handle, delimiter="\t"):
                if r["status"] in ("warn", "fail"):
                    gates[r["sample"]].append(f"{r['check']}={r['status']}")
    pipeline_rows = []
    for sample in sorted({r["sample"] for r in rows}):
        replicons = [
            r for r in rows if r["sample"] == sample and r["arm"] == "pipeline"
        ]
        info = status.get(sample, {})
        pipeline_rows.append(
            [
                sample,
                info.get("status", "-"),
                info.get("assembly_source", "-"),
                percent(sum(r["recovered"] for r in replicons), len(replicons)),
                sum(r["edit_distance"] for r in replicons if r["recovered"]),
                replicons[0]["extra_contigs"] if replicons else "-",
                ", ".join(gates.get(sample, [])) or "-",
            ]
        )
    lines += markdown_table(
        [
            "sample",
            "status",
            "source",
            "recovered",
            "edit distance",
            "extra contigs",
            "warn/fail gates",
        ],
        pipeline_rows,
    )

    lines += ["## Per-process resources (pipeline trace)", ""]
    lines += markdown_table(
        [
            "process",
            "tasks",
            "failed",
            "peak RSS GB",
            "median h",
            "max h",
            "max CPUs used",
        ],
        process_resources(read_trace(results / "pipeline_info" / "trace.txt")),
    )
    path.write_text("\n".join(lines))


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--benchmark", type=Path, required=True, help="dev/benchmark.nf --outdir"
    )
    parser.add_argument(
        "--results", type=Path, required=True, help="the pipeline's --outdir"
    )
    parser.add_argument("--arms", type=Path, default=DEV / "benchmark_arms.tsv")
    parser.add_argument(
        "--output",
        type=Path,
        default=DEV / "assembler_benchmark",
        help="prefix for .csv and .md",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    rows = score(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with Path(f"{args.output}.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    write_markdown(rows, load_arms(args.arms), args.results, Path(f"{args.output}.md"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
