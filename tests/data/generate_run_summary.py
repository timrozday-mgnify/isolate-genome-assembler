#!/usr/bin/env python3
"""Build the report's fixture, ``tests/data/run_summary/``, through collect_metrics.py.

Three fictitious samples: ``iso_pass`` (clean, with a reference), ``iso_warn`` (a missing
plasmid and a little human) and ``iso_fail`` (mixed culture, unresolved consensus). Every
input is written in its tool's format, so this is also collect_metrics.py's end-to-end
test. Rerun after changing either script:

    python3 tests/data/generate_run_summary.py
"""

from __future__ import annotations

import base64
import gzip
import importlib.util
import json
import random
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "tests" / "data" / "run_summary"
THRESHOLDS = {
    "read_depth": {"direction": "below", "warn": 50, "fail": 25},
    "dominant_species_abundance": {"direction": "below", "fail": 0.9},
    "secondary_species_abundance": {"direction": "above", "warn": 0.01, "fail": 0.05},
    "human_read_fraction": {"direction": "above", "warn": 0.001},
    "consensus_unresolved": {"direction": "above", "warn": 1},
    "plasmids_missing": {"direction": "above", "warn": 1},
    "merqury_qv": {"direction": "below", "warn": 50, "fail": 40},
    "checkm2_completeness": {"direction": "below", "warn": 98, "fail": 90},
    "checkm2_contamination": {"direction": "above", "warn": 1, "fail": 5},
}
# A 1x1 transparent PNG, standing in for Bandage and Merqury images.
PNG = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg=="
)
ECOLI = "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli"
SAMPLES = {
    "iso_pass": {
        "depth": 120.0,
        "species": [("s__Escherichia coli", 99.6)],
        "human": 0.0,
        "resolved": True,
        "plasmids": [
            ("iso_pass_plasmid_1", 60_000, 1.1),
            ("iso_pass_plasmid_2", 4_000, 12.0),
        ],
        "missing": 0,
        "qv": 62.1,
        "checkm2": (99.9, 0.2),
        "reference": True,
    },
    "iso_warn": {
        "depth": 45.0,
        "species": [("s__Escherichia coli", 99.1)],
        "human": 0.004,
        "resolved": True,
        "plasmids": [("iso_warn_plasmid_1", 90_000, 0.8)],
        "missing": 1,
        "qv": 55.3,
        "checkm2": (99.5, 0.6),
        "reference": False,
    },
    "iso_fail": {
        "depth": 80.0,
        "species": [("s__Escherichia coli", 70.2), ("s__Klebsiella pneumoniae", 28.9)],
        "human": 0.0,
        "resolved": False,
        "plasmids": [],
        "missing": 0,
        "qv": 38.0,
        "checkm2": (97.0, 22.5),
        "reference": False,
    },
}
CHROMOSOME = 200_000


def write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def tsv(path: Path, rows: list[dict], columns: str = "sample") -> Path:
    """Rows as a TSV; with no rows, a header of the given columns only."""
    header = "\t".join(rows[0]) if rows else columns
    lines = [header] + ["\t".join(str(v) for v in row.values()) for row in rows]
    return write(path, "\n".join(lines) + "\n")


def load(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "bin" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sample_inputs(
    sample: str, s: dict, work: Path, rng: random.Random
) -> list[tuple[str, Path]]:
    d = work / sample
    files: list[tuple[str, Path]] = []

    def add(kind: str, path: Path) -> None:
        files.append((kind, path))

    replicons = [(f"{sample}_chromosome", CHROMOSOME, 1.0, "chromosome")] + [
        (name, length, ratio, "plasmid") for name, length, ratio in s["plasmids"]
    ]
    bases = int(s["depth"] * CHROMOSOME)
    add(
        "read-qc",
        tsv(
            d / "read_qc.tsv",
            [
                {
                    "sample": sample,
                    "reads_in": bases // 15_000,
                    "reads_kept": bases // 15_000 - 3,
                    "reads_below_min_qv": 3,
                    "bases": bases,
                    "read_n50": 15_800,
                    "read_length_mean": 15_000,
                    "read_length_min": 1_200,
                    "read_length_max": 41_000,
                    "q20_fraction": 0.99,
                    "q30_fraction": 0.91,
                    "gc_mean": 50.6,
                    "gc_sd": 2.1,
                    "duplicate_id_fraction": 0,
                    "duplicate_sequence_fraction": 0,
                    "adapter_fraction": 0.0004,
                    "genome_size_declared": "",
                    "genome_size_autocycler": CHROMOSOME + 3_000,
                    "genome_size_genomescope": CHROMOSOME - 2_000,
                    "genome_size_used": CHROMOSOME + 3_000,
                    "genome_size_source": "autocycler",
                    "genome_size_disagreement": 0.025,
                    "heterozygosity": 0.0001,
                    "depth": s["depth"],
                }
            ],
        ),
    )
    add(
        "gc-hist",
        tsv(
            d / "gc_hist.tsv",
            [
                {
                    "gc_percent_bin": gc,
                    "read_count": int(1000 * 2.718 ** (-((gc - 50.5) ** 2) / 4)),
                }
                for gc in range(44, 58)
            ],
        ),
    )

    total = sum(a for _, a in s["species"])
    species = [
        {
            "clade": f"{ECOLI.rsplit(';', 2)[0].replace(';', '|')}|g__{n.split()[0][3:]}|{n}",
            "species": n,
            "relative_abundance": a,
            "sequence_abundance": a,
        }
        for n, a in s["species"]
    ]
    add(
        "contamination",
        write(
            d / "contamination_summary.json",
            json.dumps(
                {
                    "sample": sample,
                    "species": species,
                    "dominant_species": species[0]["species"],
                    "dominant_abundance": species[0]["sequence_abundance"] / 100,
                    "secondary_species": species[1]["species"]
                    if len(species) > 1
                    else "",
                    "secondary_abundance": species[1]["sequence_abundance"] / 100
                    if len(species) > 1
                    else 0.0,
                    "unknown_fraction": round(1 - total / 100, 4),
                    "sylph_profile_rows": len(species),
                    "human_query_hits": int(s["human"] > 0),
                    "human_query_max_ani": 99.8 if s["human"] else 0.0,
                    "human_fraction": s["human"],
                    "expected_taxon": "g__Escherichia",
                    "expected_taxon_matches": True,
                    "remove_human_mode": "auto",
                    "remove_human_applied": s["human"] > 0,
                }
            ),
        ),
    )

    source = "autocycler" if s["resolved"] else "fallback_hifiasm"
    add(
        "assembly-source",
        tsv(
            d / "assembly_source.tsv",
            [
                {
                    "sample": sample,
                    "assembly_source": source,
                    "fully_resolved": str(s["resolved"]).lower(),
                    "reason": "consensus assembly fully resolved"
                    if s["resolved"]
                    else "consensus graph has 4 unitigs",
                }
            ],
        ),
    )
    add(
        "contigs",
        tsv(
            d / "contigs.tsv",
            [
                {
                    "sample": sample,
                    "contig": name,
                    "original_contig": f"contig_{i + 1}",
                    "replicon_type": kind,
                    "length": length,
                    "depth": round(s["depth"] * ratio, 1),
                    "circular": "true",
                    "flags": "",
                    "in_final_assembly": "true",
                }
                for i, (name, length, ratio, kind) in enumerate(replicons)
            ],
        ),
    )
    add(
        "contig-depth",
        tsv(
            d / "contig_depth.tsv",
            [
                {
                    "sample": sample,
                    "contig": name,
                    "replicon_type": kind,
                    "length": length,
                    "median_depth": round(s["depth"] * ratio, 1),
                    "depth_ratio": ratio,
                }
                for name, length, ratio, kind in replicons
            ],
        ),
    )
    # Every candidate assembly scored against the others; the winner is what an
    # unresolved consensus falls back to.
    scored = [
        ("autocycler", 1 + len(s["plasmids"]), 1 + len(s["plasmids"]), 0, 59.8),
        ("hifiasm", 1 + len(s["plasmids"]), 1 + len(s["plasmids"]), 0, 58.2),
        ("flye", 2 + len(s["plasmids"]), 1, 1, 55.1),
        ("raven", 8, 0, 3, 43.0),
    ]
    add(
        "full-assemblies",
        tsv(
            d / "full_assemblies.tsv",
            [
                {
                    "sample": sample,
                    "assembler": assembler,
                    "contigs": contigs,
                    "circular_contigs": circular,
                    "n50": CHROMOSOME,
                    "total_length": total,
                    "size_ratio": round(total / CHROMOSOME, 4),
                    "merqury_qv": qv,
                    "merqury_completeness": 99.1,
                    "unmapped_read_fraction": 0.004,
                    "clipping_confirmed": clipping,
                    "filtered": "",
                    "rank": rank,
                    "selected": str(rank == 1).lower(),
                }
                for rank, (assembler, contigs, circular, clipping, qv) in enumerate(
                    scored if s["resolved"] else scored[1:], start=1
                )
                for total in [
                    CHROMOSOME + sum(length for _, length, _ in s["plasmids"])
                ]
            ],
        ),
    )
    add(
        "plasmid-audit",
        tsv(
            d / "plasmid_audit.tsv",
            [
                {
                    "sample": sample,
                    "plassembler_contig": i + 1,
                    "length": length,
                    "copy_number": ratio,
                    "plsdb_hit": "NZ_CP000001.1",
                    "status": "recovered",
                    "matched_contig": name,
                    "identity": 99.99,
                    "coverage": 100.0,
                }
                for i, (name, length, ratio) in enumerate(s["plasmids"])
            ]
            + [
                {
                    "sample": sample,
                    "plassembler_contig": 9,
                    "length": 3_100,
                    "copy_number": 30.0,
                    "plsdb_hit": "None",
                    "status": "missing",
                    "matched_contig": "",
                    "identity": "",
                    "coverage": "",
                }
            ]
            * s["missing"],
        ),
    )

    add(
        "mapping",
        tsv(
            d / "mapping.tsv",
            [
                {
                    "sample": sample,
                    "reads": 1600,
                    "unmapped_reads": 8,
                    "unmapped_read_fraction": 0.005,
                    "bases": bases,
                    "unmapped_bases": 40_000,
                    "unmapped_base_fraction": 0.0017,
                }
            ],
        ),
    )
    add(
        "coverage-regions",
        tsv(
            d / "coverage_regions.tsv",
            [
                {
                    "sample": sample,
                    "contig": f"{sample}_chromosome",
                    "kind": "high",
                    "start": 120_000,
                    "end": 126_000,
                    "mean_depth": s["depth"] * 2.1,
                    "depth_ratio": 2.1,
                }
            ]
            if sample == "iso_fail"
            else [],
        ),
    )
    add("clipping", tsv(d / "clipping.tsv", []))
    add(
        "inspector",
        write(
            d / "summary_statistics",
            "Structural error\t0\nSmall-scale assembly error /per Mbp\t1.5\nQV\t58.2\n",
        ),
    )
    high = 12 if sample == "iso_fail" else 0
    add(
        "variants",
        tsv(
            d / "variant_summary.tsv",
            [
                {
                    "sample": sample,
                    "high_af": high,
                    "high_af_homopolymer_indels": high // 3,
                    "mixed_af": 3,
                    "mixed_af_homopolymer_indels": 1,
                }
            ],
        ),
    )
    add(
        "variant-sites",
        tsv(
            d / "variants.tsv",
            [
                {
                    "sample": sample,
                    "contig": f"{sample}_chromosome",
                    "position": rng.randrange(CHROMOSOME),
                    "ref": "A",
                    "alt": "AT" if i % 3 == 0 else "G",
                    "depth": 100,
                    "alt_depth": 60 if i < high else 30,
                    "af": 0.6 if i < high else 0.3,
                    "af_bin": "high" if i < high else "mixed",
                    "type": "indel" if i % 3 == 0 else "snv",
                    "homopolymer_length": 9 if i % 3 == 0 else 1,
                    "homopolymer_indel": str(i % 3 == 0).lower(),
                }
                for i in range(high + 3)
            ],
        ),
    )
    add(
        "merqury",
        write(
            d / f"{sample}.qv",
            f"{sample}\t12\t4800000\t{s['qv']}\t{10 ** (-s['qv'] / 10):.3g}\n",
        ),
    )
    add(
        "merqury-completeness",
        write(
            d / f"{sample}.completeness.stats",
            f"{sample}\tall\t4700000\t4800000\t97.9\n",
        ),
    )
    ratios = [
        min(1.2, max(0.2, rng.gauss(1.0, 0.08 if sample != "iso_fail" else 0.2)))
        for _ in range(300)
    ]
    truncated = sum(r < 0.9 for r in ratios)
    add(
        "ideel",
        tsv(
            d / "ideel_summary.tsv",
            [
                {
                    "sample": sample,
                    "proteins": 320,
                    "proteins_with_hit": 300,
                    "truncated": truncated,
                    "truncated_fraction": round(truncated / 300, 4),
                }
            ],
        ),
    )
    add(
        "ideel-ratios",
        tsv(
            d / "ideel.tsv",
            [
                {
                    "sample": sample,
                    "protein": f"P{i:04}",
                    "hit": f"sp|Q{i:05}|X",
                    "ratio": round(r, 3),
                }
                for i, r in enumerate(ratios)
            ],
        ),
    )
    add(
        "rrna-depth",
        tsv(
            d / "rrna_depth.tsv",
            [
                {
                    "sample": sample,
                    "contig": f"{sample}_chromosome",
                    "start": 10_000 * i,
                    "end": 10_000 * i + 1_540,
                    "product": "16S ribosomal RNA",
                    "depth": s["depth"],
                    "depth_ratio": 1.0,
                }
                for i in range(1, 4)
            ],
        ),
    )
    add(
        "bakta",
        write(
            d / "bakta.txt",
            (
                "Sequence(s):\nLength: 264000\nCount: 3\nGC: 50.6\nN50: 200000\nN ratio: 0.0\n"
                "coding density: 87.9\n\nAnnotation:\ntRNAs: 86\ntmRNAs: 1\nrRNAs: 22\nncRNAs: 70\n"
                "CRISPR arrays: 1\nCDSs: 250\npseudogenes: 4\nhypotheticals: 20\n"
            ),
        ),
    )
    completeness, contamination = s["checkm2"]
    add(
        "checkm2",
        tsv(
            d / "quality_report.tsv",
            [
                {
                    "Name": sample,
                    "Completeness": completeness,
                    "Contamination": contamination,
                    "Completeness_Model_Used": "Neural Network (Specific Model)",
                    "Translation_Table_Used": 11,
                    "Coding_Density": 0.879,
                    "Contig_N50": CHROMOSOME,
                    "Genome_Size": 264_000,
                    "GC_Content": 0.51,
                }
            ],
        ),
    )
    add(
        "busco",
        tsv(
            d / "batch_summary.txt",
            [
                {
                    "Input_file": f"{sample}.fasta",
                    "Dataset": "bacteria_odb12",
                    "Complete": 99.2,
                    "Single": 98.4 if sample != "iso_fail" else 90.1,
                    "Duplicated": 0.8 if sample != "iso_fail" else 9.1,
                    "Fragmented": 0.4,
                    "Missing": 0.4,
                    "n_markers": 124,
                }
            ],
        ),
    )
    add(
        "gtdbtk",
        tsv(
            d / "gtdbtk.bac120.summary.tsv",
            [
                {
                    "user_genome": sample,
                    "classification": ECOLI,
                    "closest_genome_reference": "GCF_000005845.2",
                    "closest_genome_ani": 98.9,
                    "closest_genome_af": 0.93,
                    "classification_method": "ani_screen",
                    "warnings": "N/A",
                }
            ],
        ),
    )

    windows = [
        {
            "contig": name,
            "start": w,
            "end": min(w + 1000, length),
            "depth": round(
                s["depth"]
                * ratio
                * (
                    2.1
                    if sample == "iso_fail" and 120_000 <= w < 126_000
                    else rng.uniform(0.85, 1.15)
                ),
                2,
            ),
        }
        for name, length, ratio, _ in replicons
        for w in range(0, length, 1000)
    ]
    bed = d / "regions.bed.gz"
    with gzip.open(bed, "wt") as handle:
        handle.writelines(
            f"{w['contig']}\t{w['start']}\t{w['end']}\t{w['depth']}\n" for w in windows
        )
    add("coverage", bed)

    add(
        "assembly-attempts",
        tsv(
            d / "assembly_attempts.tsv",
            [
                {
                    "sample": sample,
                    "assembler": assembler,
                    "subset": f"{n:02}",
                    "status": "failed" if (assembler, n) == ("canu", 3) else "ok",
                    "contigs": 3,
                }
                for assembler in ("flye", "hifiasm", "raven", "canu")
                for n in range(1, 5)
            ],
        ),
    )
    add(
        "autocycler",
        tsv(
            d / "autocycler_table.tsv",
            [
                {
                    "name": sample,
                    "input_read_count": 1600,
                    "pass_cluster_count": len(replicons),
                    "fail_cluster_count": 1,
                    "overall_clustering_score": 0.97 if s["resolved"] else 0.61,
                    "consensus_assembly_bases": 264_000,
                    "consensus_assembly_unitigs": len(replicons)
                    if s["resolved"]
                    else 7,
                    "consensus_assembly_fully_resolved": str(s["resolved"]).lower(),
                }
            ],
        ),
    )
    add(
        "assembler-contribution",
        tsv(
            d / "assembler_contribution.tsv",
            [
                {
                    "sample": sample,
                    "assembler": assembler,
                    "cluster": c + 1,
                    "subsets_total": 4,
                    "subsets_in_cluster": 4 if assembler != "raven" or c == 0 else 2,
                    "contigs_in_cluster": 4,
                    "contigs_trimmed_out": 0,
                }
                for assembler in ("flye", "hifiasm", "raven", "canu")
                for c in range(len(replicons))
            ],
        ),
    )
    (d / "bandage.png").write_bytes(PNG)
    add("image-bandage", d / "bandage.png")
    # One dotplot per QC-pass cluster, so the numbering of repeated image kinds is
    # exercised as well as the report section.
    for cluster in range(1, len(replicons) + 1):
        plot = d / f"{sample}_cluster_{cluster:03d}.png"
        plot.write_bytes(PNG)
        add("image-dotplot", plot)

    if s["reference"]:
        add(
            "reference-skani",
            tsv(
                d / "reference_skani.tsv",
                [
                    {
                        "Ref_file": "reference.fasta",
                        "Query_file": f"{sample}.fasta",
                        "ANI": 99.97,
                        "Align_fraction_ref": 98.1,
                        "Align_fraction_query": 99.2,
                    }
                ],
            ),
        )
        blocks = [
            (0, 90_000, 0, 90_000, "+"),
            (90_000, 140_000, 140_000, 90_000, "-"),
            (140_000, 200_000, 140_000, 200_000, "+"),
        ]
        add(
            "dotplot",
            write(
                d / "reference.paf",
                "".join(
                    f"{sample}_chromosome\t{CHROMOSOME}\t{qs}\t{qe}\t{strand}\tref_chromosome\t{CHROMOSOME}\t"
                    f"{min(ts, te)}\t{max(ts, te)}\t{qe - qs}\t{qe - qs}\t60\n"
                    for qs, qe, ts, te, strand in blocks
                ),
            ),
        )
    return files


def qc_json(sample: str, s: dict, work: Path) -> Path:
    qc_gates = load("qc_gates")
    measurements = {
        "read_depth": s["depth"],
        "dominant_species_abundance": s["species"][0][1] / 100,
        "secondary_species_abundance": s["species"][1][1] / 100
        if len(s["species"]) > 1
        else 0.0,
        "human_read_fraction": s["human"],
        "consensus_unresolved": 0.0 if s["resolved"] else 1.0,
        "plasmids_missing": float(s["missing"]),
        "merqury_qv": s["qv"],
        "checkm2_completeness": s["checkm2"][0],
        "checkm2_contamination": s["checkm2"][1],
    }
    result = qc_gates.evaluate(measurements, THRESHOLDS)
    return write(
        work / sample / f"{sample}.qc.json",
        json.dumps({"sample": sample} | result, indent=2),
    )


def generate(outdir: Path) -> None:
    rng = random.Random(1)
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        lines = []
        for sample, s in SAMPLES.items():
            files = sample_inputs(sample, s, work, rng) + [
                ("qc", qc_json(sample, s, work))
            ]
            lines += [f"{sample}\t{kind}\t{path}" for kind, path in files]
        manifest = write(work / "manifest.tsv", "\n".join(lines) + "\n")
        versions = write(
            work / "versions.tsv", "FLYE\tflye\t2.9.6\nQUARTO_REPORT\tquarto\t1.9.37\n"
        )
        run_info = write(
            work / "run_info.json",
            json.dumps(
                {
                    "pipeline": "timrozday/isolate-genome-assembler",
                    "pipeline_version": "0.0.0-dev",
                    "commit": "",
                    "revision": "",
                    "nextflow_version": "25.10.0",
                    "run_name": "fixture_run",
                    "start": "2026-09-21T12:00:00Z",
                    "profile": "test",
                    "container_engine": "docker",
                    "params": {
                        "assemblers": "flye,hifiasm,raven,canu",
                        "subsample_count": 4,
                        "gtdbtk_db": "/dbs/gtdbtk_r226",
                        "checkm2_db": "/dbs/checkm2.dmnd",
                    },
                },
                indent=2,
            ),
        )
        if outdir.exists():
            shutil.rmtree(outdir)
        load("collect_metrics").main(
            [
                "--manifest",
                str(manifest),
                "--software-versions",
                str(versions),
                "--run-info",
                str(run_info),
                "--outdir",
                str(outdir),
            ]
        )


if __name__ == "__main__":
    generate(Path(sys.argv[1]) if len(sys.argv) > 1 else FIXTURE)
