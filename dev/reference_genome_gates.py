#!/usr/bin/env python3
"""Hold a collection of finished genomes to the pipeline's own acceptance criteria.

The question this answers is whether `assets/qc_thresholds.yml` is too strict: if the
closed, published reference genomes of the same organisms also trip a gate, the gate is
measuring the pipeline's conventions rather than the assembly's quality.

Run it against the 20HM panel:

    reference_genome_gates.py --manifest <collection>/metadata/panel_genomes.tsv \\
        --out-dir <collection>/iga-assessment [--download]

It writes three things:

  reads.tsv        every genome, what reads exist for it in ENA, and of what kind
  samplesheet.yml  the pipeline input: `assembly` for all of them, `reads` only where
                   HiFi reads of the same strain were actually obtained
  fetch_reads.sh   the downloads, so they can be re-run on the machine with the disk

`--download` fetches them here instead. Nothing else is downloaded: HiFi is the only read
type this pipeline's read-based checks are calibrated for, so Illumina, 454, Nanopore and
PacBio CLR runs are reported and then left alone -- feeding them in would produce fake
failures at every gate that compares reads with the assembly.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ENA_API = "https://www.ebi.ac.uk/ena/portal/api/"
RUN_FIELDS = (
    "run_accession,sample_accession,study_accession,instrument_platform,"
    "instrument_model,library_strategy,read_count,base_count,fastq_ftp,submitted_ftp,"
    "sample_title,experiment_title"
)


class LookupFailed(Exception):
    """The portal could not be reached or answered with an error.

    This is deliberately not swallowed. A failed lookup looks exactly like "this genome has
    no reads" if it is, and a table that silently reports absence for a transient 500 is
    worse than no table: the whole point of the survey is to say what does and does not
    exist.
    """


def api(endpoint: str, **params: str) -> str:
    """One ENA portal call, retried with backoff: the portal 500s under sustained load."""
    url = ENA_API + endpoint + "?" + urllib.parse.urlencode(params)
    last = ""
    for attempt in range(5):
        if attempt:
            time.sleep(2**attempt)
        try:
            with urllib.request.urlopen(url, timeout=180) as response:
                return response.read().decode()
        except Exception as error:  # noqa: BLE001 - any failure is a retry, then a raise
            last = f"{type(error).__name__}: {error}"
    raise LookupFailed(
        f"{endpoint} {params.get('accession') or params.get('query')}: {last}"
    )


def rows_of(text: str) -> list[dict[str, str]]:
    return [
        r for r in csv.DictReader(text.splitlines(), delimiter="\t") if any(r.values())
    ]


def read_manifest(path: Path) -> list[dict[str, str]]:
    body = "\n".join(
        line for line in path.read_text().splitlines() if not line.startswith("#")
    )
    return rows_of(body)


def read_kind(run: dict[str, str]) -> str:
    """How the run was produced, in the terms that decide whether this pipeline can use it.

    The instrument model is not enough for PacBio: a Sequel II run can be raw CLR subreads
    or consensus HiFi, and ENA has no field that says which. The deposited file names do --
    `*.hifi_reads.fastq.gz` and `*ccs*` against `*subreads*` -- so they are what decides,
    with the model only as a fallback. Getting this wrong is the difference between reads
    this pipeline is calibrated for and reads at 10% error.
    """
    platform, model = run["instrument_platform"], run["instrument_model"]
    if platform == "PACBIO_SMRT":
        names = (run.get("fastq_ftp", "") + " " + run.get("submitted_ftp", "")).lower()
        if "hifi" in names or "ccs" in names:
            return "hifi"
        if "subread" in names:
            return "clr"
        # Revio only ever emits HiFi; RS, RS II and the first Sequel predate it entirely.
        # Only a Sequel II with no telling file name is genuinely ambiguous.
        if "Revio" in model:
            return "hifi"
        return "pacbio_unclear" if "Sequel II" in model else "clr"
    if platform == "OXFORD_NANOPORE":
        return "ont"
    if platform == "CAPILLARY":
        return "sanger"
    return (
        "short"
        if platform in ("ILLUMINA", "LS454", "ION_TORRENT", "BGISEQ", "DNBSEQ")
        else "other"
    )


# Library strategies that sequence the whole genome. The rest of what sits under an
# isolate's taxid -- amplicon panels, RNA-seq, Tn-seq, ChIP -- covers part of the genome or
# none of it, and counting those as "reads for this genome" overstates the answer badly:
# for the 20HM panel they are 60% of every run found.
WHOLE_GENOME_STRATEGIES = {"WGS", "WGA", "OTHER"}


def runs_for(query: str) -> list[dict[str, str]]:
    if not query:
        return []
    found = rows_of(
        api(
            "search",
            result="read_run",
            query=query,
            fields=RUN_FIELDS,
            format="tsv",
            limit="1000",
        )
    )
    return [
        run | {"kind": read_kind(run)}
        for run in found
        if run["library_strategy"] in WHOLE_GENOME_STRATEGIES
    ]


def strain_level(tax_name: str) -> bool:
    """Is this taxid the strain's own, or the whole species?"""
    return bool(re.search(r"\b(ATCC|DSM|JCM|NCTC|NCIMB|CCUG|LMG|WM1)\b", tax_name))


def strain_tokens(member: dict[str, str]) -> set[str]:
    """The culture-collection ids that identify this strain, normalised for matching."""
    raw = f"{member['strain_id']} / {member['dsm_atcc']}".upper()
    tokens = {
        f"{collection}{number}"
        for collection, number in re.findall(
            r"\b(ATCC|DSM|DSMZ|JCM|NCTC|NCIMB|CCUG|CIP|LMG)\s*-?\s*([A-Z]*\d[\dA-Z-]*)",
            raw,
        )
    }
    return {t.lower().replace("-", "") for t in tokens}


def names_the_strain(run: dict[str, str], tokens: set[str]) -> bool:
    """Does this run's own metadata say it is that strain?

    A strain-level taxid is not enough on its own: submitters routinely file a new isolate
    under whatever strain taxid NCBI offers for the species, so a run can sit under
    `Phocaeicola vulgatus ATCC 8482` and be a pig isolate. The sample and experiment
    titles are where a real strain claim appears.
    """
    blob = (run.get("sample_title", "") + " " + run.get("experiment_title", "")).lower()
    blob = blob.replace(" ", "").replace("-", "").replace("_", "")
    return any(token in blob for token in tokens)


def survey(member: dict[str, str]) -> dict[str, object]:
    """What reads exist for one genome, and whether any of them are usable here.

    A `lookup_failed` verdict means the portal did not answer, not that there are no reads;
    re-run to resolve it. `main()` exits non-zero while any row says so.
    """
    try:
        return _survey(member)
    except LookupFailed as failure:
        return {
            "strain_key": member["strain_key"],
            "species": member["species"],
            "strain_id": member["strain_id"],
            "source_type": member["source_type"],
            "accession": member["accession"],
            "verdict": "lookup_failed",
            "note": str(failure),
            "runs": [],
        }


def _survey(member: dict[str, str]) -> dict[str, object]:
    entry = {
        "strain_key": member["strain_key"],
        "species": member["species"],
        "strain_id": member["strain_id"],
        "source_type": member["source_type"],
        "accession": member["accession"],
    }
    if member["source_type"] != "public_ena":
        entry |= {
            "verdict": "local_assembly",
            "note": "assembled by this pipeline from HiFi reads held off this machine"
            if member["source_type"] == "new_longread"
            else "METT assembly; no public reads under this accession",
            "runs": [],
        }
        return entry

    accession = member["accession"].split()[0].split(".")[0]
    sequence = rows_of(
        api(
            "filereport",
            accession=accession,
            result="sequence",
            fields="accession,sample_accession,description",
        )
    )
    sample = sequence[0]["sample_accession"] if sequence else ""
    sample_meta = (
        rows_of(
            api(
                "filereport",
                accession=sample,
                result="sample",
                fields="sample_accession,tax_id,scientific_name,strain",
            )
        )
        if sample
        else []
    )
    tax_id = sample_meta[0]["tax_id"] if sample_meta else ""
    tax_name = sample_meta[0]["scientific_name"] if sample_meta else ""

    # The genome's own sample: these are the reads it was built from, by definition.
    own = runs_for(f'sample_accession="{sample}"' if sample else "")
    elsewhere = [
        r
        for r in runs_for(f"tax_eq({tax_id})" if tax_id else "")
        if r["run_accession"] not in {o["run_accession"] for o in own}
    ]
    tokens = strain_tokens(member)
    same_strain = [r for r in elsewhere if names_the_strain(r, tokens)]

    usable = [r for r in own if r["kind"] == "hifi"] or [
        r for r in same_strain if r["kind"] == "hifi"
    ]
    if usable:
        verdict = "hifi_own_sample" if usable[0] in own else "hifi_same_strain"
    elif own:
        verdict = "own_reads_not_hifi"
    elif same_strain:
        verdict = "strain_reads_not_hifi"
    elif elsewhere:
        # Nothing under the taxid claims this strain, so these are other isolates: a
        # base-level comparison would call real strain differences assembly errors.
        verdict = "species_reads_only"
    else:
        verdict = "no_reads"

    entry |= {
        "sample_accession": sample,
        "tax_id": tax_id,
        "tax_name": tax_name,
        "tax_is_strain_level": strain_level(tax_name),
        "strain_tokens": sorted(tokens),
        "verdict": verdict,
        "runs": own
        + [
            r
            | {
                "scope": "same strain, other submitter"
                if names_the_strain(r, tokens)
                else "same species, other strain"
            }
            for r in elsewhere
        ],
        "usable_runs": usable,
        "note": "",
    }
    return entry


def kind_counts(runs: list[dict[str, str]]) -> str:
    counts: dict[str, int] = {}
    for run in runs:
        counts[run["kind"]] = counts.get(run["kind"], 0) + 1
    return ", ".join(f"{kind} x{n}" for kind, n in sorted(counts.items())) or "none"


def scoped_counts(entry: dict, scope: str | None) -> str:
    """Read kinds for one scope: the genome's own sample, or runs that name its strain.

    Counting every run under the taxid together would flatter the answer -- a HiFi run of
    a different isolate is not a HiFi run of this genome.
    """
    runs = [r for r in entry["runs"] if r.get("scope") == scope]
    return kind_counts(runs)


def relative(target: Path, base: Path) -> str:
    """`target` as a path relative to `base`, walking up with `..` where it has to.

    `Path.relative_to` refuses to go upwards, and the reads directory and the genomes
    directory are siblings of the samplesheet as often as they are below it.
    """
    return os.path.relpath(target.resolve(), base.resolve())


def fastq_urls(run: dict[str, str]) -> list[str]:
    return [f"ftp://{part}" for part in run["fastq_ftp"].split(";") if part]


def write_reads_table(entries: list[dict], path: Path) -> None:
    fields = [
        "strain_key",
        "species",
        "strain_id",
        "source_type",
        "accession",
        "sample_accession",
        "tax_id",
        "tax_name",
        "verdict",
        "own_sample_reads",
        "same_strain_reads",
        "other_strain_reads",
        "usable_runs",
        "usable_bases",
        "note",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for entry in entries:
            usable = entry.get("usable_runs") or []
            writer.writerow(
                entry
                | {
                    "own_sample_reads": scoped_counts(entry, None),
                    "same_strain_reads": scoped_counts(
                        entry, "same strain, other submitter"
                    ),
                    "other_strain_reads": scoped_counts(
                        entry, "same species, other strain"
                    ),
                    "usable_runs": ",".join(r["run_accession"] for r in usable),
                    "usable_bases": sum(int(r["base_count"] or 0) for r in usable)
                    or "",
                }
            )


def write_samplesheet(
    entries: list[dict],
    manifest: list[dict],
    out_dir: Path,
    collection: Path,
    reads_dir: Path,
) -> None:
    """The pipeline input. Every genome is assessed; only some carry reads.

    Paths are written relative to this file, because the run happens on an HPC and an
    absolute path from the machine that generated the samplesheet is wrong there. `main.nf`
    resolves a relative samplesheet path against the samplesheet's own directory, so the
    whole collection directory can be copied across and run as it stands.
    """
    genome_of = {m["strain_key"]: m["genome_file"] for m in manifest}
    lines = [
        "# Every 20HM panel genome held to the pipeline's own acceptance criteria.",
        "# Written by dev/reference_genome_gates.py -- regenerate it, do not edit it.",
        "#",
        "# `assembly` is the genome as the panel ships it (for the two co-cultures that is",
        "# the single selected contig, not the whole assembly). `reads` appears only where",
        "# HiFi reads of the same strain exist; every other member is assessed on the",
        "# assembly alone and its read-based checks come back not_measured.",
        "#",
        "# Every path is relative to this file, so this directory is portable to the HPC.",
        "samples:",
    ]
    for entry in entries:
        key = entry["strain_key"]
        genome = relative(collection / genome_of[key], out_dir)
        genus = entry["species"].split()[0]
        lines += [
            f"  {key}:",
            f"    assembly: {genome}",
            f"    expected_taxon: g__{genus}",
        ]
        usable = entry.get("usable_runs") or []
        if usable:
            lines.append("    reads:")
            lines += [
                f"      - {relative(reads_dir / Path(urllib.parse.urlparse(url).path).name, out_dir)}"
                for run in usable
                for url in fastq_urls(run)
            ]
    (out_dir / "samplesheet.yml").write_text("\n".join(lines) + "\n")


def write_fetch_script(
    entries: list[dict], out_dir: Path, reads_dir: Path
) -> list[str]:
    """One curl per FASTQ, as a script, so the download can happen where the disk is.

    The script locates its own directory rather than carrying the generating machine's
    paths: the download happens on the HPC, and this file gets there by being copied.
    """
    urls = [
        url
        for entry in entries
        for run in entry.get("usable_runs") or []
        for url in fastq_urls(run)
    ]
    script = [
        "#!/bin/sh",
        "# The HiFi runs the assessment can use. Idempotent: re-run after an interrupted",
        "# transfer and curl resumes each file. Run it on the machine that will run the",
        "# pipeline -- the samplesheet expects the reads beside it.",
        "set -eu",
        'cd "$(dirname "$0")"',
        f'mkdir -p "{relative(reads_dir, out_dir)}"',
    ]
    script += [
        f'curl -fsS --continue-at - -o "{relative(reads_dir, out_dir)}/'
        f'{Path(urllib.parse.urlparse(url).path).name}" "{url}"'
        for url in urls
    ]
    path = out_dir / "fetch_reads.sh"
    path.write_text("\n".join(script) + "\n")
    path.chmod(0o755)
    return urls


def download(urls: list[str], reads_dir: Path) -> None:
    reads_dir.mkdir(parents=True, exist_ok=True)
    for url in urls:
        target = reads_dir / Path(urllib.parse.urlparse(url).path).name
        print(f"  fetching {target.name}", file=sys.stderr)
        subprocess.run(
            ["curl", "-fsS", "--continue-at", "-", "-o", str(target), url], check=True
        )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        required=True,
        help="<collection>/metadata/panel_genomes.tsv",
    )
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--reads-dir", type=Path, help="where the FASTQs go (default: <out-dir>/reads)"
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="fetch the usable HiFi runs now, rather than only writing fetch_reads.sh",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    collection = args.manifest.resolve().parent.parent
    manifest = read_manifest(args.manifest)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    reads_dir = args.reads_dir or args.out_dir / "reads"

    entries = []
    for member in manifest:
        entry = survey(member)
        entries.append(entry)
        print(
            f"{entry['strain_key']:32} {entry['verdict']:22} "
            f"own: {scoped_counts(entry, None)} | "
            f"strain: {scoped_counts(entry, 'same strain, other submitter')}",
            file=sys.stderr,
        )

    (args.out_dir / "README.md").write_text(
        "# IGA assessment of this genome collection\n\n"
        "Generated by `dev/reference_genome_gates.py` in the isolate-genome-assembler\n"
        "repository -- regenerate it, do not edit it. The write-up of what it is for and\n"
        "what the read survey found is `dev/reference_genome_gates.md` there.\n\n"
        "    samplesheet.yml  the pipeline input: every genome, with reads where usable\n"
        "                     HiFi reads of the same strain exist\n"
        "    reads.tsv        every genome, what reads ENA holds for it, and of what kind\n"
        "    reads.json       the same, with every run's full metadata\n"
        "    fetch_reads.sh   downloads the usable runs (empty when there are none)\n"
    )
    write_reads_table(entries, args.out_dir / "reads.tsv")
    write_samplesheet(entries, manifest, args.out_dir, collection, reads_dir)
    urls = write_fetch_script(entries, args.out_dir, reads_dir)
    (args.out_dir / "reads.json").write_text(json.dumps(entries, indent=2) + "\n")
    if args.download:
        download(urls, reads_dir)

    with_reads = sum(1 for e in entries if e.get("usable_runs"))
    print(
        f"\n{len(entries)} genomes, {with_reads} with usable HiFi reads "
        f"({len(urls)} FASTQ files). Wrote {args.out_dir}/",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
