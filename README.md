# isolate-genome-assembler

Nextflow (DSL2) pipeline for complete, checked prokaryotic isolate genomes from PacBio HiFi
reads: read QC, sylph contamination screen (GTDB + human), Autocycler consensus of Flye,
hifiasm, Raven, Canu, miniasm, metaMDBG and Plassembler, finishing, assembly checks with
pass/warn/fail gates, and a Quarto report. Built for SLURM + Singularity/Apptainer.

**Status: Phase 1 complete — read QC and contamination screen.** Assembly and everything
after it are not implemented yet. See [the implementation plan](docs/implementation_plan.md).

## Quick start

Build the databases once per site:

```bash
nextflow run main.nf --prepare_databases --database_dir /shared/databases -profile singularity
```

Create a YAML samplesheet using [the example](assets/samplesheet.example.yml), then run:

```bash
nextflow run main.nf --input samples.yml --outdir results -profile slurm,singularity \
    --sylph_gtdb_db /shared/databases/gtdb-r226-c200-dbv1.syldb \
    --sylph_gtdb_taxonomy /shared/databases/gtdb_r226_metadata.tsv.gz \
    --sylph_human_db /shared/databases/sylph_human_db.syldb \
    --human_reference /shared/databases/chm13v2.0.fa.gz \
    --publish_outputs
```

Always use `-resume`, and run the head job inside `sbatch` or `tmux`.

## Databases

`--prepare_databases` downloads the sylph GTDB r226 database and its sylph-tax metadata,
downloads CHM13 v2.0 and GRCh38, and sketches those two into the human sylph database. It
writes everything to `--database_dir`. The four paths above are then required whenever
`--input` is given; the run stops before any task if one is missing.

## Parameters

Every parameter is declared with an explanatory comment in
[`nextflow.config`](nextflow.config). The ones that change results today:

| Parameter | Default | Effect |
|---|---|---|
| `--min_read_qv` | `20` | Reads below this mean accuracy are dropped and counted |
| `--remove_adapter_reads` | `true` | Use the HiFiAdapterFilt-filtered read set downstream |
| `--kmer_size` | `21` | K-mer length for the KMC + GenomeScope2 genome size estimate |
| `--remove_human` | `auto` | `auto` removes human reads only when human is detected |
| `--contam_max_human` | `0.001` | Human read fraction at which `auto` removes |
| `--publish_outputs` | `false` | Publication is opt-in; nothing is copied by default |

## Outputs

With `--publish_outputs`, `reads/<id>/` holds `read_qc.tsv` and the per-tool measurements
behind it (seqkit stats, NanoPlot, GC histogram, duplicates, adapters, GenomeScope2), and
`contamination/<id>/` holds the sylph profile and taxonomy, the human query, the human read
fraction and `contamination_summary.json`. Thresholds are not applied yet: the gates and the
report arrive in Phases 4 and 5.

## Testing

```bash
nextflow run main.nf -preview
nf-test test --tag stub --profile docker
pytest
```

The stub tests need a container engine: several modules report their version with an `eval`
output, which runs the tool even in a stub task.

The pipeline is HiFi-only: it will not accept Illumina reads and will not include short-read
or Medaka polishing.
