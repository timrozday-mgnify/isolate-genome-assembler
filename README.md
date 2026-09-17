# isolate-genome-assembler

Nextflow (DSL2) pipeline for complete, checked prokaryotic isolate genomes from PacBio HiFi
reads: read QC, sylph contamination screen (GTDB + human), Autocycler consensus of Flye,
hifiasm, Raven, Canu, miniasm, metaMDBG and Plassembler, finishing, assembly checks with
pass/warn/fail gates, and a Quarto report. Built for SLURM + Singularity/Apptainer.

**Status: Phase 0 complete — repository skeleton.** The current implementation validates the YAML
samplesheet and deliberately performs no analysis yet. See [the implementation plan](docs/implementation_plan.md).

## Quick start

Create a YAML samplesheet using [the example](assets/samplesheet.example.yml), then validate
it before submitting work to a cluster:

```bash
nextflow run main.nf --input samples.yml -preview
```

Run the empty skeleton check with:

```bash
nextflow run main.nf -preview
nf-test test --tag stub
pytest
```

The eventual pipeline is HiFi-only: it will not accept Illumina reads and will not include
short-read or Medaka polishing. Phase 0 provides the SLURM, Docker, Singularity, and
Apptainer profile scaffolding; all process output publication remains opt-in.
