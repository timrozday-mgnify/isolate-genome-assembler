# isolate-genome-assembler

Nextflow (DSL2) pipeline for complete, checked prokaryotic isolate genomes from PacBio HiFi
reads: read QC, sylph contamination screen (GTDB + human), Autocycler consensus of Flye,
hifiasm, Raven, Canu, miniasm, metaMDBG and Plassembler, finishing, assembly checks with
pass/warn/fail gates, and a Quarto report. Built for SLURM + Singularity/Apptainer.

**Status: Phase 5 complete — read QC, contamination screen, assembly, consensus,
finishing, assembly checks, pass/warn/fail gates and the report.** Benchmarked defaults
(Phase 6) are still to come. See [the implementation plan](docs/implementation_plan.md).

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
    --plassembler_db /shared/databases/plassembler_db \
    --checkm2_db /shared/databases/checkm2_db/uniref100.KO.1.dmnd \
    --bakta_db /shared/databases/db-light \
    --busco_db /shared/databases/busco_downloads \
    --gtdbtk_db /shared/databases/gtdbtk_db \
    --ideel_db /shared/databases/ideel_db.dmnd
```

Always use `-resume`, and run the head job inside `sbatch` or `tmux`.

## Databases

`--prepare_databases` downloads the sylph GTDB r226 database and its sylph-tax metadata,
downloads CHM13 v2.0 and GRCh38, sketches those two into the human sylph database, runs
`plassembler download`, fetches the CheckM2 database, the Bakta light database, the BUSCO
`--busco_lineage` dataset and the GTDB-Tk r232 package, and builds a
DIAMOND database of UniProt Swiss-Prot for the IDEEL test. It writes everything to
`--database_dir`. The ten paths above are then required whenever `--input` is given; the
run stops before any task if one is missing.
A database whose param is already set is skipped by `--prepare_databases`, so pass the
paths you have to build only the rest.

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
| `--assemblers` | seven tools | Comma list; also accepts `myloasm` and `lja` |
| `--full_read_assemblers` | `--assemblers` minus Plassembler | Assemblers also run on the whole read set, scored against each other and the consensus |
| `--subsample_count` | `4` | Independent subsampled read sets per sample |
| `--subsample_min_depth` | `25` | Minimum depth a subset may have |
| `--publish_outputs` | `true` | Copy process outputs to `--outdir`; `false` leaves them in `work/` |
| `--publish_input_assemblies` | `false` | Also publish the per-subset input assemblies |
| `--min_contig_len` | `1000` | Contigs shorter than this are flagged |
| `--chromosome_min_len` | `1000000` | A contig this long is called a chromosome |
| `--drop_flagged_contigs` | `false` | Move flagged contigs out of the final FASTA |
| `--qc_thresholds` | `assets/qc_thresholds.yml` | Pass/warn/fail thresholds for every check |
| `--coverage_low_ratio`, `--coverage_high_ratio` | `0.5`, `2.0` | Depth-window flags, relative to the contig median |
| `--clip_min_reads`, `--clip_min_fraction`, `--clip_min_length` | `5`, `0.2`, `500` | What counts as a clipping pile-up |
| `--circularise_min_identity` | `0.95` | Identity a contig's start-end self-match must reach before the duplicated span is cut |
| `--score_size_tolerance`, `--score_min_qv` | `0.1`, `40` | How far a candidate assembly's length may sit from the genome size estimate, and the QV below which it is not a candidate |
| `--assembly_selection` | `always_score` | `always_score` finishes the highest-scoring candidate, `score` keeps a resolved consensus, `flye` is the pre-scoring behaviour |
| `--clip_extend_args` | `-z 5000,5000 --end-bonus 100` | minimap2 options for the permissive second alignment that separates a misjoin from a local difference |
| `--variant_min_depth` | `10` | Minimum depth for the allele-frequency scan |
| `--homopolymer_min_len` | `8` | Homopolymer indels at least this long are tallied separately |
| `--ideel_min_ratio` | `0.9` | A protein below this fraction of its best hit is truncated |
| `--busco_lineage` | `bacteria_odb12` | BUSCO dataset, read offline from `--busco_db` |
| `--per_sample_reports` | `true` | Also render `<id>.report.html` for each sample |

## Pipeline flow

```
samplesheet (YAML): id, reads, optional reference / autocycler_dir
        |
        v
[1] READ QC ................ PREPARE_READS -> HiFiAdapterFilt -> seqkit stats, NanoPlot,
                             GC + duplicate profile, KMC -> GenomeScope2, Autocycler
                             genome size  =>  read_qc.tsv, genome size estimate
        |
        v
[2] CONTAMINATION .......... sylph profile (GTDB) + sylph-tax -> taxonomy
                             sylph query (human db), minimap2 vs CHM13 -> human fraction
                             =>  contamination_summary.json
        |
        v
[3] REMOVE_HUMAN ........... drops human reads when --remove_human says so (always runs;
                             the decision file, not the DAG, decides)
        |
        +-------------------------------+
        |                               |
        v                               v
[4] ASSEMBLY                        full read set
    Autocycler subsample            each --full_read_assemblers tool assembles all reads
    -> N subsets                    (Flye always, whether or not it is selected)
    -> every --assemblers tool             |
       assembles every subset              |
    -> normalise headers                   |
    -> Autocycler cluster/trim/            |
       resolve/combine -> consensus        |
    -> Plassembler (plasmids)              |
        |                                  |
        +----------------+-----------------+
                         v
[5] SCORING ................ every candidate (each full assembly + the consensus):
                             contig ends -> circularise, seqkit stats, meryl+Merqury QV,
                             minimap2 -> clipping pile-ups  =>  full_assemblies.tsv
                             SELECT_ASSEMBLY picks the winner per --assembly_selection
                         |
                         v
[6] FINISHING .............. contig ends -> dnaapler rotate -> classify replicons
                             (<id>_chromosome, <id>_plasmid_N); flags, never drops
                             Plassembler on full reads + skani  =>  plasmid_audit.tsv
                         |
                         v
[7] CHECKS ................. minimap2/samtools mapping -> mosdepth windows, clipping
                             pile-ups, Flye on unmapped reads; Inspector; bcftools
                             pileup -> variant scan; meryl + Merqury QV; Bakta ->
                             DIAMOND vs Swiss-Prot (IDEEL); CheckM2; BUSCO; GTDB-Tk;
                             skani vs --reference; Bandage + assembler contribution
                         |
                         v
[8] GATES + REPORT ......... qc_gates.py vs qc_thresholds.yml -> <id>.qc.json
                             collect_metrics.py -> report/run_summary/*.tsv
                             Quarto -> isolate_assembly_report.html (+ per-sample)
```

`--prepare_databases` is a separate entry point that runs none of the above: it only
downloads and builds the databases listed under [Databases](#databases).

## Outputs

`reads/<id>/` holds `read_qc.tsv` and the per-tool measurements
behind it (seqkit stats, NanoPlot, GC histogram, duplicates, adapters, GenomeScope2), and
`contamination/<id>/` holds the sylph profile and taxonomy, the human query, the human read
fraction and `contamination_summary.json`. `assemblies/<id>/` holds the selected assembly
under `final/`, Autocycler's own `autocycler_out/` directory, the cluster dotplots,
`assembly_attempts.tsv` and `assembly_source.tsv`, with the per-subset input assemblies
behind `--publish_input_assemblies`. `assemblies/<id>/full/` holds each assembler's
whole-read-set assembly and `assemblies/<id>/scoring/` holds `full_assemblies.tsv`, which
ranks them against each other and against the consensus, plus the per-candidate
circularity tables. `checks/<id>/` holds every check's table (below),
`annotation/<id>/` the Bakta annotation, and `qc/<id>.qc.json` the sample's verdict.
`report/` holds the report (below). `pipeline_info/` always holds Nextflow's trace,
timeline, report and DAG, and `software_versions.tsv`.

## Report

`bin/collect_metrics.py` gathers every sample's measurements, verdicts and images into
one `report/run_summary/` bundle of tidy TSVs, and the Quarto template in
[`assets/report/`](assets/report/) reads only that bundle. `report/isolate_assembly_report.html`
is a single self-contained file (safe to email): run overview and parameters, a status
board, the samples × checks gate matrix with every warning and failure listed, a tabbed
section per sample (Reads, Contamination, Assembly, Accuracy, Genes, Taxonomy, Plasmids),
generated methods text with citations, and software versions. `report/samples/<id>.report.html`
is the same report for one sample.

The report renders in `ghcr.io/timrozday-mgnify/isolate-genome-assembler-report`
([`containers/report/Dockerfile`](containers/report/Dockerfile)): Quarto, papermill,
pandas and plotly. [`build-images.yml`](.github/workflows/build-images.yml) publishes it
from `main`. To re-render by hand, e.g. after editing the template:

```bash
quarto render assets/report/isolate_assembly_report.qmd -P summary_dir:$PWD/results/report/run_summary
```

## Assembly

Each sample's reads are split into `--subsample_count` independent subsets, and every
enabled assembler assembles every subset in its own pinned container, with the flags
Autocycler v0.7.0's `src/helper.rs` uses. Headers are then normalised to
`<assembler>_<subset>_<n>`, keeping each tool's circularity and depth tags, so a consensus
contig can be traced back to the assemblies that voted for it.

An assembler that crashes on one subset is retried and then given up on: the other 27 input
assemblies still make a consensus, and the failure is recorded in `assembly_attempts.tsv`.
If Autocycler's consensus is not fully resolved, the sample falls back to a Flye assembly of
the full read set (`assembly_source=fallback_flye`), which runs for every sample anyway.

Re-clustering by hand needs no reassembly: point a sample's `autocycler_dir` key at a
previous run's `autocycler_out/` and pass `autocycler_cluster_args` (for example
`--manual 12,34`), and only clustering onwards is repeated.

## Finishing

The selected assembly is checked for leftover circular end-overlap, rotated with
`dnaapler all` (chromosome at *dnaA*, plasmids at *repA*, phages at *terL*), then
classified and renamed to `<id>_chromosome` and `<id>_plasmid_1..n`.

Nothing is dropped silently. Linear contigs, contigs under `--min_contig_len`, contigs
below `--min_contig_depth_ratio` of chromosome depth, contigs that could not be rotated and
contigs whose ends still overlap are all **flagged** in `<id>.contigs.tsv` and kept.
`--drop_flagged_contigs` is what moves them to `<id>.removed_contigs.fasta`.

Plassembler then runs a second time on the full read set, independently of the consensus.
Each plasmid it reports is matched against the finished contigs with skani, and
`plasmid_audit.tsv` marks it recovered or missing. Treat a missing plasmid as a prompt to
look, not a verdict — and note that HiFi library prep under-represents plasmids below
~10 kb in the first place.

## Checks and gates

The finished assembly is checked against the reads and against reference databases. Each
check writes a small table under `checks/<id>/`; none of them decides pass or fail itself.

| Check | Tool | Output |
|---|---|---|
| Unmapped reads | minimap2 `map-hifi` | `mapping.tsv`; the unmapped reads are assembled with Flye, and a circular contig among them is a possible missing replicon |
| Uneven depth | mosdepth, 1 kb windows | `coverage_regions.tsv` (low/high regions), `contig_depth.tsv` (depth relative to the chromosome) |
| Candidate misjoins | minimap2, pysam | `clipping.tsv`: places where many reads are clipped at once. Each is re-measured in a second, permissive alignment (`--clip_extend_args`) that extends through what it can: `verdict` is `confirmed` when the pile-up survives that, `resolved` when the reads carry on matching and the clip marked a local difference. `tail_target` is where the clipped tails align instead, when they agree on a place: the join the assembly should have made. Tails that land nowhere are sequence missing from the assembly. Only `confirmed` pile-ups reach the QC gate |
| Structural and small errors | Inspector | `summary_statistics`, error BEDs, QV |
| Per-base accuracy | bcftools mpileup | `variants.tsv`: AF ≥ 0.5 (likely error) and 0.2–0.5 (mixed strain or collapsed repeat), homopolymer indels apart |
| K-mer QV | meryl + Merqury | QV, completeness, spectra-cn plot |
| Frameshifts | Bakta + DIAMOND vs Swiss-Prot | `ideel.tsv`: protein/hit length ratios; `rrna_depth.tsv` |
| Completeness | CheckM2, BUSCO | completeness, contamination, duplicated BUSCOs |
| Identity | GTDB-Tk, skani vs the samplesheet `reference` | classification, ANI |
| Consensus | Autocycler, Bandage | `assembler_contribution.tsv`, graph image |

`bin/qc_gates.py` then compares every measurement with
[`assets/qc_thresholds.yml`](assets/qc_thresholds.yml) and writes `qc/<id>.qc.json`: one
`{value, threshold, status, message}` entry per check and an overall status, the worst
entry. A check that could not be measured is `not_measured` and does not change the overall
status. To change a threshold, copy the file, edit it, and pass `--qc_thresholds`.

## Tools

Every tool runs in a pinned container (Bioconda/Galaxy or Seqera Wave), listed here with
the version the pipeline pins and what it is used for. The exact pins live in each module;
`pipeline_info/software_versions.tsv` records what actually ran.

### Reads

| Tool | Version | Purpose |
|---|---|---|
| HiFiAdapterFilt | 3.0.0 | Find and remove reads carrying PacBio adapter sequence |
| seqkit | 2.13.0 | Read and assembly length/count/N50 statistics |
| NanoPlot | 1.47.0 | Read length and quality distributions, plots for the report |
| KMC | 3.2.4 | K-mer count histogram (`--kmer_size`) feeding GenomeScope2 |
| GenomeScope2 | 2.1.0 | Genome size, heterozygosity and repeat estimate from the k-mer spectrum |
| Autocycler (`subsample`, `helper genome_size`) | 0.7.0 | Independent read subsets and a second genome size estimate |

### Contamination

| Tool | Version | Purpose |
|---|---|---|
| sylph | 0.9.0 | Profile reads against GTDB, and query them against the human sketch |
| sylph-tax | 1.9.1 | Turn sylph's genome hits into a taxonomic profile |
| minimap2 | 2.30 | Align reads to CHM13 to measure (and select) human reads |

### Assembly

| Tool | Version | Purpose |
|---|---|---|
| Flye | 2.9.6 | HiFi assembler; also the fallback assembly and the unmapped-read assembly |
| hifiasm | 0.25.0 | HiFi assembler |
| Raven | 1.8.3 | Long-read assembler |
| Canu | 2.3 | Long-read assembler |
| miniasm | 0.3 | Long-read assembler (overlaps from minimap2) |
| minipolish | 0.2.1 | Polish the miniasm graph, which is otherwise unpolished |
| metaMDBG | 1.4 | Metagenome-oriented HiFi assembler, used here for its contiguity |
| myloasm | 0.7.0 | Optional assembler (`--assemblers myloasm`) |
| LJA | 0.2 | Optional assembler (`--assemblers lja`) |
| Plassembler | 1.8.5 | Dedicated plasmid assembly, on subsets and again on the full read set |
| Autocycler (`cluster`/`trim`/`resolve`/`combine`) | 0.7.0 | Consensus assembly across all input assemblies |

### Selection and finishing

| Tool | Version | Purpose |
|---|---|---|
| minimap2 | 2.30 | Self-alignment for circular end-overlap, and read alignment for scoring |
| samtools | 1.24 | Sort, index and stream the alignments the checks read |
| meryl | 1.4.1 | K-mer database of the reads, shared by every candidate's QV |
| Merqury | 1.3 | Reference-free QV and k-mer completeness per candidate |
| dnaapler | 1.4.0 | Rotate circular contigs to *dnaA* / *repA* / *terL* |
| skani | 0.2.2 | ANI: plasmid audit matches, and comparison to the samplesheet `reference` |

### Checks

| Tool | Version | Purpose |
|---|---|---|
| mosdepth | (Wave `htslib_mosdepth_gzip`) | Per-window depth for the uneven-coverage check |
| pysam | 0.24.1 | Read the BAM for clipping pile-ups and depth-window tables |
| Inspector | 1.3.1 | Structural and small-scale assembly errors, plus its own QV |
| bcftools | 1.23.1 | Pileup and call for the allele-frequency and homopolymer scan |
| Bakta | (Wave `bakta_diamond`) | Annotation: CDS, rRNA and tRNA for the gene checks |
| DIAMOND | 2.2.1 | `blastp` of Bakta proteins vs Swiss-Prot for the IDEEL frameshift test |
| CheckM2 | 1.1.0 | Completeness and contamination |
| BUSCO | 6.1.0 | Single-copy marker completeness and duplication (`--busco_lineage`) |
| GTDB-Tk | 2.7.2 | Taxonomic classification of the finished assembly |
| Bandage | 0.9.0 | Image of the Autocycler consensus graph |

### Infrastructure

| Tool | Version | Purpose |
|---|---|---|
| Nextflow (DSL2) | ≥ 25.0.0 | Workflow engine; SLURM + Singularity/Apptainer profiles |
| Python | 3.12 | Every script in [`bin/`](bin/) and the small local modules |
| Quarto, papermill, pandas, plotly | report image | Render the self-contained HTML report |
| GNU wget | 1.18 | Database downloads under `--prepare_databases` |
| nf-test, pytest, ruff, pre-commit | — | Stub tests, Python tests, lint and formatting |

## Testing

```bash
nextflow run main.nf -preview
nf-test test --tag stub --profile docker
pytest
```

The stub tests need a container engine: several modules report their version with an `eval`
output, which runs the tool even in a stub task. Until the report image is published, build
it locally first:

```bash
docker build -t ghcr.io/timrozday-mgnify/isolate-genome-assembler-report:latest containers/report
```

`tests/test_report.py` renders the template against the fixture bundle
`tests/data/run_summary/` (one pass, one warn, one fail sample) and skips outside the report
image. The fixture is generated by `tests/data/generate_run_summary.py` through
`collect_metrics.py`; rerun it after changing either script, or `tests/test_collect_metrics.py`
fails.

The pipeline is HiFi-only: it will not accept Illumina reads and will not include short-read
or Medaka polishing.
