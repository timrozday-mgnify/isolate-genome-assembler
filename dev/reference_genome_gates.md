# Are the acceptance criteria too strict? Holding reference genomes to them

Status: pipeline support implemented and the run set up, 2026-10-16. The run itself is
still to be executed; this file records the method, what reads exist, and what the read
survey already settles.

## The question

`assets/qc_thresholds.yml` decides whether an assembly passes. Its values are guesses from
[the implementation plan](../docs/implementation_plan.md), calibrated against nothing but
the pipeline's own output, and the 20HM run left every one of the eight isolates failing
`chromosome_not_circular` while three of them were single-contig chromosomes with intact
operons. A threshold that fails good assemblies is worse than no threshold.

The check on a threshold is to apply it to genomes that are known-good by an independent
route: the closed, published reference genomes of the same organisms. If a published
closed chromosome trips a gate, the gate is measuring this pipeline's conventions, not the
assembly's quality.

## What was needed in the pipeline

The samplesheet gained `assembly`, and `reads` became optional for a sample that has one.
See [the plan](../docs/implementation_plan.md#assessing-an-existing-genome-added-2026-10-16)
for the mechanics. One fix fell out of it immediately, before any genome was assessed:

**`chromosome_not_circular` was counting a missing header tag as a failure.** The check
reads the `circular=` tag the assemblers write. A submitted FASTA has no such tag, so
every closed reference genome would have failed the gate on a missing header field —
including the panel's own sources, which are complete chromosomes. `circular` is now
tri-state and an absent tag is `unknown`, which leaves the check `not_measured`. This also
explains the 20HM result above: the tag was there and said `false`, written by a tool
downstream of the assembler, so that gate was never about circularity on those isolates
either.

## Running it

```bash
python3 dev/reference_genome_gates.py \
    --manifest ~/Documents/mimicc/20hm-genome-collection/metadata/panel_genomes.tsv \
    --out-dir  ~/Documents/mimicc/20hm-genome-collection/iga-assessment
```

That writes `samplesheet.yml`, `reads.tsv`, `reads.json` and `fetch_reads.sh` into
`iga-assessment/`. Then, on a machine with the databases:

```bash
nextflow run main.nf -profile slurm,singularity -resume \
    --input  <collection>/iga-assessment/samplesheet.yml \
    --outdir <collection>/iga-assessment/results \
    --checkm2_db ... --bakta_db ... --busco_db ... --gtdbtk_db ... --ideel_db ...
```

The `--sylph_*`, `--human_reference` and `--plassembler_db` paths are not required for this
samplesheet: no sample in it has reads, so the stages that use them never run.

## What reads exist: none that can be used

Of the 24 panel members, **18 are public genomes and none of them has HiFi reads of its own
strain in ENA.** The other 6 are local (4 assembled by this pipeline, 2 from METT) and
their reads are not on this machine.

| verdict | n | members |
|---|---:|---|
| `own_reads_not_hifi` | 8 | its own sample's reads exist, but as Illumina, 454, Sanger, Nanopore or PacBio CLR |
| `strain_reads_not_hifi` | 6 | nothing from its own sample; other submitters' runs name the strain, none HiFi |
| `species_reads_only` | 2 | runs exist under the species, but none claims this strain |
| `no_reads` | 0 | |
| `local_assembly` | 6 | the 4 new long-read isolates and the 2 METT assemblies |

The per-member detail, with run accessions and read kinds split by scope, is in
`iga-assessment/reads.tsv`. `fetch_reads.sh` is empty of downloads, for the same reason.

### Two near-misses, and why they are no

Both are worth recording because both look like usable HiFi at a glance, and feeding
either in would have produced a page of confident, wrong failures.

- **`ERR14761318`** (*Phocaeicola vulgatus*, Sequel II, 85k reads at 9.7 kb, submitted as
  `*.hifi_reads.fastq.gz`) is genuinely HiFi and sits under taxid 435590, which is
  strain-level — `Phocaeicola vulgatus ATCC 8482`, the panel's own strain. It is not that
  strain: the sample title is "commensal bacteria isolated from fecal suckling pig".
  Submitters routinely file a new isolate under whatever strain taxid NCBI offers for the
  species, so a strain-level taxid is not a strain claim.
- **`SRR15171256`** really is `ATCC 29149`, the panel's *Mediterraneibacter gnavus*
  strain, and really is a Sequel II run — and its only file is
  `SRR15171256_subreads.fastq.gz`. Subreads are raw CLR at ~10% error. ENA has no field
  that distinguishes CLR from HiFi on a Sequel II, so the deposited file name is what
  decides it.

`read_kind()` in the script does both of these checks, and reports `pacbio_unclear` rather
than guessing when a Sequel II run's file names say neither.

## What this means

**The read-based half of the acceptance criteria cannot be calibrated against public
reference genomes.** `read_depth`, `unmapped_read_fraction`, `coverage_regions`,
`clipping_pileups`, `inspector_structural_errors`, `variants_high_af`, `merqury_qv`,
`plasmids_missing`, `possible_missing_replicons` and the whole contamination screen have no
independent yardstick here. They will come back `not_measured`, which is the honest answer
and not a finding.

**The assembly-side gates can be, and that is where the suspicion lies.** What the run will
produce for all 24: `checkm2_completeness`, `checkm2_contamination`, `busco_duplicated`,
`ideel_truncated_fraction`, `taxon_disagreement` and the replicon layout. Those six are the
gates to argue about, and `checkm2_contamination` especially — the 20HM run flagged
isolates 02 and 07 at 78.8% and 95.1% contamination when Autocycler had in fact assembled
each organism of a co-culture into its own complete chromosome.

Three follow-ups this does not cover:

1. **Calibrating the read-based gates needs HiFi data this pipeline did not produce.** The
   nearest source is a public HiFi assembly *with* its reads — any organism will do, since
   the question is about thresholds and not these strains. Worth a separate search.
2. **The 8 IGA isolates' own HiFi reads are the best available control** and are not on
   this machine. Re-running them with `assembly` pointed at the delivered FASTA and `reads`
   at the originals would exercise every gate, on data where the answer is known.
3. **`PENDING_ENA`.** Once the four new assemblies are submitted, their reads go with them
   and this survey answers differently.

## Non-HiFi reads: available, but they answer a different question

Asked after the first survey: if no HiFi exists, what about the Illumina, 454, Nanopore and
PacBio CLR reads that do, and can the checks be run on those?

### What is actually there

Two corrections were needed before the numbers meant anything. First, **60% of the runs
under these isolates' taxids are not whole-genome**: of 2530 runs found, 890 are AMPLICON
and 551 RNA-seq, with Tn-seq, ChIP-seq and clone libraries behind them. Counting those as
"reads for this genome" inflated *Bacteroides fragilis* from 1039x to 12097x and
*Phocaeicola vulgatus* from 7x to 4247x. `runs_for()` now keeps only `WGS`, `WGA` and
`OTHER`. Second, only runs from the genome's own sample or from a run whose metadata names
the strain count at all.

What survives, per read type, over the 18 public genomes:

| read type | genomes | typical depth | caveat |
|---|---:|---|---|
| Illumina / 454 | 13 | 33x–1670x | 3 of them are 454, which is 454-era error profile |
| Nanopore | 6 | 51x–1071x | all MinION, 2018–2021 accessions, so R9-era |
| PacBio CLR | 5 | 0x–838x | RS, RS II and one Sequel II subreads set |
| Sanger | 4 | 0x–10x | unusable depth |
| nothing attributable | 3 | | `clostridium_saccharolyticum`, `methanobrevibacter_smithii`, `salinibacter_ruber` |

So the reads exist for 15 of 18, and for 13 of them at serious Illumina depth.

### Which checks would survive the platform change

Three tiers, and the tiers are about read accuracy and length, not about "non-HiFi":

**Transfers, and is good evidence.** `merqury_qv` — Merqury was built for Illumina k-mers
against an assembly, and on a 400x NovaSeq set it is a more trustworthy QV than HiFi gives.
`variants_high_af` — bcftools pileup on Illumina is the gold standard for SNVs. The
contamination screen too: sylph takes Illumina natively.

**Transfers only after pipeline work, and degrades in the places that matter here.**
`unmapped_read_fraction`, `coverage_regions` and `human_read_fraction` need the minimap2
preset switched per platform (`sr`, `map-ont`, `map-pb` instead of `map-hifi`) — a
`read_platform` samplesheet key threaded through three modules. The semantics survive, but
short reads cannot map uniquely across rRNA operons and IS elements, which for a 16S
reference panel is precisely the region under scrutiny: coverage anomalies there would be
mapping artefacts, not assembly defects.

**Does not transfer.** `clipping_pileups` is vacuous on 150 bp reads — a short read cannot
carry the kilobase soft clip the check looks for. `inspector_structural_errors` is a
long-read assembly evaluator and will not take Illumina. `possible_missing_replicons` runs
Flye, and `plasmids_missing` runs Plassembler; both are long-read tools.
`genome_size_disagreement` compares two estimates, one of which is `autocycler genome size`
running Raven on long reads — with Illumina only there is no second estimate to disagree
with.

For the CLR and R9-era ONT sets specifically, at 5–15% error: k-mer QV collapses because
most k-mers are erroneous, and pileup variant calls are unreliable. Those sets are worse
than absent — they would fail every accuracy gate for reasons that have nothing to do with
the assembly.

### Why this does not answer the threshold question

Even where a check runs cleanly on Illumina, **a threshold calibrated on Illumina reads is
not the threshold for HiFi reads.** `merqury_qv: {warn: 50, fail: 40}` asks what QV a good
assembly reaches *when measured with HiFi k-mers*. Measuring the same assembly with
Illumina k-mers gives a different number for the same quality, so the result cannot be
carried across. The same holds for `read_depth`, where 50x means something different per
platform, and for every coverage-anomaly count.

That is a real finding rather than a dead end, because it splits the work in two:

1. **Validating the panel** — "are these assemblies actually good?" — can use the Illumina
   reads, and `merqury_qv` plus `variants_high_af` on the 13 genomes that have them would
   be strong independent evidence. Worth doing on its own terms.
2. **Calibrating the HiFi thresholds** needs HiFi reads paired with a known-good assembly,
   and the organism is irrelevant to that question.

### The calibration set does exist

ENA holds complete bacterial genomes assembled from Revio HiFi reads with the reads
deposited alongside. Of 181 distinct samples with a Revio HiFi WGS run, 60 were checked and
**22 have an ENA assembly**, most at `complete genome` level — *Streptomyces* enrichment
cultures, *Pelagimonas varians*, *Sulfitobacter* sp. SH22, *Paenibacillus* sp. SEL1 at
16 Gbp of HiFi. These are same-strain, same-platform, independently finished genomes: the
dataset the threshold question needs, and the pipeline change in this commit is already
enough to run them (`assembly` plus `reads`, no platform handling required).

Searching that set properly, rather than the 20HM organisms, is the next step for the
threshold question. The 20HM run stays worth doing for the six assembly-side gates.

## Next: assess these genomes with the reads that do exist

The read-mapping break evidence is the part worth having even without HiFi, and it needs
platform-aware checks plus one new module for the short-read signal the pipeline does not
measure at all. That is specified in
[the multi-platform read checks plan](../docs/multiplatform_read_checks_plan.md), which also
records what is fetchable per platform and that both the download and the run happen on the
HPC.
