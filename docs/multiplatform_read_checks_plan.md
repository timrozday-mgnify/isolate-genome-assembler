# Multi-platform read checks: finding breaks in finished genomes

Status: **plan only, nothing implemented.** Written 2026-10-16. The `assembly` samplesheet
key and the assembly-only assessment path it builds on are implemented and tested; see
[the implementation plan](implementation_plan.md#assessing-an-existing-genome-added-2026-10-16).

Read [the read survey](../dev/reference_genome_gates.md) first: it establishes that no 20HM
panel genome has HiFi reads, which is why this plan exists.

## Goal

Assess a finished genome against whatever reads exist for it, with **read-mapping evidence
of structural breaks as the primary output**, and the other checks as secondary. The
pipeline is HiFi-only today; every read-based check either hardcodes a HiFi preset or
silently misbehaves on anything else.

Not in scope: reassembling these genomes, adding a short-read SV caller (delly, manta,
lumpy), hybrid samples that mix platforms in one assembly check, and read subsampling.

## Where it runs

**The pipeline and the download both run on the HPC. This plan assumes no local copy of the
reads and no machine-specific paths anywhere.** Concretely:

- The fetcher is a script the user runs on an HPC login node. Python 3 standard library
  only, no new conda environment, `curl` for transfers.
- Every path in the generated samplesheet is **relative to the samplesheet's own
  directory** (`../genomes/x.fa.gz`, `reads/SRR123_1.fastq.gz`). `main.nf` already resolves
  relative samplesheet paths that way, so the whole collection directory stays portable:
  copy it to the HPC, fetch there, run there.
- The fetcher must be resumable (`curl --continue-at -`), verify `fastq_md5` from ENA, and
  be safe to re-run — a login node job gets killed and restarted as a matter of course.
- If the login node has no outbound access, the fallback is documented rather than
  automated: fetch elsewhere, `rsync` the `reads/` directory across. ~15 GB at full depth.

## What a "break" looks like, per platform

This is the table the rest of the plan follows. A break is a position where the reads
disagree with the assembly's joins, and the evidence available depends entirely on read
length and accuracy.

| evidence | HiFi | ONT (R9) | PacBio CLR | Illumina PE | 454 SE |
|---|---|---|---|---|---|
| soft-clip pile-up | strong | strong | noisy | weak | weak |
| clipped-tail destination (supplementary alignments) | strong | strong | noisy | none | none |
| coverage dropout / spike | strong | strong | strong | strong | strong |
| **discordant pairs** (mate elsewhere, wrong orientation, insert outliers) | n/a | n/a | n/a | **strong** | n/a |
| Inspector structural errors | yes (`-d hifi`) | yes (`-d nanopore`) | yes (`-d clr`) | not supported | not supported |
| unmapped reads reassembled | yes (Flye) | yes (Flye) | yes (Flye) | no Flye mode | no |

Two conclusions drive the design:

1. **For long reads the checks already exist** and only need their presets and thresholds
   changed. `clipping_pileups.py` is already the right shape — its `--min-clip` is a
   parameter, not a constant.
2. **For Illumina the primary break signal is missing entirely.** A 150 bp read cannot
   carry the 500 bp soft clip `clip_min_length` looks for. What a short-read library gives
   instead is the read pair: a misjoin shows up as pairs whose mates land on another
   contig, in the wrong orientation, or at an impossible insert distance. Nothing in the
   pipeline measures that. It is the one genuinely new module this plan adds.

**Cross-platform corroboration is the strongest evidence of all**, and 8 of the 13 panel
genomes with reads have two or more platforms. A position flagged by ONT clip pile-ups
*and* Illumina discordant pairs is a break; one flagged by only the noisy CLR set is a
candidate. This falls out for free from design decision 2 below.

## Design decisions

### 1. One platform table, read by every module

`read_platform` becomes a samplesheet key with values `hifi | ont | clr | illumina | 454`.
Explicit, never inferred: read-length heuristics cannot tell CLR from HiFi, which is the
distinction that matters most, and the fetcher knows the answer from ENA anyway.

The temptation is a `platform == 'x'` conditional in each of ten modules. Instead **one
map in `conf/platforms.config`**, keyed by platform, giving every per-platform value:

| field | hifi | ont | clr | illumina | 454 |
|---|---|---|---|---|---|
| `minimap2_preset` | `map-hifi` | `map-ont` | `map-pb` | `sr` | `sr` |
| `inspector_datatype` | `hifi` | `nanopore` | `clr` | — | — |
| `flye_mode` | `--pacbio-hifi` | `--nano-raw` | `--pacbio-raw` | — | — |
| `clip_min_length` | 500 | 300 | 300 | 20 | 40 |
| `paired` | no | no | no | yes | no |
| `min_read_qv` | 20 | 7 | 0 | 20 | 15 |
| `kmer_accurate` (meryl/Merqury, GenomeScope2) | yes | no | no | yes | no |
| `long_read` (Plassembler, Autocycler genome size, Flye) | yes | yes | yes | no | no |
| `adapter_filter` (HiFiAdapterFilt) | yes | no | no | no | no |

Modules read `meta.platform_<field>`, populated once in `main.nf` from this map when the
samplesheet is parsed. A new platform is one row, not a sweep through the modules.

`min_read_qv` is the trap to call out: the current global default of 20 would **delete
essentially every ONT R9 and CLR read** (mean QV ~10 and ~8) and leave an empty FASTQ that
every downstream check would then report as `not_measured`. It has to come from the table.

### 2. One platform per sample; a genome with two platforms is two samples

No hybrid handling. A genome with ONT and Illumina becomes
`roseburia_intestinalis.ont` and `roseburia_intestinalis.illumina`, both pointing `assembly`
at the same FASTA. This needs no code: the samplesheet already allows it, and the ids
differ. The report then carries one row per (genome, platform), which is exactly the
comparison the cross-platform corroboration above wants.

The alternative — one sample carrying several read sets — multiplies every channel in
`CHECKS` by the number of platforms, to produce the same numbers in a worse shape.

### 3. New module: `PAIR_ANOMALIES` (+ `bin/pair_anomalies.py`)

Illumina only (`when: meta.platform_paired`). Reads the same BAM `MAP_READS` already
produces, and writes one row per candidate position, binned as `clipping_pileups.py` bins:

- mate mapped to a different contig
- mate orientation not FR
- insert size beyond the run's own 99.9th percentile, estimated from the sample itself
  rather than assumed
- mate unmapped while this read maps

Each as a count and as a fraction of local properly-paired depth, so the measure is
depth-independent. Output shape mirrors `clipping.tsv` so the report can draw both tracks
against one coordinate axis. Gate: `discordant_pair_pileups`, thresholded like
`clipping_pileups`.

Deliberately **not** a structural-variant caller. The question is "does the assembly join
things the reads say it should not", and a per-bin anomaly fraction answers it. Adding
manta or delly means a new container, a new database of expectations, and VCF parsing, to
reach the same table.

### 4. Candidate breaks are merged in the report, not in a gate

Each method keeps its own measurement and its own gate. The report gains a **breaks table**:
one row per candidate position, columns for which methods support it and at what strength,
across every platform for that genome. No new combined gate and no scoring function — the
arithmetic of "two weak signals plus one strong" has no principled answer, and the table
lets a human see the agreement directly.

### 5. Provenance of the reads travels with the sample

7 of the 13 best available runs are **not from the genome's own sample** — they are another
submitter's culture of the same strain. A real passage difference, a lost plasmid or a
genuine rearrangement in their isolate would show up as a break in this genome. So
`read_provenance` (`own_sample | same_strain`) becomes a samplesheet key, carried to the
report and shown beside every break. Without it the output invites exactly the wrong
conclusion.

### 6. Per-platform thresholds

A clip pile-up fraction that means a misjoin in HiFi means nothing in CLR. `assets/qc_thresholds.yml`
gains an optional `platforms:` section whose per-platform entries are merged over the base
values; `qc_gates.py` takes `--platform` and does the merge. The base file keeps today's
HiFi numbers, so nothing changes for an existing run.

Starting values are guesses again, and the first 20HM run is what calibrates them. Checks a
platform cannot measure are simply absent from its channels and come back `not_measured`,
which is the existing behaviour and needs no new mechanism.

## What each check does per platform

| check | hifi | ont | clr | illumina | 454 |
|---|---|---|---|---|---|
| `read_depth` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `unmapped_read_fraction` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `coverage_regions` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `clipping_pileups` | ✓ | ✓ | ✓ (noisy) | weak | weak |
| `discordant_pair_pileups` (new) | — | — | — | ✓ | — |
| `inspector_structural_errors` | ✓ | ✓ | ✓ | — | — |
| `variants_high_af` | ✓ | ✓ (noisy) | — | ✓ | ✓ (homopolymers) |
| `merqury_qv` | ✓ | — | — | ✓ | — |
| `genome_size_disagreement` | ✓ | ✓ | ✓ | — | — |
| `possible_missing_replicons` | ✓ | ✓ | ✓ | — | — |
| `plasmids_missing` | ✓ | ✓ | ✓ | — | — |
| contamination screen (sylph) | ✓ | ✓ | ✓ | ✓ | ✓ |
| assembly-only checks (CheckM2, BUSCO, GTDB-Tk, IDEEL) | ✓ | ✓ | ✓ | ✓ | ✓ |

Illumina loses `possible_missing_replicons` because Flye has no short-read mode, and
`plasmids_missing` because Plassembler is long-read-first. Both are noted rather than
replaced; adding SPAdes for one check is not worth a container.

## The fetcher

Extends `dev/reference_genome_gates.py`, which already surveys which runs exist and
classifies them. New behaviour:

- `--platforms hifi,ont,clr,illumina` — which read types to fetch, best run per
  (genome, platform) by depth.
- `--min-depth 30` — skip a run too shallow to say anything. This drops *Phocaeicola
  vulgatus*'s 5x Illumina set and *Streptococcus salivarius*'s CLR run, which has no FASTQ
  on the portal at all.
- Paired runs: ENA's `fastq_ftp` is semicolon-separated; `_1`/`_2` go to `reads` and
  `reads_2`, a third unpaired file is ignored.
- `fastq_md5` checked after each transfer; a mismatch deletes the file and fails loudly,
  rather than leaving a truncated FASTQ for minimap2 to choke on.
- Writes `samplesheet.yml` with one sample per (genome, platform), relative paths, and
  `read_platform` / `read_provenance` set.
- `fetch_reads.sh` stays the artefact for a machine that cannot run Python with network
  access, and gains the same md5 checks.

What is actually available, as of the 2026-10-16 survey (best run per genome and platform,
whole-genome libraries only, own sample or a run naming the strain):

| platform | genomes | FASTQ files | total | median depth |
|---|---:|---:|---:|---:|
| Illumina / 454 | 12 | 19 | 19.1 Gbp | 374x |
| Nanopore (MinION, R9-era) | 5 | 5 | 6.2 Gbp | 179x |
| PacBio CLR | 5 | 4 | 8.8 Gbp | 243x |

13 of the 18 public genomes have at least one usable run; 8 have two or more platforms.
The 6 local members (4 IGA isolates, 2 METT) are unaffected — their reads are not public.

## Phases

| phase | work | done when |
|---|---|---|
| **A** | `conf/platforms.config`; `read_platform` + `read_provenance` + `reads_2` samplesheet keys; `meta.platform_*` populated in `main.nf`; `MAP_READS`, `INSPECTOR`, `FLYE_UNMAPPED`, `PREPARE_READS` take their values from it; `when:` conditions on `HIFIADAPTERFILT`, `KMC_HISTOGRAM`/`GENOMESCOPE2`, `AUTOCYCLER_GENOME_SIZE`, `MERYL_COUNT`/`MERQURY`, `PLASSEMBLER` | an ONT and an Illumina sample both complete in a stub run; pytest for the platform map; nf-test asserting the long-read-only processes are absent for an Illumina sample |
| **B** | `bin/pair_anomalies.py` + `PAIR_ANOMALIES` module; `discordant_pair_pileups` in `qc_gates.py` and the thresholds file | pytest over a synthetic BAM with a planted misjoin; the gate trips on it and not on a clean one |
| **C** | per-platform threshold overrides in `qc_thresholds.yml` and `qc_gates.py --platform` | pytest showing a CLR sample and a HiFi sample getting different verdicts from the same measurement |
| **D** | fetcher: `--platforms`, `--min-depth`, paired output, md5, relative paths, per-platform samplesheet rows | a samplesheet that `nextflow run -preview` accepts, with no absolute paths in it |
| **E** | breaks table in the report: candidate positions by method and platform | renders for a genome with two platforms |
| **F** | run the 20HM set on the HPC; record what tripped and recalibrate | `dev/reference_genome_gates.md` carries the results and the revised thresholds |

Phases A and D are independent and can run in parallel; B depends on A for `meta.platform_paired`;
C and E depend on B for the new measurement; F depends on all of them.

## Risks

- **ONT R9 at 5–10% error will make clip pile-ups noisy.** `clip_min_length` 300 and a
  raised `clip_min_fraction` are guesses. If phase F shows every ONT sample tripping, the
  honest outcome may be that ONT contributes coverage and Inspector calls only.
- **CLR may contribute nothing but coverage.** It is included because it is cheap once the
  platform table exists, and excluded from any conclusion if phase F says so.
- **454 at 21–57x on three genomes** is thin, single-end and homopolymer-prone. Coverage
  only, and `variants_high_af` should probably be disabled for it in phase C.
- **The `sr` preset changes what "unmapped" means.** Short reads fail to map uniquely across
  rRNA operons and IS elements — exactly the regions a 16S reference panel cares about. The
  breaks table must not treat a multi-mapping dropout as a break, so phase B needs to
  distinguish zero coverage from zero *unique* coverage.
- **Same-strain-not-same-sample reads (7 of 13)** can disagree with the assembly for real
  biological reasons. Design decision 5 surfaces it; it cannot be solved.
