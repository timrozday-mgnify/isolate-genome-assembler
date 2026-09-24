# Full-read assemblies and assembly scoring

Status: all four steps landed 2026-09-23; the real validation is the Phase 6 run. Extends
[implementation_plan.md](implementation_plan.md) as Phase 7.

Decisions taken 2026-09-23: run **all** enabled assemblers on the full read set (cost
accepted); **exclude Plassembler**; `always_score` is **on by default**.

## Goal

Run **every enabled assembler on the full read set**, not just Flye, measure the quality
of each resulting assembly, and use those measurements three ways:

1. **Selection.** When the Autocycler consensus is not fully resolved, fall back to the
   *best* full-read assembly rather than always to Flye.
2. **Reporting.** Show the user what each assembler produced for their sample and why one
   was chosen, instead of a single unexplained `fallback_flye`.
3. **Benchmarking.** Answer the question Phase 6 cannot currently ask: for this data type,
   how does each assembler do on its own, and does the consensus actually beat the best
   single assembler?

The third is the real prize. The 20hm run fell back to Flye on five of eight samples, and
we have no evidence that Flye was the right fallback for any of them.

## What already exists and is reused

The pipeline already runs one assembler on the full read set, so the shape is proven:

- `FLYE as FLYE_FULL` in [assembly.nf](../subworkflows/local/assembly.nf), fed
  `[meta, 'full', reads]` — the same tuple shape as a subset, with `'full'` in place of
  the subset id.
- `NORMALISE_HEADERS as NORMALISE_HEADERS_FULL`, already handling the `full` subset label.
- `PLASSEMBLER as PLASSEMBLER_FULL` in [finishing.nf](../subworkflows/local/finishing.nf),
  showing the aliasing pattern across subworkflows.
- `SEQKIT_STATS` (already used in `READ_QC`) for contig counts, N50 and total length.
- `MERYL_COUNT` + `MERQURY_MERQURY` in `CHECKS` for reference-free QV and k-mer
  completeness.
- `MAP_READS` and `CLIPPING_PILEUPS` for read support and misjoin evidence.

Almost nothing here is new machinery; the work is generalising `_FULL` from one assembler
to all of them, and adding one scoring step.

## Design decisions

### 1. Full-read assemblies must **not** enter the Autocycler consensus

This is the decision most likely to be got wrong, so it is first.

Autocycler's consensus is only meaningful because its inputs are **independent**: each
subset is a different draw of reads, so a random assembly error appears in one input and
is outvoted. A full-read assembly shares reads with every subset. Feeding it to
`autocycler compress` would add a vote correlated with all 28 others, inflating bridge
support for whatever that assembler got wrong and biasing `resolve` toward it.

Full-read assemblies therefore go into their own channel and their own published
directory. `ch_assemblies` — the channel that reaches `AUTOCYCLER_CONSENSUS` — keeps
carrying subset assemblies only.

### 2. Scoring is reference-free and cheap

Scoring must work on a real sample with no truth genome, and must not cost more than the
assemblies themselves. Four signals, all from tools already in the pipeline:

| Signal | Source | What it catches |
|---|---|---|
| `contigs`, `circular_contigs`, `n50`, `total_length` | `SEQKIT_STATS` + FASTA headers | Fragmentation; a finished genome is few circular contigs |
| `size_ratio` = `total_length` / `genome_size_used` | `read_qc.tsv` | Collapsed repeats (too small), duplicated haplotypes (too large) |
| `merqury_qv`, `merqury_completeness` | `MERYL_COUNT` + `MERQURY_MERQURY` | Base accuracy, and whether the assembly holds the reads' k-mers |
| `unmapped_read_fraction`, `clipping_pileups` (confirmed only) | `MAP_READS` + `CLIPPING_PILEUPS` | Missing sequence, and misjoins |

CheckM2 and BUSCO are deliberately **excluded** from scoring: they need large databases,
cost more than the rest combined, and measure gene content rather than assembly structure.
They still run once, on the selected assembly, in `CHECKS`.

`MERYL_COUNT` must move earlier so the meryl database is built once per sample and shared
between scoring and `CHECKS` (see [Component changes](#component-changes)).

### 3. The ranking rule is ordered and deterministic, not a weighted sum

A weighted score invites endless tuning and hides why a choice was made. Instead: hard
filters, then an ordered comparison, with a final alphabetical tie-break so the same
inputs always give the same answer.

**Hard filters** — an assembly failing any of these is ranked last, never selected unless
nothing else survives:

- `size_ratio` outside `[1 - params.score_size_tolerance, 1 + params.score_size_tolerance]`
  (default tolerance `0.1`)
- `merqury_qv` below `params.score_min_qv` (default `40`, matching the existing
  `merqury_qv` fail threshold in `assets/qc_thresholds.yml`)
- empty or unreadable FASTA

**Ordered comparison**, first difference wins:

1. Fewest **confirmed** clipping pile-ups — a misjoin is the worst thing an assembly can
   have, and the permissive-realignment verdict means this now counts real ones only
2. Most circular contigs, **as determined by `CIRCULARISE` (decision 6), never by the
   assembler's own tag** — Raven, metaMDBG and LJA emit no circularity at all
   ([normalise_headers.py:239](../bin/normalise_headers.py)), so ranking on the tag would
   zero those three regardless of how good their assemblies are
3. Fewest contigs
4. Highest `merqury_completeness`
5. Highest `merqury_qv`
6. Assembler name, alphabetically (tie-break only)

Rationale for the order: structural correctness before structural completeness before base
accuracy. An assembly with one misjoin and QV 60 is worse than a slightly noisier one that
is correctly joined, because polishing can fix the second and nothing downstream fixes the
first.

Every input to the comparison is written to `full_assemblies.tsv`, so a user can always
see why one assembly beat another.

### 4. The consensus still wins when it is resolved

Selection order is unchanged in spirit: a fully resolved Autocycler consensus is still
preferred over any single assembler, because it is a consensus. Scoring only decides the
fallback. `params.assembly_selection` allows overriding this:

- `score` (default) — consensus if resolved, else best-scoring full-read assembly
- `flye` — current behaviour, consensus if resolved, else Flye; kept so an existing run
  can be reproduced exactly
- `always_score` (**default**) — score the consensus alongside the full-read assemblies
  and take the winner outright.

`always_score` being the default means the consensus has to earn its place on every
sample rather than by assumption. The asymmetry this creates is real and is dealt with in
[decision 7](#7-there-is-no-single-assembly-equivalent-of-resolve); it is not a reason to
avoid the comparison, but it must be disclosed in the report next to the result.

### 5. Cost is opt-outable, not free

Full-read assembly of 4 subsets' worth of reads costs roughly what all 4 subsets cost for
that assembler. Running all seven roughly **doubles** stage-3 compute.

**Decided: the cost is accepted.** `params.full_read_assemblers` defaults to
`params.assemblers` minus Plassembler, so every real assembler runs on the full read set.
The param still exists so a user can narrow it (`flye,hifiasm,raven,metamdbg` is the cheap
set) without touching `--assemblers` and changing the consensus inputs. Canu is the one to
watch: it is already the slowest per subset.

**Plassembler is excluded** and this is enforced, not merely defaulted: it is a
plasmid-only assembler, so as a whole-genome candidate it would always be filtered on
`size_ratio` and would only add noise to the table. It already runs on the full read set in
`FINISHING` for the plasmid audit, so including it here would also duplicate that work.
`score_assemblies.py` drops it with a stated reason if a user forces it in.

### 6. `CIRCULARISE`: the trim equivalent, applied uniformly

`autocycler trim` cannot be reused. It takes `--cluster_dir <dir containing
1_untrimmed.gfa>` (verified against the pinned 0.7.0 image), so it is cluster- and
graph-shaped, not FASTA-shaped. Feeding it a single assembly would mean a contrived
one-assembly `compress` + `cluster` first, and its `--mad` length-outlier logic is
meaningless with one cluster member. Not worth it.

What `trim` does that matters here is one thing: **remove the duplicated sequence a
circular contig carries at its ends.** The standard way to do that for HiFi is minimap2
self-alignment of the contig ends plus a trim, and the pipeline already does the alignment
half — it just throws the answer away.

#### Reuse, not new alignment code

`CONTIG_ENDS` already emits exactly the PAF this needs: it cuts the first and last
`contig_end_window` bases (default 10 kb) of every contig into `<contig>_start` and
`<contig>_end` records and self-aligns them with `minimap2 -x asm5`, keeping only
cross-record hits. `classify_replicons.py` then reads that PAF for *existence* of a
`_start` vs `_end` hit and discards the coordinates
([classify_replicons.py:51](../bin/classify_replicons.py)).

So `CIRCULARISE` is `CONTIG_ENDS` aliased into `scoring.nf` plus one new script that reads
the coordinates the current parser drops. No new alignment, no new container, no new
dependency.

#### The trim, precisely

For each `<contig>_start` → `<contig>_end` PAF record on the same contig, with `W` the
window actually used (`min(contig_end_window, contig_length)`) and PAF fields
`qstart/qend/strand/tstart/tend/matches/block_len`:

A record is a **wrap** only when all of:

- `strand == '+'` — a reverse hit is an inverted repeat, not a wrap
- `qstart <= slop` — the match starts at the very beginning of the contig
- `tend >= W - slop` — and runs to the very end of it
- `matches / block_len >= circularise_min_identity` (default 0.95)

with `slop` a small tolerance (propose 50 bp, a named constant, not a param).

Then the duplicated span is `L = qend - qstart`, and the contig is rewritten as
`seq[:len(seq) - L]` with `circular=true`. Trimming the tail rather than the head keeps
the contig's start coordinate stable, which matters because `dnaapler` rotates afterwards
and `CONTIG_ENDS` in `FINISHING` runs before rotation for the same reason.

A record that fails the anchoring conditions but passes identity is an **internal repeat
near the ends** and must be left alone. That is the tandem-repeat case in the tests.

#### Prerequisite bug: the two windows can overlap each other

`CONTIG_ENDS` cuts `_start = seq[0:W]` and `_end = seq[-W:]`. When `contig_length < 2 * W`
those two windows **overlap inside the contig**, so they self-align for a trivial reason and
`self_overlapping()` records a false `end_overlap`. When `contig_length <= W` they are the
same sequence outright.

This is not hypothetical. **All 14 `end_overlap` flags in the 20hm run are false positives**,
and the arithmetic is exact — the shared region begins at `L - W`, which is where the
alignment starts:

| contig | length | predicted `qstart` (`L - W`) | observed `qstart` |
|---|---|---|---|
| flye_full_1 | 10930 | 930 | 932 |
| flye_full_6 | 11934 | 1934 | 1934 |
| flye_full_7 | 11389 | 1389 | 1392 |
| flye_full_12 | 10161 | 161 | 179 |
| flye_full_21 | 10056 | 56 | 70 |
| flye_full_32 | 11813 | 1813 | 1819 |
| flye_full_79 | 10764 | 764 | 779 |

The remaining seven (isolate04 `flye_full_5/9/11/14/53`, isolate08 `flye_full_1/6`, all
2885–9166 bp) are shorter than `W`, so their two windows are byte-identical and the PAF
shows `qstart == tstart` and `qend == tend`.

**Fix, in `CONTIG_ENDS`, before anything else in this plan:**

```awk
window = (len < 2 * window) ? int(len / 2) : window
```

so the windows are always disjoint. A genuine wrap longer than `contig_length / 2` is not
a wrap, it is a tandem duplication of the whole contig, so nothing real is lost.

This is worth shipping on its own: it is a one-line fix to a flag that is currently wrong
on every contig it fires for, and `classify_replicons.py` needs no change.

Note that the anchoring conditions above would already have rejected all 14 — none has
`qstart ≈ 0` and `tend ≈ W` — so `circularise.py` is safe against this class even before
the fix. The fix matters because the flag itself is user-visible and wrong.

#### When the overlap fills the window

If a real overlap is longer than `W`, the alignment fills the window on both sides
(`qstart ≈ 0`, `qend ≈ W`, `tend ≈ W`) and the true overlap is unknowable from this PAF.
Trimming `W` would silently under-trim and leave a contig that still looks circular but is
not. So when `qend - qstart >= W - slop`, `circularise.py` **does not trim**: it writes
`circular=false`, records `overlap=W+`, and flags `overlap_exceeds_window`. The fix is a
larger `--contig_end_window` and a re-run, which is why the flag names the cause.

#### Calibrate before trusting it

The anchoring rule above is derived from geometry, and the only real PAFs available
(`assemblies/*/finishing/*.end_overlaps.paf` in the 20hm run) contain **no** genuine wrap
to confirm it against — every record is the window-overlap artefact. So `circularise.py`
trims only the unambiguous signature and flags everything else rather than guessing.

Before Step 3 is allowed to select on circularity, run `circularise.py` over the 20hm
assemblies and confirm it trims nothing. A run that trims anything there has found either a
genuine missed wrap or a bug, and either way wants looking at.

#### Why it also fixes the scoring

Circularity is currently whatever each assembler chose to report, and **Raven, metaMDBG and
LJA report none at all** — the comment at
[normalise_headers.py:239](../bin/normalise_headers.py) says so explicitly, while Canu even
gets its overlap trimmed there. Ranking on the assembler's tag would zero those three
regardless of quality. Deriving circularity from the sequence makes it assembler-agnostic,
which is a precondition for ranking on it at all (rung 2 of decision 3).

### 7. There is no single-assembly equivalent of `resolve`

`autocycler resolve` reconciles disagreement **between** inputs: it finds anchor unitigs
present once in every input, builds bridges from what each input says lies between them,
and culls the conflicting ones. A single assembly has no disagreement to reconcile, so
there is nothing to port. Repeat resolution within one assembly is the assembler's own job
and has already happened by the time we see the FASTA.

So the fair comparison is **consensus after trim+resolve vs full-read assembly after
`CIRCULARISE`**, and the missing `resolve` is an inherent asymmetry to disclose, not to
close.

It can, however, be **measured**. The clipping analysis added in the same release reports
confirmed pile-ups with `tail_target` — the joins the reads support that the assembly did
not make. That is precisely the work `resolve` would have had to do. Ranking rung 1 is
therefore already the right handicap: an assembly that needed resolving and did not get it
carries confirmed pile-ups with agreeing tails, and loses on the first comparison. The
report should say this next to the `always_score` result rather than presenting the
comparison as neutral.

## Component changes

### Params (`nextflow.config`)

```groovy
// Assemblers also run on the full read set, scored, and used as the fallback when the
// consensus does not resolve. Defaults to --assemblers; narrow it when the full-read
// arm's cost matters more than covering every assembler.
full_read_assemblers = null            // null -> params.assemblers minus plassembler
assembly_selection   = 'always_score'  // always_score | score | flye
score_size_tolerance = 0.1
score_min_qv         = 40
// CIRCULARISE reuses the existing contig_end_window; this is the identity a start-end
// self-match must reach before the duplicated span is cut.
circularise_min_identity = 0.95
```

### `subworkflows/local/assembly.nf`

- Add `as *_FULL` aliases for every assembler, mirroring the existing `FLYE_FULL` line.
  The miniasm arm needs three (`MINIASM_OVERLAP_FULL`, `MINIASM_FULL`, `MINIPOLISH_FULL`).
- Build `ch_full_in = ch_reads.map { meta, reads -> [meta, 'full', reads] }` and gate each
  `_FULL` call on membership of `full_read_assemblers`, exactly as the subset calls gate on
  `assemblers`.
- Mix the `_FULL` outputs into `ch_native_full`, through `NORMALISE_HEADERS_FULL`, into a
  new emit `full_assemblies` — **not** into `ch_assemblies`.
- `FLYE_FULL` stays as it is and is simply one member of the new set when `flye` is in
  `full_read_assemblers`. When it is not, Flye still runs on full reads as the
  last-resort fallback, so `select_assembly.py` always has something.

### New: `subworkflows/local/scoring.nf`

Takes `[meta, assembler, assembly]` for every full-read assembly plus the consensus, and
the shared meryl DB. For each candidate: `CIRCULARISE` first (decision 6), then
`SEQKIT_STATS`, `MERQURY_MERQURY`, `MAP_READS`, `CLIPPING_PILEUPS` on the circularised
FASTA. Emits one `[meta, assembler, metrics.tsv]` per candidate, plus the circularised
assembly, which is what selection and `FINISHING` then use.

The consensus goes through `CIRCULARISE` too. It should be a no-op — `autocycler trim`
already removed its start-end overlaps — and if it is not, that is worth knowing.

`MERYL_COUNT` moves here from `CHECKS`, and its DB is emitted so `CHECKS` joins it rather
than rebuilding. This is the only change to an existing subworkflow's internals.

### New: `bin/score_assemblies.py` + `modules/local/score_assemblies`

Reads every candidate's metrics for one sample, applies the filters and ordered comparison
from decision 3, writes `<id>.full_assemblies.tsv`: one row per candidate with every input
signal plus `filtered` (why it was excluded, blank if not), `rank`, and `selected`.

Pure Python over TSVs, so it is unit-testable with no containers — same shape as
`qc_gates.py`.

### New: `bin/circularise.py` + `modules/local/circularise`

Inputs: the candidate FASTA and the PAF from `CONTIG_ENDS` (aliased as
`CONTIG_ENDS_CANDIDATE` in `scoring.nf` — see decision 6; no new alignment is written).
Outputs: the trimmed FASTA and `<id>.<assembler>.circularity.tsv`, one row per contig:

    contig  length_before  overlap  length_after  identity  circular  flag

where `flag` is blank, `overlap_exceeds_window`, or `internal_repeat`. Pure Python over a
PAF, so it is unit-testable with no container, same shape as `qc_gates.py`.

`CONTIG_ENDS` keeps its existing job in `FINISHING` unchanged — it runs there before
rotation, and rotation would move a duplicated end into the middle where neither check
would see it.

`classify_replicons.py`'s `self_overlapping()` can then be simplified to read
`circularity.tsv` instead of re-parsing the PAF, or left alone. Not required by this plan;
noted so the duplication is deliberate rather than forgotten.

### `bin/select_assembly.py`

- Gains `--scores <full_assemblies.tsv>` and `--candidates <assembler>=<fasta> ...`.
- The `fallback_flye` label becomes `fallback_<assembler>`, so
  `qc_gates.consensus_unresolved` (which tests `== "autocycler"`) keeps working unchanged.
- `reason` gains the ranking evidence, e.g.
  `consensus assembly not fully resolved; best full-read assembly: hifiasm (1 circular contig, QV 58.2, 0 clipping pile-ups)`.
- Under `assembly_selection = 'flye'` the current code path is taken verbatim.

Four existing tests assert the literal `fallback_flye`
([test_select_assembly.py:59,69](../tests/test_select_assembly.py),
[test_qc_gates.py:113](../tests/test_qc_gates.py),
[test_collect_metrics.py:48](../tests/test_collect_metrics.py)) and need updating to the
new label.

### Reporting

- `collect_metrics.py`: new kind `full-assemblies` → `full_assemblies.tsv`, concatenated
  as an ordinary headed TSV (no parser needed).
- `assets/report/isolate_assembly_report.qmd`, `assembly_tab`:
  - The hard-coded fallback sentence
    ([qmd:238](../assets/report/isolate_assembly_report.qmd)) must stop saying "The Flye
    assembly of all reads" and name the assembler actually chosen.
  - New **Full-read assemblies** section: the `full_assemblies.tsv` table with the selected
    row highlighted, plus a grouped bar chart of contigs / circular contigs / QV per
    assembler so the comparison is visible at a glance.
- Status board (`sample_rows` in `collect_metrics.py`): `assembly_source` already flows
  through and will show `fallback_hifiasm` etc. with no code change.

**Fold in while here:** five run_summary tables are computed but never rendered —
`coverage_regions`, `unmapped_assembly`, `plasmid_audit`, `assembler_contribution`,
`merqury_completeness`. `unmapped_assembly` in particular is what explains a sample like
isolate02. These are independent of this plan and can ship separately, but the report work
is in the same file.

### Benchmark (`dev/`)

- `benchmark_arms.tsv` gains a `full_read_assemblers` column.
- `assembler_benchmark.py`:
  - Score **every** full-read assembly against truth, not only the delivered one. This is
    the new primary table: per assembler, replicons recovered, edit distance, cost.
  - New comparison table: per sample, consensus vs best full-read assembly, so the
    "is Autocycler worth it?" question has an answer per arm.
  - `ASSEMBLER_PROCESSES` maps bare process names (`FLYE`, `CANU`, …) to assemblers for
    cost attribution. The `_FULL` aliases produce task names like `FLYE_FULL`, which will
    **not** match — the mapping needs `_FULL` entries, and the cost table needs a
    subset-vs-full split so the arm cost stays interpretable.
  - The docstring's "otherwise the full-read Flye fallback, as the pipeline would" is now
    wrong and must track `assembly_selection`.

## Testing

| Level | Test |
|---|---|
| pytest | `score_assemblies.py`: each hard filter excludes; each comparison rung decides when the ones above tie; alphabetical tie-break is stable |
| pytest | `select_assembly.py`: picks the top-ranked candidate; falls back to Flye when every candidate is filtered; `assembly_selection='flye'` reproduces current behaviour |
| pytest | updated `fallback_flye` → `fallback_<assembler>` assertions in the three existing test files |
| nf-test (stub) | with `full_read_assemblers = 'flye,raven'`, both full assemblies are produced, `full_assemblies.tsv` has two rows, and the consensus input still contains subset assemblies only |
| nf-test (stub) | an unresolved consensus selects a non-Flye assembler when that one scores better |
| pytest | `classify_replicons.self_overlapping()`: a contig shorter than `2 * contig_end_window` no longer self-flags — regression for the 14 false positives in the 20hm run |
| pytest | `circularise.py`, one case per branch of decision 6: a clean start-end overlap is cut once and marked circular; identity below `circularise_min_identity` is left alone; `qstart` well past 0, or `tend` well short of `W`, is an internal repeat and is left alone; a `-` strand hit (inverted repeat) is left alone; a window-filling match is **not** trimmed and is flagged `overlap_exceeds_window`; a short contig uses `W = length / 2` and produces no spurious wrap |
| render | report renders with the new section from a fixture; `tests/data/generate_run_summary.py` gains `full_assemblies` and `circularity` tables |

The stub caveat from the clipping work applies again: stub mode does not run real
assemblers, so these tests pin wiring and selection logic, not assembly quality. Real
validation is the Phase 6 benchmark run.

## Risks and open questions

1. ~~**Cost.**~~ **Resolved 2026-09-23:** accepted. All assemblers run on the full read
   set. Expect stage 3 to roughly double.
2. **`always_score` compares a consensus against single assemblies, which is not
   apples-to-apples.** The consensus has had `resolve` applied and single assemblies
   structurally cannot ([decision 7](#7-there-is-no-single-assembly-equivalent-of-resolve)).
   `CIRCULARISE` closes the `trim` half of the gap; the `resolve` half stays open by
   construction. Mitigation is disclosure plus ranking rung 1, which charges an assembly
   for exactly the joins `resolve` would have made. **The report must state this beside the
   comparison** rather than presenting it as neutral.
3. **A circular contig is not necessarily correct.** Ranking rung 2 rewards circularity, and
   a chimeric join can be circular — `CIRCULARISE` will happily close a contig whose ends
   overlap for the wrong reason. Rung 1 (confirmed clipping pile-ups) is the guard, which is
   why it sits above it. Worth re-examining after the benchmark run.
4. ~~**Plassembler**~~ **Resolved 2026-09-23:** excluded, enforced in
   `score_assemblies.py`, reasons in [decision 5](#5-cost-is-opt-outable-not-free).
6. **`CIRCULARISE` is new code cutting sequence.** Trimming the wrong span silently
   shortens a genome. Guards: a minimum identity (`circularise_min_identity`, default 0.95),
   the match must run to the contig end, the trimmed span is recorded per contig in
   `circularity.tsv`, and `length_before`/`length_after` are both kept so an unexpected cut
   is visible rather than inferred. Unit tests cover a clean overlap, a partial match that
   must not be cut, and a tandem repeat at one end that must not be mistaken for a wrap.
5. **Does the fallback assembly deserve finishing?** It goes through `FINISHING` exactly as
   the consensus does today, so no change — but a fallback chosen by score is more likely to
   be a good assembly than Flye-by-default was, which may change what `dnaapler` and
   `classify_replicons` see. No action, just something to watch in the first real run.

## Delivery

Three steps, each independently shippable and testable:

**Step 0 — fix the end-window overlap. Landed 2026-09-23.** One-line `CONTIG_ENDS` change so the two windows
are disjoint (decision 6). Independent of everything else here, and fixes a flag that is
currently wrong on all 14 contigs it fires for in the 20hm run. Ships on its own.

**Step 1 — produce the assemblies. Landed 2026-09-23.** `_FULL` aliases for every assembler except
Plassembler, `full_read_assemblers` param, new `full_assemblies` emit, published under
`assemblies/<id>/full/<assembler>.fasta`. Nothing consumes them yet. Stub test asserts they
exist and that consensus inputs are unchanged.

**Step 2 — circularise and score. Landed 2026-09-23.** `circularise.py` and its module, `scoring.nf`,
`MERYL_COUNT` moved and shared, `score_assemblies.py`, `full_assemblies.tsv`, report
section. Selection still uses Flye. The table is useful on its own: it tells you what the
fallback *would* have been, and `CIRCULARISE` on the existing Flye fallback is a fix in its
own right — untrimmed end overlaps are inflating that assembly's length today.

**Step 3 — select on score. Landed 2026-09-23.** `select_assembly.py` changes, `assembly_selection` param
defaulting to `always_score`, label change and its test updates, benchmark changes.

Steps 1–2 change nothing a user receives except the Flye fallback's end overlaps; Step 3 is
where selection changes. Splitting there keeps the risky change isolated and reviewable.
