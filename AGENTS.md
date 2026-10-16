# Repository instructions

Before implementing or changing any part of the pipeline, read
[the implementation plan](docs/implementation_plan.md). It is the approved specification
and records which phases are done.

Two further plans cover work in progress:

- [full read assemblies](docs/full_read_assemblies_plan.md)
- [multi-platform read checks](docs/multiplatform_read_checks_plan.md) -- finding breaks in
  finished genomes with non-HiFi reads. **Plan only, nothing implemented**; it is the one
  place where the HiFi-only constraint below is scheduled to change, and it does not change
  until that plan says a phase is done.

Key constraints to preserve:

- HiFi-only input: no short-read or Medaka polishing.
- Assemblers run as separate processes on pinned containers; Autocycler steps 3–7 run as
  one process. Keep the Flye fallback and the `autocycler_dir` re-entry point.
- Every check writes data consumed by `bin/qc_gates.py`; thresholds live in
  `assets/qc_thresholds.yml`, not in code. The Quarto report only reads `run_summary/`.
- Follow the repository standards table in the plan (YAML samplesheet, publish opt-in,
  closure `ext.args`, stubs, nf-test + pytest, pre-commit).
