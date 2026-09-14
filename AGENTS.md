# Instructions for agents working in this repository

Read `PLAN.md` first: it states the question, the six masks, the simulator and the metrics.

## Data rules (hard)

- Never read, copy, quote or summarise anything from `Contexto/` or from any path outside this
  repository into a tracked file. That material is private. Parameters and definitions used in
  the code must come from public sources (papers, public data releases) and be cited where used.
- `data/` is not versioned. Public data is obtained only with `python scripts/fetch_public_data.py`,
  which verifies SHA-256 digests. Do not commit data files.

## Environment

```
conda env create -f environment.yml
conda activate skmask
python scripts/fetch_public_data.py
pytest
```

## Provenance (every result)

- A figure or number that appears in the deliverables is produced by one script in `analysis/`.
- The script writes, next to its output, a JSON sidecar: script path, git commit, input files
  with SHA-256, parameters, random seed, date.
- Add an entry to `PROVENANCE.md`: output, producing script, inputs, what was written from
  scratch versus called from a library, choices that had a defensible alternative, and **how the
  result was checked**. A result with no recorded check is not finished.

## Code

- Python, NumPy/SciPy. Simulator and masks live in `src/skmask/`; tests in `tests/`.
- Every adaptive mask must pass a null test (defect absent: masked fraction consistent with its
  alpha) and an injection test (defect present: recovered) before it is used on real data.
- Fix random seeds in analysis scripts.

## Language

Code, comments, documentation and deliverables are in English.
