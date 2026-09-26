# Verification logs

The claims in `PROVENANCE.md` about reproduction are not assertions about runs that happened
somewhere invisible: these are the logs those runs printed. They are raw, including the failures.
Total 272 KB, no data, no secrets.

Times in these logs are local time on the machine the project was produced on.

| File | What it is | Which claim it backs |
|---|---|---|
| `verify_clone_short.log` | The first clean-clone pass, 2026-09-15. Ends `2 of 3 result files match` | The R1 discrepancy that was first blamed on a scipy version. The per-field relative differences are here |
| `verify_clone_r123.log` | The second clone, with Python pinned to 3.11 | The correction: the same numbers came back, so the version was not the cause; the BLAS library was |
| `rerun_B_seed_*.log` | The five production runs of R4, one per seed | The oracle constants chosen per sensor, the timings, and the `edge check: passed` line each run prints. `rerun_B_seed_20260917_retry.log` is the re-run after the crash that does not reproduce; `rerun_B_seed_20260921_early.log` is the third process, which is why that seed's sidecar names a different commit |
| `clone4_seed_*.log` | The five R4 seeds re-run in a clean clone | "all five seeds … identical to the committed file, largest relative difference 0.0e+00" |
| `rerun_all.log` | The full re-run after the truth-leak fix, 71 KB | The two-stale-seed incident: it opens with the `git rm` refusal that left seeds 20260918 and 20260919 on disk, which is how a figure came to say "7 independent seeds" |
| `r5.log` | The one-defect-at-a-time run | The segmentation fault that killed the first attempt (`attempt 1 exit 139`), and the 20 runs of the second |

## What is not here

- The wrappers that drove these runs, and the diagnostic that investigated the seed-20260917 crash,
  were scratch files and were not kept. The crash is described in `PROVENANCE.md` from its output;
  the diagnostic itself no longer exists, and that is a gap in the evidence rather than a detail.
- The clean clones themselves (each carries its own virtual environment) are not committed. What a
  reader needs to repeat the check is in `README.md`: clone, `uv sync`, re-run, `compare_results.py`.
