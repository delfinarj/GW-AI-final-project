# What is in `results/`

One folder per analysis. Every file has a `<name>.provenance.json` beside it recording the script,
the git commit that was checked out when the run started, the SHA-256 of each input, the parameters
and the seed. The pages are built from these files, so changing a number here and rebuilding changes
the pages; no number is typed into a page by hand.

| Folder | Analysis | Contents |
|---|---|---|
| `release_rate/` | R1 | the single-electron rate of the public SENSEI SNOLAB release, and its figure |
| `release_mask_bits/` | R2 | which bit of the release's `mask` branch is which mask, tested by geometric signature |
| `null_false_positive_rates/` | R3 | how often each adaptive mask fires with every defect switched off |
| `muon_mask_null_rate/` | R7 | how often the muon mask fires with no tracks at all |
| `cross_defect_false_positives/` | R5 | what each mask fires on when one defect at a time is switched on |
| `compare_masks/` | R4 | the transplant comparison, one folder per seed |
| `adaptive_on_public/` | R6 | the adaptive masks run on the real release |
| `report_numbers.json` | &mdash; | the headline figures of R4, written by `analysis/build_report.py` so that the five-page summary reads them instead of recomputing them |

## The seeds of `compare_masks/`

Five seeds, each a folder:

- **20260915, 20260916, 20260917** were used while the masks were being developed. Anything seen in
  them could have influenced the code, so on their own they would not show that the method was not
  tuned to them.
- **20260920 and 20260921** were run once, after the code was frozen, and were never used to change
  anything. They are the held-out seeds, and the page reports them separately.
- Seeds **20260918 and 20260919** appear in the git history and are not here. They were made with an
  earlier version in which the pixel threshold read the simulator's truth; when that was fixed they
  were superseded, but a failed `git rm` left them on disk and a figure briefly reported "7
  independent seeds". They were removed in commit `11cc397`. `PROVENANCE.md` tells the story.

`relative_fom_summary.json` is derived: `analysis/figure_compare_masks.py` reads the five seed files
and writes it alongside the figure, so the figure and the summary cannot disagree.
