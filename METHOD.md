# How this project was made

Every line of code, every analysis and every figure in this repository was written and run by AI
agents. This file records how that was organised, what the agents got wrong, and what was done to
make the result trustworthy rather than merely fast. It exists because a reader can check the
*results* from `PROVENANCE.md` and the generated page, but cannot otherwise check the *method* — and
the method is the part of this work that was actually novel for us.

## The division of labour

**The people** chose the problem and the question, set the scope (six masks, three sensor presets,
simulator plus public data), decided what would count as an answer, supplied the private background
reading and the rule that it must never enter the repository, decided the deadlines and what was
worth two hours of computer time, and are responsible for the result being correct.

**The implementing agent** (Claude Opus 5, in the Claude Code CLI) wrote the simulator, the six masks
in both forms, every analysis script, the tests, the figures, the generated page and the prose of
`PROVENANCE.md`; ran everything; and recorded provenance as it went.

**The reviewing agents** were separate agents given no access to the conversation that produced the
work. They saw only the repository. Their job was to attack it.

## The rules, written before the work

`AGENTS.md` was written first and holds the standing rules an agent working here must follow. The
ones that did the most work:

- **No number on any presented page is typed by hand.** `analysis/build_report.py` and
  `analysis/build_summary.py` read the result JSONs and write the pages, so a re-run changes the
  pages consistently and a sentence cannot drift away from the number it describes.
- **A result with no recorded check is not finished.** Every entry in `PROVENANCE.md` has to say how
  the result was checked, not only where it came from.
- **Every output carries a sidecar** (`<output>.provenance.json`) with the script, the git commit at
  the start of the run, the SHA-256 of every input, the parameters and the seed.
- **Nothing from `Contexto/` may be read into a tracked file.** The private material (a licenciatura
  thesis, a colleague's doctoral thesis, an internal mask specification) was used as background
  reading only. Where detail was missing, the order was: look for it on the public web, then in the
  references of those documents, then ask the authors. Every number in `src/skmask/presets.py` cites
  a public source or is marked as a scenario value.

## The loop that was actually used

1. **State the question and the expectation first.** Before an analysis was run, what it was expected
   to show was written into `PLAN.md` and committed. This is why the page can say that the trail
   calibration on the real sensor did *not* do what was expected: the prediction is in the history,
   before the result.
2. **Freeze the code, then run.** The commit is recorded in the sidecar at the start of the run, not
   at the end, so a long run cannot silently claim code that was written while it was running.
3. **Record the result and its check** in `PROVENANCE.md`, including the choices that had a
   defensible alternative and what the result does *not* show.
4. **Hand the repository to a reviewing agent with no memory of the work** and ask it to find what is
   wrong. Apply what it finds, or record why not.
5. **Reproduce from a clean clone** — a fresh checkout that shares nothing with the working copy but
   the commit and the pinned environment — and compare numerically at a relative tolerance of
   10⁻⁹.

## What the reviews found

Four independent reviews were run, each by an agent that had not seen the work being produced.

| Review | What it was asked | What it found |
|---|---|---|
| First, on the simulator and the first cross-sensor comparison | find anything wrong | The pixel threshold of every "adaptive" mask was set from the simulator's *true* mean charge — truth leaking into a procedure whose whole claim is that it uses only data. Also: the headline compared masks against an oracle without ever asking whether masking helped at all |
| Second, on the corrected analyses | check the statistics and the claims | Exclusions applied to medians instead of per seed; harm decided by rounding; oracle optima sitting on the edges of their parameter grids; Python version unpinned |
| Third, on the one-defect-at-a-time test and the real-data analysis | check for overstatement | The halo mask measured distance in superpixels on a sensor that bins 32 rows into one, so a radius of 15 covered the whole height of the frame and hid the loudest column from the hot-column calibration; the page reported a sum of overlapping masks as if it were their union; the page gave only the flattering direction of the agreement with the published mask; a verdict applied to sixty cells with no correction for their number |
| Fourth, an audit against the assignment | what is missing | The development material lived outside the repository; there was no account of the method (this file); `PLAN.md` had gone stale; sidecar hashes depended on the operating system's line endings |

Findings that were **accepted and not fixed** are recorded too, in `PROVENANCE.md`: that R5's
pre-registration cannot be demonstrated from the commit order alone, and that the two images of a
run are not independent trials for two of the masks.

## The mistakes the agents made

These are in `PROVENANCE.md` and on the report page as "what the checks caught". They are collected
here because in a course about AI-assisted research the failure modes are the point:

- **Truth leaking into a procedure that claimed not to use it.** The most serious. Caught by a
  review, not by a test, because the test suite had been written by the same agent that wrote the bug.
- **A simulator that redrew its hot columns in every image**, which would have made every
  stack-based calibration meaningless. Caught by an injection test.
- **Geometry assumed to be isotropic** on a sensor whose pixels are not. Caught by a review.
- **Statistics applied to the wrong unit**: a family of sixty tests treated as one.
- **A reproduction command that did nothing.** The repository ships a finished result file, and the
  analysis saw it and skipped the work, so the instruction in `README.md` silently no-opped in a
  clean clone. Caught by actually doing the clean-clone run rather than trusting it.
- **Commit messages that claimed more than the commit held** ("five seeds" in a commit carrying
  four), and one wrong diagnosis written confidently into `PROVENANCE.md` (a numerical difference
  blamed on a scipy version, when a pinned clone showed the cause was the BLAS library). Both were
  corrected in place, with the correction left visible.

The pattern is consistent: the agent's own tests and prose are not a check on the agent. What caught
things were injection tests written against simulated truth, independent agents with no context, and
running the thing from scratch somewhere else.

## What is uncomfortable about this method, stated plainly

- The agents wrote the prose that describes their own work, including this file. The reviews are the
  counterweight, and they were themselves agents.
- Pre-registration is only as strong as the commit order that demonstrates it. For one analysis (R5)
  that order does not prove what the text claims, and `PROVENANCE.md` says so.
- An agent that is asked to find problems will find *something*. The reviews were therefore asked to
  state plainly when a section was clean, and their findings were checked against the data before
  being acted on; at least one was downgraded from "wrong" to "worth saying" on inspection.

## Reproducing the method, not just the results

`README.md` reproduces the numbers. To reproduce the *method*, the prompts matter more than the code:
each reviewing agent was given the repository path, the constraint never to open `Contexto/`, a list
of specific things to attack, and the instruction to separate "definitely wrong" from "arguable" and
to say plainly when something was fine. The clean-clone check was run by an agent told to behave as
someone who had never spoken to the authors and to follow only `README.md`.
