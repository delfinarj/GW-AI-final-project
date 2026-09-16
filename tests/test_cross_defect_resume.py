"""A run of R5 that is interrupted and resumed must land where an uninterrupted one lands.

PROVENANCE.md claims this, because the first attempt at R5 was killed by a segmentation fault and the
analysis was made resumable rather than restarted. The claim rests on two things: each run draws from
a generator seeded by (seed, run index), so it does not depend on the runs before it, and the output
file carries the counts it was built from. Both are exercised here with a stand-in for the simulation,
which is what makes the test fast; the simulation itself is covered by the other tests.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ANALYSIS = Path(__file__).resolve().parents[1] / "analysis" / "cross_defect_false_positives.py"


@pytest.fixture
def module(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location("cross_defect_under_test", ANALYSIS)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    monkeypatch.setattr(mod, "OUT_DIR", tmp_path)
    monkeypatch.setattr(mod, "write_sidecar", lambda *a, **k: None)
    monkeypatch.setattr(mod, "configurations",
                        lambda: [("sensor_a", "cti", None), ("sensor_b", "halo", None)])

    def fake_run(sensor, rng):
        """Consume the generator the way the real one does, so the draws have to line up too."""
        draws = rng.random(4)
        return {"cti": (bool(draws[0] < 0.5), float(draws[0])),
                "hot_columns": (bool(draws[1] < 0.5), float(draws[1])),
                "halo": (bool(draws[2] < 0.5), float(draws[2])) if draws[3] < 0.9 else None,
                "serial": [(bool(d < 0.5), float(d)) for d in rng.random(mod.IMAGES_PER_RUN)],
                "low_energy_clusters": [(bool(d < 0.5), float(d)) for d in rng.random(mod.IMAGES_PER_RUN)]}

    monkeypatch.setattr(mod, "one_run", fake_run)
    return mod


def result(mod):
    return json.loads((mod.OUT_DIR / "cross_defect.json").read_text(encoding="utf-8"))


def test_resumed_run_matches_uninterrupted_run(module):
    module.main(4)
    uninterrupted = result(module)

    (module.OUT_DIR / "cross_defect.json").unlink()
    module.main(2)
    interrupted = result(module)
    assert interrupted["n_runs_completed"] == 2

    module.main(4)          # resumes from the file left behind
    resumed = result(module)

    # the numbers must match; the record of how the file was produced must not pretend they do
    without_record = lambda d: {k: v for k, v in d.items() if k != "code_commits"}
    assert without_record(resumed) == without_record(uninterrupted)
    assert [c["from_run"] for c in uninterrupted["code_commits"]] == [0]
    assert [c["from_run"] for c in resumed["code_commits"]] == [0, 2]


def test_a_finished_run_is_not_repeated(module, capsys):
    module.main(2)
    before = result(module)
    module.main(2)
    capsys.readouterr()
    assert result(module) == before


def test_counts_grow_with_runs(module):
    module.main(1)
    one = result(module)
    (module.OUT_DIR / "cross_defect.json").unlink()
    module.main(3)
    three = result(module)
    for key, cell in three["counts_this_was_built_from"].items():
        assert cell["trials"]["serial"] == 3 * mod_images(one, key, "serial")


def mod_images(one_run_result, key, mask):
    return one_run_result["counts_this_was_built_from"][key]["trials"][mask]


def test_a_file_from_other_settings_does_not_poison_a_run(module):
    (module.OUT_DIR / "cross_defect.json").write_text(
        json.dumps({"seed": 1, "alpha": 0.5, "images_per_run": 7, "n_runs_completed": 9,
                    "counts_this_was_built_from": {}}), encoding="utf-8")
    module.main(2)
    assert result(module)["n_runs_completed"] == 2
    assert np.isfinite(result(module)["alpha"])


def test_a_different_configuration_order_starts_over(module, capsys):
    """Within a run the configurations draw from one generator in sequence, so order is part of the run."""
    module.main(2)
    capsys.readouterr()
    module.configurations = lambda: [("sensor_b", "halo", None), ("sensor_a", "cti", None)]
    module.main(4)
    assert "starting over" in capsys.readouterr().out
    assert result(module)["n_runs_completed"] == 4
