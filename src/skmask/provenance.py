"""Sidecar records: every output file gets a JSON next to it saying where it came from."""
import datetime
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_state():
    def run(*args):
        return subprocess.run(["git", "-C", str(REPO_ROOT), *args],
                              capture_output=True, text=True).stdout.strip()
    # dirty = uncommitted changes to tracked files; new untracked files cannot affect this output
    return {"commit": run("rev-parse", "HEAD"),
            "dirty": bool(run("status", "--porcelain", "--untracked-files=no"))}


def _relative(path):
    path = Path(path).resolve()
    try:
        return str(path.relative_to(REPO_ROOT)).replace("\\", "/")
    except ValueError:
        return str(path)


def write_sidecar(output, script, inputs=(), parameters=None, seed=None, results=None, notes=""):
    """Write `<output>.provenance.json` and return its path."""
    import numpy
    import scipy

    record = {
        "output": _relative(output),
        "output_sha256": sha256(output) if Path(output).exists() else None,
        "script": _relative(script),
        "git": git_state(),
        "inputs": [{"path": _relative(p), "sha256": sha256(p)} for p in inputs],
        "parameters": parameters or {},
        "seed": seed,
        "results": results or {},
        "notes": notes,
        "created": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "environment": {"python": sys.version.split()[0], "platform": platform.platform(),
                        "numpy": numpy.__version__, "scipy": scipy.__version__},
    }
    sidecar = Path(str(output) + ".provenance.json")
    sidecar.write_text(json.dumps(record, indent=2, default=float), encoding="utf-8")
    return sidecar
