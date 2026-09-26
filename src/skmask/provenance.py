"""Sidecar records: every output file gets a JSON next to it saying where it came from."""
import datetime
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


# Text files are hashed with their line endings normalised to LF. Git checks the same file out with
# CRLF on Windows and LF elsewhere, so hashing the bytes on disk would make every recorded hash of a
# text input depend on the operating system of whoever cloned the repository, and a reader on Linux
# would find that none of them matches. Binary inputs (the ROOT data, the figures, the PDFs) are
# hashed byte for byte, which is what they are compared on.
TEXT_SUFFIXES = {".json", ".md", ".py", ".html", ".css", ".txt", ".yml", ".yaml", ".toml", ".cfg", ".C"}


def sha256(path):
    path = Path(path)
    digest = hashlib.sha256()
    if path.suffix in TEXT_SUFFIXES:
        digest.update(path.read_bytes().replace(b"\r\n", b"\n"))
        return digest.hexdigest()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_state():
    def run(*args):
        return subprocess.run(["git", "-C", str(REPO_ROOT), *args],
                              capture_output=True, text=True).stdout.strip()
    # dirty = uncommitted changes to tracked files that can change an output: code, tests, docs,
    # environment. Changes under results/ and report/ are outputs, not inputs: a run that follows another
    # in the same session always finds the earlier outputs rewritten, which says nothing about the code.
    # "XY path" per line; the whole output is stripped, so split on whitespace instead of slicing columns
    changed = [line.split(None, 1)[-1] for line in run("status", "--porcelain", "--untracked-files=no").splitlines()
               if line.strip()]
    code_changes = [f for f in changed if not f.startswith(("results/", "report/"))]
    return {"commit": run("rev-parse", "HEAD"), "dirty": bool(code_changes),
            "changed_outputs": len(changed) - len(code_changes)}


# The code that produces a result is the code present when the run starts. Capturing the state only
# when the sidecar is written would record whatever was committed or edited during a long run.
STATE_AT_START = git_state()


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
        "git": STATE_AT_START,
        "git_at_write": git_state(),
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
