"""Check every output against its provenance sidecar.

For each `<file>.provenance.json` under `results/` and `report/`, this verifies that the file next to
it still hashes to what the sidecar recorded, and that every input the sidecar lists still hashes to
what it recorded. It is the check that an outside reader can run without re-running any analysis, and
the one that fails if a result is edited by hand.

Two hashing conventions exist in this repository and the sidecar says which one applies:

- sidecars written from 2026-09-25 carry `hash_convention` and hash text with its line endings
  normalised to LF, so the value does not depend on the operating system the repository was cloned on;
- older sidecars have no such field and hash the bytes on disk, which on the machine that produced
  them means CRLF for text. Those are reported as `raw-byte sidecar` rather than as failures, since
  re-writing them would claim runs that did not happen.

Exit code 1 if any file that should match does not.

Run:  python scripts/check_sidecars.py
"""
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from skmask.provenance import TEXT_SUFFIXES, sha256 as normalised_sha256   # noqa: E402


def raw_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def matches(path, recorded, normalised):
    """True if the file hashes to what the sidecar recorded under the convention it declares."""
    if not Path(path).exists():
        return None
    if normalised:
        return normalised_sha256(path) == recorded
    # an old sidecar: the bytes may since have been renormalised in git, so accept either form
    return raw_sha256(path) == recorded or normalised_sha256(path) == recorded


def main():
    failures, old_style, checked = [], 0, 0
    for sidecar in sorted(list((ROOT / "results").rglob("*.provenance.json"))
                          + list((ROOT / "report").rglob("*.provenance.json"))):
        record = json.loads(sidecar.read_text(encoding="utf-8"))
        normalised = "hash_convention" in record
        old_style += not normalised
        output = sidecar.parent / sidecar.name[:-len(".provenance.json")]
        for path, recorded, what in ([(output, record.get("output_sha256"), "output")]
                                     + [(ROOT / i["path"], i["sha256"], "input") for i in record.get("inputs", [])]):
            if recorded is None:
                continue
            checked += 1
            verdict = matches(path, recorded, normalised)
            if verdict is None:
                failures.append(f"missing {what}: {path.relative_to(ROOT)} (from {sidecar.name})")
            elif not verdict:
                failures.append(f"{what} changed since its sidecar: {path.relative_to(ROOT)} (from {sidecar.name})")

    print(f"checked {checked} hashes over {len(list((ROOT / 'results').rglob('*.provenance.json')))} "
          f"result sidecars and {len(list((ROOT / 'report').rglob('*.provenance.json')))} page sidecars")
    if old_style:
        print(f"{old_style} sidecars predate the line-ending convention and were checked against either form")
    for failure in failures:
        print(f"  FAIL {failure}")
    print("every file matches its sidecar" if not failures else f"{len(failures)} mismatches")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
