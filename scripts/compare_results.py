"""Compare a reproduced results/ tree with a reference one (for example, a clean clone against the
committed results).

Every *.json of the reference except the *.provenance.json sidecars (which hold dates, commits and
paths that legitimately differ) must exist in the reproduction with the same keys, the same strings
and numbers equal within a relative tolerance. Figures are compared by SHA-256 for information only,
since image encoders may embed metadata.

Usage:  python scripts/compare_results.py <reference_results_dir> <reproduced_results_dir> [--rel 1e-9]
Exit code 0 if every JSON matches, 1 otherwise.
"""
import argparse
import hashlib
import json
import math
import sys
from pathlib import Path


def compare(a, b, path, rel, problems, worst):
    if isinstance(a, dict) and isinstance(b, dict):
        for key in sorted(set(a) | set(b)):
            if key not in a or key not in b:
                problems.append(f"{path}/{key}: present in only one file")
            else:
                compare(a[key], b[key], f"{path}/{key}", rel, problems, worst)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            problems.append(f"{path}: length {len(a)} vs {len(b)}")
        for i, (x, y) in enumerate(zip(a, b)):
            compare(x, y, f"{path}[{i}]", rel, problems, worst)
    elif isinstance(a, bool) or isinstance(b, bool) or a is None or b is None or isinstance(a, str):
        if a != b:
            problems.append(f"{path}: {a!r} vs {b!r}")
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if math.isnan(a) and math.isnan(b):
            return
        scale = max(abs(a), abs(b), 1e-300)
        diff = abs(a - b) / scale
        worst[0] = max(worst[0], diff)
        if diff > rel:
            problems.append(f"{path}: {a!r} vs {b!r} (relative difference {diff:.2e})")
    elif a != b:
        problems.append(f"{path}: {a!r} vs {b!r}")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reference")
    parser.add_argument("reproduced")
    parser.add_argument("--rel", type=float, default=1e-9, help="relative tolerance for numbers")
    parser.add_argument("--only", nargs="*", default=None, help="restrict to these sub-paths of results/")
    args = parser.parse_args()
    ref, new = Path(args.reference), Path(args.reproduced)

    files = sorted(p for p in ref.rglob("*.json") if not p.name.endswith(".provenance.json"))
    if args.only:
        files = [p for p in files if any(str(p.relative_to(ref)).replace("\\", "/").startswith(o) for o in args.only)]
    failed = 0
    for f in files:
        rel_path = f.relative_to(ref)
        other = new / rel_path
        if not other.exists():
            print(f"MISSING   {rel_path}")
            failed += 1
            continue
        problems, worst = [], [0.0]
        compare(json.loads(f.read_text(encoding="utf-8")), json.loads(other.read_text(encoding="utf-8")),
                "", args.rel, problems, worst)
        status = "MATCH    " if not problems else "DIFFER   "
        failed += bool(problems)
        print(f"{status} {rel_path}  (largest relative difference {worst[0]:.1e})")
        for p in problems[:10]:
            print(f"            {p}")
        if len(problems) > 10:
            print(f"            ... {len(problems) - 10} more")

    for png in sorted(ref.rglob("*.png")):
        if args.only and not any(str(png.relative_to(ref)).replace("\\", "/").startswith(o) for o in args.only):
            continue
        other = new / png.relative_to(ref)
        if other.exists():
            same = "identical bytes" if sha256(png) == sha256(other) else "bytes differ (information only)"
            print(f"FIGURE    {png.relative_to(ref)}: {same}")

    print(f"\n{len(files) - failed} of {len(files)} result files match")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
