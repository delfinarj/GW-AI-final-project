"""Download the public SENSEI SNOLAB data release and verify it byte for byte.

Source: https://github.com/sensei-skipper/DataReleases (folder Snolab_binned_1), CC-BY-4.0,
the data behind arXiv:2410.18716. The SHA-256 digests below were recorded when the files were
first downloaded (2026-09-14); a mismatch means the upstream release changed and every result
built on it must be re-checked.

Usage:  python scripts/fetch_public_data.py
"""
import hashlib
import sys
import urllib.request
from pathlib import Path

BASE_URL = "https://raw.githubusercontent.com/sensei-skipper/DataReleases/main/Snolab_binned_1"
DEST = Path(__file__).resolve().parents[1] / "data" / "public" / "sensei_snolab_binned_1"

FILES = {
    "RELEASE_hits_blinded_EXP0_13.root":
        "e38824424e019e2e5096786e8fa38834da7c2bf61a6087214fe0ec11df02a585",
    "RELEASE_hits_blinded_EXP7200_13.root":
        "21d7294f4238dd83450b9aec56178599ff4eabf62729aeb16a94621b6b58e712",
    "RELEASE_hits_blinded_EXP21600_13.root":
        "9864c9e342b34077759560ee4240faab42f57a8c960f8e832ee517cd71d0f7e6",
    "RELEASE_hits_blinded_EXP72000_13.root":
        "79a09c60ecf8e05ca2f742b7816985cf5a31b3ce51b9f925c96f0e47fe4624fb",
    "plotRate.C":
        "11cd68984866cf1450e8ba2fa305b9662d64ccbd32b87b695cb91c440c85876c",
}


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    failures = 0
    for name, expected in FILES.items():
        path = DEST / name
        if not path.exists() or sha256(path) != expected:
            print(f"downloading {name}")
            urllib.request.urlretrieve(f"{BASE_URL}/{name}", path)
        actual = sha256(path)
        status = "ok" if actual == expected else "CHECKSUM MISMATCH"
        failures += actual != expected
        print(f"{status:>17}  {name}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
