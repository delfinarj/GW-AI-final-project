"""Print the report page to PDF with a headless Chromium browser (Chrome or Edge).

The PDF is the same document as the HTML page, so no number is copied by hand between them.

Usage:  python scripts/make_pdf.py [report/index.html] [report/report.pdf]
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CANDIDATES = [
    shutil.which("chrome"), shutil.which("google-chrome"), shutil.which("chromium"), shutil.which("msedge"),
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
]


def find_browser():
    for candidate in CANDIDATES:
        if candidate and os.path.exists(candidate):
            return candidate
    sys.exit("No Chrome, Chromium or Edge found; open the HTML page and print it to PDF instead.")


def main():
    html = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "report" / "index.html"
    pdf = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "report" / "report.pdf"
    browser = find_browser()
    subprocess.run([browser, "--headless=new", "--disable-gpu", "--no-pdf-header-footer",
                    f"--print-to-pdf={pdf.resolve()}", html.resolve().as_uri()], check=True)

    sys.path.insert(0, str(ROOT / "src"))
    from skmask.provenance import write_sidecar

    write_sidecar(pdf, __file__, inputs=[html], parameters={"browser": browser},
                  notes="the report page printed to PDF; same document, no number retyped")
    print(f"wrote {pdf} with {browser}")


if __name__ == "__main__":
    main()
