"""The English and Spanish summaries must differ in prose and in nothing else.

Both are generated from `skmask.summary_numbers`, so they cannot read different values; what this
catches is the other way a translation goes wrong, a number typed into one of the two by hand, and it
catches a stale build: if one was regenerated after a result changed and the other was not, the two
files stop carrying the same numbers.
"""
import re
import html
from collections import Counter
from pathlib import Path

import pytest

REPORT = Path(__file__).resolve().parents[1] / "report"
ENGLISH, SPANISH = REPORT / "summary.html", REPORT / "resumen.html"


def numbers(path):
    text = html.unescape(re.sub(r"<[^>]+>", " ", path.read_text(encoding="utf-8")))
    body = text.split("</style>")[-1]          # the stylesheet's own numbers are not content
    return Counter(re.findall(r"\d+\.\d+|\b\d+\b", body))


@pytest.mark.skipif(not (ENGLISH.exists() and SPANISH.exists()),
                    reason="the summaries have not been built in this checkout")
def test_both_summaries_carry_the_same_numbers():
    english, spanish = numbers(ENGLISH), numbers(SPANISH)
    only_english = sorted((english - spanish).elements())
    only_spanish = sorted((spanish - english).elements())
    assert not only_english and not only_spanish, (
        f"numbers in one summary and not the other - English only: {only_english}, "
        f"Spanish only: {only_spanish}")
