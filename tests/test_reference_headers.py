"""Headers checked against independently derived fields (see tests/reference/)."""

from __future__ import annotations

import io
import json
from pathlib import Path

import pytest

from fastatacular import read_fasta, write_fasta

_REF = json.loads((Path(__file__).parent / "reference" / "header_reference.json").read_text(encoding="utf-8"))


def _parse(header: str):
    (entry,) = read_fasta(io.StringIO(f"{header}\nMKV\n"))
    return entry


@pytest.mark.parametrize("case", _REF["cases"], ids=lambda c: f"{c['source']}:{c['header'][1:30]}")
def test_header_fields_match_reference(case):
    entry = _parse(case["header"])
    assert {k: getattr(entry, k) for k in case["expected"]} == case["expected"]


@pytest.mark.parametrize("case", _REF["cases"], ids=lambda c: f"{c['source']}:{c['header'][1:30]}")
def test_reference_header_writes_back_byte_exact(case):
    buf = io.StringIO()
    write_fasta([_parse(case["header"])], buf)
    assert buf.getvalue() == f"{case['header']}\nMKV\n"


@pytest.mark.parametrize("case", _REF["known_differences"], ids=lambda c: c["header"][1:30])
def test_known_differences(case):
    entry = _parse(case["header"])
    actual = {k: getattr(entry, k) for k in case["expected"]}
    if actual == case["expected"]:
        pytest.fail("known difference is gone; move this case to 'cases'")
    pytest.xfail(case["reason"])
