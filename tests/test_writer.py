"""Tests for the FASTA writer."""

from __future__ import annotations

import dataclasses
import io

import pytest

from fastatacular import FastaWriteError, SequenceEntry, read_fasta, write_fasta


def _write_str(entries, **kw) -> str:
    buf = io.StringIO()
    write_fasta(entries, buf, **kw)
    return buf.getvalue()


def test_write_minimal_entry():
    entry = SequenceEntry(identifier="foo", sequence="ACDEFG")
    out = _write_str([entry])
    assert out == ">foo\nACDEFG\n"


def test_write_with_pname_and_keys():
    entry = SequenceEntry(
        identifier="sp|P12345|EX_HUMAN",
        sequence="MKTIIALSYIFCLVFA",
        pname="Example protein",
        os_name="Homo sapiens",
        ncbi_tax_id=9606,
        gname="EXMP",
        pe=1,
        sv=2,
    )
    out = _write_str([entry])
    assert out.startswith(">sp|P12345|EX_HUMAN Example protein OS=Homo sapiens OX=9606 GN=EXMP PE=1 SV=2\n")
    assert "MKTIIALSYIFCLVFA\n" in out


def test_write_wraps_at_line_width():
    entry = SequenceEntry(identifier="x", sequence="A" * 130)
    out = _write_str([entry])
    seq_lines = out.splitlines()[1:]
    assert all(len(line) <= 60 for line in seq_lines)
    assert sum(len(line) for line in seq_lines) == 130


def test_write_custom_line_width():
    entry = SequenceEntry(identifier="x", sequence="ACDEFGHIJK")
    out = _write_str([entry], line_width=4)
    assert out == ">x\nACDE\nFGHI\nJK\n"


def test_write_no_wrap_when_line_width_zero():
    entry = SequenceEntry(identifier="x", sequence="A" * 200)
    out = _write_str([entry], line_width=0)
    assert out == ">x\n" + "A" * 200 + "\n"


def test_raw_header_round_trips_exactly():
    entry = SequenceEntry(
        identifier="anything",
        sequence="AAA",
        raw_header="sp|P00001|FOO Original header text OS=Mus OX=10090",
    )
    out = _write_str([entry])
    assert out.startswith(">sp|P00001|FOO Original header text OS=Mus OX=10090\n")


def test_write_empty_identifier_raises():
    entry = SequenceEntry(identifier="", sequence="A")
    with pytest.raises(FastaWriteError):
        _write_str([entry])


def test_write_empty_sequence_raises():
    entry = SequenceEntry(identifier="x", sequence="")
    with pytest.raises(FastaWriteError):
        _write_str([entry])


def test_extra_keys_are_emitted():
    entry = SequenceEntry(
        identifier="x",
        sequence="A",
        pname="name",
        extra={"FOO": "bar"},
    )
    out = _write_str([entry])
    assert out.startswith(">x name FOO=bar\n")


def test_rebuild_from_parsed_description_does_not_duplicate_keys():
    (parsed,) = read_fasta(io.StringIO(">x OS=Homo sapiens GN=A\nA\n"))
    assert parsed.pname is None
    entry = dataclasses.replace(parsed, raw_header="")
    assert _write_str([entry]) == ">x OS=Homo sapiens GN=A\nA\n"


def test_rebuild_from_description_uses_structured_fields_and_keeps_name():
    entry = SequenceEntry(
        identifier="x",
        sequence="A",
        description="Some protein OS=Mus musculus XY=1",
        os_name="Homo sapiens",
        gname="G",
    )
    assert _write_str([entry]) == ">x Some protein OS=Homo sapiens GN=G XY=1\nA\n"
