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
    (entry,) = read_fasta(io.StringIO(">sp|P00001|FOO Original header text OS=Mus OX=10090\nAAA\n"))
    out = _write_str([entry])
    assert out.startswith(">sp|P00001|FOO Original header text OS=Mus OX=10090\n")


def test_raw_header_that_disagrees_with_fields_is_ignored():
    entry = SequenceEntry(
        identifier="anything",
        sequence="AAA",
        raw_header="sp|P00001|FOO Original header text OS=Mus OX=10090",
    )
    assert _write_str([entry]) == ">anything\nAAA\n"


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


_UNIPROT = ">sp|P12345|EX_HUMAN Example protein OS=Homo sapiens OX=9606 GN=ABC PE=1 SV=2\nACDE\n"


def test_unedited_parsed_entry_writes_raw_header_byte_exact():
    src = ">sp|P12345|EX_HUMAN  Example   protein\tOS=Homo sapiens  GN=ABC \nACDE\n"
    (parsed,) = read_fasta(io.StringIO(src))
    assert _write_str([parsed]) == src


def test_edited_structured_field_wins_over_raw_header():
    (parsed,) = read_fasta(io.StringIO(_UNIPROT))
    edited = dataclasses.replace(parsed, gname="XYZ")
    assert _write_str([edited]).splitlines()[0] == (
        ">sp|P12345|EX_HUMAN Example protein OS=Homo sapiens OX=9606 GN=XYZ PE=1 SV=2"
    )


def test_edited_identifier_wins_over_raw_header():
    (parsed,) = read_fasta(io.StringIO(">x some protein GN=A\nACDE\n"))
    edited = dataclasses.replace(parsed, identifier="DECOY_x")
    assert _write_str([edited]).splitlines()[0] == ">DECOY_x some protein GN=A"


def test_edited_extra_wins_over_raw_header():
    (parsed,) = read_fasta(io.StringIO(">x name FOO=1\nACDE\n"))
    edited = dataclasses.replace(parsed, extra={"FOO": "2"})
    assert _write_str([edited]).splitlines()[0] == ">x name FOO=2"


def test_cleared_field_stays_cleared_when_entry_has_no_name():
    # Found by Hypothesis: with no protein name the rebuild fell back to the
    # parsed (stale) description and re-added the cleared key from it.
    (parsed,) = read_fasta(io.StringIO(">x OS=Homo sapiens GN=A\nACDE\n"))
    assert _write_str([dataclasses.replace(parsed, gname=None)]).splitlines()[0] == ">x OS=Homo sapiens"
    (parsed,) = read_fasta(io.StringIO(">UPI0000000005 status=active\nACDE\n"))
    assert _write_str([dataclasses.replace(parsed, extra={})]).splitlines()[0] == ">UPI0000000005"


def test_cleared_pname_stays_cleared():
    (parsed,) = read_fasta(io.StringIO(">x Old name GN=A\nACDE\n"))
    assert _write_str([dataclasses.replace(parsed, pname=None)]).splitlines()[0] == ">x GN=A"


def test_edited_description_is_still_used_as_fallback():
    (parsed,) = read_fasta(io.StringIO(">x OS=Homo sapiens\nACDE\n"))
    edited = dataclasses.replace(parsed, description="New name OS=Homo sapiens")
    assert _write_str([edited]).splitlines()[0] == ">x New name OS=Homo sapiens"


@pytest.mark.parametrize(
    "change",
    [
        {"pname": "two\nlines"},
        {"os_name": "Homo\rsapiens"},
        {"extra": {"K": "a\nb"}},
        {"identifier": "a\nb"},
        {"identifier": "a b"},
        {"raw_header": "x\ny", "identifier": "x\ny"},
        {"sequence": "AC\n>DE"},
        {"sequence": "AC DE"},
        {"sequence": ">ACDE"},
        {"sequence": ";ACDE"},
    ],
)
def test_values_that_would_corrupt_the_file_raise(change):
    # Found by Hypothesis: a line break in a field was written as-is, so the
    # file read back with a different header or an extra entry.
    entry = dataclasses.replace(SequenceEntry(identifier="x", sequence="ACDE"), **change)
    with pytest.raises(FastaWriteError):
        _write_str([entry])


def test_gt_or_semicolon_inside_a_sequence_line_is_written():
    entry = SequenceEntry(identifier="x", sequence="AC>D;E")
    assert _write_str([entry]) == ">x\nAC>D;E\n"
    with pytest.raises(FastaWriteError):
        _write_str([entry], line_width=2)
