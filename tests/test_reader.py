"""Tests for the FASTA reader."""

from __future__ import annotations

import io

import pytest

from fastatacular import FastaParseError, FastaReader, read_fasta


def _read_str(text: str):
    return read_fasta(io.StringIO(text))


def test_plain_single_entry():
    entries = _read_str(">myid simple description\nACDEFG\nHIJKLM\n")
    assert len(entries) == 1
    e = entries[0]
    assert e.identifier == "myid"
    assert e.sequence == "ACDEFGHIJKLM"
    assert e.description == "simple description"
    assert e.pname == "simple description"
    assert e.prefix is None
    assert e.accession is None
    assert e.gname is None
    assert e.extra == {}


def test_multiple_entries():
    text = ">a one\nAAA\n>b two\nCCC\nGGG\n>c three\nTTT\n"
    entries = _read_str(text)
    assert [e.identifier for e in entries] == ["a", "b", "c"]
    assert [e.sequence for e in entries] == ["AAA", "CCCGGG", "TTT"]


def test_uniprot_style_header():
    header = ">sp|P12345|EX_HUMAN Example protein OS=Homo sapiens OX=9606 GN=EXMP PE=1 SV=2\nMKTIIALSYIFCLVFA\n"
    [e] = _read_str(header)
    assert e.prefix == "sp"
    assert e.accession == "P12345"
    assert e.entry_name == "EX_HUMAN"
    assert e.pname == "Example protein"
    assert e.os_name == "Homo sapiens"
    assert e.ncbi_tax_id == 9606
    assert e.gname == "EXMP"
    assert e.pe == 1
    assert e.sv == 2
    assert e.sequence == "MKTIIALSYIFCLVFA"


def test_ncbi_style_pipe_id():
    [e] = _read_str(">gi|12345|ref|NP_000001.1| some description\nACDEFG\n")
    assert e.identifier == "gi|12345|ref|NP_000001.1|"
    assert e.prefix == "gi"
    assert e.accession == "12345"
    assert e.entry_name is None
    assert e.description == "some description"


def test_extra_keys_captured():
    [e] = _read_str(">id name part FOO=bar BAZ=qux extra trailing\nA\n")
    assert e.extra == {"FOO": "bar", "BAZ": "qux extra trailing"}
    assert e.pname == "name part"


def test_skips_blank_and_comment_lines():
    text = "\n; this is a comment\n>id description\nACDE\n\nFGHI\n"
    [e] = _read_str(text)
    assert e.sequence == "ACDEFGHI"


def test_sequence_before_header_raises():
    with pytest.raises(FastaParseError) as exc:
        _read_str("ACDEFG\n>id description\nACDE\n")
    assert exc.value.line == 1


def test_empty_sequence_raises():
    with pytest.raises(FastaParseError):
        _read_str(">id description\n>next other\nACDE\n")


def test_empty_header_raises():
    with pytest.raises(FastaParseError):
        _read_str("> \nACDE\n")


def test_lazy_reader_context_manager(tmp_path):
    p = tmp_path / "x.fasta"
    p.write_text(">a x\nAAA\n>b y\nCCC\n")
    with FastaReader(p) as reader:
        ids = [e.identifier for e in reader]
    assert ids == ["a", "b"]


def test_lazy_reader_requires_context_manager():
    reader = FastaReader(io.StringIO(">a x\nAAA\n"))
    with pytest.raises(RuntimeError):
        next(iter(reader))


def test_header_with_no_description():
    [e] = _read_str(">justanid\nACDE\n")
    assert e.identifier == "justanid"
    assert e.description is None
    assert e.pname is None
    assert e.sequence == "ACDE"


def test_invalid_ox_falls_back_to_extra():
    [e] = _read_str(">id name OX=notanumber\nACDE\n")
    assert e.ncbi_tax_id is None
    assert e.extra == {"OX": "notanumber"}
