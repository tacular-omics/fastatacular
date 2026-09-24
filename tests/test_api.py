"""1.0 API contract: error hierarchy, hashing policy, reader behaviour."""

from __future__ import annotations

import io
from collections.abc import Hashable

import pytest

import fastatacular
from fastatacular import (
    FastaError,
    FastaParseError,
    FastaReader,
    FastaWriteError,
    SequenceEntry,
    read_fasta,
    write_fasta,
)


def test_fasta_error_is_exported_value_error_base():
    assert "FastaError" in fastatacular.__all__
    assert issubclass(FastaError, ValueError)
    assert issubclass(FastaParseError, FastaError)
    assert issubclass(FastaWriteError, FastaError)


def test_parse_and_write_errors_are_caught_by_base():
    with pytest.raises(FastaError):
        read_fasta(io.StringIO("ACDE\n"))
    with pytest.raises(FastaError):
        write_fasta([SequenceEntry(identifier="", sequence="A")], io.StringIO())


def test_errors_carry_a_hint_note():
    with pytest.raises(FastaParseError) as info:
        read_fasta(io.StringIO("ACDE\n"))
    assert info.value.hint
    assert any(note.startswith("hint: ") for note in info.value.__notes__)
    with pytest.raises(FastaWriteError) as winfo:
        write_fasta([SequenceEntry(identifier="x", sequence="")], io.StringIO())
    assert winfo.value.hint


def test_sequence_entry_is_explicitly_unhashable():
    entry = SequenceEntry(identifier="x", sequence="A")
    assert SequenceEntry.__hash__ is None
    assert not isinstance(entry, Hashable)
    with pytest.raises(TypeError):
        hash(entry)
    with pytest.raises(TypeError):
        {entry}  # noqa: B018
    # Equality still compares by value.
    assert entry == SequenceEntry(identifier="x", sequence="A")


def test_write_validates_every_entry_before_writing(tmp_path):
    good = SequenceEntry(identifier="ok", sequence="ACDE")
    bad = SequenceEntry(identifier="bad", sequence="")
    buf = io.StringIO()
    with pytest.raises(FastaWriteError):
        write_fasta([good, bad], buf)
    assert buf.getvalue() == ""

    path = tmp_path / "out.fasta"
    with pytest.raises(FastaWriteError):
        write_fasta([good, bad], path)
    assert not path.exists()


def test_write_error_names_the_entry_index():
    good = SequenceEntry(identifier="ok", sequence="ACDE")
    bad = SequenceEntry(identifier="bad", sequence="")
    with pytest.raises(FastaWriteError, match=r"^Entry 2: "):
        write_fasta([good, good, bad], io.StringIO())
    with pytest.raises(FastaWriteError) as info:
        write_fasta([good, good, bad], io.StringIO())
    assert info.value.index == 2


def test_write_accepts_a_generator():
    entries = (SequenceEntry(identifier=f"e{i}", sequence="AC") for i in range(3))
    buf = io.StringIO()
    write_fasta(entries, buf)
    assert buf.getvalue() == ">e0\nAC\n>e1\nAC\n>e2\nAC\n"


def test_reader_enter_returns_self_and_second_iteration_is_empty():
    with FastaReader(io.StringIO(">a\nAC\n>b\nDE\n")) as reader:
        assert isinstance(reader, FastaReader)
        assert [e.identifier for e in reader] == ["a", "b"]
        assert list(reader) == []
