"""The 1.1 fast paths give the same result as the reference implementations.

- ``_split_kv`` against the ``_KV_PATTERN`` regex it replaces.
- The chunked path reader (``_iter_path_entries``, used for file paths) against the
  line-by-line reader (``_iter_entries``, used for text handles), at many chunk sizes.
"""

from __future__ import annotations

import dataclasses
import io
from pathlib import Path

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from fastatacular import FastaParseError, FastaReader, read_fasta
from fastatacular import _parser as parser


def _reference_kv(text: str) -> tuple[int | None, list[tuple[str, str]]]:
    matches = list(parser._KV_PATTERN.finditer(text))
    return (matches[0].start() if matches else None), [(m["key"], m["val"].strip()) for m in matches]


# Text dense in '=' , key characters and every kind of whitespace (incl. \n, \x1c, NBSP).
_kv_alphabet = st.sampled_from(list("ab_Z09=( )|.-") + ["\t", "\n", "\r", "\x0b", "\x1c", "\xa0", " ", "é", "Ω"])
_kv_text = st.lists(_kv_alphabet, max_size=40).map("".join) | st.text(max_size=60)


@given(_kv_text)
@settings(max_examples=2000)
def test_split_kv_matches_reference_regex(text: str) -> None:
    assert parser._split_kv(text) == _reference_kv(text)


@pytest.mark.parametrize(
    "text",
    [
        "Example protein OS=Homo sapiens OX=9606 GN=EXMP PE=1 SV=2",
        "Protein(EC=2.7.1) OS=Homo sapiens",
        "pH=7 sensor",
        "a=b=c d\t=e f\tG=h",
        "x K=v",
        "1ab=c _k=v",
        "a=\nb=c",
        "no keys here",
        "",
    ],
)
def test_split_kv_examples(text: str) -> None:
    assert parser._split_kv(text) == _reference_kv(text)


# ------------------------------------------------------------ path reader vs line reader

_line = st.one_of(
    st.from_regex(r"\A>[A-Za-z0-9|_]{0,8}( [A-Za-z =]{0,20})?\Z"),
    st.from_regex(r"\A[ACDEFGHIKLMNPQRSTVWY]{1,12}\Z"),
    st.from_regex(r"\A[ ACDK*\t-]{0,8}\Z"),
    st.sampled_from([";comment", "# PEFF 1.0", "# //", "#x", "", "  ", ">", " >x", "﻿>b"]),
)


def _outcome(fn):  # type: ignore[no-untyped-def]
    try:
        return [dataclasses.astuple(e) for e in fn()]
    except FastaParseError as e:
        return ("error", str(e), e.line, e.context)


@given(
    st.lists(_line, max_size=12),
    st.booleans(),
    st.booleans(),
    st.integers(1, 9),
)
@settings(suppress_health_check=[HealthCheck.function_scoped_fixture], max_examples=400)
def test_path_reader_matches_line_reader(
    tmp_path: Path, lines: list[str], trailing_newline: bool, bom: bool, chunk: int
) -> None:
    text = "\n".join(lines) + ("\n" if trailing_newline and lines else "")
    if bom:
        text = "﻿" + text
    path = tmp_path / "x.fasta"
    path.write_text(text, encoding="utf-8")
    expected = _outcome(lambda: read_fasta(io.StringIO(text.removeprefix("﻿"))))
    old_chunk = parser._CHUNK_CHARS
    parser._CHUNK_CHARS = chunk
    try:
        got = _outcome(lambda: read_fasta(path))
    finally:
        parser._CHUNK_CHARS = old_chunk
    assert got == expected


def test_path_reader_long_record_across_many_chunks(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(parser, "_CHUNK_CHARS", 7)
    seq = "ACDEFGHIKLMNPQRSTVWY" * 50
    body = "\n".join(seq[i : i + 60] for i in range(0, len(seq), 60))
    path = tmp_path / "long.fasta"
    path.write_text(f">a\n{body}\n>b desc\nMK\n")
    a, b = read_fasta(path)
    assert (a.identifier, a.sequence, b.identifier, b.sequence) == ("a", seq, "b", "MK")


def test_path_reader_error_line_numbers(tmp_path: Path) -> None:
    path = tmp_path / "bad.fasta"
    path.write_text(";c\n\n# x\nACD\n>a\nM\n")
    with pytest.raises(FastaParseError) as info:
        read_fasta(path)
    assert (info.value.line, info.value.context) == (4, "ACD")
    path.write_text(">a\nM\n\n>b\n;only a comment\n>c\nK\n")
    with pytest.raises(FastaParseError, match="no sequence data") as info:
        read_fasta(path)
    assert info.value.line == 4


def test_path_reader_streams_through_fasta_reader(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(parser, "_CHUNK_CHARS", 3)
    path = tmp_path / "s.fasta"
    path.write_text(">a\nMK\n>b\nPE\n>c\nRR")
    with FastaReader(path) as reader:
        it = iter(reader)
        assert next(it).identifier == "a"
        assert [e.sequence for e in it] == ["PE", "RR"]


def test_empty_file(tmp_path: Path) -> None:
    path = tmp_path / "empty.fasta"
    path.write_text("")
    assert read_fasta(path) == []
