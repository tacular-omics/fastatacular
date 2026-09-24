"""Property-based tests (Hypothesis) for header fields, sequences and file framing."""

from __future__ import annotations

import dataclasses
import io
import re
from typing import Any

import pytest
from hypothesis import assume, given
from hypothesis import strategies as st

from fastatacular import FastaParseError, FastaWriteError, SequenceEntry, read_fasta, write_fasta

HEADER_FIELDS = (
    "identifier",
    "prefix",
    "accession",
    "entry_name",
    "pname",
    "gname",
    "os_name",
    "ncbi_tax_id",
    "pe",
    "sv",
    "extra",
)
TYPED_KEYS = {"OS", "OX", "GN", "PE", "SV"}
_KEY_RE = re.compile(r"(?:^|\s)[A-Za-z_][A-Za-z0-9_]*=")

# ---------------------------------------------------------------- strategies

# Printable text without line breaks: what a single header line can hold.
_line_chars = st.characters(
    blacklist_categories=("Cs", "Cc", "Zl", "Zp"),
    blacklist_characters="\x85﻿",
)
_word = st.text(_line_chars.filter(lambda c: not c.isspace() and c not in "|="), min_size=1, max_size=8)
_phrase = st.lists(_word, min_size=1, max_size=5).map(" ".join)


def _free_text(s: str) -> bool:
    """Text usable as a name/value: no ``KEY=`` token, no edge whitespace."""
    return s == s.strip() and not _KEY_RE.search(s)


name_text = _phrase.filter(_free_text)
accession = st.from_regex(r"\A[A-Z][0-9][A-Z0-9]{3}[0-9](-[0-9]{1,2})?\Z")
entry_name = st.from_regex(r"\A[A-Z0-9]{1,10}_[A-Z0-9]{1,5}\Z")
key = st.from_regex(r"\A[A-Za-z_][A-Za-z0-9_]{0,6}\Z").filter(lambda k: k not in TYPED_KEYS)
seq = st.text(alphabet="ACDEFGHIKLMNPQRSTVWYXBZUO*-", min_size=1, max_size=150)


@st.composite
def identifiers(draw: st.DrawFn) -> dict[str, Any]:
    kind = draw(st.sampled_from(["uniprot", "ncbi_gi", "plain"]))
    if kind == "uniprot":
        db, acc, name = draw(st.sampled_from(["sp", "tr"])), draw(accession), draw(entry_name)
        return {"identifier": f"{db}|{acc}|{name}", "prefix": db, "accession": acc, "entry_name": name}
    if kind == "ncbi_gi":
        gi = str(draw(st.integers(1, 10**9)))
        return {"identifier": f"gi|{gi}|ref|NP_{gi}.1|", "prefix": "gi", "accession": gi, "entry_name": None}
    ident = draw(st.from_regex(r"\A[A-Za-z0-9_.:-]{1,20}\Z"))
    return {"identifier": ident, "prefix": None, "accession": None, "entry_name": None}


@st.composite
def entries(draw: st.DrawFn) -> SequenceEntry:
    """Entries built from structured fields only (no raw header, no description)."""
    ident = draw(identifiers())
    return SequenceEntry(
        sequence=draw(seq),
        pname=draw(st.none() | name_text),
        os_name=draw(st.none() | name_text),
        ncbi_tax_id=draw(st.none() | st.integers(0, 10**7)),
        gname=draw(st.none() | name_text),
        pe=draw(st.none() | st.integers(1, 5)),
        sv=draw(st.none() | st.integers(1, 99)),
        extra=draw(st.dictionaries(key, name_text, max_size=3)),
        **ident,
    )


def _write(entries_: list[SequenceEntry], **kw: int) -> str:
    buf = io.StringIO()
    write_fasta(entries_, buf, **kw)
    return buf.getvalue()


def _fields(entry: SequenceEntry) -> dict[str, object]:
    return {name: getattr(entry, name) for name in (*HEADER_FIELDS, "sequence")}


# ------------------------------------------------ 1. parse(write(e)) == e


@given(st.lists(entries(), min_size=1, max_size=4), st.integers(0, 80))
def test_parse_write_roundtrip(ents: list[SequenceEntry], width: int) -> None:
    back = read_fasta(io.StringIO(_write(ents, line_width=width)))
    assert [_fields(e) for e in back] == [_fields(e) for e in ents]


@given(entries())
def test_write_parse_write_is_stable(entry: SequenceEntry) -> None:
    once = _write([entry])
    assert _write(read_fasta(io.StringIO(once))) == once


_FIELD_STRATEGIES: dict[str, st.SearchStrategy[object]] = {
    "pname": st.none() | name_text,
    "os_name": st.none() | name_text,
    "gname": st.none() | name_text,
    "ncbi_tax_id": st.none() | st.integers(0, 10**7),
    "pe": st.none() | st.integers(1, 5),
    "sv": st.none() | st.integers(1, 99),
    "extra": st.dictionaries(key, name_text, max_size=3),
}


@given(entries(), st.data())
def test_edited_fields_survive(entry: SequenceEntry, data: st.DataObject) -> None:
    """Parse, edit (set, change or clear) one field, write, parse: the edit is kept."""
    parsed = read_fasta(io.StringIO(_write([entry])))[0]
    name = data.draw(st.sampled_from(["pname", "os_name", "ncbi_tax_id", "gname", "pe", "sv", "extra"]))
    new = data.draw(_FIELD_STRATEGIES[name])
    edited = dataclasses.replace(parsed, **{name: new})
    back = read_fasta(io.StringIO(_write([edited])))[0]
    assert getattr(back, name) == new
    for other in HEADER_FIELDS:
        assert getattr(back, other) == getattr(edited, other), other


# ------------------------------------------- 2. unedited raw headers byte-exact

raw_header = st.text(_line_chars, min_size=1, max_size=120).filter(lambda s: s.strip() != "")


@given(st.lists(raw_header, min_size=1, max_size=4), seq)
def test_unedited_raw_header_is_byte_exact(headers: list[str], sequence: str) -> None:
    text = "".join(f">{h}\n{sequence}\n" for h in headers)
    ents = read_fasta(io.StringIO(text))
    assert [e.raw_header for e in ents] == headers
    assert _write(ents, line_width=0) == text


# ------------------------------------------- 3. BOM, CRLF and blank lines

blank = st.sampled_from(["", " ", "\t", "  \t "])


@given(
    st.lists(entries(), min_size=1, max_size=3),
    st.booleans(),
    st.sampled_from(["\n", "\r\n"]),
    st.lists(blank, max_size=3),
    st.integers(0, 3),
)
def test_bom_crlf_blank_lines(ents: list[SequenceEntry], bom: bool, eol: str, blanks: list[str], every: int) -> None:
    lines = _write(ents, line_width=7).splitlines()
    out: list[str] = []
    for i, line in enumerate(lines):
        out.append(line)
        if every and i % every == 0:
            out.extend(blanks)
    text = ("﻿" if bom else "") + eol.join([*blanks, *out]) + eol
    via_handle = read_fasta(io.StringIO(text, newline=""))
    assert [_fields(e) for e in via_handle] == [_fields(e) for e in ents]
    # Header text never keeps a stray carriage return or BOM.
    assert all("\r" not in e.raw_header and "﻿" not in e.raw_header for e in via_handle)


@given(st.lists(entries(), min_size=1, max_size=3), st.booleans(), st.sampled_from(["\n", "\r\n"]))
def test_bom_crlf_from_path(
    tmp_path_factory: pytest.TempPathFactory, ents: list[SequenceEntry], bom: bool, eol: str
) -> None:
    path = tmp_path_factory.mktemp("p") / "x.fasta"
    text = _write(ents).replace("\n", eol)
    path.write_bytes((b"\xef\xbb\xbf" if bom else b"") + text.encode("utf-8"))
    assert [_fields(e) for e in read_fasta(path)] == [_fields(e) for e in ents]


# --------------------------------------------- malformed input: typed errors


@given(st.text(max_size=200))
def test_arbitrary_text_parses_or_raises_parse_error(text: str) -> None:
    try:
        read_fasta(io.StringIO(text))
    except FastaParseError:
        pass


@given(entries(), st.sampled_from(["pname", "os_name", "gname", "identifier", "sequence"]), st.data())
def test_unwritable_entries_raise_write_error(entry: SequenceEntry, name: str, data: st.DataObject) -> None:
    """A value that would corrupt the file raises FastaWriteError instead of writing it."""
    bad = data.draw(st.sampled_from(["\n", "\r", "a\nb", "x\r\ny"]))
    edited = dataclasses.replace(entry, **{name: bad})
    assume(edited != entry)
    with pytest.raises(FastaWriteError):
        _write([edited])
