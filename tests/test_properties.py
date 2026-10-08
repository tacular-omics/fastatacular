"""Property-based tests (Hypothesis) for header fields, sequences and file framing."""

from __future__ import annotations

import bz2
import dataclasses
import gzip
import io
import lzma
import re
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
from hypothesis import assume, given
from hypothesis import strategies as st

from fastatacular import (
    Compression,
    FastaIndex,
    FastaParseError,
    FastaReader,
    FastaWriteError,
    SequenceEntry,
    is_decoy,
    make_decoys,
    read_fasta,
    write_fasta,
)

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
        return {"identifier": f"gi|{gi}|ref|NP_{gi}.1|", "prefix": "gi", "accession": f"NP_{gi}.1", "entry_name": None}
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


# ------------------------------- 4. write -> read across compression and sinks
#
# Every combination of compression and destination kind: write_fasta then read_fasta
# gives back the same entries (full equality, raw_header and description included),
# and FastaIndex agrees with the streaming reader on every key of a plain file.

_COMPRESS: dict[str, Callable[[bytes], bytes]] = {"gzip": gzip.compress, "bz2": bz2.compress, "xz": lzma.compress}
_MAGIC = {"gzip": b"\x1f\x8b", "bz2": b"BZh", "xz": b"\xfd7zXZ\x00"}
_SUFFIX = {None: ".fasta", "gzip": ".fasta.gz", "bz2": ".fa.bz2", "xz": ".fasta.xz"}

# Residues as found in real files: upper and lower case, ambiguity codes, gaps, '*'.
_residue_text = st.text(alphabet="ACDEFGHIKLMNPQRSTVWYXBZUOacdefghiklmnpqrstvwyxu-", min_size=1, max_size=60)
any_seq = st.builds(
    lambda core, repeat, stop: core * repeat + ("*" if stop else ""),
    _residue_text,
    st.sampled_from([1, 1, 1, 2, 50, 120]),  # some sequences span many lines
    st.booleans(),
)
# Header text with tabs and any printable character, as long as it holds a word.
odd_header = st.text(_line_chars | st.just("\t"), min_size=1, max_size=120).filter(lambda s: s.strip() != "")


@st.composite
def file_entries(draw: st.DrawFn) -> SequenceEntry:
    """An entry as read from a file: a structured header or arbitrary header text."""
    sequence = draw(any_seq)
    if draw(st.booleans()):
        return dataclasses.replace(draw(entries()), sequence=sequence)
    return read_fasta(io.StringIO(f">{draw(odd_header)}\n{sequence}\n"))[0]


def _canonical(ents: list[SequenceEntry]) -> list[SequenceEntry]:
    """Entries as read back from a text buffer (sets raw_header/description)."""
    return read_fasta(io.StringIO(_write(ents)))


def _check_index(path: Path, ents: list[SequenceEntry]) -> None:
    """FastaIndex lookups equal the streaming reader's first entry for every key."""
    for key in ("identifier", "accession"):
        first: dict[str, SequenceEntry] = {}
        for e in ents:
            first.setdefault(e.identifier if key == "identifier" else e.accession or e.identifier, e)
        index = FastaIndex(path, key=key, duplicates="first")  # type: ignore[arg-type]
        assert list(index) == list(first)
        for k, e in first.items():
            assert index[k] == e


_SINK_CASES = [
    *[("path", c) for c in (None, "infer", "gzip", "bz2", "xz")],
    *[("binary", c) for c in (None, "infer", "gzip", "bz2", "xz")],
    *[("text", c) for c in (None, "infer")],  # a text handle cannot carry compressed bytes
]


@pytest.mark.parametrize(("sink", "compression"), _SINK_CASES)
@given(ents=st.lists(file_entries(), min_size=1, max_size=4), width=st.integers(0, 80), data=st.data())
def test_write_read_roundtrip_every_compression_and_sink(
    tmp_path_factory: pytest.TempPathFactory,
    sink: str,
    compression: Compression,
    ents: list[SequenceEntry],
    width: int,
    data: st.DataObject,
) -> None:
    ents = _canonical(ents)
    path = tmp_path_factory.mktemp("rt") / "x.fasta"
    # "infer" on a path compresses by suffix; on a handle it means plain text.
    kind = compression
    if compression == "infer":
        kind = data.draw(st.sampled_from([None, "gzip", "bz2", "xz"])) if sink == "path" else None
        path = path.with_name("x" + _SUFFIX[kind])
    if sink == "path":
        write_fasta(ents, path, line_width=width, compression=compression)
        back = read_fasta(path, compression=compression)
    elif sink == "binary":
        with path.open("wb") as fh:
            write_fasta(ents, fh, line_width=width, compression=compression)
            assert not fh.closed
        with path.open("rb") as fh:
            back = read_fasta(fh, compression=compression)
    else:
        # Default newline translation: \r\n on disk on Windows, read back as \n.
        with path.open("w", encoding="utf-8") as fh:
            write_fasta(ents, fh, line_width=width, compression=compression)
        with path.open(encoding="utf-8") as fh:
            back = read_fasta(fh, compression=compression)
    assert back == ents
    raw = path.read_bytes()
    if kind is None:
        assert raw.startswith(b">")
        _check_index(path, ents)
        if sink != "text":
            assert b"\r" not in raw  # the writer emits \n on every OS
    else:
        assert raw.startswith(_MAGIC[kind])
        assert read_fasta(path) == ents  # magic-byte detection


@given(
    ents=st.lists(file_entries(), min_size=1, max_size=4),
    bom=st.booleans(),
    eol=st.sampled_from(["\n", "\r\n"]),
    final_eol=st.booleans(),
    kind=st.sampled_from([None, "gzip", "bz2", "xz"]),
)
def test_crlf_bom_input_from_every_source(
    tmp_path_factory: pytest.TempPathFactory,
    ents: list[SequenceEntry],
    bom: bool,
    eol: str,
    final_eol: bool,
    kind: str | None,
) -> None:
    """Windows line endings and a BOM on input, plain or compressed, from a path or a handle."""
    ents = _canonical(ents)
    text = _write(ents, line_width=13).replace("\n", eol)
    if not final_eol:
        text = text.removesuffix(eol)
    raw = (b"\xef\xbb\xbf" if bom else b"") + text.encode("utf-8")
    path = tmp_path_factory.mktemp("in") / "x.fasta"
    path.write_bytes(raw if kind is None else _COMPRESS[kind](raw))
    assert read_fasta(path) == ents
    assert read_fasta(path, compression="infer" if kind is None else kind) == ents  # type: ignore[arg-type]
    with path.open("rb") as fh:
        assert read_fasta(fh, compression=kind) == ents  # type: ignore[arg-type]
    if kind is None:
        with path.open(encoding="utf-8", newline="") as fh:  # \r\n reaches the parser untranslated
            assert read_fasta(fh) == ents
        with FastaReader(path) as reader:
            assert list(reader) == ents
        _check_index(path, ents)


# ------------------------------------- 5. header templates of each database style
#
# Headers generated from each database's template, with the decoy and contaminant
# tags search engines put in front of the identifier, must give exactly the fields
# the template was filled with; rebuilding the header from those fields must read
# back to the same fields.

# Tags in front of the identifier. All of them are also recognized in front of ``gi|``.
TAGS = [
    "",
    "DECOY_",
    "decoy_",
    "REV_",
    "rev_",
    "Reverse_",
    "CON_",
    "CONTAM_",
    "contam_",
    "rev-0-",
    "REV-12-",
    "DECOY-3-",
]
# The UniProtKB accession format (https://www.uniprot.org/help/accession_numbers).
uniprot_acc = st.from_regex(
    r"\A(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]{1,2})?\Z"
)
refseq_acc = st.from_regex(r"\A(?:NP|XP|YP|WP|AP)_[0-9]{6,9}\.[0-9]{1,2}\Z")
genbank_acc = st.from_regex(r"\A[A-Z]{3}[0-9]{5}\.[0-9]\Z")
pdb_id = st.from_regex(r"\A[1-9][A-Z0-9]{3}\Z")
pdb_chain = st.from_regex(r"\A[A-Za-z0-9]{1,2}\Z")
number = st.integers(1, 10**9).map(str)


@st.composite
def uniprot_kv(draw: st.DrawFn) -> tuple[str, dict[str, Any]]:
    """`` OS=... OX=... GN=... PE=... SV=...`` with a random subset of keys, and the fields."""
    fields = {
        "os_name": draw(st.none() | name_text),
        "ncbi_tax_id": draw(st.none() | st.integers(1, 10**7)),
        "gname": draw(st.none() | st.from_regex(r"\A[A-Za-z0-9-]{1,10}\Z")),
        "pe": draw(st.none() | st.integers(1, 5)),
        "sv": draw(st.none() | st.integers(1, 99)),
    }
    keys = [("OS", "os_name"), ("OX", "ncbi_tax_id"), ("GN", "gname"), ("PE", "pe"), ("SV", "sv")]
    text = "".join(f" {k}={fields[f]}" for k, f in keys if fields[f] is not None)
    return text, fields


@st.composite
def templated_headers(draw: st.DrawFn) -> tuple[str, dict[str, Any]]:
    """``(header line, expected fields)`` for one database style and decoy/contaminant tag."""
    tag = draw(st.sampled_from(TAGS))
    style = draw(
        st.sampled_from(
            ["uniprot", "uniprot_archived", "gi", "gi_pdb", "ncbi_pipe", "pdb", "refseq", "uniref", "uniparc"]
        )
    )
    expected: dict[str, Any] = {"entry_name": None, "extra": {}}
    name = draw(name_text)
    if style == "uniprot":
        db, acc = draw(st.sampled_from(["sp", "tr"])), draw(uniprot_acc)
        entry = draw(entry_name)
        kv, fields = draw(uniprot_kv())
        ident = f"{tag}{db}|{acc}|{entry}"
        expected |= {"prefix": tag + db, "accession": acc, "entry_name": entry, "pname": name, **fields}
        return f">{ident} {name}{kv}", {"identifier": ident, **expected}
    if style == "uniprot_archived":
        db, acc, sv = draw(st.sampled_from(["sp", "tr"])), draw(uniprot_acc), draw(st.integers(1, 9))
        ident = f"{tag}{db}|{acc}"
        pname = "archived from Release 18.0 01-MAY-1991"
        expected |= {"prefix": tag + db, "accession": acc, "pname": pname, "sv": sv}
        return f">{ident} {pname} SV={sv}", {"identifier": ident, **expected}
    if style == "gi":
        db = draw(st.sampled_from(["ref", "gb", "emb", "dbj", "sp"]))
        acc = draw(refseq_acc if db == "ref" else genbank_acc)
        trailing = draw(st.sampled_from(["|", ""]))
        ident = f"{tag}gi|{draw(number)}|{db}|{acc}{trailing}"
        pname = f"{name} [{draw(name_text)}]"
        expected |= {"prefix": tag + "gi", "accession": acc, "pname": pname}
        return f">{ident} {pname}", {"identifier": ident, **expected}
    if style == "gi_pdb":
        pdb, chain = draw(pdb_id), draw(st.just("") | pdb_chain)
        ident = f"{tag}gi|{draw(number)}|pdb|{pdb}|{chain}"
        expected |= {"prefix": tag + "gi", "accession": f"{pdb}_{chain}" if chain else pdb, "pname": name}
        return f">{ident} {name}", {"identifier": ident, **expected}
    if style == "ncbi_pipe":
        db = draw(st.sampled_from(["ref", "gb", "emb", "dbj"]))
        acc = draw(refseq_acc if db == "ref" else genbank_acc)
        ident = f"{tag}{db}|{acc}|"
        expected |= {"prefix": tag + db, "accession": acc, "pname": name}
        return f">{ident} {name}", {"identifier": ident, **expected}
    if style == "pdb":
        # NCBI ``pdb|ENTRY|CHAIN``: three fields, so the chain is the entry_name.
        pdb, chain = draw(pdb_id), draw(pdb_chain)
        ident = f"{tag}pdb|{pdb}|{chain}"
        expected |= {"prefix": tag + "pdb", "accession": pdb, "entry_name": chain, "pname": name}
        return f">{ident} {name}", {"identifier": ident, **expected}
    # No pipes: the whole first word is the identifier, with no prefix or accession.
    expected |= {"prefix": None, "accession": None}
    if style == "refseq":
        ident = tag + draw(refseq_acc)
        pname = f"{name} [{draw(name_text)}]"
        return f">{ident} {pname}", {"identifier": ident, **expected, "pname": pname}
    if style == "uniref":
        ident = f"{tag}UniRef{draw(st.sampled_from(['50', '90', '100']))}_{draw(uniprot_acc)}"
        extra = {"n": draw(number), "Tax": draw(name_text), "TaxID": draw(number), "RepID": draw(entry_name)}
        kv = "".join(f" {k}={v}" for k, v in extra.items())
        return f">{ident} {name}{kv}", {"identifier": ident, **expected, "pname": name, "extra": extra}
    ident = tag + draw(st.from_regex(r"\AUPI[0-9A-F]{10}\Z"))
    kv, fields = draw(uniprot_kv())
    return f">{ident} {name}{kv}", {"identifier": ident, **expected, "pname": name, **fields}


@given(templated_headers(), seq)
def test_templated_headers_parse_exactly(case: tuple[str, dict[str, Any]], sequence: str) -> None:
    header, expected = case
    [entry] = read_fasta(io.StringIO(f"{header}\n{sequence}\n"))
    got = {name: getattr(entry, name) for name in HEADER_FIELDS}
    want = dict.fromkeys(HEADER_FIELDS) | expected
    assert got == want
    assert entry.raw_header == header[1:]
    assert entry.description == header[1:].split(maxsplit=1)[1]
    # Unedited: written verbatim. Rebuilt from the fields alone: same fields back, and stable.
    assert _write([entry], line_width=0) == f"{header}\n{sequence}\n"
    rebuilt = _write([dataclasses.replace(entry, raw_header="", description=None)])
    [again] = read_fasta(io.StringIO(rebuilt))
    assert _fields(again) == _fields(entry)
    assert _write([dataclasses.replace(again, raw_header="", description=None)]) == rebuilt


@given(templated_headers(), st.sampled_from(["DECOY_", "rev_", "REV_", "CON_", "rev-1-"]))
def test_decoy_headers_keep_accession_and_fields(case: tuple[str, dict[str, Any]], tag: str) -> None:
    """make_decoys puts the tag in front of the identifier; accession and description fields stay."""
    header, expected = case
    [target] = read_fasta(io.StringIO(f"{header}\nMPEPTIDEKAAR\n"))
    [decoy] = make_decoys([target], method="reverse", prefix=tag)
    assert decoy.identifier == tag + target.identifier
    assert is_decoy(decoy, prefix=tag)
    for name in HEADER_FIELDS:
        if name not in ("identifier", "prefix"):
            assert getattr(decoy, name) == getattr(target, name), name
    assert decoy.prefix == (None if target.prefix is None else tag + target.prefix)
