"""FASTA file parser — converts text into model objects."""

from __future__ import annotations

import re
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from types import TracebackType
from typing import IO, Self

from fastatacular._models import SequenceEntry
from fastatacular.errors import FastaParseError

# Matches ``KEY=value`` pairs in UniProt-style headers. A key starts the
# description or follows whitespace, so ``Protein(EC=2.7.1)`` and
# ``[organism=Homo sapiens]`` are text, not keys. Value runs up to the next
# ``KEY=`` token or end-of-string, then trailing whitespace is trimmed.
_KV_PATTERN = re.compile(r"(?:^|(?<=\s))(?P<key>[A-Za-z_][A-Za-z0-9_]*)=(?P<val>.*?)(?=\s+[A-Za-z_][A-Za-z0-9_]*=|$)")

# UniProt FASTA identifier: ``db|ACCESSION|ENTRY_NAME`` (e.g. ``sp|P12345|EX_HUMAN``).
# The prefix is any run without ``|``, so decoy/contaminant tags such as
# ``Reverse_sp|...``, ``rev_sp|...`` and ``DECOY-0-sp|...`` keep their accession.
_UNIPROT_ID = re.compile(r"^(?P<prefix>[^|\s]+)\|(?P<accession>[^|]+)\|(?P<entry_name>[^|\s]+)$")

# NCBI-ish ``db|ID`` or ``db|ID|...`` identifier — accept the leading two fields.
_PIPE_ID = re.compile(r"^(?P<prefix>[^|\s]+)\|(?P<accession>[^|\s]+)(?:\|.*)?$")


@dataclass(slots=True)
class _ParsedHeader:
    """Mutable scratch space built up while parsing a description line."""

    identifier: str
    raw_header: str
    prefix: str | None = None
    accession: str | None = None
    entry_name: str | None = None
    description: str | None = None
    pname: str | None = None
    gname: str | None = None
    os_name: str | None = None
    ncbi_tax_id: int | None = None
    pe: int | None = None
    sv: int | None = None
    extra: dict[str, str] = field(default_factory=dict)


def _parse_header_line(line: str, line_no: int) -> _ParsedHeader:
    """Parse a single header line into structured fields."""
    if not line.startswith(">"):
        raise FastaParseError("Header line must start with '>'", line=line_no, context=line)

    raw = line[1:].rstrip("\r\n")
    stripped = raw.strip()
    if not stripped:
        raise FastaParseError(
            "Empty FASTA header", line=line_no, context=line, hint="Put an identifier right after '>'"
        )

    # Split on the first run of any whitespace (space or tab), as UniProt, BLAST,
    # Biopython and samtools do.
    identifier, *tail = stripped.split(maxsplit=1)
    rest = tail[0] if tail else ""

    header = _ParsedHeader(identifier=identifier, raw_header=raw)

    if m := _UNIPROT_ID.match(identifier):
        header.prefix = m["prefix"]
        header.accession = m["accession"]
        header.entry_name = m["entry_name"]
    elif m := _PIPE_ID.match(identifier):
        header.prefix = m["prefix"]
        header.accession = m["accession"]

    if not rest:
        return header

    first_match_start: int | None = None
    for m in _KV_PATTERN.finditer(rest):
        key = m["key"]
        value = m["val"].strip()
        if first_match_start is None:
            first_match_start = m.start()
        if key == "GN":
            header.gname = value
        elif key == "OS":
            header.os_name = value
        elif key == "OX":
            try:
                header.ncbi_tax_id = int(value)
            except ValueError:
                header.extra[key] = value
        elif key == "PE":
            try:
                header.pe = int(value)
            except ValueError:
                header.extra[key] = value
        elif key == "SV":
            try:
                header.sv = int(value)
            except ValueError:
                header.extra[key] = value
        else:
            header.extra[key] = value

    header.description = rest
    if first_match_start is None:
        header.pname = rest
    elif first_match_start > 0:
        header.pname = rest[:first_match_start].rstrip() or None

    return header


def _build_entry(header: _ParsedHeader, seq_chunks: list[str], header_line_no: int) -> SequenceEntry:
    sequence = "".join(chunk for chunk in seq_chunks if chunk)
    if not sequence:
        raise FastaParseError(
            f"Entry {header.identifier!r} has no sequence data",
            line=header_line_no,
            context=header.raw_header,
            hint="Add sequence lines after the header, or remove the header",
        )
    return SequenceEntry(
        identifier=header.identifier,
        sequence=sequence,
        prefix=header.prefix,
        accession=header.accession,
        entry_name=header.entry_name,
        description=header.description,
        pname=header.pname,
        gname=header.gname,
        os_name=header.os_name,
        ncbi_tax_id=header.ncbi_tax_id,
        pe=header.pe,
        sv=header.sv,
        extra=header.extra,
        raw_header=header.raw_header,
    )


def _iter_entries(fh: IO[str]) -> Iterator[SequenceEntry]:
    """Yield ``SequenceEntry`` objects from an open text-mode file."""
    header: _ParsedHeader | None = None
    seq_chunks: list[str] = []
    header_line_no: int = 0

    for line_no, line in enumerate(fh, start=1):
        if line_no == 1 and line.startswith("\ufeff"):
            # UTF-8 byte-order mark read through a plain ``utf-8`` text handle.
            line = line[1:]
        if not line or line.isspace():
            continue
        if line.startswith(";"):
            # Comment line (NCBI / legacy FASTA convention) — skip.
            continue
        if line.startswith(">"):
            if header is not None:
                yield _build_entry(header, seq_chunks, header_line_no)
            header = _parse_header_line(line, line_no)
            header_line_no = line_no
            seq_chunks = []
        else:
            if header is None:
                raise FastaParseError(
                    "Sequence data appears before any '>' header",
                    line=line_no,
                    context=line.rstrip("\n"),
                    hint="A FASTA file must start with a '>' header line",
                )
            seq_chunks.append("".join(line.split()))

    if header is not None:
        yield _build_entry(header, seq_chunks, header_line_no)


class FastaReader:
    """Iterate over a FASTA file lazily without loading the entire file.

    Use it as a context manager. A reader is single-pass: iterating it again
    continues from where the previous iteration stopped, so after a full pass
    a second ``for`` loop yields nothing. Open a new reader (or use
    ``read_fasta``) to read the file again.
    """

    def __init__(self, source: str | Path | IO[str]) -> None:
        self._source = source
        self._fh: IO[str] | None = None
        self._owns_fh = False

    def __enter__(self) -> Self:
        if isinstance(self._source, (str, Path)):
            self._fh = Path(self._source).open(encoding="utf-8-sig")
            self._owns_fh = True
        else:
            self._fh = self._source
            self._owns_fh = False
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        if self._owns_fh and self._fh is not None:
            self._fh.close()
        self._fh = None

    def __iter__(self) -> Iterator[SequenceEntry]:
        if self._fh is None:
            raise RuntimeError("FastaReader must be used as a context manager (`with FastaReader(...) as r:`)")
        return _iter_entries(self._fh)


def read_fasta(source: str | Path | IO[str]) -> list[SequenceEntry]:
    """Read an entire FASTA file into a list of ``SequenceEntry`` objects."""
    if isinstance(source, (str, Path)):
        with Path(source).open(encoding="utf-8-sig") as fh:
            return list(_iter_entries(fh))
    return list(_iter_entries(source))


__all__ = ["FastaReader", "read_fasta"]
