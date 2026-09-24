"""FASTA file writer — serializes models back to FASTA format."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import IO

from fastatacular._models import SequenceEntry
from fastatacular._parser import _KV_PATTERN, _parse_header_line, _ParsedHeader
from fastatacular.errors import FastaParseError, FastaWriteError

_SEQ_LINE_WIDTH = 60

_HEADER_FIELDS = (
    "identifier",
    "prefix",
    "accession",
    "entry_name",
    "description",
    "pname",
    "gname",
    "os_name",
    "ncbi_tax_id",
    "pe",
    "sv",
    "extra",
)


def _parse_raw_header(entry: SequenceEntry) -> _ParsedHeader | None:
    """Parse ``entry.raw_header``, or return None if it is empty or unparseable."""
    if not entry.raw_header:
        return None
    try:
        return _parse_header_line(">" + entry.raw_header, 0)
    except FastaParseError:
        return None


def _build_header_line(entry: SequenceEntry) -> str:
    """Reconstruct a FASTA description line for ``entry``.

    Priority: if ``raw_header`` is set and still parses to the entry's current
    header fields (identifier, description, ``gname``, ``extra``, ...), write it
    verbatim, so unedited entries round-trip byte-exact. If any of those fields
    was changed (e.g. with ``dataclasses.replace``), or ``raw_header`` is empty,
    rebuild the header from the structured fields so the edit is kept.

    When falling back to ``description`` (``pname`` unset), its ``KEY=value``
    text is not copied verbatim: the structured fields are written instead, and
    only keys that no structured field or ``extra`` covers are kept from it.
    A ``description`` that is still the one parsed from ``raw_header`` is not
    used at all: it is stale, and reading it would undo a cleared field.
    """
    parsed = _parse_raw_header(entry)
    if parsed is not None and all(getattr(parsed, name) == getattr(entry, name) for name in _HEADER_FIELDS):
        return f">{entry.raw_header}"

    description = entry.description
    if parsed is not None and description == parsed.description:
        description = None

    parts: list[str] = [entry.identifier]
    desc_extra: dict[str, str] = {}
    if entry.pname:
        parts.append(entry.pname)
    elif description:
        matches = list(_KV_PATTERN.finditer(description))
        name = description[: matches[0].start()].rstrip() if matches else description
        if name:
            parts.append(name)
        desc_extra = {m["key"]: m["val"].strip() for m in matches}

    if entry.os_name is not None:
        parts.append(f"OS={entry.os_name}")
    if entry.ncbi_tax_id is not None:
        parts.append(f"OX={entry.ncbi_tax_id}")
    if entry.gname is not None:
        parts.append(f"GN={entry.gname}")
    if entry.pe is not None:
        parts.append(f"PE={entry.pe}")
    if entry.sv is not None:
        parts.append(f"SV={entry.sv}")
    for k, v in entry.extra.items():
        parts.append(f"{k}={v}")
    typed = {
        "OS": entry.os_name,
        "OX": entry.ncbi_tax_id,
        "GN": entry.gname,
        "PE": entry.pe,
        "SV": entry.sv,
    }
    for k, v in desc_extra.items():
        if typed.get(k) is None and k not in entry.extra:
            parts.append(f"{k}={v}")

    return ">" + " ".join(parts)


def _write_entry(entry: SequenceEntry, out: IO[str], line_width: int) -> None:
    if not entry.identifier:
        raise FastaWriteError("SequenceEntry has an empty identifier")
    if not entry.sequence:
        raise FastaWriteError(f"SequenceEntry {entry.identifier!r} has an empty sequence")
    if any(c.isspace() for c in entry.identifier):
        raise FastaWriteError(f"SequenceEntry identifier {entry.identifier!r} contains whitespace")
    if any(c.isspace() for c in entry.sequence):
        raise FastaWriteError(f"SequenceEntry {entry.identifier!r} sequence contains whitespace")

    header = _build_header_line(entry)
    if "\n" in header or "\r" in header:
        raise FastaWriteError(f"SequenceEntry {entry.identifier!r} header contains a line break")

    seq = entry.sequence
    step = line_width if line_width > 0 else len(seq)
    lines = [seq[i : i + step] for i in range(0, len(seq), step)]
    if any(line[0] in ">;" for line in lines):
        # The reader would take such a line as a header or a comment.
        raise FastaWriteError(f"SequenceEntry {entry.identifier!r} sequence would start a line with '>' or ';'")
    out.write(header + "\n")
    for line in lines:
        out.write(line + "\n")


def write_fasta(
    entries: Iterable[SequenceEntry],
    dest: str | Path | IO[str],
    *,
    line_width: int = _SEQ_LINE_WIDTH,
) -> None:
    """Write a sequence of ``SequenceEntry`` objects to FASTA.

    ``dest`` may be a path or an already-opened text-mode file object.
    ``line_width`` controls sequence wrapping; pass ``0`` (or any value ``<= 0``)
    to emit each sequence on a single line.
    """
    if isinstance(dest, (str, Path)):
        with Path(dest).open("w", encoding="utf-8") as fh:
            for entry in entries:
                _write_entry(entry, fh, line_width)
    else:
        for entry in entries:
            _write_entry(entry, dest, line_width)


__all__ = ["write_fasta"]
