"""FASTA file writer — serializes models back to FASTA format."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import IO

from fastatacular._models import SequenceEntry
from fastatacular._parser import _KV_PATTERN
from fastatacular.errors import FastaWriteError

_SEQ_LINE_WIDTH = 60


def _build_header_line(entry: SequenceEntry) -> str:
    """Reconstruct a FASTA description line for ``entry``.

    Priority: if ``raw_header`` was preserved during parsing, round-trip it
    exactly. Otherwise rebuild from the structured fields.

    When falling back to ``description`` (``pname`` unset), its ``KEY=value``
    text is not copied verbatim: the structured fields are written instead, and
    only keys that no structured field or ``extra`` covers are kept from it.
    """
    if entry.raw_header:
        return f">{entry.raw_header}"

    parts: list[str] = [entry.identifier]
    desc_extra: dict[str, str] = {}
    if entry.pname:
        parts.append(entry.pname)
    elif entry.description:
        matches = list(_KV_PATTERN.finditer(entry.description))
        name = entry.description[: matches[0].start()].rstrip() if matches else entry.description
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

    out.write(_build_header_line(entry) + "\n")

    seq = entry.sequence
    if line_width <= 0:
        out.write(seq + "\n")
        return
    for i in range(0, len(seq), line_width):
        out.write(seq[i : i + line_width] + "\n")


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
