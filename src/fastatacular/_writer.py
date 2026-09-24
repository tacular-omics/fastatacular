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


def _build_header_line(entry: SequenceEntry) -> tuple[str, bool]:
    """Reconstruct a FASTA description line for ``entry``.

    Returns the line and whether it was rebuilt from the fields (``False``
    when ``raw_header`` was written verbatim).

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
        return f">{entry.raw_header}", False

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

    return ">" + " ".join(parts), True


_TYPED_FIELDS = (
    ("os_name", "OS"),
    ("ncbi_tax_id", "OX"),
    ("gname", "GN"),
    ("pe", "PE"),
    ("sv", "SV"),
)


def _check_rebuilt_header(entry: SequenceEntry, header: str, index: int) -> None:
    """Raise ``FastaWriteError`` if a rebuilt ``header`` would not read back to ``entry``.

    Free text holding a ``KEY=`` token (``os_name="Homo sapiens GN=X"``), an
    ``extra`` key that is not a single word, or an ``extra`` key that a typed
    field owns (``extra={"OS": ...}``) would be split up differently on read.
    Only fields the entry sets are checked: keys taken from a ``description``
    fallback are expected to appear.
    """
    parsed = _parse_header_line(header, 0)
    if entry.pname and parsed.pname != entry.pname:
        raise FastaWriteError(
            f"SequenceEntry {entry.identifier!r} pname {entry.pname!r} would read back as {parsed.pname!r}",
            index=index,
            hint="Remove 'KEY=' text and leading/trailing whitespace from pname",
        )
    for name, key in _TYPED_FIELDS:
        value = getattr(entry, name)
        if value is not None and getattr(parsed, name) != value:
            raise FastaWriteError(
                f"SequenceEntry {entry.identifier!r} {name} {value!r} would read back as {getattr(parsed, name)!r}",
                index=index,
                hint=f"Remove 'KEY=' text and leading/trailing whitespace from {name} (written as {key}=)",
            )
    for key, value in entry.extra.items():
        if key not in parsed.extra or parsed.extra[key] != value:
            raise FastaWriteError(
                f"SequenceEntry {entry.identifier!r} extra[{key!r}] = {value!r} would not read back as written",
                index=index,
                hint=(
                    "extra keys must be one word matching [A-Za-z_][A-Za-z0-9_]* and not OS/OX/GN/PE/SV "
                    "(use the typed field); values must not contain 'KEY=' text or edge whitespace"
                ),
            )


def _prepare_entry(entry: SequenceEntry, index: int, line_width: int) -> tuple[str, int]:
    """Validate ``entry`` and return its header line and sequence line step.

    Raises ``FastaWriteError`` (with ``index``) if the entry cannot be written
    as FASTA that reads back to the same entry.
    """
    if not entry.identifier:
        raise FastaWriteError(
            "SequenceEntry has an empty identifier", index=index, hint="Set identifier to a non-empty token"
        )
    if not entry.sequence:
        raise FastaWriteError(
            f"SequenceEntry {entry.identifier!r} has an empty sequence",
            index=index,
            hint="A FASTA entry needs at least one residue",
        )
    if any(c.isspace() for c in entry.identifier):
        raise FastaWriteError(
            f"SequenceEntry identifier {entry.identifier!r} contains whitespace",
            index=index,
            hint="The identifier ends at the first whitespace; put the rest in pname or description",
        )
    if any(c.isspace() for c in entry.sequence):
        raise FastaWriteError(
            f"SequenceEntry {entry.identifier!r} sequence contains whitespace",
            index=index,
            hint="Pass the residues only; the writer wraps the sequence itself",
        )

    header, rebuilt = _build_header_line(entry)
    if "\n" in header or "\r" in header:
        raise FastaWriteError(
            f"SequenceEntry {entry.identifier!r} header contains a line break",
            index=index,
            hint="A FASTA header is one line; remove the \\n or \\r from the header fields",
        )
    if rebuilt:
        _check_rebuilt_header(entry, header, index)

    seq = entry.sequence
    step = line_width if line_width > 0 else len(seq)
    if any(seq[i] in ">;" for i in range(0, len(seq), step)):
        # The reader would take such a line as a header or a comment.
        raise FastaWriteError(
            f"SequenceEntry {entry.identifier!r} sequence would start a line with '>' or ';'",
            index=index,
            hint="Change line_width, or remove '>' / ';' from the sequence",
        )
    return header, step


def _write_prepared(items: list[tuple[SequenceEntry, str, int]], out: IO[str]) -> None:
    for entry, header, step in items:
        seq = entry.sequence
        out.write(header + "\n")
        for i in range(0, len(seq), step):
            out.write(seq[i : i + step] + "\n")


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

    Every entry is validated before anything is written: if one is unwritable,
    ``FastaWriteError`` is raised (its ``index`` names the entry), a path
    ``dest`` is not created or truncated, and nothing is written to a handle.
    """
    if not isinstance(line_width, int) or isinstance(line_width, bool):
        raise FastaWriteError(
            f"line_width must be an int, got {line_width!r}",
            hint="Pass an int; 0 or less writes each sequence on one line",
        )
    items = [(entry, *_prepare_entry(entry, i, line_width)) for i, entry in enumerate(entries)]
    if isinstance(dest, (str, Path)):
        with Path(dest).open("w", encoding="utf-8") as fh:
            _write_prepared(items, fh)
    else:
        _write_prepared(items, dest)


__all__ = ["write_fasta"]
