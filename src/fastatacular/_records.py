"""Flat ``dict`` records for building data frames (pandas, polars) without a dependency."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import IO, TYPE_CHECKING, cast

if TYPE_CHECKING:
    from fastatacular._models import SequenceEntry

# Column order of every record; the keys never change between entries or releases.
RECORD_KEYS: tuple[str, ...] = (
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
    "description",
    "extra",
    "raw_header",
    "length",
    "sequence",
)


def entry_to_record(entry: SequenceEntry) -> dict[str, str | int | None]:
    extra = " ".join(f"{k}={v}" for k, v in entry.extra.items()) or None
    return {
        "identifier": entry.identifier,
        "prefix": entry.prefix,
        "accession": entry.accession,
        "entry_name": entry.entry_name,
        "pname": entry.pname,
        "gname": entry.gname,
        "os_name": entry.os_name,
        "ncbi_tax_id": entry.ncbi_tax_id,
        "pe": entry.pe,
        "sv": entry.sv,
        "description": entry.description,
        "extra": extra,
        "raw_header": entry.raw_header,
        "length": len(entry.sequence),
        "sequence": entry.sequence,
    }


def to_records(source: str | Path | IO[str] | Iterable[SequenceEntry]) -> list[dict[str, str | int | None]]:
    """Return one flat ``dict`` per FASTA entry, ready for ``pandas.DataFrame(records)``.

    ``source`` is a path (``str`` or ``Path``, may be compressed), an open text handle,
    or an iterable of :class:`SequenceEntry` (e.g. a :class:`FastaReader` or a list).
    Every record has the same keys, in this order:

    ============== ========= =======================================================
    key            type      value
    ============== ========= =======================================================
    identifier     str       text after ``>`` up to the first whitespace
    prefix         str|None  ``sp`` in ``sp|P12345|NAME``
    accession      str|None  ``P12345``
    entry_name     str|None  ``NAME`` (three-field ``db|ACC|NAME`` ids only)
    pname          str|None  protein name (description before the first ``KEY=``)
    gname          str|None  ``GN=``
    os_name        str|None  ``OS=``
    ncbi_tax_id    int|None  ``OX=``
    pe             int|None  ``PE=``
    sv             int|None  ``SV=``
    description    str|None  everything after the identifier
    extra          str|None  other ``KEY=value`` pairs, space-separated
    raw_header     str       the header line without ``>``
    length         int       ``len(sequence)``
    sequence       str       the residues
    ============== ========= =======================================================

    The package does not use or require pandas or polars; the records are plain dicts.
    """
    from fastatacular._parser import FastaReader

    if isinstance(source, (str, Path)) or hasattr(source, "read"):
        with FastaReader(cast("str | Path | IO[str]", source)) as reader:
            return [entry_to_record(e) for e in reader]
    return [entry_to_record(e) for e in cast("Iterable[SequenceEntry]", source)]
