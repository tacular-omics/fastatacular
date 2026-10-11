"""Flat ``dict`` records for building data frames (pandas, polars) without a dependency."""

from __future__ import annotations

import numbers
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import IO, TYPE_CHECKING, Any, cast

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


_STR_FIELDS = ("identifier", "prefix", "accession", "entry_name", "pname", "gname", "os_name", "description")
_INT_FIELDS = ("ncbi_tax_id", "pe", "sv")
# Fields that make up a rebuilt description; with none of them, ``description`` is used verbatim.
_DESC_FIELDS = ("pname", "gname", "os_name", "ncbi_tax_id", "pe", "sv", "extra")


def _absent(value: object) -> bool:
    """``None``, float NaN and pandas/polars missing markers (``NA``, ``NaT``) count as absent."""
    if value is None:
        return True
    try:
        return bool(value != value)
    except TypeError:  # pandas.NA: ``NA != NA`` is NA, whose truth value raises
        return True


def _bad(index: int, key: str, value: object, what: str) -> Exception:
    from fastatacular.errors import FastaWriteError

    return FastaWriteError(f"{key} must be {what}, got {value!r}", index=index)


def _as_int(value: object, index: int, key: str) -> int:
    """An int from an int, an integral float (``9606.0`` from a frame column with gaps) or digits."""
    if isinstance(value, bool):
        raise _bad(index, key, value, "an int")
    if isinstance(value, numbers.Integral):
        return int(value)
    if isinstance(value, numbers.Real) and (number := float(value)).is_integer():
        return int(number)
    if isinstance(value, str):
        try:
            return int(value.strip())
        except ValueError:
            pass
    raise _bad(index, key, value, "an int")


def _as_str(value: object, index: int, key: str) -> str:
    if isinstance(value, str):
        return value
    if isinstance(value, numbers.Real) and not isinstance(value, bool):
        # A column of numeric ids reads back as numbers; 12.0 means "12".
        return str(_as_int(value, index, key))
    raise _bad(index, key, value, "a str")


def _as_extra(value: object, index: int) -> dict[str, str]:
    from fastatacular._parser import _split_kv

    if isinstance(value, Mapping):
        return {str(k): str(v) for k, v in value.items()}
    if not isinstance(value, str):
        raise _bad(index, "extra", value, "space-separated KEY=value pairs or a dict")
    text = value.strip()
    if not text:
        return {}
    first_key, pairs = _split_kv(text)
    if first_key != 0:
        raise _bad(index, "extra", value, "space-separated KEY=value pairs")
    return dict(pairs)


def _record_entry(record: Mapping[str, object], index: int) -> SequenceEntry:
    from fastatacular._models import SequenceEntry
    from fastatacular._parser import _build_entry, _parse_header_line, _split_identifier
    from fastatacular._writer import _build_header_line, _prepare_entry
    from fastatacular.errors import FastaWriteError

    if not isinstance(record, Mapping):
        raise _bad(index, "the record", record, "a mapping")
    values: dict[str, object] = {k: v for k, v in record.items() if k in RECORD_KEYS and not _absent(v)}
    fields: dict[str, Any] = {}
    for key in _STR_FIELDS:
        if key in values:
            fields[key] = _as_str(values[key], index, key)
    for key in _INT_FIELDS:
        if key in values:
            fields[key] = _as_int(values[key], index, key)
    extra = _as_extra(values["extra"], index) if "extra" in values else {}
    sequence = _as_str(values["sequence"], index, "sequence") if "sequence" in values else ""

    # Identifier: the column while it still splits into prefix|accession|entry_name, else
    # rebuilt from those columns, so editing them changes the output.
    identifier = fields.get("identifier")
    parts = (fields.get("prefix"), fields.get("accession"), fields.get("entry_name"))
    if parts != (None, None, None) and (identifier is None or (_split_identifier(identifier) or (None,) * 3) != parts):
        if parts[1] is None:
            raise FastaWriteError(
                f"prefix/entry_name {parts[0]!r}/{parts[2]!r} without an accession "
                f"do not match identifier {identifier!r}",
                index=index,
                hint="Set accession, or clear prefix/accession/entry_name to use the identifier column",
            )
        identifier = "|".join(p for p in parts if p is not None)
        if _split_identifier(identifier) != parts:
            raise FastaWriteError(
                f"prefix/accession/entry_name {parts!r} do not form an identifier "
                f"that reads back the same ({identifier!r})",
                index=index,
                hint="Set the identifier column and clear prefix/accession/entry_name",
            )
    if not identifier:
        raise FastaWriteError("record has no identifier", index=index, hint="Set identifier, or prefix and accession")

    structured = any(fields.get(k) is not None for k in _DESC_FIELDS) or bool(extra)
    entry = SequenceEntry(
        identifier=identifier,
        sequence=sequence,
        pname=fields.get("pname"),
        gname=fields.get("gname"),
        os_name=fields.get("os_name"),
        ncbi_tax_id=fields.get("ncbi_tax_id"),
        pe=fields.get("pe"),
        sv=fields.get("sv"),
        extra=extra,
        description=None if structured else fields.get("description"),
    )
    _prepare_entry(entry, index, 0)  # identifier, sequence and header checks, with the index
    description = entry.description
    if structured or not description:
        header, _ = _build_header_line(entry)
    else:
        header = f">{identifier} {description}"
    parsed = _parse_header_line(header, 0)
    if description and not structured and parsed.description != description:
        raise FastaWriteError(
            f"description {description!r} would read back as {parsed.description!r}",
            index=index,
            hint="Remove leading/trailing whitespace from description",
        )
    return _build_entry(parsed, [sequence])


def from_records(records: Iterable[Mapping[str, object]]) -> list[SequenceEntry]:
    """Build :class:`SequenceEntry` objects from flat records, the inverse of :func:`to_records`.

    Writes a data frame back to FASTA::

        write_fasta(from_records(df.to_dict("records")), "out.fasta")

    Each record is a mapping with keys from :data:`RECORD_KEYS`. Missing keys and missing
    values (``None``, float NaN, ``pandas.NA``) count as absent; other keys are ignored.

    The header is rebuilt from the columns, so an edited column changes the output:

    * identifier: ``prefix|accession|entry_name`` (the parts that are set) when any of
      those columns is set and the ``identifier`` column does not already split into
      exactly them; otherwise the ``identifier`` column.
    * description: ``pname``, then ``OS=``, ``OX=``, ``GN=``, ``PE=``, ``SV=`` (UniProt
      order), then the ``extra`` pairs. A record with none of those but a
      ``description`` uses the description verbatim.

    ``raw_header`` and ``length`` are ignored. ``ncbi_tax_id``, ``pe`` and ``sv`` accept
    ints, integral floats (``9606.0`` from a frame column with gaps writes ``OX=9606``)
    and digit strings. ``extra`` is space-separated ``KEY=value`` text or a dict. A missing
    ``sequence`` gives ``""``. The returned entries are what reading the written FASTA
    gives back (``raw_header`` is the rebuilt header).

    Raises:
        FastaWriteError: A record with a bad value (non-integer ``pe``, no identifier,
            text that would not read back as written, ...); ``index`` is its position.
    """
    return [_record_entry(record, i) for i, record in enumerate(records)]
