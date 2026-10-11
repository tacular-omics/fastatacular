"""Flat ``dict`` records for building data frames (pandas, polars) without a dependency."""

from __future__ import annotations

import numbers
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import IO, TYPE_CHECKING, Any, cast

if TYPE_CHECKING:
    import re

    from fastatacular._models import SequenceEntry
    from fastatacular._parser import _ParsedHeader

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
_ID_PARTS = ("prefix", "accession", "entry_name")
# Fields that make up a description besides its free text.
_DESC_FIELDS = ("pname", "gname", "os_name", "ncbi_tax_id", "pe", "sv", "extra")
# Typed fields and their header keys, in the order a rebuilt or appended header writes them.
_TYPED_KEYS = (("os_name", "OS"), ("ncbi_tax_id", "OX"), ("gname", "GN"), ("pe", "PE"), ("sv", "SV"))
# Header fields compared against the parse of ``raw_header`` to find edited columns.
_HEADER_KEYS = ("identifier", *_ID_PARTS, *_DESC_FIELDS, "description")


def _absent(value: object) -> bool:
    """``None``, ``""``, float NaN and pandas/polars missing markers (``NA``, ``NaT``) count as absent."""
    if value is None or (isinstance(value, str) and not value):
        return True
    try:
        return bool(value != value)
    except TypeError:  # pandas.NA: ``NA != NA`` is NA, whose truth value raises
        return True


def _error(message: str, index: int, hint: str | None = None) -> Exception:
    from fastatacular.errors import FastaWriteError

    return FastaWriteError(message, index=index, hint=hint)


def _bad(index: int, key: str, value: object, what: str) -> Exception:
    return _error(f"{key} must be {what}, got {value!r}", index)


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


def _as_extra(value: object, index: int) -> dict[str, str] | None:
    from fastatacular._parser import _split_kv

    if isinstance(value, Mapping):
        return {str(k): str(v) for k, v in value.items()} or None
    if not isinstance(value, str):
        raise _bad(index, "extra", value, "space-separated KEY=value pairs or a dict")
    first_key, pairs = _split_kv(value.strip())
    if first_key != 0:
        raise _bad(index, "extra", value, "space-separated KEY=value pairs")
    extra: dict[str, str] = {}
    for key, val in pairs:
        if key in extra:
            raise _bad(index, "extra", value, "pairs with distinct keys")
        extra[key] = val
    return extra


def _clean(record: Mapping[str, object], index: int) -> dict[str, Any]:
    """The record's known columns, converted; a missing value is ``None``, a missing key is left out."""
    if not isinstance(record, Mapping):
        raise _bad(index, "the record", record, "a mapping")
    values: dict[str, Any] = {}
    for key in RECORD_KEYS:
        if key not in record or key == "length":
            continue
        value = record[key]
        if _absent(value):
            values[key] = None
        elif key in _INT_FIELDS:
            values[key] = _as_int(value, index, key)
        elif key == "extra":
            values[key] = _as_extra(value, index)
        else:
            values[key] = _as_str(value, index, key)
    return values


def _fields_of(parsed: _ParsedHeader) -> dict[str, Any]:
    """Header fields of a parsed line, normalised like :func:`_clean` (``""`` is ``None``)."""
    fields = {key: getattr(parsed, key) for key in _HEADER_KEYS if key != "extra"}
    fields = {key: None if value == "" else value for key, value in fields.items()}
    fields["extra"] = dict(parsed.extra) or None
    return fields


def _parse_description(description: str | None) -> dict[str, Any]:
    from fastatacular._parser import _parse_header_line

    return _fields_of(_parse_header_line(f">x {description or ''}", 0))


def _parts_of(identifier: str) -> tuple[str | None, ...]:
    from fastatacular._parser import _split_identifier

    return _split_identifier(identifier) or (None, None, None)


def _check_identifier(identifier: str | None, index: int) -> str:
    if not identifier:
        raise _error("record has no identifier", index, "Set identifier, or prefix and accession")
    if any(c.isspace() for c in identifier):
        raise _error(
            f"identifier {identifier!r} contains whitespace",
            index,
            "The identifier ends at the first whitespace; put the rest in pname or description",
        )
    return identifier


def _conflict(index: int, what: str, columns: list[str], hint: str) -> Exception:
    return _error(f"{what} and {', '.join(columns)} were both set and disagree", index, hint)


def _substitute_parts(identifier: str, want: tuple[str | None, ...]) -> str | None:
    """``identifier`` with its prefix/accession/entry_name tokens replaced by ``want``.

    Inverts :func:`_split_identifier` for each identifier style, so every other pipe field
    and the layout stay: ``gi|N|ref|ACC|`` keeps ``gi|N|ref|``, ``pdb|1MBA|A`` takes an
    accession ``ENTRY_CHAIN`` back apart. ``None`` if ``want`` cannot be written in this style.
    """
    from fastatacular._parser import _EMPTY_MIDDLE_DBS, _EMPTY_MIDDLE_ID, _GI_ID, _PDB_ID, _PIPE_ID, _UNIPROT_ID

    prefix, accession, name = want
    if prefix is None or accession is None:
        return None
    if _parts_of(identifier) == (None, None, None):  # plain id: build db|ACC[|NAME]
        return "|".join(p for p in want if p is not None)
    edits: list[tuple[int, int, str]] = []

    def entry_chain(m: re.Match[str], entry: str, chain: str | None) -> None:
        # accession ENTRY_CHAIN back into the entry and chain fields
        if chain is not None and m[chain]:
            head, sep, tail = accession.rpartition("_")
            edits.extend([(*m.span(entry), head if sep else accession), (*m.span(chain), tail if sep else "")])
        else:
            edits.append((*m.span(entry), accession))

    if m := _PDB_ID.match(identifier):
        if name is not None:
            return None
        edits.append((*m.span("prefix"), prefix))
        entry_chain(m, "entry", "chain")
    elif m := _UNIPROT_ID.match(identifier):
        edits += [(*m.span("prefix"), prefix), (*m.span("accession"), accession), (*m.span("entry_name"), name or "")]
    elif m := _EMPTY_MIDDLE_ID.match(identifier):
        if name is not None:
            return None
        edits += [(*m.span("prefix"), prefix), (*m.span("accession"), accession)]
    elif (m := _GI_ID.match(identifier)) and (m["accession"] or (m["db"].lower() in _EMPTY_MIDDLE_DBS and m["name"])):
        edits.append((*m.span("prefix"), prefix))
        if not m["accession"]:  # gi|N|pir||ENTRY: the entry is the fifth field
            if name is not None:
                return None
            edits.append((*m.span("name"), accession))
        elif m["db"].lower() == "pdb":
            if name is not None:
                return None
            entry_chain(m, "accession", "name" if m["name"] is not None else None)
        else:
            edits.append((*m.span("accession"), accession))
            if m["name"] is not None:
                edits.append((*m.span("name"), name or ""))
            elif name is not None:
                edits.append((m.end("accession"), m.end("accession"), "|" + name))
    else:  # db|ID or db|ID|...
        m = cast("re.Match[str]", _PIPE_ID.match(identifier))
        edits += [(*m.span("prefix"), prefix), (*m.span("accession"), accession)]
        if name is not None:
            if m.end("accession") != len(identifier):
                return None
            edits.append((len(identifier), len(identifier), "|" + name))
    for start, end, text in sorted(edits, reverse=True):
        identifier = identifier[:start] + text + identifier[end:]
    return identifier


def _edited_identifier(values: dict[str, Any], old: dict[str, Any], changed: set[str], index: int) -> str:
    parts_changed = [k for k in _ID_PARTS if k in changed]
    want = tuple(values[k] if k in values else old[k] for k in _ID_PARTS)
    if "identifier" in changed:
        identifier = _check_identifier(values["identifier"], index)
        if parts_changed and _parts_of(identifier) != want:
            raise _conflict(
                index, f"identifier {identifier!r}", parts_changed, "Edit the identifier or its parts, not both"
            )
        return identifier
    if not parts_changed:
        return old["identifier"]
    identifier = _substitute_parts(old["identifier"], want)
    if identifier is None or _parts_of(identifier) != want:
        raise _error(
            f"prefix/accession/entry_name {want!r} cannot be written into identifier {old['identifier']!r}",
            index,
            "Edit the identifier column instead",
        )
    return identifier


def _key_values(fields: dict[str, Any], index: int) -> dict[str, str]:
    """The ``KEY=value`` pairs a description needs for ``fields`` (typed keys, then ``extra``)."""
    pairs = {key: str(fields[name]) for name, key in _TYPED_KEYS if fields[name] is not None}
    for key, value in (fields["extra"] or {}).items():
        if key in pairs:
            raise _error(f"{key}= is set both in extra and in its typed column", index, f"Set {key}= in one place")
        pairs[key] = value
    return pairs


def _edit_description(text: str, old: dict[str, Any], new: dict[str, Any], index: int) -> str | None:
    """``text`` with the changed fields edited in place; order, spacing and other text untouched.

    The value of a changed ``KEY=`` (its first occurrence) is replaced, a removed key loses
    every occurrence, a new key is appended, and a changed ``pname`` replaces the text
    before the first key.
    """
    from fastatacular._parser import _KV_PATTERN

    matches = list(_KV_PATTERN.finditer(text))
    first = matches[0].start() if matches else len(text)
    name = text[:first].rstrip() if new["pname"] == old["pname"] else (new["pname"] or "")
    out = [name]
    prev = len(text[:first].rstrip())
    old_kv, new_kv = _key_values(old, index), _key_values(new, index)
    done: set[str] = set()
    for m in matches:
        key, gap = m["key"], text[prev : m.start()]
        prev = m.end()
        token = text[m.start() : m.end()]
        if key in old_kv and key not in new_kv:
            continue  # removed
        if key not in done and key in new_kv and new_kv[key] != old_kv.get(key):
            raw = m["val"]
            lead, core = raw[: len(raw) - len(raw.lstrip())], raw.strip()
            token = f"{key}={lead}{new_kv[key]}{raw[len(lead) + len(core) :]}"
        done.add(key)
        out.append((gap or " ") + token)
    out.append(text[prev:])
    out += [f" {key}={value}" for key, value in new_kv.items() if key not in done and key not in old_kv]
    edited = "".join(out).lstrip()
    got = _parse_description(edited)
    bad = [k for k in _DESC_FIELDS if got[k] != new[k]]
    if bad:
        raise _error(
            f"{', '.join(bad)} would not read back as written in description {edited!r}",
            index,
            "Remove 'KEY=' text and leading/trailing whitespace from the field",
        )
    return edited or None


def _edited_header(values: dict[str, Any], raw: str, index: int) -> tuple[str, str | None]:
    """Header for a record with ``raw_header``: verbatim unless a column differs from its parse."""
    from fastatacular._parser import _parse_header_line
    from fastatacular.errors import FastaParseError

    try:
        old = _fields_of(_parse_header_line(">" + raw, 0))
    except FastaParseError as e:
        raise _error(f"raw_header {raw!r} is not a FASTA header", index, "Clear raw_header") from e
    changed = {k for k in _HEADER_KEYS if k in values and values[k] != old[k]}
    identifier = _edited_identifier(values, old, changed, index)

    struct = [k for k in _DESC_FIELDS if k in changed]
    if "description" in changed:
        description = values["description"]
        if struct:
            got = _parse_description(description)
            bad = [k for k in struct if got[k] != values[k]]
            if bad:
                raise _conflict(index, "description", bad, "Edit the description or the structured columns, not both")
    elif struct:
        new = {k: values[k] if k in values else old[k] for k in _DESC_FIELDS}
        description = _edit_description(old["description"] or "", old, new, index)
    else:
        description = old["description"]

    lead = len(raw) - len(raw.lstrip())
    after = raw[lead + len(old["identifier"]) :]
    if description == old["description"]:
        return f">{raw[:lead]}{identifier}{after}", description
    sep = after[: len(after) - len(after.lstrip())] if old["description"] else " "
    return f">{raw[:lead]}{identifier}{sep + description if description else ''}", description


def _built_header(values: dict[str, Any], index: int) -> tuple[str, str | None]:
    """Header for a hand-built record (no ``raw_header``)."""
    from fastatacular._models import SequenceEntry
    from fastatacular._writer import _build_header_line, _check_rebuilt_header

    identifier = values.get("identifier")
    parts = tuple(values.get(k) for k in _ID_PARTS)
    if parts != (None, None, None):
        if identifier is None:
            if parts[1] is None:
                raise _error("prefix/entry_name without an accession", index, "Set accession or identifier")
            identifier = "|".join(p for p in parts if p is not None)
            if _parts_of(identifier) != parts:
                raise _error(
                    f"prefix/accession/entry_name {parts!r} do not form an identifier "
                    f"that reads back the same ({identifier!r})",
                    index,
                    "Set the identifier column and clear prefix/accession/entry_name",
                )
        elif _parts_of(identifier) != parts:
            bad = [
                k for k, a, b in zip(_ID_PARTS, parts, _parts_of(identifier), strict=True) if a is not None and a != b
            ]
            if bad:
                raise _conflict(index, f"identifier {identifier!r}", bad, "Set the identifier or its parts, not both")
    identifier = _check_identifier(identifier, index)

    description = values.get("description")
    struct = [k for k in _DESC_FIELDS if values.get(k) is not None]
    if description is not None:
        if struct:
            got = _parse_description(description)
            bad = [k for k in struct if got[k] != values[k]]
            if bad:
                raise _conflict(index, "description", bad, "Set the description or the structured columns, not both")
        return f">{identifier} {description}", description
    if not struct:
        return f">{identifier}", None
    entry = SequenceEntry(
        identifier=identifier,
        sequence="",
        pname=values.get("pname"),
        gname=values.get("gname"),
        os_name=values.get("os_name"),
        ncbi_tax_id=values.get("ncbi_tax_id"),
        pe=values.get("pe"),
        sv=values.get("sv"),
        extra=values.get("extra") or {},
    )
    header, _ = _build_header_line(entry)
    _check_rebuilt_header(entry, header, index)
    return header, header.partition(" ")[2] or None


def _record_entry(record: Mapping[str, object], index: int) -> SequenceEntry:
    from fastatacular._parser import _build_entry, _parse_header_line
    from fastatacular._writer import _prepare_entry

    values = _clean(record, index)
    raw = values.get("raw_header")
    header, description = _edited_header(values, raw, index) if raw is not None else _built_header(values, index)
    if "\n" in header or "\r" in header:
        raise _error(f"header {header!r} contains a line break", index, "A FASTA header is one line")
    parsed = _parse_header_line(header, 0)
    if parsed.description != description:
        raise _error(
            f"description {description!r} would read back as {parsed.description!r}",
            index,
            "Remove leading/trailing whitespace from description",
        )
    entry = _build_entry(parsed, [values.get("sequence") or ""])
    _prepare_entry(entry, index, 0)  # identifier, sequence and header checks, with the index
    return entry


def from_records(records: Iterable[Mapping[str, object]]) -> list[SequenceEntry]:
    """Build :class:`SequenceEntry` objects from flat records, the inverse of :func:`to_records`.

    Writes a data frame back to FASTA::

        write_fasta(from_records(df.to_dict("records")), "out.fasta")

    Each record is a mapping with keys from :data:`RECORD_KEYS`; other keys and ``length``
    are ignored. ``None``, ``""``, float NaN and ``pandas.NA`` are missing values.

    A record with ``raw_header`` (as from :func:`to_records`) is compared with the parse
    of that header, so only edited columns change the output:

    * Nothing edited: ``raw_header`` is written verbatim (byte-identical round trip).
    * ``identifier`` edited: used verbatim. Otherwise an edited ``prefix``, ``accession``
      or ``entry_name`` replaces only that token in the original identifier, keeping its
      other pipe fields (``gi|N|ref|ACC|``, ``tr|ACC|NAME|more``; a ``pdb`` accession
      ``ENTRY_CHAIN`` goes back into ``ENTRY|CHAIN``).
    * ``description`` edited: used verbatim. Otherwise edited ``pname``, ``gname``,
      ``os_name``, ``ncbi_tax_id``, ``pe``, ``sv`` or ``extra`` edit the original
      description in place: a changed ``KEY=`` gets the new value, a removed one is
      deleted, a new one is appended, ``pname`` replaces the text before the first key;
      order, spacing and other text stay.
    * A key left out of the record keeps its original value. After an edit the
      identifier and description keep their original separator.

    A record without ``raw_header`` is built from its columns: the ``identifier`` (or
    ``prefix|accession|entry_name``), then ``description`` verbatim, or, without one,
    ``pname``, ``OS=``, ``OX=``, ``GN=``, ``PE=``, ``SV=`` and the ``extra`` pairs.

    ``ncbi_tax_id``, ``pe`` and ``sv`` accept ints, integral floats (``9606.0`` from a frame
    column with gaps) and digit strings. ``extra`` is space-separated ``KEY=value`` text or
    a dict. A missing ``sequence`` gives ``""``. The returned entries are what reading the
    written FASTA gives back.

    Raises:
        FastaWriteError: A bad value (non-integer ``pe``, no identifier, text that would
            not read back as written, ...), or an identifier and its parts, or a
            description and the structured columns, that were both set and disagree;
            ``index`` is the record's position.
    """
    return [_record_entry(record, i) for i, record in enumerate(records)]
