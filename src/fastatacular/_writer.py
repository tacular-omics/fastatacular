"""FASTA file writer — serializes models back to FASTA format."""

from __future__ import annotations

import io
from collections.abc import Iterable
from pathlib import Path
from typing import IO, TYPE_CHECKING

from fastatacular._compression import _KINDS, Compression, _check_compression, _DetachingText, _is_text_handle
from fastatacular._models import SequenceEntry
from fastatacular._parser import _MODULES, _parse_header_line, _ParsedHeader, _split_kv
from fastatacular.errors import FastaError, FastaParseError, FastaWriteError

if TYPE_CHECKING:
    import gzip

_SEQ_LINE_WIDTH = 60

# Fields a header line can carry. ``prefix``, ``accession`` and ``entry_name`` are not
# here: they come from the identifier and a rebuilt header cannot write them either, so
# they do not decide between ``raw_header`` and a rebuild. A decoy with a custom prefix
# (``XXX_pdb|1MBA|A``, accession ``1MBA_A``) is thus still written verbatim.
_HEADER_FIELDS = (
    "identifier",
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
        first_key, pairs = _split_kv(description)
        name = description[:first_key].rstrip() if first_key is not None else description
        if name:
            parts.append(name)
        desc_extra = dict(pairs)

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


_SUFFIX_KINDS = {".gz": "gz", ".bz2": "bz2", ".xz": "xz"}


def _write_kind(path: Path | None, compression: Compression) -> str | None:
    """The compression to write with: by suffix (``"infer"``; none for a handle), forced, or none."""
    if compression is None:
        return None
    if compression == "infer":
        return None if path is None else _SUFFIX_KINDS.get(path.suffix.lower())
    return _KINDS[compression]


def _compressor(kind: str, target: Path | IO[bytes], name: object) -> IO[str]:
    """A UTF-8 text writer compressing ``kind`` into a path or a binary handle.

    Closing it finishes the compressed stream; a handle ``target`` is left open.
    """

    try:
        if kind == "gz":
            import gzip

            # mtime=0 and filename="" so the same entries always give the same bytes
            # (gzip.open stamps the current time and the file name into the header).
            if isinstance(target, Path):
                binary = _owning_gzip(target)
            else:
                binary = gzip.GzipFile(filename="", fileobj=target, mode="wb", mtime=0)
        elif kind == "bz2":
            import bz2

            binary = bz2.BZ2File(target, mode="wb")
        else:
            import lzma

            binary = lzma.LZMAFile(target, mode="wb")
    except ImportError as e:
        raise FastaWriteError(
            f"Cannot write {name}: this Python has no {_MODULES[kind]} module",
            hint=f"Write uncompressed output, or use a Python built with {_MODULES[kind]} support",
        ) from e
    return io.TextIOWrapper(binary, encoding="utf-8", newline="\n")


def _owning_gzip(path: Path) -> gzip.GzipFile:
    """A reproducible gzip writer to ``path`` (no name or time in the header) that closes the file."""
    import gzip

    raw = path.open("wb")

    class _Gzip(gzip.GzipFile):
        def close(self) -> None:
            try:
                super().close()
            finally:
                raw.close()

    try:
        return _Gzip(filename="", fileobj=raw, mode="wb", mtime=0)
    except BaseException:
        raw.close()
        raise


def _open_for_write(path: Path, compression: Compression = "infer") -> IO[str]:
    """Open ``path`` for UTF-8 text output, compressed per ``compression`` (see :func:`write_fasta`)."""
    kind = _write_kind(path, compression)
    if kind is None:
        return path.open("w", encoding="utf-8", newline="\n")
    return _compressor(kind, path, path)


def _write_items(
    items: list[tuple[SequenceEntry, str, int]], dest: str | Path | IO[str] | IO[bytes], compression: Compression
) -> None:
    if isinstance(dest, (str, Path)):
        with _open_for_write(Path(dest), compression) as fh:
            _write_prepared(items, fh)
        return
    kind = _write_kind(None, compression)
    if kind is None:
        if _is_text_handle(dest):
            _write_prepared(items, dest)  # ty: ignore[invalid-argument-type]
        else:
            # Plain UTF-8 to a binary handle; detach (not close) so the handle stays open.
            with _DetachingText(dest, encoding="utf-8", newline="\n") as fh:  # ty: ignore[invalid-argument-type]
                _write_prepared(items, fh)
        return
    with _compressor(kind, dest, "the handle") as fh:  # ty: ignore[invalid-argument-type]
        _write_prepared(items, fh)


def _check_dest(dest: object, compression: Compression) -> None:
    """Reject a text handle for compressed output before anything is written."""
    if isinstance(dest, (str, Path)) or _write_kind(None, compression) is None:
        return
    if _is_text_handle(dest):
        err = FastaError(f"compression={compression!r} needs a binary handle, got a text-mode one")
        err.add_note("hint: open the file with mode 'wb', or pass a path")
        raise err


def write_fasta(
    entries: Iterable[SequenceEntry],
    dest: str | Path | IO[str] | IO[bytes],
    *,
    line_width: int = _SEQ_LINE_WIDTH,
    compression: Compression = "infer",
) -> None:
    """Write a sequence of ``SequenceEntry`` objects to FASTA.

    ``dest`` may be a path or an already-opened text or binary file object.
    ``line_width`` controls sequence wrapping; pass ``0`` (or any value ``<= 0``)
    to emit each sequence on a single line.

    Every entry is validated before anything is written: if one is unwritable,
    ``FastaWriteError`` is raised (its ``index`` names the entry), a path
    ``dest`` is not created or truncated, and nothing is written to a handle.

    Args:
        compression: ``"infer"`` (default) compresses a path ending in ``.gz``,
            ``.bz2`` or ``.xz`` (any case) with gzip, bzip2 or xz and writes plain
            text otherwise and to a handle (UTF-8 bytes to a binary handle).
            ``"gzip"``, ``"bz2"`` or ``"xz"`` forces that format whatever the
            suffix; a handle must then be binary (``"wb"``). ``None`` always writes
            plain text. A handle is never closed, even on error. gzip output has
            mtime 0 and no file name in its header, so the same entries give the
            same bytes.

    Raises:
        FastaError: An unknown ``compression``, or a text handle with an explicit one.
        FastaWriteError: An entry that cannot be written.
    """
    compression = _check_compression(compression)
    _check_dest(dest, compression)
    if not isinstance(line_width, int) or isinstance(line_width, bool):
        raise FastaWriteError(
            f"line_width must be an int, got {line_width!r}",
            hint="Pass an int; 0 or less writes each sequence on one line",
        )
    items = [(entry, *_prepare_entry(entry, i, line_width)) for i, entry in enumerate(entries)]
    _write_items(items, dest, compression)


__all__ = ["write_fasta"]
