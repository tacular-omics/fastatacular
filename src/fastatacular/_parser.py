"""FASTA file parser — converts text into model objects."""

from __future__ import annotations

import io
import re
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from types import TracebackType
from typing import IO, Any, Self

from fastatacular._models import SequenceEntry
from fastatacular.errors import FastaError, FastaParseError

# Matches ``KEY=value`` pairs in UniProt-style headers. A key starts the
# description or follows whitespace, so ``Protein(EC=2.7.1)`` and
# ``[organism=Homo sapiens]`` are text, not keys. Value runs up to the next
# ``KEY=`` token or end-of-string, then trailing whitespace is trimmed.
_KV_PATTERN = re.compile(r"(?:^|(?<=\s))(?P<key>[A-Za-z_][A-Za-z0-9_]*)=(?P<val>.*?)(?=\s+[A-Za-z_][A-Za-z0-9_]*=|$)")

# Characters of a ``KEY`` in ``KEY=value`` (``[A-Za-z0-9_]``; it may not start with a digit).
_KEY_CHARS = frozenset("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_")


def _split_kv(text: str) -> tuple[int | None, list[tuple[str, str]]]:
    """Return the start of the first ``KEY=`` token and the ``(key, value)`` pairs.

    Same result as ``_KV_PATTERN.finditer`` (keys start the text or follow
    whitespace; a value runs to the next key and is stripped), but found by
    looking back from each ``=`` instead of trying the lazy pattern at every
    character.
    """
    eq = text.find("=")
    if eq == -1:
        return None, []
    if "\n" in text:
        # ``.`` in ``_KV_PATTERN`` stops at a line break; keep its exact behaviour.
        matches = list(_KV_PATTERN.finditer(text))
        return (matches[0].start() if matches else None), [(m["key"], m["val"].strip()) for m in matches]
    keys: list[tuple[int, int]] = []  # (key start, index of its '=')
    key_chars = _KEY_CHARS
    while eq != -1:
        # Common case: the token before '=' runs back to a space (or the start).
        start = text.rfind(" ", 0, eq) + 1
        token = text[start:eq]
        if token.isascii() and token.isidentifier():
            keys.append((start, eq))
        elif token.split() != [token]:
            # Other whitespace (tab, no-break space, ...) inside: walk back over key characters.
            start = eq
            while start > 0 and text[start - 1] in key_chars:
                start -= 1
            if start < eq and not text[start].isdigit() and (start == 0 or text[start - 1].isspace()):
                keys.append((start, eq))
        eq = text.find("=", eq + 1)
    if not keys:
        return None, []
    pairs: list[tuple[str, str]] = []
    for i in range(len(keys) - 1):
        start, eq = keys[i]
        pairs.append((text[start:eq], text[eq + 1 : keys[i + 1][0]].strip()))
    start, eq = keys[-1]
    pairs.append((text[start:eq], text[eq + 1 :].strip()))
    return keys[0][0], pairs


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

    first_match_start, pairs = _split_kv(rest)
    for key, value in pairs:
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


# Leading bytes of the compressed formats read transparently. The magic bytes alone
# decide: a plain-text file named ``x.fasta.gz`` is read as plain text.
_MAGIC = ((b"\x1f\x8b", "gz"), (b"BZh", "bz2"), (b"\xfd7zXZ\x00", "xz"))
_MODULES = {"gz": "gzip", "bz2": "bz2", "xz": "lzma"}


def _compressed_errors() -> tuple[type[Exception], ...]:
    """What a corrupt or truncated compressed stream raises while it is read.

    Built on use: ``lzma`` is optional in CPython builds (``_lzma`` may be missing).
    """
    errors: list[type[Exception]] = [UnicodeDecodeError, EOFError, OSError]
    try:
        import lzma
    except ImportError:
        pass
    else:
        errors.append(lzma.LZMAError)
    return tuple(errors)


def _kind(head: bytes) -> str | None:
    """Return ``"gz"``, ``"bz2"``, ``"xz"`` or ``None`` for the first bytes of a file."""
    for magic, kind in _MAGIC:
        if head.startswith(magic):
            return kind
    return None


def _compression(path: Path) -> str | None:
    """Return the compression of a regular file from its magic bytes (``None`` if plain)."""
    with path.open("rb") as raw:
        return _kind(raw.read(6))


class _TextOverRaw(io.TextIOWrapper):
    """Text over a decompressor; closing it also closes the underlying file.

    ``gzip``/``bz2``/``lzma`` never close a file object they were given.
    """

    def __init__(self, buffer: io.BufferedIOBase, raw: io.BufferedReader) -> None:
        super().__init__(buffer, encoding="utf-8-sig")  # ty: ignore[invalid-argument-type]
        self._raw_file = raw

    def close(self) -> None:
        try:
            super().close()
        finally:
            self._raw_file.close()


def _decompressor(kind: str, raw: io.BufferedReader) -> io.BufferedIOBase:
    try:
        if kind == "gz":
            import gzip

            return gzip.GzipFile(fileobj=raw, mode="rb")
        if kind == "bz2":
            import bz2

            return bz2.BZ2File(raw, mode="rb")
        import lzma

        return lzma.LZMAFile(raw, mode="rb")
    except ImportError as e:
        err = FastaError(f"Cannot read {kind}-compressed input: this Python has no {_MODULES[kind]} module")
        err.add_note(f"hint: decompress the file first, or use a Python built with {_MODULES[kind]} support")
        raise err from e


class _Prefixed(io.RawIOBase):
    """Bytes already read from ``raw`` (``head``), then the rest of ``raw``."""

    def __init__(self, head: bytes, raw: io.BufferedReader) -> None:
        super().__init__()
        self._head = head
        self._raw = raw

    def readable(self) -> bool:
        return True

    def readinto(self, buffer: Any) -> int:
        if self._head:
            n = min(len(buffer), len(self._head))
            buffer[:n] = self._head[:n]
            self._head = self._head[n:]
            return n
        return self._raw.readinto(buffer)

    def close(self) -> None:
        try:
            self._raw.close()
        finally:
            super().close()


def _with_head(raw: io.BufferedReader, n: int = 6) -> tuple[io.BufferedReader, bytes]:
    """Return a reader positioned at the start of ``raw`` and its first ``n`` bytes.

    ``peek`` returns only what one read delivered, which on a pipe can be shorter than
    ``n`` (a writer that sends one byte first). Then read until ``n`` bytes or the end
    of input and put them back in front of the stream.
    """
    head = raw.peek(n)[:n]
    if len(head) >= n:
        return raw, head
    head = raw.read(n)  # blocks until n bytes or EOF
    return io.BufferedReader(_Prefixed(head, raw)), head


def _open_path(source: str | Path) -> tuple[IO[str], bool]:
    """Open a FASTA path as UTF-8 text, decompressing gzip/bzip2/xz input.

    The file is opened once and its first bytes are peeked, not read, so pipes,
    FIFOs, ``/dev/stdin`` and process substitution work. Returns the handle and
    whether it is compressed.
    """
    raw = open(source, "rb")  # noqa: SIM115 - closed by the returned handle
    try:
        raw, head = _with_head(raw)
        kind = _kind(head)
        if kind is None:
            return io.TextIOWrapper(raw, encoding="utf-8-sig"), False
        return _TextOverRaw(_decompressor(kind, raw), raw), True
    except BaseException:
        raw.close()
        raise


def _iter_entries(fh: IO[str]) -> Iterator[SequenceEntry]:
    """Yield ``SequenceEntry`` objects from an open text-mode file, line by line.

    Undecodable input raises ``FastaParseError`` chained to the original error.
    """
    header: _ParsedHeader | None = None
    seq_chunks: list[str] = []
    header_line_no: int = 0
    line_no = 0

    try:
        for line_no, line in enumerate(fh, start=1):
            first = line[:1]
            if first == ">":
                if header is not None:
                    yield _build_entry(header, seq_chunks, header_line_no)
                header = _parse_header_line(line, line_no)
                header_line_no = line_no
                seq_chunks = []
                continue
            if first == ";":
                # Comment line (NCBI / legacy FASTA convention) — skip.
                continue
            if first == "\ufeff" and line_no == 1:
                # UTF-8 byte-order mark read through a plain ``utf-8`` text handle.
                line = line[1:]
                first = line[:1]
                if first == ">":
                    header = _parse_header_line(line, line_no)
                    header_line_no = line_no
                    continue
                if first == ";":
                    continue
            if first == "#" and header is None:
                # PEFF file header (``# PEFF 1.0``, ``# //``) before the first entry — skip.
                continue
            chunk = line.rstrip()
            if not chunk:
                # Blank or whitespace-only line.
                continue
            if header is None:
                raise FastaParseError(
                    "Sequence data appears before any '>' header",
                    line=line_no,
                    context=line.rstrip("\n"),
                    hint="A FASTA file must start with a '>' header line",
                )
            seq_chunks.append(chunk if chunk.isalpha() else "".join(chunk.split()))
    except UnicodeDecodeError as e:
        raise FastaParseError(
            f"Cannot read the input after line {line_no}: {e}",
            hint="FASTA input must be UTF-8 text, optionally gzip, bzip2 or xz compressed",
        ) from e

    if header is not None:
        yield _build_entry(header, seq_chunks, header_line_no)


# Characters read per ``read()`` call by the whole-record path reader.
_CHUNK_CHARS = 1 << 20


def _iter_path_entries(fh: IO[str], *, compressed: bool = False) -> Iterator[SequenceEntry]:
    """Yield ``SequenceEntry`` objects from a handle opened by ``_open_path``.

    Same result as ``_iter_entries``, but reads large chunks and splits them into
    records at ``\\n>`` instead of handling every line in Python. Only used on
    handles opened with universal newlines (``read()`` and line iteration then
    agree on where lines end); a record whose sequence is not plain letters is
    handled line by line exactly as ``_iter_entries`` does.
    """
    errors: tuple[type[Exception], ...] = _compressed_errors() if compressed else (UnicodeDecodeError,)
    # A virtual "\n" in front makes a header on line 1 split like any other, and
    # makes line indices of the text equal to 1-based file line numbers.
    pending: list[str] = ["\n"]
    header_line = 0  # line number of the next record's header; 0 until the preamble is done
    try:
        chunk = fh.read(_CHUNK_CHARS)
        if chunk.startswith("\ufeff"):
            # UTF-8 byte-order mark read through a plain ``utf-8`` text handle.
            chunk = chunk[1:] or fh.read(_CHUNK_CHARS)
        while chunk:
            if "\n>" in chunk or (chunk[0] == ">" and pending[-1].endswith("\n")):
                pieces = ("".join(pending) + chunk).split("\n>")
                pending = [pieces.pop()]
                for piece in pieces:
                    if header_line == 0:
                        header_line = _check_preamble(piece)
                    else:
                        yield _record_entry(piece, header_line)
                        header_line += piece.count("\n") + 1
            else:
                pending.append(chunk)
            chunk = fh.read(_CHUNK_CHARS)
    except errors as e:
        raise FastaParseError(
            f"Cannot read the input after line {max(header_line - 1, 0)}: {e}",
            hint="FASTA input must be UTF-8 text, optionally gzip, bzip2 or xz compressed",
        ) from e
    last = "".join(pending)
    if header_line == 0:
        _check_preamble(last)
    else:
        yield _record_entry(last, header_line, final=True)


def _check_preamble(text: str) -> int:
    """Check the text before the first header; return the first header's line number.

    ``text`` starts with the virtual ``\\n``, so ``lines[k]`` is file line ``k``.
    Blank, ``;`` comment and PEFF ``#`` lines are allowed; anything else is
    sequence data before a header.
    """
    lines = text.split("\n")
    for line_no in range(1, len(lines)):
        line = lines[line_no]
        if not line or line.isspace() or line[0] in ";#":
            continue
        raise FastaParseError(
            "Sequence data appears before any '>' header",
            line=line_no,
            context=line,
            hint="A FASTA file must start with a '>' header line",
        )
    return len(lines)


def _record_entry(record: str, header_line: int, *, final: bool = False) -> SequenceEntry:
    """Build the entry for ``record`` (``header\\nsequence lines``, without the ``>``).

    ``final`` is the last record of the input: only its header line can lack a ``\\n``.
    """
    head, sep, body = record.partition("\n")
    header = _parse_header_line(">" + head + (sep or ("" if final else "\n")), header_line)
    sequence = body.replace("\n", "")
    if not sequence.isalpha():
        # Comment lines, whitespace, '*', '-', ...: the line-by-line rules.
        chunks = []
        for line in body.split("\n"):
            if not line or line.isspace() or line[0] == ";":
                continue
            chunks.append("".join(line.split()))
        sequence = "".join(chunks)
    return _build_entry(header, [sequence], header_line)


class FastaReader:
    """Iterate over a FASTA file lazily without loading the entire file.

    Use it as a context manager. A reader is single-pass: iterating it again
    continues from where the previous iteration stopped, so after a full pass
    a second ``for`` loop yields nothing. Open a new reader (or use
    ``read_fasta``) to read the file again.

    A path may be plain or gzip/bzip2/xz compressed, as for ``read_fasta``.
    """

    def __init__(self, source: str | Path | IO[str]) -> None:
        self._source = source
        self._fh: IO[str] | None = None
        self._owns_fh = False
        self._compressed = False

    def __enter__(self) -> Self:
        if isinstance(self._source, (str, Path)):
            self._fh, self._compressed = _open_path(self._source)
            self._owns_fh = True
        else:
            self._fh = self._source
            self._owns_fh = False
            self._compressed = False
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

    def to_records(self) -> list[dict[str, str | int | None]]:
        """Read the remaining entries as flat dicts (see :func:`fastatacular.to_records`)."""
        from fastatacular._records import entry_to_record

        return [entry_to_record(e) for e in self]

    def __iter__(self) -> Iterator[SequenceEntry]:
        if self._fh is None:
            raise RuntimeError("FastaReader must be used as a context manager (`with FastaReader(...) as r:`)")
        if self._owns_fh:
            return _iter_path_entries(self._fh, compressed=self._compressed)
        return _iter_entries(self._fh)


def read_fasta(source: str | Path | IO[str]) -> list[SequenceEntry]:
    """Read an entire FASTA file into a list of ``SequenceEntry`` objects.

    A path may be plain or gzip/bzip2/xz compressed (detected from the magic
    bytes; a plain file with a ``.gz`` name is read as plain text).
    """
    if isinstance(source, (str, Path)):
        fh, compressed = _open_path(source)
        with fh:
            return list(_iter_path_entries(fh, compressed=compressed))
    return list(_iter_entries(source))


__all__ = ["FastaReader", "read_fasta"]
