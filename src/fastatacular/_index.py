"""Random access to entries of an uncompressed FASTA file by accession."""

from __future__ import annotations

import io
import mmap
import os
from collections.abc import Iterator, Mapping
from pathlib import Path

from fastatacular._models import SequenceEntry
from fastatacular._parser import _PIPE_ID, _UNIPROT_ID, _check_preamble, _compression, _iter_entries
from fastatacular.errors import FastaError, FastaParseError

_BOM = b"\xef\xbb\xbf"


def _key(identifier: str) -> str:
    """The index key of an identifier: its accession, or the identifier when it has none."""
    if m := _UNIPROT_ID.match(identifier):
        return m["accession"]
    if m := _PIPE_ID.match(identifier):
        return m["accession"]
    return identifier


def _identifier(header: bytes, start: int, path: Path) -> str:
    """First word of a header line (without ``>``); ``start`` is its byte offset, for errors."""
    try:
        text = header.decode("utf-8")
    except UnicodeDecodeError as e:
        raise FastaParseError(
            f"Cannot decode the header at byte {start} of {path}: {e}",
            hint="FASTA input must be UTF-8 text",
        ) from e
    words = text.split(maxsplit=1)
    if not words:
        raise FastaParseError(
            f"Empty FASTA header at byte {start} of {path}",
            context=">" + text,
            hint="Put an identifier right after '>'",
        )
    return words[0]


class FastaIndex(Mapping[str, SequenceEntry]):
    """Random access to the entries of an uncompressed FASTA file, keyed by accession.

    Building the index reads the file once and keeps only each entry's key and byte
    range. ``index[accession]`` then reads and parses that one entry from disk::

        index = FastaIndex("human.fasta")
        entry = index["P31946"]          # a SequenceEntry
        "P31946" in index, len(index)    # no file access
        for accession in index: ...      # keys in file order

    The key is the accession (``P31946`` for ``sp|P31946|1433B_HUMAN``, the second pipe
    field as in :class:`SequenceEntry`), or the whole identifier when it has no pipe.
    A missing key raises :class:`KeyError`. The file is opened on every lookup and must
    not change while the index is used.

    Only plain (uncompressed) files can be indexed; gzip, bgzip, bzip2 and xz input
    raises :class:`FastaError`. Decompress the file first (``gunzip -k file.fasta.gz``).

    :meth:`write_fai` and :meth:`from_fai` save and load a samtools-compatible ``.fai``.

    Raises:
        FastaError: The file is compressed, or two entries have the same accession (the
            message names the first repeated one).
        FastaParseError: Undecodable or empty headers, or sequence data before the
            first header. Errors inside an entry's sequence are raised when that entry
            is read.
    """

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        _require_plain(self.path)
        self._spans: dict[str, tuple[int, int]] = {}
        self._names: dict[str, str] = {}  # key -> identifier (the .fai name)
        with self.path.open("rb") as fh:
            size = os.fstat(fh.fileno()).st_size
            if size == 0:
                return
            with mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ) as mm:
                start = 3 if mm[:3] == _BOM else 0
                h = start if mm[start : start + 1] == b">" else mm.find(b"\n>", start)
                if h != start and h != -1:
                    h += 1
                _check_bytes_preamble(mm[start : size if h == -1 else h], self.path)
                while h != -1:
                    eol = mm.find(b"\n", h)
                    header_end = size if eol == -1 else eol
                    nxt = -1 if eol == -1 else mm.find(b"\n>", eol)
                    end = size if nxt == -1 else nxt + 1
                    identifier = _identifier(mm[h + 1 : header_end], h, self.path)
                    self._add(_key(identifier), identifier, h, end)
                    h = -1 if nxt == -1 else nxt + 1

    def _add(self, key: str, identifier: str, start: int, end: int) -> None:
        if key in self._spans:
            raise FastaError(
                f"Duplicate accession {key!r} in {self.path} (entry {identifier!r} at byte {start})",
            )
        self._spans[key] = (start, end)
        self._names[key] = identifier

    # -- Mapping --------------------------------------------------------------------

    def __getitem__(self, key: str) -> SequenceEntry:
        start, end = self._spans[key]
        with self.path.open("rb") as fh:
            fh.seek(start)
            data = fh.read(end - start)
        try:
            text = data.decode("utf-8")
            entries = list(_iter_entries(io.StringIO(text, newline=None)))
        except FastaParseError as e:
            raise FastaParseError(
                f"Entry {key!r} at byte {start} of {self.path}: {e}",
                hint=e.hint,
            ) from e
        except UnicodeDecodeError as e:
            raise FastaParseError(f"Cannot decode entry {key!r} at byte {start} of {self.path}: {e}") from e
        if len(entries) != 1 or entries[0].identifier != self._names[key]:
            raise FastaError(
                f"{self.path} changed since it was indexed: no entry {self._names[key]!r} at byte {start}",
            )
        return entries[0]

    def __contains__(self, key: object) -> bool:
        return key in self._spans

    def __iter__(self) -> Iterator[str]:
        return iter(self._spans)

    def __len__(self) -> int:
        return len(self._spans)

    def __repr__(self) -> str:
        return f"FastaIndex({str(self.path)!r}, {len(self)} entries)"

    def locate(self, key: str) -> tuple[int, int]:
        """Return ``(offset, length)`` in bytes of the entry's text, from its ``>`` to the next ``>``."""
        start, end = self._spans[key]
        return start, end - start

    def identifier(self, key: str) -> str:
        """Return the identifier (first header word, the ``.fai`` name) of the entry with ``key``."""
        return self._names[key]

    # -- .fai -----------------------------------------------------------------------

    def write_fai(self, fai_path: str | Path | None = None) -> Path:
        """Write a samtools-compatible ``.fai`` next to the file (or to ``fai_path``).

        Columns: name (first header word), residue count, byte offset of the first
        residue, residues per line, bytes per line (including the line ending). Like
        samtools, every sequence line of an entry but the last must have the same length;
        blank lines are only allowed after the sequence, and ``;`` comments not at all.

        Returns:
            The path written.

        Raises:
            FastaError: An entry has irregular line lengths or characters that samtools
                would not count as residues (whitespace inside a line, comments inside
                the sequence). Nothing is written.
        """
        out = Path(fai_path) if fai_path is not None else Path(f"{self.path}.fai")
        rows: list[str] = []
        with self.path.open("rb") as fh:
            for key, (start, end) in self._spans.items():
                fh.seek(start)
                rows.append(self._fai_row(key, start, fh.read(end - start)))
        out.write_text("".join(rows), encoding="utf-8", newline="\n")
        return out

    def _fai_row(self, key: str, start: int, record: bytes) -> str:
        name = self._names[key]
        eol = record.find(b"\n")
        if eol == -1:
            raise FastaError(f"Entry {name!r} has no sequence data")
        offset = start + eol + 1
        lines = record[eol + 1 :].split(b"\n")
        while lines and not lines[-1].strip():
            lines.pop()  # blank lines after the sequence
        if not lines:
            raise FastaError(f"Entry {name!r} has no sequence data")
        widths = [len(line) + 1 for line in lines]
        bases = [len(line.rstrip(b"\r")) for line in lines]
        for line, n in zip(lines, bases, strict=True):
            body = line[:n]
            if not body or body[:1] == b";" or len(body.split()) != 1 or len(body.split()[0]) != n:
                raise _irregular(name, "a blank, comment or whitespace-containing line inside the sequence")
        linebases, linewidth = bases[0], widths[0]
        if any(b != linebases for b in bases[:-1]) or any(w != linewidth for w in widths[:-1]):
            raise _irregular(name, "sequence lines of different lengths")
        if bases[-1] > linebases:
            raise _irregular(name, "a last line longer than the others")
        return f"{name}\t{sum(bases)}\t{offset}\t{linebases}\t{linewidth}\n"

    @classmethod
    def from_fai(cls, path: str | Path, fai_path: str | Path | None = None) -> FastaIndex:
        """Load an index from a samtools ``.fai`` (default ``path + ".fai"``) instead of scanning.

        The ``.fai`` name column is the first header word; it is mapped to the accession
        the same way as when building the index. Each entry's header line is checked
        against its name while loading, so a ``.fai`` that does not belong to the file
        raises :class:`FastaError`. Residues are not re-counted: a ``.fai`` made for a
        different version of the file with the same headers and offsets is not
        detected.

        Raises:
            FastaError: The file is compressed, a ``.fai`` line is malformed, an offset
                does not point just after a header with that name, or two names map to
                the same accession.
        """
        index = cls.__new__(cls)
        index.path = Path(path)
        _require_plain(index.path)
        index._spans = {}
        index._names = {}
        fai = Path(fai_path) if fai_path is not None else Path(f"{index.path}.fai")
        rows: list[tuple[int, int, str]] = []  # (header start, sequence offset, name)
        with index.path.open("rb") as fh:
            size = os.fstat(fh.fileno()).st_size
            text = fai.read_text(encoding="utf-8")
            if size == 0 or not text.strip():
                if text.strip():
                    raise FastaError(f"{fai} lists entries but {index.path} is empty")
                return index
            with mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ) as mm:
                for line_no, line in enumerate(text.splitlines(), 1):
                    if not line.strip():
                        continue
                    fields = line.split("\t")
                    if len(fields) != 5 or not all(f.isdigit() for f in fields[1:]):
                        raise FastaError(
                            f"{fai} line {line_no}: expected 5 tab-separated columns "
                            f"(name, length, offset, linebases, linewidth), got {line!r}",
                        )
                    name, offset = fields[0], int(fields[2])
                    if not 0 < offset <= size or mm[offset - 1 : offset] != b"\n":
                        raise _stale(fai, name, offset)
                    hstart = mm.rfind(b"\n", 0, offset - 1) + 1
                    if hstart == 0 and mm[:3] == _BOM:
                        hstart = 3
                    if mm[hstart : hstart + 1] != b">":
                        raise _stale(fai, name, offset)
                    header = mm[hstart + 1 : offset - 1].rstrip(b"\r")
                    if _identifier(header, hstart, index.path) != name:
                        raise _stale(fai, name, offset)
                    rows.append((hstart, offset, name))
        rows.sort()
        for i, (hstart, _, name) in enumerate(rows):
            end = rows[i + 1][0] if i + 1 < len(rows) else size
            index._add(_key(name), name, hstart, end)
        # Keep file order for iteration (already sorted by offset).
        return index


def _require_plain(path: Path) -> None:
    if (kind := _compression(path)) is not None:
        err = FastaError(f"Cannot index {path}: it is {kind}-compressed")
        err.add_note(
            "hint: FastaIndex needs an uncompressed file; decompress it first "
            "(e.g. gunzip -k file.fasta.gz). bgzip/.gzi is not supported."
        )
        raise err


def _check_bytes_preamble(data: bytes, path: Path) -> None:
    if not data.strip():
        return
    try:
        _check_preamble("\n" + data.decode("utf-8"))
    except UnicodeDecodeError as e:
        raise FastaParseError(f"Cannot decode the start of {path}: {e}") from e


def _irregular(name: str, what: str) -> FastaError:
    err = FastaError(f"Entry {name!r} cannot be written to a .fai: {what}")
    err.add_note("hint: rewrite the file with write_fasta(read_fasta(src), dst) for uniform 60-residue lines")
    return err


def _stale(fai: Path, name: str, offset: int) -> FastaError:
    err = FastaError(f"{fai}: entry {name!r} at offset {offset} does not match the FASTA file")
    err.add_note(
        "hint: the .fai is out of date or belongs to another file; rebuild it with FastaIndex(path).write_fai()"
    )
    return err
