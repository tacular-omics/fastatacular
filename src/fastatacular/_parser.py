"""FASTA file parser — converts text into model objects."""

from __future__ import annotations

import io
import re
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from types import TracebackType
from typing import IO, Any, Self, cast

from fastatacular._compression import (
    _KINDS,
    _NAMES,
    Compression,
    _check_compression,
    _DetachingText,
    _is_text_handle,
)
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

# Decoy/contaminant prefix allowed in front of ``gi|``, ``pdb|``, ``pir||`` and ``prf||``:
# any run without ``|`` or whitespace that ends in a character other than a letter or digit
# (``DECOY_``, ``rev_``, ``XXX_``, ``REV__``, ``DECOY-0-``, ``decoy.``). Prefixes stack: the
# decoy of a contaminant is ``DECOY_CON_gi|...``. The final separator keeps ``fungi|...``
# and ``xgi|...`` from being read as gi ids.
_TAGS = r"(?:[^|\s]*[^|\sA-Za-z0-9])?"

# Legacy NCBI ``gi|NUMBER|db|ACCESSION|NAME`` identifier (``gi|4557757|ref|NP_000240.1|``):
# the accession is the fourth field, not the GenInfo number. The prefix is ``gi`` itself,
# optionally after a decoy/contaminant prefix (``_TAGS``).
# The fifth field is the name NCBI puts after the accession (``gi|1|sp|P12345|NAME_HUMAN``),
# kept as the entry name, as for the bare ``sp|P12345|NAME_HUMAN`` form.
# For ``pdb`` the fifth field is the chain: ``gi|229552|pdb|1MBA|A`` gives ``1MBA_A``
# (NCBI's entry_chain form), so chains of one entry keep distinct accessions.
_GI_ID = re.compile(
    rf"^(?P<prefix>{_TAGS}gi)"
    r"\|(?P<gi>[0-9]+)\|(?P<db>[A-Za-z]+)\|(?P<accession>[^|\s]*)"
    r"(?:\|(?P<name>[^|\s]*)(?:\|.*)?)?$"
)

# Bare NCBI ``pdb|ENTRY|CHAIN`` identifier: the accession is ``ENTRY_CHAIN``, the same as
# for ``gi|N|pdb|ENTRY|CHAIN``, with no entry name. Checked before the UniProt rule, which
# would otherwise take the chain for the entry name and give every chain the same accession.
_PDB_ID = re.compile(rf"^(?P<prefix>{_TAGS}pdb)\|(?P<entry>[^|\s]+)\|(?P<chain>[^|\s]*)(?:\|.*)?$")

# NCBI ``pir||ENTRY`` and ``prf||NAME`` identifiers: the middle field is empty by
# definition, so the accession is the third field. ``gi|N|pir||ENTRY`` is handled by ``_GI_ID``.
_EMPTY_MIDDLE_DBS = frozenset({"pir", "prf"})
_EMPTY_MIDDLE_ID = re.compile(rf"^(?P<prefix>{_TAGS}(?:pir|prf))\|\|(?P<accession>[^|\s]+)(?:\|.*)?$")

# NCBI-ish ``db|ID`` or ``db|ID|...`` identifier — accept the leading two fields.
_PIPE_ID = re.compile(r"^(?P<prefix>[^|\s]+)\|(?P<accession>[^|\s]+)(?:\|.*)?$")

# Current RCSB PDB entity header: ``ENTRY_ENTITY|Chain(s) X|MOLECULE|ORGANISM``
# (``1MBA_1|Chain A|MYOGLOBIN|Aplysia limacina``). Its second field holds a space, so the
# usual whitespace split would cut the identifier inside ``Chain A`` and the pipe rules
# would then read ``Chain`` as the accession. The identifier is the first field instead,
# and the rest of the line is the description. Anchored on ``|Chain`` right after the
# first field, so the older ``pdb|1MBA|A Chain A`` form (second field ``1MBA``) is
# untouched and keeps its ``1MBA_A`` accession.
_RCSB_ID = re.compile(r"^(?P<identifier>[^|\s]+)\|Chains?\s")

# Identifiers that are a bare accession, with no ``db|`` prefix. Each may carry a
# decoy/contaminant prefix (``_TAGS``), which is not part of the accession, and each
# names the database whose description conventions the rest of the header follows
# (see ``_bare_accession`` and ``_parse_header_line``).
#
# NCBI protein accession, RefSeq (``NP_000549.1``) or GenPept (``AAA16334.1``). The
# RefSeq prefixes are the protein ones only. The GenPept form must carry its version, so
# an ordinary word such as ``ABC12345`` is not read as an accession.
_NCBI_ACC = re.compile(
    rf"^{_TAGS}(?P<accession>(?:AP|NP|XP|YP|WP|ZP)_\d+(?:\.\d+)?|[A-Z]{{3}}\d{{5}}(?:\d{{2}})?\.\d+)$"
)

# Ensembl protein identifier of an Ensembl ``pep`` file: ``ENSP00000322421.5``, with an
# optional species code (``ENSMUSP00000020316.9``). Ensembl Genomes uses gene-model names
# instead, which have no fixed shape and are left unparsed.
_ENSEMBL_ACC = re.compile(rf"^{_TAGS}(?P<accession>ENS[A-Z]*P\d{{6,}}(?:\.\d+)?)$")

# UniRef cluster identifier: ``UniRef90_P69905`` (also ``UniRef50_`` and ``UniRef100_``).
_UNIREF_ID = re.compile(rf"^{_TAGS}(?P<accession>UniRef(?:50|90|100)_[^|\s]+)$")

_BARE_ID_KINDS = ((_NCBI_ACC, "ncbi"), (_ENSEMBL_ACC, "ensembl"), (_UNIREF_ID, "uniref"))

# ``gene_symbol:HBA1`` in an Ensembl ``pep`` description. Ensembl writes ``key:value``,
# not ``KEY=value``, so ``_split_kv`` does not see these tokens.
_ENSEMBL_GENE_SYMBOL = re.compile(r"(?:^|(?<=\s))gene_symbol:(?P<gname>\S+)")

# The trailing ``[organism]`` of an NCBI defline (``hemoglobin subunit alpha
# [Homo sapiens]``). NCBI also writes ``[key=value]`` qualifiers, which name no organism,
# so a bracket holding ``=`` does not count.
_NCBI_ORGANISM = re.compile(r"\[(?P<os_name>[^\[\]=]+)\]$")


def _split_identifier(identifier: str) -> tuple[str, str, str | None] | None:
    """Return ``(prefix, accession, entry_name)`` of a pipe-delimited identifier, or None.

    ``sp|P12345|NAME_HUMAN`` gives ``("sp", "P12345", "NAME_HUMAN")``;
    ``gi|4557757|ref|NP_000240.1|`` gives ``("gi", "NP_000240.1", None)`` (the GenInfo
    number stays in the identifier) and ``gi|1|sp|P12345|NAME_HUMAN`` gives
    ``("gi", "P12345", "NAME_HUMAN")``; ``pdb|1MBA|A`` and ``gi|229552|pdb|1MBA|A`` give
    the accession ``1MBA_A``; ``pir||S71500``, ``prf||0601246A`` and
    ``gi|7428543|pir||S71500`` give the third (or fifth) field; other ``db|ID|...`` ids
    give their first two fields.
    """
    if m := _PDB_ID.match(identifier):
        chain = m["chain"]
        return m["prefix"], f"{m['entry']}_{chain}" if chain else m["entry"], None
    if m := _UNIPROT_ID.match(identifier):
        return m["prefix"], m["accession"], m["entry_name"]
    if m := _EMPTY_MIDDLE_ID.match(identifier):
        return m["prefix"], m["accession"], None
    if (m := _GI_ID.match(identifier)) and (m["accession"] or (m["db"].lower() in _EMPTY_MIDDLE_DBS and m["name"])):
        accession, name = m["accession"], m["name"] or None
        if not accession:  # gi|N|pir||ENTRY: the entry is the fifth field
            return m["prefix"], cast(str, name), None
        if m["db"].lower() == "pdb":
            return m["prefix"], f"{accession}_{name}" if name else accession, None
        return m["prefix"], accession, name
    if m := _PIPE_ID.match(identifier):
        return m["prefix"], m["accession"], None
    return None


def _bare_accession(identifier: str) -> tuple[str, str] | None:
    """Return ``(accession, kind)`` of an identifier that is a bare accession, or None.

    ``kind`` is the database whose header conventions apply (``"ncbi"``, ``"ensembl"`` or
    ``"uniref"``); ``_parse_header_line`` uses it to read the fields those conventions put
    in the description. Only shapes that identify their database on their own are
    recognized, so an arbitrary first word stays an identifier with no accession.
    """
    for pattern, kind in _BARE_ID_KINDS:
        if m := pattern.match(identifier):
            return m["accession"], kind
    return None


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

    if rcsb := _RCSB_ID.match(stripped):
        # RCSB PDB entity header: the identifier is the first pipe field (see ``_RCSB_ID``).
        identifier = rcsb["identifier"]
        rest = stripped[rcsb.end("identifier") + 1 :]
    else:
        # Split on the first run of any whitespace (space or tab), as UniProt, BLAST,
        # Biopython and samtools do.
        identifier, *tail = stripped.split(maxsplit=1)
        rest = tail[0] if tail else ""

    header = _ParsedHeader(identifier=identifier, raw_header=raw)

    kind: str | None = None
    if (parts := _split_identifier(identifier)) is not None:
        header.prefix, header.accession, header.entry_name = parts
    elif rcsb:
        # ``1MBA_1`` is the entity's accession at the RCSB; it carries no db prefix.
        header.accession = identifier
    elif bare := _bare_accession(identifier):
        # A bare NCBI, Ensembl or UniRef accession: the whole identifier (minus any
        # decoy/contaminant tag) is the accession, and it carries no db prefix.
        header.accession, kind = bare

    if not rest:
        return header

    first_match_start, pairs = _split_kv(rest)
    seen: set[str] = set()
    for key, value in pairs:
        # A repeated key keeps its first value (``OS=Human OS=Mouse`` gives ``Human``);
        # ``description`` and ``raw_header`` still hold every occurrence.
        if key in seen:
            continue
        seen.add(key)
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

    if kind is not None:
        _read_database_conventions(header, kind)
    return header


def _read_database_conventions(header: _ParsedHeader, kind: str) -> None:
    """Fill the fields that ``kind``'s own description conventions name.

    The databases whose identifier is a bare accession do not use UniProt's
    ``OS=``/``OX=``/``GN=`` keys throughout:

    - ``ncbi``: the organism is the defline's trailing ``[Homo sapiens]``. It stays in
      ``pname`` as well, because ``pname`` is the description text before its first
      ``KEY=`` token whatever the database.
    - ``ensembl``: the gene symbol is the ``gene_symbol:`` token. ``gene:`` and
      ``transcript:`` hold Ensembl gene and transcript ids, which no field takes, and
      stay in ``description`` only.
    - ``uniref``: the organism is ``Tax=`` and its NCBI taxonomy id ``TaxID=``, which move
      out of ``extra`` the way ``OS=`` and ``OX=`` never enter it. ``n=`` (cluster size)
      and ``RepID=`` (the representative member's entry name) stay in ``extra``.

    An explicit UniProt key wins over all of these, and a ``TaxID=`` that is not an
    integer stays in ``extra``, as a non-integer ``OX=`` does.
    """
    if kind == "uniref":
        if header.os_name is None and "Tax" in header.extra:
            header.os_name = header.extra.pop("Tax")
        if header.ncbi_tax_id is None and "TaxID" in header.extra:
            try:
                header.ncbi_tax_id = int(header.extra["TaxID"])
            except ValueError:
                pass
            else:
                del header.extra["TaxID"]
    elif kind == "ensembl":
        if header.gname is None and (m := _ENSEMBL_GENE_SYMBOL.search(header.description or "")):
            header.gname = m["gname"]
    elif kind == "ncbi":
        if header.os_name is None and header.pname and (m := _NCBI_ORGANISM.search(header.pname)):
            header.os_name = m["os_name"]


def _build_entry(header: _ParsedHeader, seq_chunks: list[str]) -> SequenceEntry:
    # A header with no residue lines (next line is a header, or end of input) gives "".
    sequence = "".join(chunk for chunk in seq_chunks if chunk)
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
    """Bytes already read from ``raw`` (``head``), then the rest of ``raw``.

    ``close_raw=False`` leaves ``raw`` open on close (a caller's handle).
    """

    def __init__(self, head: bytes, raw: IO[bytes], *, close_raw: bool = True) -> None:
        super().__init__()
        self._head = head
        self._raw = raw
        self._close_raw = close_raw

    def readable(self) -> bool:
        return True

    def readinto(self, buffer: Any) -> int:
        if self._head:
            n = min(len(buffer), len(self._head))
            buffer[:n] = self._head[:n]
            self._head = self._head[n:]
            return n
        readinto = getattr(self._raw, "readinto", None)
        if readinto is not None:
            return readinto(buffer)
        data = self._raw.read(len(buffer))
        buffer[: len(data)] = data
        return len(data)

    def close(self) -> None:
        try:
            if self._close_raw:
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


def _read_kind(head: bytes, compression: Compression, name: object) -> str | None:
    """The compression to read with: sniffed (``"infer"``), forced, or none (``None``).

    A forced format must match the magic bytes, else :class:`FastaParseError`.
    """
    if compression is None:
        return None
    found = _kind(head)
    if compression == "infer":
        return found
    kind = _KINDS[compression]
    if found != kind:
        what = "plain text (no compression magic bytes)" if found is None else f"{_NAMES[found]}-compressed"
        raise FastaParseError(
            f"compression={compression!r} but {name} looks {what}",
            hint="pass compression='infer' to detect the format from the magic bytes",
        )
    return kind


def _open_path(source: str | Path, compression: Compression = "infer") -> tuple[IO[str], bool]:
    """Open a FASTA path as UTF-8 text, decompressing gzip/bzip2/xz input.

    The file is opened once and its first bytes are peeked, not read, so pipes,
    FIFOs, ``/dev/stdin`` and process substitution work. Returns the handle and
    whether it is compressed. ``compression`` is as for :func:`read_fasta`.
    """
    raw = open(source, "rb")  # noqa: SIM115 - closed by the returned handle
    try:
        raw, head = _with_head(raw)
        kind = _read_kind(head, compression, source)
        if kind is None:
            return io.TextIOWrapper(raw, encoding="utf-8-sig"), False
        return _TextOverRaw(_decompressor(kind, raw), raw), True
    except BaseException:
        raw.close()
        raise


def _open_stream(handle: IO[Any], compression: Compression) -> tuple[IO[str], bool] | None:
    """Wrap a caller's binary handle in a UTF-8 text reader; return it and whether it decompresses.

    Returns ``None`` for a text handle with ``"infer"`` or ``None`` (used as is). A binary
    handle with ``"infer"`` or ``None`` is decoded as plain UTF-8 (no sniffing); an explicit
    ``compression`` decompresses it. Closing the returned reader leaves ``handle`` open.
    """
    if compression is None or compression == "infer":
        if _is_text_handle(handle):
            return None
        return _DetachingText(handle, encoding="utf-8-sig"), False
    if _is_text_handle(handle):
        err = FastaError(f"compression={compression!r} needs a binary handle, got a text-mode one")
        err.add_note("hint: open the file with mode 'rb', or pass a path")
        raise err
    head = handle.read(6)
    kind = _read_kind(head, compression, "the input")
    assert kind is not None
    buffered = io.BufferedReader(_Prefixed(head, handle, close_raw=False))
    try:
        return io.TextIOWrapper(_decompressor(kind, buffered), encoding="utf-8-sig"), True  # ty: ignore[invalid-argument-type]
    except BaseException:
        buffered.close()
        raise


def _iter_entries(fh: IO[str]) -> Iterator[SequenceEntry]:
    """Yield ``SequenceEntry`` objects from an open text-mode file, line by line.

    Undecodable input raises ``FastaParseError`` chained to the original error.
    """
    header: _ParsedHeader | None = None
    seq_chunks: list[str] = []
    line_no = 0

    try:
        for line_no, line in enumerate(fh, start=1):
            first = line[:1]
            if first == "\ufeff" and line_no == 1:
                # UTF-8 byte-order mark read through a plain ``utf-8`` text handle.
                line = line[1:]
                first = line[:1]
            if first.isspace():
                # Leading whitespace is ignored on every line, so an indented
                # ``>`` still starts a header (as in FastaFrames).
                line = line.lstrip()
                first = line[:1]
            if first == ">":
                if header is not None:
                    yield _build_entry(header, seq_chunks)
                header = _parse_header_line(line, line_no)
                seq_chunks = []
                continue
            if first == ";":
                # Comment line (NCBI / legacy FASTA convention) — skip.
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
        yield _build_entry(header, seq_chunks)


# Characters read per ``read()`` call by the whole-record path reader.
_CHUNK_CHARS = 1 << 20
# A header line with leading whitespace, and whitespace after the last line break.
_INDENTED_HEADER = re.compile(r"\n[^\S\n]+>")
_TRAILING_INDENT = re.compile(r"\n[^\S\n]+\Z")


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
        tail = ""  # indentation at the end of the previous chunk, maybe before a '>'
        while chunk:
            # An indented header (``  >b``) splits like ``>b``: drop the indentation.
            # Indentation at the end of a chunk is held back until the next one shows
            # whether a '>' follows it.
            chunk = tail + chunk
            line_start = pending[-1].endswith("\n")
            if line_start:
                chunk = "\n" + chunk
            chunk = _INDENTED_HEADER.sub("\n>", chunk)
            cut = _TRAILING_INDENT.search(chunk)
            if cut is not None:
                tail, chunk = chunk[cut.start() + 1 :], chunk[: cut.start() + 1]
            else:
                tail = ""
            if line_start:
                chunk = chunk[1:]
            if not chunk:
                more = fh.read(_CHUNK_CHARS)
                if not more:
                    break
                chunk = more
                continue
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
    pending.append(tail)
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
        line = line.lstrip()
        if not line or line[0] in ";#":
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
            line = line.lstrip()
            if not line or line[0] == ";":
                continue
            chunks.append("".join(line.split()))
        sequence = "".join(chunks)
    return _build_entry(header, [sequence])


class FastaReader:
    """Iterate over a FASTA file lazily without loading the entire file.

    Use it as a context manager. A reader is single-pass: iterating it again
    continues from where the previous iteration stopped, so after a full pass
    a second ``for`` loop yields nothing. Open a new reader (or use
    ``read_fasta``) to read the file again.

    A path may be plain or gzip/bzip2/xz compressed; ``compression`` is as for
    ``read_fasta``.
    """

    def __init__(self, source: str | Path | IO[str] | IO[bytes], *, compression: Compression = "infer") -> None:
        self._source = source
        self._compression = _check_compression(compression)
        self._fh: IO[str] | None = None
        self._owns_fh = False
        self._compressed = False

    def __enter__(self) -> Self:
        if isinstance(self._source, (str, Path)):
            self._fh, self._compressed = _open_path(self._source, self._compression)
            self._owns_fh = True
        elif (wrapped := _open_stream(self._source, self._compression)) is not None:
            self._fh, self._compressed = wrapped
            self._owns_fh = True
        else:
            self._fh = cast("IO[str]", self._source)
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


def read_fasta(source: str | Path | IO[str] | IO[bytes], *, compression: Compression = "infer") -> list[SequenceEntry]:
    """Read an entire FASTA file into a list of ``SequenceEntry`` objects.

    ``source`` is a path or an open file object.

    Args:
        compression: ``"infer"`` (default) detects gzip, bzip2 or xz from a path's
            magic bytes (a plain file with a ``.gz`` name is read as plain text) and
            reads a handle as text. ``"gzip"``, ``"bz2"`` or ``"xz"`` forces that
            format: a path or **binary** handle whose bytes are not in it raises
            :class:`FastaParseError`, and a text handle raises :class:`FastaError`. ``None`` reads plain text
            even if the bytes look compressed. A binary handle with ``"infer"`` or
            ``None`` is decoded as plain UTF-8 (no sniffing); a handle is never closed.

    Raises:
        FastaError: An unknown ``compression``, or a text handle with an explicit one.
        FastaParseError: Malformed or undecodable input, or input whose magic bytes
            do not match an explicit ``compression``.
    """
    compression = _check_compression(compression)
    if isinstance(source, (str, Path)):
        fh, compressed = _open_path(source, compression)
        with fh:
            return list(_iter_path_entries(fh, compressed=compressed))
    if (wrapped := _open_stream(source, compression)) is not None:
        fh, compressed = wrapped
        with fh:
            return list(_iter_path_entries(fh, compressed=compressed))
    return list(_iter_entries(source))  # ty: ignore[invalid-argument-type]


__all__ = ["Compression", "FastaReader", "read_fasta"]
