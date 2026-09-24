"""Compressed input, undecodable input and PEFF file headers."""

from __future__ import annotations

import bz2
import gzip
import io
import lzma
import os
import subprocess
import sys
import threading
from collections.abc import Callable
from pathlib import Path

import pytest

from fastatacular import FastaError, FastaParseError, FastaReader, read_fasta

TEXT = ">sp|P12345|EX_HUMAN Example OS=Homo sapiens OX=9606 GN=EXMP PE=1 SV=2\nMKTIIALSYI\nFCLVFA\n>b\nPEPTIDE\n"

COMPRESSORS: dict[str, Callable[[bytes], bytes]] = {
    ".gz": gzip.compress,
    ".bz2": bz2.compress,
    ".xz": lzma.compress,
}


@pytest.fixture
def plain(tmp_path: Path) -> list:  # type: ignore[type-arg]
    path = tmp_path / "plain.fasta"
    path.write_text(TEXT)
    return read_fasta(path)


@pytest.mark.parametrize("suffix", list(COMPRESSORS))
def test_compressed_path_reads_like_plain(tmp_path: Path, plain: list, suffix: str) -> None:  # type: ignore[type-arg]
    path = tmp_path / f"x.fasta{suffix}"
    path.write_bytes(COMPRESSORS[suffix](TEXT.encode()))
    assert read_fasta(path) == plain
    assert read_fasta(str(path)) == plain
    with FastaReader(path) as reader:
        assert list(reader) == plain


@pytest.mark.parametrize("suffix", list(COMPRESSORS))
def test_compression_detected_from_magic_bytes(tmp_path: Path, plain: list, suffix: str) -> None:  # type: ignore[type-arg]
    path = tmp_path / "no_suffix.fasta"
    path.write_bytes(COMPRESSORS[suffix](TEXT.encode()))
    assert read_fasta(path) == plain


def test_compressed_with_bom(tmp_path: Path, plain: list) -> None:  # type: ignore[type-arg]
    path = tmp_path / "bom.fasta.gz"
    path.write_bytes(gzip.compress(("﻿" + TEXT).encode()))
    assert read_fasta(path) == plain


def test_gz_suffix_on_plain_text_reads_as_plain(tmp_path: Path, plain: list) -> None:  # type: ignore[type-arg]
    # The magic bytes decide, not the file name.
    path = tmp_path / "not_really.fasta.gz"
    path.write_text(TEXT, encoding="utf-8")
    assert read_fasta(path) == plain


@pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="needs os.mkfifo")
@pytest.mark.parametrize("suffix", ["", *COMPRESSORS])
def test_fifo_is_read_once(tmp_path: Path, plain: list, suffix: str) -> None:  # type: ignore[type-arg]
    # A pipe can be read only once: the format sniff must not consume its first bytes.
    data = COMPRESSORS[suffix](TEXT.encode()) if suffix else TEXT.encode()
    fifo = tmp_path / "pipe.fasta"
    os.mkfifo(fifo)

    def feed() -> None:
        with fifo.open("wb") as w:
            w.write(data)

    result: list = []  # type: ignore[type-arg]
    writer = threading.Thread(target=feed, daemon=True)
    reader = threading.Thread(target=lambda: result.append(read_fasta(fifo)), daemon=True)
    writer.start()
    reader.start()
    reader.join(timeout=10)
    if reader.is_alive():  # a second open() of the FIFO blocks: give it EOF, then fail
        os.close(os.open(fifo, os.O_WRONLY | os.O_NONBLOCK))
        pytest.fail("reading the FIFO blocked (was it opened twice?)")
    writer.join(timeout=10)
    assert result == [plain]


@pytest.mark.skipif(sys.platform == "win32", reason="/dev/fd is POSIX only")
def test_os_pipe_path(plain: list) -> None:  # type: ignore[type-arg]
    r, w = os.pipe()
    writer = threading.Thread(target=lambda: (os.write(w, TEXT.encode()), os.close(w)))
    writer.start()
    try:
        with FastaReader(f"/dev/fd/{r}") as reader:
            assert list(reader) == plain
    finally:
        writer.join(timeout=10)
        os.close(r)


_NO_LZMA_BZ2 = """
import sys
sys.modules["_lzma"] = None
sys.modules["_bz2"] = None
import fastatacular
from fastatacular import FastaError, read_fasta
assert [e.identifier for e in read_fasta(sys.argv[1])] == ["sp|P12345|EX_HUMAN", "b"]
for path, module in ((sys.argv[2], "lzma"), (sys.argv[3], "bz2")):
    try:
        read_fasta(path)
    except FastaError as e:
        assert module in str(e), e
    else:
        raise SystemExit(f"{path} read without {module}")
print("ok")
"""


def test_import_without_lzma_and_bz2(tmp_path: Path) -> None:
    # Minimal Python builds (pyenv, slim images) can lack _lzma and _bz2.
    plain_path = tmp_path / "p.fasta"
    plain_path.write_text(TEXT, encoding="utf-8")
    xz = tmp_path / "x.fasta.xz"
    xz.write_bytes(lzma.compress(TEXT.encode()))
    bz = tmp_path / "x.fasta.bz2"
    bz.write_bytes(bz2.compress(TEXT.encode()))
    result = subprocess.run(
        [sys.executable, "-c", _NO_LZMA_BZ2, str(plain_path), str(xz), str(bz)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "ok"


@pytest.mark.parametrize("suffix", list(COMPRESSORS))
def test_truncated_compressed_raises_fasta_error(tmp_path: Path, suffix: str) -> None:
    data = COMPRESSORS[suffix]((TEXT * 50).encode())
    path = tmp_path / f"cut.fasta{suffix}"
    path.write_bytes(data[: len(data) // 2])
    with pytest.raises(FastaError, match="Cannot read the input") as info:
        read_fasta(path)
    assert info.value.__cause__ is not None


def test_non_utf8_path_raises_fasta_error(tmp_path: Path) -> None:
    path = tmp_path / "latin1.fasta"
    path.write_bytes(">x caf\xe9\nMK\n".encode("latin-1"))
    with pytest.raises(FastaParseError, match="Cannot read the input") as info:
        read_fasta(path)
    assert isinstance(info.value.__cause__, UnicodeDecodeError)
    with pytest.raises(FastaParseError), FastaReader(path) as reader:
        list(reader)


def test_non_utf8_gzip_raises_fasta_error(tmp_path: Path) -> None:
    path = tmp_path / "latin1.fasta.gz"
    path.write_bytes(gzip.compress(">x caf\xe9\nMK\n".encode("latin-1")))
    with pytest.raises(FastaParseError) as info:
        read_fasta(path)
    assert isinstance(info.value.__cause__, UnicodeDecodeError)


def test_non_utf8_text_handle_raises_fasta_error() -> None:
    handle = io.TextIOWrapper(io.BytesIO(b">a\nMK\n>b caf\xe9\nPE\n"), encoding="utf-8")
    with pytest.raises(FastaParseError, match="Cannot read the input") as info:
        read_fasta(handle)
    assert isinstance(info.value.__cause__, UnicodeDecodeError)


# ------------------------------------------------------------------ PEFF file headers

PEFF = """# PEFF 1.0
# //
# DbName=neXtProt
# Prefix=nxp
# //
>nxp:NX_P1 \\DbUniqueId=NX_P1 \\PName=Example \\Length=6
MKTIIA
>nxp:NX_P2 \\DbUniqueId=NX_P2
PEPT
"""


@pytest.mark.parametrize("as_path", [True, False])
def test_peff_file_header_is_skipped(tmp_path: Path, as_path: bool) -> None:
    path = tmp_path / "x.peff"
    path.write_text(PEFF)
    entries = read_fasta(path) if as_path else read_fasta(io.StringIO(PEFF))
    assert [(e.identifier, e.sequence) for e in entries] == [("nxp:NX_P1", "MKTIIA"), ("nxp:NX_P2", "PEPT")]
    assert entries[0].description == "\\DbUniqueId=NX_P1 \\PName=Example \\Length=6"


def test_hash_line_after_first_entry_is_sequence_data() -> None:
    # Only the file header before the first '>' is skipped; later '#' lines are unchanged.
    [entry] = read_fasta(io.StringIO(">a\nMK\n#X\n"))
    assert entry.sequence == "MK#X"
