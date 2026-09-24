"""Compressed input, undecodable input and PEFF file headers."""

from __future__ import annotations

import bz2
import gzip
import io
import lzma
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


def test_gz_suffix_on_plain_text_raises_fasta_error(tmp_path: Path) -> None:
    path = tmp_path / "not_really.fasta.gz"
    path.write_text(TEXT)
    with pytest.raises(FastaParseError, match="Cannot read the input") as info:
        read_fasta(path)
    assert isinstance(info.value.__cause__, OSError)


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
