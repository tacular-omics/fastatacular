"""The ``compression`` parameter of the readers and writers."""

from __future__ import annotations

import bz2
import gzip
import io
import lzma
from collections.abc import Callable
from pathlib import Path

import pytest

from fastatacular import (
    Compression,
    FastaError,
    FastaIndex,
    FastaParseError,
    FastaReader,
    SequenceEntry,
    read_fasta,
    train_markov_model,
    write_decoy_fasta,
    write_fasta,
)

TEXT = ">sp|P12345|EX_HUMAN Example OS=Homo sapiens OX=9606 GN=EXMP PE=1 SV=2\nMKTIIALSYI\n>b\nPEPTIDE\n"

COMPRESS: dict[str, Callable[[bytes], bytes]] = {"gzip": gzip.compress, "bz2": bz2.compress, "xz": lzma.compress}
DECOMPRESS: dict[str, Callable[[bytes], bytes]] = {
    "gzip": gzip.decompress,
    "bz2": bz2.decompress,
    "xz": lzma.decompress,
}
MAGIC = {"gzip": b"\x1f\x8b", "bz2": b"BZh", "xz": b"\xfd7zXZ\x00"}
SUFFIX = {"gzip": ".gz", "bz2": ".bz2", "xz": ".xz"}
FORMATS = list(COMPRESS)


@pytest.fixture
def entries() -> list[SequenceEntry]:
    return read_fasta(io.StringIO(TEXT))


def _format(data: bytes) -> str | None:
    return next((name for name, magic in MAGIC.items() if data.startswith(magic)), None)


def test_compression_alias_exported() -> None:
    assert Compression is not None
    import fastatacular

    assert "Compression" in fastatacular.__all__


# --- validation ---------------------------------------------------------------------


@pytest.mark.parametrize("bad", ["gz", "GZIP", "zip", "", 1, True])
def test_unknown_compression_rejected_everywhere(tmp_path: Path, entries: list[SequenceEntry], bad: object) -> None:
    path = tmp_path / "x.fasta"
    path.write_text(TEXT)
    calls: list[Callable[[], object]] = [
        lambda: read_fasta(path, compression=bad),  # type: ignore[arg-type]
        lambda: FastaReader(path, compression=bad),  # type: ignore[arg-type]
        lambda: write_fasta(entries, tmp_path / "o.fasta", compression=bad),  # type: ignore[arg-type]
        lambda: FastaIndex(path, compression=bad),  # type: ignore[arg-type]
        lambda: FastaIndex.from_fai(path, compression=bad),  # type: ignore[arg-type]
        lambda: train_markov_model(path, compression=bad),  # type: ignore[arg-type]
        lambda: write_decoy_fasta(path, tmp_path / "d.fasta", method="reverse", compression=bad),  # type: ignore[arg-type]
    ]
    for call in calls:
        with pytest.raises(FastaError, match="'infer', 'gzip', 'bz2', 'xz' or None"):
            call()
    assert not (tmp_path / "o.fasta").exists()
    assert not (tmp_path / "d.fasta").exists()


# --- writing to a path --------------------------------------------------------------


@pytest.mark.parametrize("fmt", FORMATS)
def test_write_infer_by_suffix_case_insensitive(tmp_path: Path, entries: list[SequenceEntry], fmt: str) -> None:
    path = tmp_path / f"o.fasta{SUFFIX[fmt].upper()}"
    write_fasta(entries, path)
    data = path.read_bytes()
    assert _format(data) == fmt
    assert DECOMPRESS[fmt](data).decode() == TEXT


def test_write_infer_plain_without_suffix(tmp_path: Path, entries: list[SequenceEntry]) -> None:
    path = tmp_path / "o.fasta"
    write_fasta(entries, path, compression="infer")
    assert path.read_text() == TEXT


@pytest.mark.parametrize("fmt", FORMATS)
def test_write_explicit_on_plain_suffix(tmp_path: Path, entries: list[SequenceEntry], fmt: str) -> None:
    path = tmp_path / "o.fasta"
    write_fasta(entries, path, compression=fmt)  # type: ignore[arg-type]
    data = path.read_bytes()
    assert _format(data) == fmt
    assert read_fasta(path) == entries


@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize("suffix_fmt", FORMATS)
def test_write_explicit_overrides_suffix(
    tmp_path: Path, entries: list[SequenceEntry], fmt: str, suffix_fmt: str
) -> None:
    path = tmp_path / f"o.fasta{SUFFIX[suffix_fmt]}"
    write_fasta(entries, path, compression=fmt)  # type: ignore[arg-type]
    assert _format(path.read_bytes()) == fmt
    assert read_fasta(path) == entries


@pytest.mark.parametrize("fmt", FORMATS)
def test_write_none_is_plain_even_with_suffix(tmp_path: Path, entries: list[SequenceEntry], fmt: str) -> None:
    path = tmp_path / f"o.fasta{SUFFIX[fmt]}"
    write_fasta(entries, path, compression=None)
    assert path.read_text() == TEXT


def test_write_gzip_explicit_is_reproducible(tmp_path: Path, entries: list[SequenceEntry]) -> None:
    a, b = tmp_path / "a.out", tmp_path / "b.out"
    write_fasta(entries, a, compression="gzip")
    write_fasta(entries, b, compression="gzip")
    assert a.read_bytes()[4:8] == b"\x00\x00\x00\x00"  # mtime
    buf = io.BytesIO()
    write_fasta(entries, buf, compression="gzip")
    assert buf.getvalue()[4:8] == b"\x00\x00\x00\x00"
    assert gzip.decompress(a.read_bytes()) == gzip.decompress(b.read_bytes())


# --- writing to a handle ------------------------------------------------------------


def test_write_text_handle_infer_and_none_plain(entries: list[SequenceEntry]) -> None:
    for compression in ("infer", None):
        buf = io.StringIO()
        write_fasta(entries, buf, compression=compression)
        assert buf.getvalue() == TEXT


@pytest.mark.parametrize("fmt", FORMATS)
def test_write_binary_handle_explicit(entries: list[SequenceEntry], fmt: str) -> None:
    buf = io.BytesIO()
    write_fasta(entries, buf, compression=fmt)  # type: ignore[arg-type]
    assert not buf.closed  # the caller's handle stays open
    assert DECOMPRESS[fmt](buf.getvalue()).decode() == TEXT


@pytest.mark.parametrize("fmt", FORMATS)
def test_write_binary_file_handle_explicit(tmp_path: Path, entries: list[SequenceEntry], fmt: str) -> None:
    path = tmp_path / "o.fasta"
    with path.open("wb") as fh:
        write_fasta(entries, fh, compression=fmt)  # type: ignore[arg-type]
        assert not fh.closed
    assert read_fasta(path) == entries


@pytest.mark.parametrize("fmt", FORMATS)
def test_write_text_handle_explicit_raises(entries: list[SequenceEntry], fmt: str) -> None:
    buf = io.StringIO()
    with pytest.raises(FastaError, match="binary handle") as exc:
        write_fasta(entries, buf, compression=fmt)  # type: ignore[arg-type]
    assert any("'wb'" in note for note in exc.value.__notes__)
    assert buf.getvalue() == ""


@pytest.mark.parametrize("fmt", FORMATS)
def test_write_decoy_fasta_compression(tmp_path: Path, fmt: str) -> None:
    src = tmp_path / "t.fasta"
    src.write_text(TEXT)
    dst = tmp_path / "d.fasta"
    write_decoy_fasta(src, dst, method="reverse", compression=fmt)  # type: ignore[arg-type]
    assert _format(dst.read_bytes()) == fmt
    assert len(read_fasta(dst)) == 4
    with pytest.raises(FastaError, match="binary handle"):
        write_decoy_fasta(src, io.StringIO(), method="reverse", compression=fmt)  # type: ignore[arg-type]
    plain = tmp_path / f"d.fasta{SUFFIX[fmt]}"
    write_decoy_fasta(src, plain, method="reverse", compression=None)
    assert plain.read_text().startswith(">sp|")


# --- reading a path -----------------------------------------------------------------


@pytest.mark.parametrize("fmt", FORMATS)
def test_read_explicit_matches(tmp_path: Path, entries: list[SequenceEntry], fmt: str) -> None:
    path = tmp_path / "x.fasta"  # suffix says plain; explicit wins
    path.write_bytes(COMPRESS[fmt](TEXT.encode()))
    assert read_fasta(path, compression=fmt) == entries  # type: ignore[arg-type]
    with FastaReader(path, compression=fmt) as reader:  # type: ignore[arg-type]
        assert list(reader) == entries


@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize("actual", [*FORMATS, None])
def test_read_explicit_wrong_format_raises(tmp_path: Path, fmt: str, actual: str | None) -> None:
    if actual == fmt:
        return
    path = tmp_path / f"x.fasta{SUFFIX[fmt]}"  # suffix agrees with the wrong explicit value
    path.write_bytes(COMPRESS[actual](TEXT.encode()) if actual else TEXT.encode())
    with pytest.raises(FastaParseError, match=f"compression='{fmt}'"):
        read_fasta(path, compression=fmt)  # type: ignore[arg-type]
    with pytest.raises(FastaParseError), FastaReader(path, compression=fmt):  # type: ignore[arg-type]
        pass


@pytest.mark.parametrize("fmt", FORMATS)
def test_read_none_treats_compressed_as_plain(tmp_path: Path, fmt: str) -> None:
    path = tmp_path / f"x.fasta{SUFFIX[fmt]}"
    path.write_bytes(COMPRESS[fmt](TEXT.encode()))
    assert read_fasta(path, compression="infer")  # sniffed
    with pytest.raises(FastaParseError):
        read_fasta(path, compression=None)


def test_read_none_plain_file(tmp_path: Path, entries: list[SequenceEntry]) -> None:
    path = tmp_path / "x.fasta.gz"
    path.write_text(TEXT)
    assert read_fasta(path, compression=None) == entries


# --- reading a handle ---------------------------------------------------------------


def test_read_text_handle_infer_and_none(entries: list[SequenceEntry]) -> None:
    for compression in ("infer", None):
        assert read_fasta(io.StringIO(TEXT), compression=compression) == entries
        with FastaReader(io.StringIO(TEXT), compression=compression) as reader:
            assert list(reader) == entries


@pytest.mark.parametrize("fmt", FORMATS)
def test_read_binary_handle_explicit(entries: list[SequenceEntry], fmt: str) -> None:
    buf = io.BytesIO(COMPRESS[fmt](TEXT.encode()))
    assert read_fasta(buf, compression=fmt) == entries  # type: ignore[arg-type]
    assert not buf.closed
    buf = io.BytesIO(COMPRESS[fmt](TEXT.encode()))
    with FastaReader(buf, compression=fmt) as reader:  # type: ignore[arg-type]
        assert list(reader) == entries
    assert not buf.closed


@pytest.mark.parametrize("fmt", FORMATS)
def test_read_binary_handle_wrong_format(fmt: str) -> None:
    other = next(f for f in FORMATS if f != fmt)
    with pytest.raises(FastaParseError, match=f"compression='{fmt}'"):
        read_fasta(io.BytesIO(COMPRESS[other](TEXT.encode())), compression=fmt)  # type: ignore[arg-type]
    with pytest.raises(FastaParseError, match="plain text"):
        read_fasta(io.BytesIO(TEXT.encode()), compression=fmt)  # type: ignore[arg-type]


@pytest.mark.parametrize("fmt", FORMATS)
def test_read_text_handle_explicit_raises(fmt: str) -> None:
    with pytest.raises(FastaError, match="binary handle"):
        read_fasta(io.StringIO(TEXT), compression=fmt)  # type: ignore[arg-type]


@pytest.mark.parametrize("fmt", FORMATS)
def test_handle_roundtrip(entries: list[SequenceEntry], fmt: str) -> None:
    buf = io.BytesIO()
    write_fasta(entries, buf, compression=fmt)  # type: ignore[arg-type]
    buf.seek(0)
    assert read_fasta(buf, compression=fmt) == entries  # type: ignore[arg-type]


# --- FastaIndex and train_markov_model ----------------------------------------------


@pytest.mark.parametrize("fmt", FORMATS)
def test_index_explicit_compression_refused(tmp_path: Path, fmt: str) -> None:
    path = tmp_path / "x.fasta"
    path.write_bytes(COMPRESS[fmt](TEXT.encode()))
    with pytest.raises(FastaError, match="Cannot index"):
        FastaIndex(path, compression=fmt)  # type: ignore[arg-type]
    with pytest.raises(FastaError, match="Cannot index"):
        FastaIndex(path)  # infer: sniffed and refused, as before


def test_index_none_skips_sniff(tmp_path: Path) -> None:
    path = tmp_path / "x.fasta.gz"
    path.write_text(TEXT)
    index = FastaIndex(path, compression=None)
    assert index["b"].sequence == "PEPTIDE"
    index.write_fai()
    assert FastaIndex.from_fai(path, compression=None)["b"].sequence == "PEPTIDE"


@pytest.mark.parametrize("fmt", FORMATS)
def test_train_markov_model_compression(tmp_path: Path, fmt: str) -> None:
    path = tmp_path / "x.fasta"
    path.write_bytes(COMPRESS[fmt](TEXT.encode()))
    model = train_markov_model(path, order=0, compression=fmt)  # type: ignore[arg-type]
    assert model.metadata["entries"] == 2
    with pytest.raises(FastaParseError):
        train_markov_model(path, order=0, compression=None)


# --- binary handles without an explicit format, reproducible gzip ---------------------


@pytest.mark.parametrize("compression", ["infer", None])
def test_read_binary_handle_plain(tmp_path: Path, entries: list[SequenceEntry], compression: Compression) -> None:
    buf = io.BytesIO(("\ufeff" + TEXT).encode())
    assert read_fasta(buf, compression=compression) == entries
    assert not buf.closed
    buf = io.BytesIO(TEXT.encode())
    with FastaReader(buf, compression=compression) as reader:
        assert list(reader) == entries
    assert not buf.closed
    path = tmp_path / "x.fasta"
    path.write_text(TEXT)
    with path.open("rb", buffering=0) as raw:
        assert read_fasta(raw, compression=compression) == entries
        assert not raw.closed
    with path.open("rb") as fh:
        n = write_decoy_fasta(fh, tmp_path / "d.fasta", method="reverse")
        assert not fh.closed
    assert n == 2


@pytest.mark.parametrize("fmt", FORMATS)
def test_read_binary_handle_infer_does_not_sniff(fmt: str) -> None:
    buf = io.BytesIO(COMPRESS[fmt](TEXT.encode()))
    with pytest.raises(FastaParseError):
        read_fasta(buf)
    assert not buf.closed


@pytest.mark.parametrize("compression", ["infer", None])
def test_write_binary_handle_plain(entries: list[SequenceEntry], compression: Compression) -> None:
    buf = io.BytesIO()
    write_fasta(entries, buf, compression=compression)
    assert not buf.closed
    assert buf.getvalue() == TEXT.encode()
    buf.write(b"after")  # still usable


class _FailingWriter(io.BytesIO):
    def write(self, data: object) -> int:  # type: ignore[override]
        raise OSError("disk full")


@pytest.mark.parametrize("compression", ["infer", None, *FORMATS])
def test_write_error_leaves_handle_open(entries: list[SequenceEntry], compression: Compression) -> None:
    import gc

    buf = _FailingWriter()
    with pytest.raises(OSError, match="disk full"):
        write_fasta(entries, buf, compression=compression)
    gc.collect()
    assert not buf.closed


def test_gzip_bytes_independent_of_name(tmp_path: Path, entries: list[SequenceEntry]) -> None:
    a, b = tmp_path / "a.fasta.gz", tmp_path / "other_name.out"
    write_fasta(entries, a)
    write_fasta(entries, b, compression="gzip")
    with (tmp_path / "handle_name.bin").open("wb") as fh:
        write_fasta(entries, fh, compression="gzip")
    buf = io.BytesIO()
    write_fasta(entries, buf, compression="gzip")
    assert a.read_bytes() == b.read_bytes() == (tmp_path / "handle_name.bin").read_bytes() == buf.getvalue()
    assert a.read_bytes()[3] & 0x08 == 0  # no FNAME field in the header
