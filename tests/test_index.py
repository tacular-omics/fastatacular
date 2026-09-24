"""FastaIndex: random access by identifier or accession, and samtools .fai read/write."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from fastatacular import (
    FastaError,
    FastaIndex,
    FastaKeyError,
    FastaParseError,
    SequenceEntry,
    read_fasta,
    write_decoy_fasta,
    write_fasta,
)

_REF = json.loads((Path(__file__).parent / "reference" / "header_reference.json").read_text(encoding="utf-8"))


def _reference_fasta() -> str:
    seen: set[str] = set()
    parts = []
    for i, case in enumerate(_REF["cases"]):
        acc = case["expected"].get("accession") or case["expected"]["identifier"]
        if acc in seen:
            continue
        seen.add(acc)
        parts.append(f"{case['header']}\n{'ACDEFGHIKLMNPQRSTVWY' * (i % 5 + 1)}\n")
    return "".join(parts)


# (name, text, whether samtools-compatible line geometry)
FASTAS: list[tuple[str, str, bool]] = [
    (
        "uniprot",
        ">sp|P12345|EX_HUMAN Example OS=Homo sapiens OX=9606 GN=EX PE=1 SV=2\nMKTAYIAKQR\nQISFVKSHFS\nQIS\n"
        ">tr|Q00001|Q00001_MOUSE Other OS=Mus musculus OX=10090\nACDEFGHIK\n",
        True,
    ),
    ("plain_ids", ">a\nMK\n>b desc with > inside\nPEPTIDE\n>c\nA\n", True),
    ("no_final_newline", ">x\nACGT\nAC\n>y\nMMM", True),
    ("crlf", ">sp|P1|A_HUMAN n\r\nMKVL\r\nMK\r\n>sp|P2|B_HUMAN\r\nAAAA\r\n", True),
    ("bom", "﻿>first\nMKV\n>second\nAAA\n", True),
    ("preamble", "# PEFF 1.0\n;comment\n\n>gi|12345|ref|NP_000001.1| x\nMMM\n", True),
    ("blank_after", ">a\nMK\n\n\n>b\nAA\n\n", True),
    ("irregular_lines", ">a\nMK\nMKVL\nM\n>b\nAA\n", False),
    ("inner_blank_and_comment", ">a\nMK\n\nMK\n;c\nMK\n>b\nA A\n", False),
    ("reference_headers", _reference_fasta(), True),
]


def _write(tmp_path: Path, name: str, text: str) -> Path:
    p = tmp_path / f"{name}.fasta"
    p.write_bytes(text.encode("utf-8"))
    return p


def _key(e: SequenceEntry, key: str = "identifier") -> str:
    return e.identifier if key == "identifier" else e.accession or e.identifier


@pytest.mark.parametrize("key", ["identifier", "accession"])
@pytest.mark.parametrize(("name", "text", "regular"), FASTAS, ids=[f[0] for f in FASTAS])
def test_index_matches_full_parse(tmp_path: Path, name: str, text: str, regular: bool, key: str) -> None:
    path = _write(tmp_path, name, text)
    entries = read_fasta(path)
    index = FastaIndex(path, key=key)  # type: ignore[arg-type]
    assert index.key == key
    assert len(index) == len(entries)
    assert list(index) == [_key(e, key) for e in entries]
    for e in entries:
        assert _key(e, key) in index
        assert index[_key(e, key)] == e
        assert index.identifier(_key(e, key)) == e.identifier
    assert dict(index.items()) == {_key(e, key): e for e in entries}


@pytest.mark.parametrize(("name", "text", "regular"), FASTAS, ids=[f[0] for f in FASTAS])
def test_fai_roundtrip(tmp_path: Path, name: str, text: str, regular: bool) -> None:
    path = _write(tmp_path, name, text)
    index = FastaIndex(path)
    if not regular:
        with pytest.raises(FastaError, match="cannot be written to a .fai"):
            index.write_fai()
        assert not Path(f"{path}.fai").exists()
        return
    fai = index.write_fai()
    assert fai == Path(f"{path}.fai")
    loaded = FastaIndex.from_fai(path)
    assert list(loaded) == list(index)
    assert all(loaded.locate(k) == index.locate(k) for k in index)
    assert all(loaded[k] == index[k] for k in index)
    # each row reads back the right residues, as samtools faidx would
    raw = path.read_bytes()
    for line, e in zip(fai.read_text().splitlines(), read_fasta(path), strict=True):
        name_, length, offset, linebases, linewidth = line.split("\t")
        assert name_ == e.identifier and int(length) == len(e.sequence)
        n, off, lb, lw = int(length), int(offset), int(linebases), int(linewidth)
        full = (n - 1) // lb
        chunk = raw[off : off + full * lw + (n - full * lb)]
        assert chunk.replace(b"\r", b"").replace(b"\n", b"").decode() == e.sequence


def test_samtools_fai_values(tmp_path: Path) -> None:
    path = _write(tmp_path, "s", ">one desc\nACGTA\nCGTAC\nGT\n>two\nAAAA\n")
    FastaIndex(path).write_fai()
    # the values `samtools faidx` writes for this file
    assert Path(f"{path}.fai").read_text() == "one\t12\t10\t5\t6\ntwo\t4\t30\t4\t5\n"


def test_fai_custom_path_and_written_by_write_fasta(tmp_path: Path) -> None:
    src = _write(tmp_path, "src", FASTAS[0][1])
    dst = tmp_path / "dst.fasta"
    write_fasta(read_fasta(src), dst)
    out = FastaIndex(dst).write_fai(tmp_path / "x.fai")
    assert FastaIndex.from_fai(dst, out)["sp|P12345|EX_HUMAN"] == read_fasta(src)[0]
    assert FastaIndex.from_fai(dst, out, key="accession")["P12345"] == read_fasta(src)[0]


def test_locate(tmp_path: Path) -> None:
    text = ">a\nMK\n>b\nAAA\n"
    index = FastaIndex(_write(tmp_path, "l", text))
    assert index.locate("a") == (0, 6)
    assert index.locate("b") == (6, 7)
    assert "FastaIndex(" in repr(index) and "2 entries" in repr(index)


def test_missing_key(tmp_path: Path) -> None:
    index = FastaIndex(_write(tmp_path, "m", ">a\nMK\n"))
    for lookup in (index.__getitem__, index.locate, index.identifier):
        with pytest.raises(FastaKeyError) as info:
            lookup("nope")
        assert isinstance(info.value, KeyError) and isinstance(info.value, FastaError)
        assert info.value.args[0] == "nope" and str(info.value) == "nope"
    assert "nope" not in index
    assert index.get("nope") is None


def test_bad_key_mode(tmp_path: Path) -> None:
    path = _write(tmp_path, "k", ">a\nMK\n")
    with pytest.raises(FastaError, match="key must be"):
        FastaIndex(path, key="name")  # type: ignore[arg-type]
    with pytest.raises(FastaError, match="key must be"):
        FastaIndex.from_fai(path, key="name")  # type: ignore[arg-type]


def test_duplicate_accession_names_the_first(tmp_path: Path) -> None:
    text = ">sp|P1|A_HUMAN\nMK\n>tr|P2|B\nMK\n>tr|P1|C_HUMAN\nMK\n>x|P2\nMK\n"
    path = _write(tmp_path, "d", text)
    with pytest.raises(FastaError, match=r"Duplicate accession 'P1'") as info:
        FastaIndex(path, key="accession")
    assert any('key="identifier"' in n for n in info.value.__notes__)
    assert list(FastaIndex(path)) == ["sp|P1|A_HUMAN", "tr|P2|B", "tr|P1|C_HUMAN", "x|P2"]


def test_duplicate_identifier(tmp_path: Path) -> None:
    with pytest.raises(FastaError, match=r"Duplicate identifier 'a'"):
        FastaIndex(_write(tmp_path, "d", ">a\nMK\n>a x\nMK\n"))


def test_decoy_database(tmp_path: Path) -> None:
    # Targets and their decoys share accessions; the default key tells them apart.
    src = _write(tmp_path, "t", FASTAS[0][1])
    dst = tmp_path / "td.fasta"
    write_decoy_fasta(src, dst, method="reverse")
    entries = read_fasta(dst)
    index = FastaIndex(dst)
    assert list(index) == [e.identifier for e in entries]
    assert index["DECOY_sp|P12345|EX_HUMAN"] == entries[2]
    assert index["sp|P12345|EX_HUMAN"] == entries[0]
    index.write_fai()
    assert list(FastaIndex.from_fai(dst)) == list(index)
    with pytest.raises(FastaError, match="Duplicate accession"):
        FastaIndex(dst, key="accession")


def test_reverse_prefix_database(tmp_path: Path) -> None:
    # IP2/ProLuCID style: Reverse_sp|P1|X next to sp|P1|X.
    text = ">sp|P1|A_HUMAN x\nMKV\n>Reverse_sp|P1|A_HUMAN x\nVKM\n>contaminant_K1\nAA\n"
    index = FastaIndex(_write(tmp_path, "r", text))
    assert list(index) == ["sp|P1|A_HUMAN", "Reverse_sp|P1|A_HUMAN", "contaminant_K1"]
    assert index["Reverse_sp|P1|A_HUMAN"].sequence == "VKM"
    assert "sp|P1|A_HUMAN" in index


def test_empty_file(tmp_path: Path) -> None:
    path = _write(tmp_path, "e", "")
    index = FastaIndex(path)
    assert len(index) == 0 and list(index) == []
    index.write_fai()
    assert Path(f"{path}.fai").read_text() == ""
    assert len(FastaIndex.from_fai(path)) == 0


def test_comments_only(tmp_path: Path) -> None:
    assert len(FastaIndex(_write(tmp_path, "c", "# PEFF 1.0\n; x\n\n"))) == 0


def test_sequence_before_header(tmp_path: Path) -> None:
    with pytest.raises(FastaParseError, match="before any '>' header"):
        FastaIndex(_write(tmp_path, "b", "MKV\n>a\nMK\n"))


def test_empty_header(tmp_path: Path) -> None:
    with pytest.raises(FastaParseError, match="Empty FASTA header"):
        FastaIndex(_write(tmp_path, "h", ">a\nMK\n>  \nMK\n"))


def test_undecodable_header(tmp_path: Path) -> None:
    p = tmp_path / "u.fasta"
    p.write_bytes(b">a\xff\nMK\n")
    with pytest.raises(FastaParseError, match="Cannot decode"):
        FastaIndex(p)


def test_entry_errors_are_raised_on_access(tmp_path: Path) -> None:
    index = FastaIndex(_write(tmp_path, "n", ">a\n>b\nMK\n"))
    assert index["b"].sequence == "MK"
    with pytest.raises(FastaParseError, match="Entry 'a' at byte 0.*no sequence data"):
        index["a"]
    with pytest.raises(FastaError, match="no sequence data"):
        index.write_fai()


def test_undecodable_entry(tmp_path: Path) -> None:
    p = tmp_path / "u.fasta"
    p.write_bytes(b">a\nMK\xff\n")
    with pytest.raises(FastaParseError, match="Cannot decode entry 'a'"):
        FastaIndex(p)["a"]


def test_changed_file(tmp_path: Path) -> None:
    path = _write(tmp_path, "c", ">a\nMK\n>b\nAA\n")
    index = FastaIndex(path)
    path.write_text(">z\nMK\n>b\nAA\n")
    with pytest.raises(FastaError, match="changed since it was indexed"):
        index["a"]


def test_compressed_is_rejected(tmp_path: Path) -> None:
    p = tmp_path / "x.fasta.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(">a\nMK\n")
    with pytest.raises(FastaError, match="gz-compressed") as info:
        FastaIndex(p)
    assert any("gunzip" in n for n in info.value.__notes__)
    with pytest.raises(FastaError, match="gz-compressed"):
        FastaIndex.from_fai(p)


@pytest.mark.parametrize(
    "fai",
    [
        "a\t2\t3\t2\n",  # 4 columns
        "a\tx\t3\t2\t3\n",  # not a number
        "a\t2\t4\t2\t3\n",  # offset not after a newline
        "b\t2\t3\t2\t3\n",  # wrong name
        "a\t2\t999\t2\t3\n",  # past the end
        "c\t2\t6\t2\t3\n",  # offset after a sequence line, not a header
        "a\t2\t9\t2\t3\n",  # offset of another entry
    ],
)
def test_bad_fai(tmp_path: Path, fai: str) -> None:
    path = _write(tmp_path, "f", ">a\nMK\n>c\nMK\n")
    Path(f"{path}.fai").write_text(fai)
    with pytest.raises(FastaError):
        FastaIndex.from_fai(path)


def test_fai_for_empty_fasta(tmp_path: Path) -> None:
    path = _write(tmp_path, "e", "")
    Path(f"{path}.fai").write_text("a\t2\t3\t2\t3\n")
    with pytest.raises(FastaError, match="is empty"):
        FastaIndex.from_fai(path)


def test_fai_duplicate_accession(tmp_path: Path) -> None:
    path = _write(tmp_path, "d", ">sp|P1|A\nMK\n>tr|P1|B\nMK\n")
    Path(f"{path}.fai").write_text("sp|P1|A\t2\t9\t2\t3\ntr|P1|B\t2\t21\t2\t3\n")
    assert list(FastaIndex.from_fai(path)) == ["sp|P1|A", "tr|P1|B"]
    with pytest.raises(FastaError, match="Duplicate accession 'P1'"):
        FastaIndex.from_fai(path, key="accession")


def test_fai_blank_lines_ignored_and_bom(tmp_path: Path) -> None:
    path = _write(tmp_path, "b", "﻿>a\nMK\n")
    Path(f"{path}.fai").write_text("a\t2\t6\t2\t3\n\n")
    assert FastaIndex.from_fai(path)["a"].sequence == "MK"


def test_fai_rejects_header_only_last_entry(tmp_path: Path) -> None:
    index = FastaIndex(_write(tmp_path, "h", ">a\nMK\n>b"))
    with pytest.raises(FastaError, match="no sequence data"):
        index.write_fai()


def test_fai_trailing_spaces_like_samtools(tmp_path: Path) -> None:
    # samtools counts trailing spaces in the line width but not in the residues.
    path = _write(tmp_path, "s", ">a\nMKV  \nMKV  \nMK \n")
    FastaIndex(path).write_fai()
    assert Path(f"{path}.fai").read_text(encoding="utf-8") == "a\t8\t3\t3\t6\n"
    assert FastaIndex.from_fai(path)["a"].sequence == "MKVMKVMK"


def test_fai_rejects_long_last_line(tmp_path: Path) -> None:
    index = FastaIndex(_write(tmp_path, "l", ">a\nMK\nMKV\n"))
    with pytest.raises(FastaError, match="last line longer"):
        index.write_fai()


def test_undecodable_preamble(tmp_path: Path) -> None:
    p = tmp_path / "p.fasta"
    p.write_bytes(b"; \xff\n>a\nMK\n")
    with pytest.raises(FastaParseError, match="Cannot decode the start"):
        FastaIndex(p)
