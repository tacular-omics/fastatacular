import io
import re
from pathlib import Path

import pytest

from fastatacular import (
    RECORD_KEYS,
    FastaReader,
    FastaWriteError,
    SequenceEntry,
    from_records,
    read_fasta,
    to_records,
    write_fasta,
)

FASTA = """>sp|P12345|EX_HUMAN Example protein OS=Homo sapiens OX=9606 GN=EX PE=1 SV=2 XX=1 YY=a b
MKTAYIAKQR
QISFVKSHFS
>tr|Q00001|Q00001_MOUSE Other OS=Mus musculus OX=10090
ACDEFGHIK
>plain_id
PEPTIDE
>gi|12345|ref|NP_000001.1| something
MMM
"""

README = Path(__file__).resolve().parents[1] / "README.md"


@pytest.fixture
def fasta_path(tmp_path: Path) -> Path:
    p = tmp_path / "human.fasta"
    p.write_text(FASTA)
    return p


def test_keys_are_stable(fasta_path: Path) -> None:
    records = to_records(fasta_path)
    assert len(records) == 4
    for rec in records:
        assert tuple(rec) == RECORD_KEYS


def test_roundtrip_against_entries(fasta_path: Path) -> None:
    entries = read_fasta(fasta_path)
    for rec, e in zip(to_records(fasta_path), entries, strict=True):
        for key in RECORD_KEYS:
            if key == "extra":
                expected = " ".join(f"{k}={v}" for k, v in e.extra.items()) or None
            elif key == "length":
                expected = len(e.sequence)
            else:
                expected = getattr(e, key)
            assert rec[key] == expected
        # the flat record rebuilds the same entry
        rebuilt = read_fasta(io.StringIO(f">{rec['raw_header']}\n{rec['sequence']}\n"))[0]
        assert rebuilt == e


def test_values(fasta_path: Path) -> None:
    first = to_records(fasta_path)[0]
    assert first["accession"] == "P12345"
    assert first["ncbi_tax_id"] == 9606
    assert first["extra"] == "XX=1 YY=a b"
    assert first["length"] == 20
    plain = to_records(fasta_path)[2]
    assert plain["prefix"] is None and plain["extra"] is None and plain["description"] is None


def test_sources_agree(fasta_path: Path) -> None:
    expected = to_records(fasta_path)
    assert to_records(str(fasta_path)) == expected
    assert to_records(io.StringIO(FASTA)) == expected
    assert to_records(read_fasta(fasta_path)) == expected
    assert to_records(iter(read_fasta(fasta_path))) == expected
    with FastaReader(fasta_path) as reader:
        assert reader.to_records() == expected
    assert [e.to_record() for e in read_fasta(fasta_path)] == expected


def test_empty(tmp_path: Path) -> None:
    p = tmp_path / "empty.fasta"
    p.write_text("")
    assert to_records(p) == []
    assert to_records([]) == []
    assert to_records(io.StringIO("")) == []


def test_record_from_constructed_entry() -> None:
    rec = SequenceEntry(identifier="x", sequence="AC").to_record()
    assert rec["identifier"] == "x" and rec["length"] == 2 and rec["raw_header"] == ""


def _readme_block() -> str:
    text = README.read_text(encoding="utf-8")
    section = text.split("## Tables with pandas or polars", 1)[1]
    return re.search(r"```python\n(.*?)```", section, re.S).group(1)  # type: ignore[union-attr]


def test_pandas(fasta_path: Path) -> None:
    pd = pytest.importorskip("pandas")
    df = pd.DataFrame(to_records(fasta_path))
    assert list(df.columns) == list(RECORD_KEYS)
    assert len(df) == 4


def test_polars(fasta_path: Path) -> None:
    pl = pytest.importorskip("polars")
    df = pl.DataFrame(to_records(fasta_path))
    assert df.columns == list(RECORD_KEYS)
    assert df.height == 4


def test_readme_example(fasta_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    pytest.importorskip("pandas")
    pytest.importorskip("polars")
    monkeypatch.chdir(fasta_path.parent)
    ns: dict[str, object] = {}
    exec(_readme_block(), ns)
    assert ns["human"].height == 1  # type: ignore[attr-defined]
    assert read_fasta("human_out.fasta") == read_fasta(fasta_path)[:1]


# One entry per header style; every one is written by write_fasta in this exact form.
ROUND_TRIP = """\
>sp|P31946|1433B_HUMAN 14-3-3 protein beta/alpha OS=Homo sapiens OX=9606 GN=YWHAB PE=1 SV=3
MTMDKSELVQ
>tr|A0A024R161|A0A024R161_HUMAN Guanine nucleotide-binding protein OS=Homo sapiens OX=9606 PE=3 SV=1
MGAQ
>myid some desc
MK
>gi|4557757|ref|NP_000240.1| mutL homolog 1 [Homo sapiens]
MSF
>pdb|1MBA|A Chain A, Myoglobin
GLS
>x desc OS=Mus musculus FOO=bar BAZ=a b
KK
>empty_seq Nothing here
>plain
PEP
>onlykv OS=Human
"""


def _write(entries: list[SequenceEntry]) -> str:
    out = io.StringIO()
    write_fasta(entries, out)
    return out.getvalue()


def test_from_records_round_trip() -> None:
    entries = read_fasta(io.StringIO(ROUND_TRIP))
    rebuilt = from_records(to_records(entries))
    assert rebuilt == entries
    assert _write(rebuilt) == ROUND_TRIP
    assert rebuilt[6].sequence == ""


def test_from_records_nan_and_float_ints() -> None:
    """A data frame turns None into NaN and int columns with gaps into floats."""
    nan = float("nan")
    records = [
        {k: nan if v is None else float(v) if k in ("ncbi_tax_id", "pe", "sv") else v for k, v in rec.items()}
        for rec in to_records(io.StringIO(ROUND_TRIP))
    ]
    assert _write(from_records(records)) == ROUND_TRIP


# Unedited records write raw_header byte for byte, whatever its spacing or key order.
VERBATIM = [
    "sp|Q3SA23|Q3SA23_9VIRU Polyprotein OS=Hepatitis E virus 1  OX=0000 GN=g PE=4 SV=1",
    "tr|A0A|NAME|extra Foo OS=H",
    "gi|4557757|ref|NP_000240.1| mutL homolog 1 [Homo sapiens]",
    "pdb|1MBA|A Chain A",
    "plain_id",
    "plain\tdesc after a tab",
    "sp|P1|A_HUMAN X OS=Human OS=Mouse",
    "x Name XX=1 OS=H OX=9606",
]


@pytest.mark.parametrize("header", VERBATIM)
def test_from_records_unedited_is_verbatim(header: str) -> None:
    text = f">{header}\nMK\n"
    assert _write(from_records(to_records(io.StringIO(text)))) == text


def test_repeated_key_keeps_first_value() -> None:
    (entry,) = read_fasta(io.StringIO(">sp|P1|A X OS=Human OX=9606 XX=1 OS=Mouse XX=2\nMK\n"))
    assert (entry.os_name, entry.extra) == ("Human", {"XX": "1"})


def _edit(header: str, **changes: object) -> str:
    (rec,) = to_records(io.StringIO(f">{header}\nMK\n"))
    rec.update(changes)
    return from_records([rec])[0].raw_header


@pytest.mark.parametrize(
    ("header", "changes", "expected"),
    [
        ("sp|P1|A_HUMAN X OS=H", {"identifier": "DECOY_sp|P1|A_HUMAN"}, "DECOY_sp|P1|A_HUMAN X OS=H"),
        ("gi|4557757|ref|NP_000240.1| mutL", {"accession": "NP_999.1"}, "gi|4557757|ref|NP_999.1| mutL"),
        ("tr|A0A|NAME|extra Foo", {"accession": "B0B"}, "tr|B0B|NAME|extra Foo"),
        ("pdb|1MBA|A Chain", {"prefix": "rev_pdb"}, "rev_pdb|1MBA|A Chain"),
        ("gi|229552|pdb|1MBA|A Chain", {"accession": "2XYZ_B"}, "gi|229552|pdb|2XYZ|B Chain"),
        ("plain\tdesc", {"identifier": "p2"}, "p2\tdesc"),
        ("sp|P1|A_HUMAN X OS=H GN=g", {"description": "Real desc"}, "sp|P1|A_HUMAN Real desc"),
        ("sp|P1|A Poly OS=Virus 1  OX=0000 GN=g", {"gname": "NEW"}, "sp|P1|A Poly OS=Virus 1  OX=0000 GN=NEW"),
        ("sp|P1|A Poly OS=Virus 1  OX=0000 GN=g", {"gname": None}, "sp|P1|A Poly OS=Virus 1  OX=0000"),
        ("sp|P1|A Poly OS=Virus 1  OX=0000", {"pname": "New", "pe": 1.0}, "sp|P1|A New OS=Virus 1  OX=0000 PE=1"),
        ("x Name XX=1 OS=H OX=9606", {"extra": "XX=2 YY=3"}, "x Name XX=2 OS=H OX=9606 YY=3"),
        ("sp|P1|A X OS=H OS=M", {"os_name": ""}, "sp|P1|A X"),
    ],
)
def test_from_records_edits_only_the_changed_token(header: str, changes: dict[str, object], expected: str) -> None:
    assert _edit(header, **changes) == expected


@pytest.mark.parametrize(
    ("record", "match"),
    [
        ({"identifier": "sp|P2|A", "accession": "P3"}, "identifier 'sp|P2|A' and accession"),
        ({"description": "Real desc", "gname": "ABC"}, "description and gname"),
        ({"accession": "P1_X", "entry_name": "N"}, "cannot be written"),
    ],
)
def test_from_records_conflicting_edits(record: dict[str, object], match: str) -> None:
    with pytest.raises(FastaWriteError, match=re.escape(match)):
        _edit("pdb|1MBA|A Chain" if "entry_name" in record else "sp|P1|A X GN=g", **record)


def test_from_records_without_raw_header() -> None:
    (plain,) = from_records([{"identifier": "id", "description": "free text OS=Mouse", "sequence": "MK"}])
    assert (plain.pname, plain.os_name) == ("free text", "Mouse")
    (built,) = from_records([{"prefix": "sp", "accession": "P1", "pname": "N", "gname": "G", "description": ""}])
    assert built.raw_header == "sp|P1 N GN=G"
    assert from_records([{"identifier": "id"}])[0].sequence == ""
    with pytest.raises(FastaWriteError, match="description and gname"):
        from_records([{"identifier": "x", "description": "Real desc", "gname": "ABC"}])


@pytest.mark.parametrize(
    ("record", "match"),
    [
        ({"identifier": "a", "pe": "high"}, "pe must be an int"),
        ({"identifier": "a", "ncbi_tax_id": 9606.5}, "ncbi_tax_id must be an int"),
        ({"identifier": "a", "pe": True}, "pe must be an int"),
        ({"description": "x"}, "no identifier"),
        ({"identifier": "a", "extra": "no pairs"}, "extra must be"),
        ({"identifier": "a", "pname": "x GN=Y"}, "pname"),
        ({"identifier": "a b"}, "whitespace"),
    ],
)
def test_from_records_bad_values(record: dict[str, object], match: str) -> None:
    with pytest.raises(FastaWriteError, match=match) as info:
        from_records([{"identifier": "ok"}, record])
    assert info.value.index == 1


def test_from_records_pandas_round_trip() -> None:
    pd = pytest.importorskip("pandas")
    entries = read_fasta(io.StringIO(ROUND_TRIP))
    df = pd.DataFrame(to_records(entries))
    assert from_records(df.to_dict("records")) == entries
    df = df.convert_dtypes()  # pandas.NA and nullable Int64
    assert from_records(df.to_dict("records")) == entries
