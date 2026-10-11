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


def test_from_records_columns_drive_the_header() -> None:
    (rec,) = to_records(io.StringIO(ROUND_TRIP.split(">tr|")[0]))
    rec.update(accession="Q99999", pname="Edited", gname=None, extra="KEY=v", raw_header="ignored", length=0)
    (entry,) = from_records([rec])
    assert entry.raw_header == "sp|Q99999|1433B_HUMAN Edited OS=Homo sapiens OX=9606 PE=1 SV=3 KEY=v"
    # description is only used verbatim when no structured field is set
    (plain,) = from_records([{"identifier": "id", "description": "free text OS=Mouse", "sequence": "MK"}])
    assert (plain.pname, plain.os_name) == ("free text", "Mouse")
    assert from_records([{"identifier": "id"}])[0].sequence == ""


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
