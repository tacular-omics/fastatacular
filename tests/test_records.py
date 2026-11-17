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


def test_polars_filter_round_trip(fasta_path: Path, tmp_path: Path) -> None:
    pl = pytest.importorskip("polars")
    human = pl.DataFrame(to_records(fasta_path)).filter(pl.col("ncbi_tax_id") == 9606)
    assert human.height == 1
    out = tmp_path / "human_out.fasta"
    write_fasta(from_records(human.to_dicts()), out)
    assert read_fasta(out) == read_fasta(fasta_path)[:1]


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
        (
            "NP_000549.1 hemoglobin subunit beta [Homo sapiens]",
            {"prefix": "DECOY"},
            "DECOY|NP_000549.1 hemoglobin subunit beta [Homo sapiens]",
        ),
        (
            "UniRef90_P69905 Hemoglobin subunit alpha n=2 Tax=Homo sapiens TaxID=9606",
            {"prefix": "DECOY"},
            "DECOY|UniRef90_P69905 Hemoglobin subunit alpha n=2 Tax=Homo sapiens TaxID=9606",
        ),
        (
            "ENSP00000322421.5 pep gene_symbol:HBB",
            {"prefix": "DECOY"},
            "DECOY|ENSP00000322421.5 pep gene_symbol:HBB",
        ),
    ],
)
def test_from_records_edits_only_the_changed_token(header: str, changes: dict[str, object], expected: str) -> None:
    assert _edit(header, **changes) == expected


@pytest.mark.parametrize(
    ("header", "identifier"),
    [
        ("1MBA_1|Chain A|MYOGLOBIN|Aplysia limacina (6502)", "1MBA_1"),
        ("4HHB_1|Chains A, C|Hemoglobin subunit alpha|Homo sapiens (9606)", "4HHB_1"),
    ],
)
@pytest.mark.parametrize(
    ("changes", "added"),
    [
        ({"gname": "MB"}, " GN=MB"),
        ({"os_name": "Example organism"}, " OS=Example organism"),
        ({"ncbi_tax_id": 101}, " OX=101"),
        ({"pe": 1}, " PE=1"),
        ({"sv": 2}, " SV=2"),
        ({"extra": {"XX": "value"}}, " XX=value"),
    ],
)
def test_from_records_edits_rcsb_entity_description(
    header: str, identifier: str, changes: dict[str, object], added: str
) -> None:
    """A structured-field edit appends its key and leaves the rest of the header alone."""
    [entry] = read_fasta(io.StringIO(f">{header}\nMKV\n"))
    [edited] = from_records([entry.to_record() | changes])
    assert edited.raw_header == f"{header}{added}"
    assert (edited.identifier, edited.accession) == (identifier, identifier)
    for field, value in changes.items():
        assert getattr(edited, field) == value

    output = io.StringIO()
    write_fasta([edited], output)
    assert output.getvalue() == f">{header}{added}\nMKV\n"


@pytest.mark.parametrize(
    ("header", "expected"),
    [
        ("1MBA_1|Chain A|MYOGLOBIN|Aplysia limacina (6502)", "1MBA_1|Chain A|NEW"),
        ("4HHB_1|Chains A, C|Hemoglobin subunit alpha|Homo sapiens (9606)", "4HHB_1|Chain A|NEW"),
    ],
)
def test_from_records_replaces_rcsb_entity_description(header: str, expected: str) -> None:
    """Replacing the description keeps the pipe, so the header still reads back as one entity."""
    [entry] = read_fasta(io.StringIO(f">{header}\nMKV\n"))
    [edited] = from_records([entry.to_record() | {"description": "Chain A|NEW"}])
    assert edited.raw_header == expected
    assert (edited.identifier, edited.accession, edited.description) == (
        expected.split("|")[0],
        expected.split("|")[0],
        "Chain A|NEW",
    )

    output = io.StringIO()
    write_fasta([edited], output)
    assert output.getvalue() == f">{expected}\nMKV\n"


@pytest.mark.parametrize("changes", [{"prefix": "DECOY"}, {"pname": "MYOGLOBIN"}])
def test_from_records_rejects_edits_that_break_rcsb_entity_shape(changes: dict[str, str]) -> None:
    with pytest.raises(FastaWriteError):
        _edit("1MBA_1|Chain A|MYOGLOBIN|Aplysia limacina (6502)", **changes)


@pytest.mark.parametrize(
    ("identifier", "description", "accession"),
    [
        ("1MBA_A", "Chain A, MYOGLOBIN", None),
        ("1MBA_1", "Chain A, MYOGLOBIN", None),
        ("P1", "Chains of command", None),
        ("P1", "Chain A|MYOGLOBIN|Aplysia limacina", None),
        ("NP_111111.1", "Chain A, MYOGLOBIN", "NP_111111.1"),
        ("ENSP00000322421.5", "Chains A, B|MYOGLOBIN", "ENSP00000322421.5"),
    ],
)
def test_from_records_without_raw_header_keeps_the_space_before_a_chain_description(
    identifier: str, description: str, accession: str | None
) -> None:
    """A description that merely starts with 'Chain ' does not make a pipe-delimited id.

    Only an RCSB PDB entity header puts the identifier before a pipe, and a record built
    from columns says it is one by setting ``accession`` to the entity id. Without that,
    1.5.0 wrote ``>ID|Chain ...`` and read the identifier back as an accession.
    """
    [entry] = from_records([{"identifier": identifier, "description": description, "sequence": "MK"}])
    assert (entry.identifier, entry.description, entry.accession) == (identifier, description, accession)

    output = io.StringIO()
    write_fasta([entry], output)
    assert output.getvalue() == f">{identifier} {description}\nMK\n"
    assert read_fasta(io.StringIO(output.getvalue())) == [entry]


@pytest.mark.parametrize(
    ("identifier", "description", "header"),
    [
        # Only the pipe reads back an accession for these identifiers, so it is written.
        ("1MBA_1", "Chain A|MYOGLOBIN|Aplysia limacina", ">1MBA_1|Chain A|MYOGLOBIN|Aplysia limacina"),
        ("4HHB_1", "Chains A, C|Hemoglobin subunit alpha", ">4HHB_1|Chains A, C|Hemoglobin subunit alpha"),
        ("P1", "Chain A, MYOGLOBIN", ">P1|Chain A, MYOGLOBIN"),
        # A bare accession already reads back from the space-separated header, so it keeps it.
        ("NP_111111.1", "Chain A, MYOGLOBIN", ">NP_111111.1 Chain A, MYOGLOBIN"),
        ("UniRef90_P69905", "Chain A, MYOGLOBIN", ">UniRef90_P69905 Chain A, MYOGLOBIN"),
    ],
)
def test_from_records_writes_the_entity_pipe_only_when_the_accession_needs_it(
    identifier: str, description: str, header: str
) -> None:
    """With ``accession`` set to the identifier, the delimiter is whatever reads it back."""
    record = {"identifier": identifier, "accession": identifier, "description": description, "sequence": "MK"}
    [entry] = from_records([record])
    assert (entry.identifier, entry.accession, entry.description) == (identifier, identifier, description)

    output = io.StringIO()
    write_fasta([entry], output)
    assert output.getvalue() == f"{header}\nMK\n"


@pytest.mark.parametrize(
    "header",
    [
        "1MBA_1|Chain A|MYOGLOBIN|Aplysia limacina (6502)",
        "4HHB_1|Chains A, C|Hemoglobin subunit alpha|Homo sapiens (9606)",
    ],
)
def test_from_records_rebuilds_an_rcsb_entity_header_without_raw_header(header: str) -> None:
    """to_records minus raw_header still rebuilds the entity header byte for byte."""
    [entry] = read_fasta(io.StringIO(f">{header}\nMKV\n"))
    record = entry.to_record()
    del record["raw_header"]
    [rebuilt] = from_records([record])
    assert rebuilt == entry

    output = io.StringIO()
    write_fasta([rebuilt], output)
    assert output.getvalue() == f">{header}\nMKV\n"


def test_from_records_bare_accession_prefix_edit() -> None:
    """Setting prefix on an NCBI, UniRef, or Ensembl entry with raw_header builds db|ACC."""
    for header in (
        "NP_000549.1 hemoglobin subunit beta [Homo sapiens]",
        "UniRef90_P69905 Hemoglobin subunit alpha n=2 Tax=Homo sapiens TaxID=9606",
        "ENSP00000322421.5 pep gene_symbol:HBB",
    ):
        [entry] = read_fasta(io.StringIO(f">{header}\nMKV\n"))
        rec = entry.to_record()
        rec["prefix"] = "DECOY"
        [out] = from_records([rec])
        assert out.identifier == f"DECOY|{entry.accession}"
        assert out.prefix == "DECOY"
        assert out.accession == entry.accession
        assert out.raw_header == f"DECOY|{header}"


@pytest.mark.parametrize(
    ("record", "match"),
    [
        ({"identifier": "sp|P2|A", "accession": "P3"}, "identifier 'sp|P2|A' and accession"),
        ({"description": "Real desc", "gname": "ABC"}, "description and gname"),
        ({"accession": "P1_X", "entry_name": "N"}, "cannot be written"),
        ({"accession": "NP_999.1"}, "cannot be written"),
    ],
)
def test_from_records_conflicting_edits(record: dict[str, object], match: str) -> None:
    header = (
        "NP_000549.1 hemoglobin"
        if record.get("accession") == "NP_999.1"
        else ("pdb|1MBA|A Chain" if "entry_name" in record else "sp|P1|A X GN=g")
    )
    with pytest.raises(FastaWriteError, match=re.escape(match)):
        _edit(header, **record)


def test_from_records_without_raw_header() -> None:
    (plain,) = from_records([{"identifier": "id", "description": "free text OS=Mouse", "sequence": "MK"}])
    assert (plain.pname, plain.os_name) == ("free text", "Mouse")
    (built,) = from_records([{"prefix": "sp", "accession": "P1", "pname": "N", "gname": "G", "description": ""}])
    assert built.raw_header == "sp|P1 N GN=G"
    assert from_records([{"identifier": "id"}])[0].sequence == ""
    with pytest.raises(FastaWriteError, match="description and gname"):
        from_records([{"identifier": "x", "description": "Real desc", "gname": "ABC"}])
    with pytest.raises(FastaWriteError, match="identifier 'NP_000549.1' and accession were both set and disagree"):
        from_records([{"identifier": "NP_000549.1", "accession": "NP_111111.2", "sequence": "MKV"}])
    with pytest.raises(FastaWriteError, match="identifier 'plain_id' and accession were both set and disagree"):
        from_records([{"identifier": "plain_id", "accession": "some_acc"}])

    # Bare accessions without raw_header roundtrip to the same identifier and accession (issue #59).
    (repro,) = from_records([{"identifier": "NP_000549.1", "sequence": "MKV", "accession": "NP_000549.1"}])
    assert (repro.identifier, repro.accession) == ("NP_000549.1", "NP_000549.1")

    headers = [
        ">sp|P31946|1433B_HUMAN 14-3-3 protein beta/alpha OS=Homo sapiens OX=9606 GN=YWHAB PE=1 SV=3\nMKV\n",
        ">pdb|1MBA|A Chain A, MYOGLOBIN\nMKV\n",
        ">gi|12345|ref|NP_000001.1| some protein\nMKV\n",
        ">plainid some description\nMKV\n",
        ">NP_000549.1 hemoglobin subunit beta [Homo sapiens]\nMKV\n",
        ">UniRef90_P69905 Hemoglobin subunit alpha n=2 Tax=Homo sapiens TaxID=9606\nMKV\n",
        ">ENSP00000322421.5 pep gene_symbol:HBB\nMKV\n",
        ">1MBA_1|Chain A|MYOGLOBIN|Aplysia limacina\nMKV\n",
    ]
    for text in headers:
        [entry] = read_fasta(io.StringIO(text))
        rec = entry.to_record()
        del rec["raw_header"]
        [out] = from_records([rec])
        assert (out.identifier, out.accession) == (entry.identifier, entry.accession)


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


def test_from_records_pandas_frame_hint() -> None:
    pd = pytest.importorskip("pandas")
    df = pd.DataFrame(to_records(io.StringIO(ROUND_TRIP)))
    with pytest.raises(FastaWriteError, match="pandas DataFrame") as info:
        from_records(df)
    assert info.value.hint is not None
    assert 'df.to_dict("records")' in info.value.hint
    assert len(str(info.value)) < 200


def test_from_records_polars_frame_hint() -> None:
    pl = pytest.importorskip("polars")
    df = pl.DataFrame(to_records(io.StringIO(ROUND_TRIP)))
    with pytest.raises(FastaWriteError, match="polars DataFrame") as info:
        from_records(df)
    assert info.value.hint is not None
    assert "df.iter_rows(named=True)" in info.value.hint
    assert len(str(info.value)) < 200
    assert from_records(df.iter_rows(named=True)) == read_fasta(io.StringIO(ROUND_TRIP))


def test_from_records_non_mapping_repr_is_capped() -> None:
    pl = pytest.importorskip("polars")
    series = pl.Series("identifier", [f"id{i}" for i in range(500)])
    with pytest.raises(FastaWriteError, match="the record must be a mapping") as info:
        from_records([series])
    assert info.value.index == 0
    assert len(str(info.value)) < 150
    assert "..." in str(info.value)
    with pytest.raises(FastaWriteError) as info:
        from_records(["x" * 5000])
    assert len(str(info.value)) < 150


def test_from_records_list_of_dicts_unchanged() -> None:
    [entry] = from_records([{"identifier": "a", "sequence": "MK"}])
    assert (entry.identifier, entry.sequence) == ("a", "MK")
