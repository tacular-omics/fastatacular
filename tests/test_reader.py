"""Tests for the FASTA reader."""

from __future__ import annotations

import io

import pytest

from fastatacular import FastaParseError, FastaReader, read_fasta


def _read_str(text: str):
    return read_fasta(io.StringIO(text))


def test_plain_single_entry():
    entries = _read_str(">myid simple description\nACDEFG\nHIJKLM\n")
    assert len(entries) == 1
    e = entries[0]
    assert e.identifier == "myid"
    assert e.sequence == "ACDEFGHIJKLM"
    assert e.description == "simple description"
    assert e.pname == "simple description"
    assert e.prefix is None
    assert e.accession is None
    assert e.gname is None
    assert e.extra == {}


def test_multiple_entries():
    text = ">a one\nAAA\n>b two\nCCC\nGGG\n>c three\nTTT\n"
    entries = _read_str(text)
    assert [e.identifier for e in entries] == ["a", "b", "c"]
    assert [e.sequence for e in entries] == ["AAA", "CCCGGG", "TTT"]


def test_uniprot_style_header():
    header = ">sp|P12345|EX_HUMAN Example protein OS=Homo sapiens OX=9606 GN=EXMP PE=1 SV=2\nMKTIIALSYIFCLVFA\n"
    [e] = _read_str(header)
    assert e.prefix == "sp"
    assert e.accession == "P12345"
    assert e.entry_name == "EX_HUMAN"
    assert e.pname == "Example protein"
    assert e.os_name == "Homo sapiens"
    assert e.ncbi_tax_id == 9606
    assert e.gname == "EXMP"
    assert e.pe == 1
    assert e.sv == 2
    assert e.sequence == "MKTIIALSYIFCLVFA"


def test_ncbi_style_pipe_id():
    [e] = _read_str(">gi|12345|ref|NP_000001.1| some description\nACDEFG\n")
    assert e.identifier == "gi|12345|ref|NP_000001.1|"
    assert e.prefix == "gi"
    assert e.accession == "NP_000001.1"
    assert e.entry_name is None
    assert e.description == "some description"


@pytest.mark.parametrize(
    ("identifier", "prefix", "accession", "entry_name"),
    [
        ("gi|4557757|ref|NP_000240.1|", "gi", "NP_000240.1", None),
        ("gi|4557757|ref|NP_000240.1", "gi", "NP_000240.1", None),
        ("gi|2914639|gb|AAC04434.1|", "gi", "AAC04434.1", None),
        ("gi|1234|emb|CAA12345.1|", "gi", "CAA12345.1", None),
        ("gi|1234|pdb|1ABC|A", "gi", "1ABC_A", None),
        ("gi|1234|pdb|1ABC|", "gi", "1ABC", None),
        ("gi|1234|pdb|1ABC", "gi", "1ABC", None),
        ("rev_gi|4557757|ref|NP_000240.1|", "rev_gi", "NP_000240.1", None),
        ("CONTAM_gi|4557757|ref|NP_000240.1|", "CONTAM_gi", "NP_000240.1", None),
        ("DECOY-0-gi|4557757|ref|NP_000240.1|", "DECOY-0-gi", "NP_000240.1", None),
        # Tags stack: make_decoys on a tagged contaminant gives DECOY_CON_gi|...
        ("DECOY_CON_gi|4557757|ref|NP_000240.1|", "DECOY_CON_gi", "NP_000240.1", None),
        ("rev-1-DECOY_gi|4557757|ref|NP_000240.1|", "rev-1-DECOY_gi", "NP_000240.1", None),
        # Only ``gi`` (optionally decoy/contaminant-tagged) is a GenInfo prefix.
        ("fungi|5|ab|Q|", "fungi", "5", None),
        ("xgi|5|ab|Q|", "xgi", "5", None),
        ("DECOY_gi|4557757|ref|NP_000240.1|", "DECOY_gi", "NP_000240.1", None),
        # The fifth field is the entry name (Swiss-Prot) or locus, as in the bare form.
        ("gi|1|sp|P12345|NAME_HUMAN", "gi", "P12345", "NAME_HUMAN"),
        ("DECOY_gi|1|sp|P12345|NAME_HUMAN", "DECOY_gi", "P12345", "NAME_HUMAN"),
        ("gi|1|sp|P12345|NAME_HUMAN|", "gi", "P12345", "NAME_HUMAN"),
        ("gi|2914639|gb|AAC04434.1|LOCUS1", "gi", "AAC04434.1", "LOCUS1"),
        # Bare NCBI pdb|ENTRY|CHAIN: accession ENTRY_CHAIN, no entry name, like gi|N|pdb|...
        ("pdb|1MBA|A", "pdb", "1MBA_A", None),
        ("pdb|1MBA|", "pdb", "1MBA", None),
        ("pdb|1MBA", "pdb", "1MBA", None),
        ("rev_pdb|1MBA|A", "rev_pdb", "1MBA_A", None),
        ("DECOY_CON_pdb|1MBA|B", "DECOY_CON_pdb", "1MBA_B", None),
        ("CONTAM_pdb|1MBA|A", "CONTAM_pdb", "1MBA_A", None),
        # NCBI pir||ENTRY and prf||NAME: empty middle field, the accession is the third.
        ("pir||S71500", "pir", "S71500", None),
        ("prf||0601246A", "prf", "0601246A", None),
        ("pir||S71500|", "pir", "S71500", None),
        ("DECOY_pir||S71500", "DECOY_pir", "S71500", None),
        ("rev_prf||0601246A", "rev_prf", "0601246A", None),
        ("gi|7428543|pir||S71500", "gi", "S71500", None),
        # Custom decoy prefixes ending in a non-alphanumeric character work like the built-in tags.
        ("XXX_gi|1|ref|NP_1.1|", "XXX_gi", "NP_1.1", None),
        ("XXX_gi|1|sp|P1|A_B", "XXX_gi", "P1", "A_B"),
        ("XXX_gi|229552|pdb|1MBA|A", "XXX_gi", "1MBA_A", None),
        ("XXX_pir||S71500", "XXX_pir", "S71500", None),
        ("REV__prf||0601246A", "REV__prf", "0601246A", None),
        ("REV__pdb|1MBA|A", "REV__pdb", "1MBA_A", None),
        ("decoy.gi|7428543|pir||S71500", "decoy.gi", "S71500", None),
        ("gi|7428543|prf||0601246A", "gi", "0601246A", None),
        ("DECOY_gi|7428543|pir||S71500", "DECOY_gi", "S71500", None),
        # An empty middle field in other dbs is not a known form: first two fields as before.
        ("gi|5|ref||X", "gi", "5", None),
        ("pir||", None, None, None),
        ("ref|NP_000240.1|", "ref", "NP_000240.1", None),
        ("sp|P12345|NAME_HUMAN", "sp", "P12345", "NAME_HUMAN"),
        ("tr|A0A024R161|A0A024R161_HUMAN", "tr", "A0A024R161", "A0A024R161_HUMAN"),
    ],
)
def test_ncbi_gi_accession(identifier: str, prefix: str, accession: str, entry_name: str | None) -> None:
    [e] = _read_str(f">{identifier} desc\nACDEFG\n")
    assert (e.identifier, e.prefix, e.accession, e.entry_name) == (identifier, prefix, accession, entry_name)


def test_rcsb_entity_header():
    """The header from issue #36: the identifier is the entity id, not ``1MBA_1|Chain``."""
    [e] = _read_str(">1MBA_1|Chain A|MYOGLOBIN|Aplysia limacina\nSLSAAEADLA\n")
    assert e.identifier == "1MBA_1"
    assert e.accession == "1MBA_1"
    assert e.prefix is None
    assert e.entry_name is None
    assert e.description == "Chain A|MYOGLOBIN|Aplysia limacina"
    assert e.pname == "Chain A|MYOGLOBIN|Aplysia limacina"
    assert e.raw_header == "1MBA_1|Chain A|MYOGLOBIN|Aplysia limacina"
    assert e.sequence == "SLSAAEADLA"


@pytest.mark.parametrize(
    ("header", "identifier", "description"),
    [
        # One chain, several chains, and an author-assigned chain label.
        ("1MBA_1|Chain A|MYOGLOBIN|Aplysia limacina", "1MBA_1", "Chain A|MYOGLOBIN|Aplysia limacina"),
        ("6M0J_1|Chains A, B|Spike|SARS-CoV-2", "6M0J_1", "Chains A, B|Spike|SARS-CoV-2"),
        ("7ABC_2|Chain A[auth B]|Nanobody|Lama glama", "7ABC_2", "Chain A[auth B]|Nanobody|Lama glama"),
        # Only the field right after the first one counts, and it must be ``Chain``/``Chains``.
        ("1MBA_1|Chaining|MYOGLOBIN", "1MBA_1|Chaining|MYOGLOBIN", None),
        ("1MBA_1|A|Chain A|MYOGLOBIN", "1MBA_1|A|Chain", "A|MYOGLOBIN"),
    ],
)
def test_rcsb_identifier_boundary(header: str, identifier: str, description: str | None) -> None:
    [e] = _read_str(f">{header}\nACDEFG\n")
    assert (e.identifier, e.description) == (identifier, description)


@pytest.mark.parametrize(
    ("header", "identifier", "accession"),
    [
        # The older RCSB/NCBI form keeps ENTRY_CHAIN: its second field is not ``Chain``.
        ("pdb|1MBA|A Chain A,  Myoglobin", "pdb|1MBA|A", "1MBA_A"),
        ("gi|229552|pdb|1MBA|A Chain A", "gi|229552|pdb|1MBA|A", "1MBA_A"),
        # A description that happens to start with "Chain" is still a description.
        ("sp|P02185|MYG_PHYCD Chain A of myoglobin", "sp|P02185|MYG_PHYCD", "P02185"),
    ],
)
def test_pdb_pipe_form_not_affected(header: str, identifier: str, accession: str) -> None:
    [e] = _read_str(f">{header}\nACDEFG\n")
    assert (e.identifier, e.accession) == (identifier, accession)


def test_rcsb_entity_header_from_a_path(tmp_path):
    """The chunked path reader agrees with the line-by-line handle reader (two readers)."""
    text = ">1MBA_1|Chain A|MYOGLOBIN|Aplysia limacina\nSLSAAEADLA\n>6M0J_1|Chains A, B|Spike|SARS-CoV-2\nACDE\n"
    p = tmp_path / "rcsb.fasta"
    p.write_text(text)
    assert read_fasta(p) == _read_str(text)
    assert [e.accession for e in read_fasta(p)] == ["1MBA_1", "6M0J_1"]


# --------------------------------------------------- bare accessions (issue #37)


def test_ncbi_defline_fields():
    """The NCBI RefSeq/GenPept defline from issue #37: bare accession, organism in brackets."""
    [e] = _read_str(">NP_000549.1 hemoglobin subunit alpha [Homo sapiens]\nMVLSPADKTN\n")
    assert e.identifier == "NP_000549.1"
    assert e.accession == "NP_000549.1"
    assert e.prefix is None
    assert e.entry_name is None
    assert e.os_name == "Homo sapiens"
    assert e.ncbi_tax_id is None  # an NCBI defline carries no taxonomy id
    # ``pname`` is the description before its first ``KEY=`` token, brackets and all.
    assert e.pname == "hemoglobin subunit alpha [Homo sapiens]"
    assert e.description == "hemoglobin subunit alpha [Homo sapiens]"
    assert e.sequence == "MVLSPADKTN"


@pytest.mark.parametrize(
    ("header", "accession", "os_name"),
    [
        # Every RefSeq protein prefix, with and without the version.
        ("NP_000549.1 x [Homo sapiens]", "NP_000549.1", "Homo sapiens"),
        ("XP_016880112.1 x [Homo sapiens]", "XP_016880112.1", "Homo sapiens"),
        (
            "YP_009724390.1 x [Severe acute respiratory syndrome coronavirus 2]",
            "YP_009724390.1",
            "Severe acute respiratory syndrome coronavirus 2",
        ),
        ("WP_000184067.1 MULTISPECIES: MbtH family protein [Bacillus]", "WP_000184067.1", "Bacillus"),
        ("AP_000123 x [Homo sapiens]", "AP_000123", "Homo sapiens"),
        # GenPept: three letters, five or seven digits, and the version is required.
        ("AAA59174.1 insulin receptor precursor [Homo sapiens]", "AAA59174.1", "Homo sapiens"),
        ("ABC1234567.2 x [Homo sapiens]", "ABC1234567.2", "Homo sapiens"),
        ("AAA59174 x [Homo sapiens]", None, None),
        # A decoy or contaminant tag keeps the target's accession, as for ``gi|`` ids.
        ("DECOY_NP_000549.1 x [Homo sapiens]", "NP_000549.1", "Homo sapiens"),
        ("rev-0-AAA59174.1 x [Homo sapiens]", "AAA59174.1", "Homo sapiens"),
        # Any other first word is an identifier only, so no organism is read either.
        ("sequence_1 x [Homo sapiens]", None, None),
        ("NC_000001.11 x [Homo sapiens]", None, None),
        # ``[key=value]`` is an NCBI qualifier, not an organism, and only a trailing
        # bracket is one.
        ("NP_000549.1 x [organism=Homo sapiens]", "NP_000549.1", None),
        ("NP_000549.1 x [Homo sapiens] isoform a", "NP_000549.1", None),
        # An explicit ``OS=`` wins over the brackets.
        ("NP_000549.1 x [Homo sapiens] OS=Mus musculus", "NP_000549.1", "Mus musculus"),
        # A pipe-delimited NCBI identifier keeps its own rules: the brackets stay text.
        ("gi|4504345|ref|NP_000549.1| x [Homo sapiens]", "NP_000549.1", None),
    ],
)
def test_ncbi_bare_accession_boundary(header: str, accession: str | None, os_name: str | None) -> None:
    [e] = _read_str(f">{header}\nACDEFG\n")
    assert (e.accession, e.os_name) == (accession, os_name)
    assert e.prefix is None or e.prefix.endswith("gi")


def test_ensembl_pep_fields():
    """The Ensembl ``pep`` header from issue #37: ``gene_symbol:`` is the gene name."""
    header = (
        ">ENSP00000322421.5 pep chromosome:GRCh38:16:176680:177522:1 gene:ENSG00000206172.8 "
        "transcript:ENST00000320868.9 gene_biotype:protein_coding gene_symbol:HBA1"
    )
    [e] = _read_str(f"{header}\nMVLSPADKTN\n")
    assert e.accession == "ENSP00000322421.5"
    assert e.gname == "HBA1"
    # Ensembl writes ``key:value``, so the whole description is the name text, and the
    # gene and transcript ids have no field of their own.
    assert e.pname == e.description == header[1:].split(maxsplit=1)[1]
    assert (e.os_name, e.ncbi_tax_id, e.extra) == (None, None, {})


@pytest.mark.parametrize(
    ("header", "accession", "gname"),
    [
        ("ENSP00000322421.5 pep gene_symbol:HBA1", "ENSP00000322421.5", "HBA1"),
        ("ENSP00000322421 pep gene_symbol:HBA1", "ENSP00000322421", "HBA1"),
        ("ENSMUSP00000020316.9 pep gene_symbol:Hba-a1", "ENSMUSP00000020316.9", "Hba-a1"),
        ("REV_ENSP00000322421.5 pep gene_symbol:HBA1", "ENSP00000322421.5", "HBA1"),
        # No gene symbol in the header, and a gene id is not a gene name.
        ("ENSP00000322421.5 pep gene:ENSG00000206172.8", "ENSP00000322421.5", None),
        # An explicit ``GN=`` wins; a word ending in the token is not the token.
        ("ENSP00000322421.5 pep gene_symbol:HBA1 GN=HBA2", "ENSP00000322421.5", "HBA2"),
        ("ENSP00000322421.5 pep xgene_symbol:HBA1", "ENSP00000322421.5", None),
        # Transcript and gene ids are not protein ids, and Ensembl Genomes gene-model
        # names have no fixed shape.
        ("ENST00000320868.9 cdna gene_symbol:HBA1", None, None),
        ("AT1G01010.1 pep gene_symbol:NAC001", None, None),
    ],
)
def test_ensembl_accession_boundary(header: str, accession: str | None, gname: str | None) -> None:
    [e] = _read_str(f">{header}\nACDEFG\n")
    assert (e.accession, e.gname) == (accession, gname)


def test_uniref_fields():
    """The UniRef header from issue #37: ``Tax=``/``TaxID=`` are the organism fields."""
    header = ">UniRef90_P69905 Hemoglobin subunit alpha n=12 Tax=Primates TaxID=9443 RepID=HBA_HUMAN"
    [e] = _read_str(f"{header}\nMVLSPADKTN\n")
    assert e.accession == "UniRef90_P69905"
    assert e.os_name == "Primates"
    assert e.ncbi_tax_id == 9443
    assert e.pname == "Hemoglobin subunit alpha"
    # The cluster size and the representative member have no field of their own.
    assert e.extra == {"n": "12", "RepID": "HBA_HUMAN"}
    assert e.gname is None


@pytest.mark.parametrize(
    ("header", "accession", "os_name", "ncbi_tax_id", "extra"),
    [
        ("UniRef50_P69905 x Tax=Primates TaxID=9443", "UniRef50_P69905", "Primates", 9443, {}),
        ("UniRef100_UPI0000000005 x Tax=Primates", "UniRef100_UPI0000000005", "Primates", None, {}),
        ("DECOY_UniRef90_P69905 x Tax=Primates", "UniRef90_P69905", "Primates", None, {}),
        # An explicit UniProt key wins, as for a repeated key, and the key it leaves
        # unused stays in ``extra``.
        (
            "UniRef90_P69905 x OS=Homo sapiens OX=9606 Tax=Primates TaxID=9443",
            "UniRef90_P69905",
            "Homo sapiens",
            9606,
            {"Tax": "Primates", "TaxID": "9443"},
        ),
        # A non-integer taxonomy id stays in ``extra``, as a non-integer ``OX=`` does.
        ("UniRef90_P69905 x Tax=Primates TaxID=9443.0", "UniRef90_P69905", "Primates", None, {"TaxID": "9443.0"}),
        # Only the three cluster identities exist, and ``Tax=`` is UniRef's own key.
        ("UniRef70_P69905 x Tax=Primates", None, None, None, {"Tax": "Primates"}),
        ("P69905 x Tax=Primates", None, None, None, {"Tax": "Primates"}),
    ],
)
def test_uniref_boundary(
    header: str, accession: str | None, os_name: str | None, ncbi_tax_id: int | None, extra: dict[str, str]
) -> None:
    [e] = _read_str(f">{header}\nACDEFG\n")
    assert (e.accession, e.os_name, e.ncbi_tax_id, e.extra) == (accession, os_name, ncbi_tax_id, extra)


def test_bare_accession_headers_from_a_path(tmp_path):
    """The chunked path reader agrees with the line-by-line handle reader (two readers)."""
    text = (
        ">NP_000549.1 hemoglobin subunit alpha [Homo sapiens]\nMVLS\n"
        ">ENSP00000322421.5 pep gene:ENSG00000206172.8 gene_symbol:HBA1\nMVLS\n"
        ">UniRef90_P69905 Hemoglobin subunit alpha n=12 Tax=Primates TaxID=9443 RepID=HBA_HUMAN\nMVLS\n"
    )
    p = tmp_path / "bare.fasta"
    p.write_text(text)
    assert read_fasta(p) == _read_str(text)
    assert [(e.accession, e.os_name, e.gname) for e in read_fasta(p)] == [
        ("NP_000549.1", "Homo sapiens", None),
        ("ENSP00000322421.5", None, "HBA1"),
        ("UniRef90_P69905", "Primates", None),
    ]


def test_extra_keys_captured():
    [e] = _read_str(">id name part FOO=bar BAZ=qux extra trailing\nA\n")
    assert e.extra == {"FOO": "bar", "BAZ": "qux extra trailing"}
    assert e.pname == "name part"


def test_skips_blank_and_comment_lines():
    text = "\n; this is a comment\n>id description\nACDE\n\nFGHI\n"
    [e] = _read_str(text)
    assert e.sequence == "ACDEFGHI"


def test_sequence_before_header_raises():
    with pytest.raises(FastaParseError) as exc:
        _read_str("ACDEFG\n>id description\nACDE\n")
    assert exc.value.line == 1


def test_empty_sequence_parses():
    entries = _read_str(">id description\n>next other\nACDE\n>last\n")
    assert [(e.identifier, e.sequence) for e in entries] == [("id", ""), ("next", "ACDE"), ("last", "")]


def test_empty_header_raises():
    with pytest.raises(FastaParseError):
        _read_str("> \nACDE\n")


def test_lazy_reader_context_manager(tmp_path):
    p = tmp_path / "x.fasta"
    p.write_text(">a x\nAAA\n>b y\nCCC\n")
    with FastaReader(p) as reader:
        ids = [e.identifier for e in reader]
    assert ids == ["a", "b"]


def test_lazy_reader_requires_context_manager():
    reader = FastaReader(io.StringIO(">a x\nAAA\n"))
    with pytest.raises(RuntimeError):
        next(iter(reader))


def test_header_with_no_description():
    [e] = _read_str(">justanid\nACDE\n")
    assert e.identifier == "justanid"
    assert e.description is None
    assert e.pname is None
    assert e.sequence == "ACDE"


def test_invalid_ox_falls_back_to_extra():
    [e] = _read_str(">id name OX=notanumber\nACDE\n")
    assert e.ncbi_tax_id is None
    assert e.extra == {"OX": "notanumber"}


def test_tab_separates_identifier_from_description():
    (entry,) = _read_str(">sp|P12345|EX_HUMAN\tExample protein OS=Homo sapiens\nMK\n")
    assert entry.identifier == "sp|P12345|EX_HUMAN"
    assert entry.entry_name == "EX_HUMAN"
    assert entry.pname == "Example protein"
    assert entry.os_name == "Homo sapiens"
    assert entry.raw_header == "sp|P12345|EX_HUMAN\tExample protein OS=Homo sapiens"


def test_internal_whitespace_in_sequence_lines_is_removed():
    (entry,) = _read_str(">x\nA C\tD\nE  F\n")
    assert entry.sequence == "ACDEF"


def test_utf8_bom_file_is_accepted(tmp_path):
    path = tmp_path / "bom.fasta"
    path.write_bytes(b"\xef\xbb\xbf>x desc\nACDE\n")
    (e,) = read_fasta(path)
    assert e.identifier == "x"
    assert e.raw_header == "x desc"
    with FastaReader(path) as reader:
        assert [r.identifier for r in reader] == ["x"]


def test_utf8_bom_in_text_handle_is_accepted():
    (e,) = _read_str("﻿>x desc\nACDE\n")
    assert e.identifier == "x"


def test_whitespace_only_lines_before_first_header_are_skipped():
    entries = _read_str("   \n\t\r\n>x\nAC\n  \nDE\n")
    assert [(e.identifier, e.sequence) for e in entries] == [("x", "ACDE")]


@pytest.mark.parametrize(
    ("header", "pname", "extra", "os_name"),
    [
        # ``KEY=`` inside a word is not a key: a key starts the description or follows whitespace.
        (">x Protein(EC=2.7.1) OS=Homo sapiens", "Protein(EC=2.7.1)", {}, "Homo sapiens"),
        (
            ">NP_000001.1 insulin [organism=Homo sapiens] [gene=INS]",
            "insulin [organism=Homo sapiens] [gene=INS]",
            {},
            None,
        ),
        (">x a-b=c d", "a-b=c d", {}, None),
        # A key at the start of the description, or after whitespace, is still a key.
        (">x pH=7 sensor", None, {"pH": "7 sensor"}, None),
        (">x name FOO=bar", "name", {"FOO": "bar"}, None),
    ],
)
def test_keys_must_start_a_word(header, pname, extra, os_name):
    [e] = _read_str(f"{header}\nA\n")
    assert e.pname == pname
    assert e.extra == extra
    assert e.os_name == os_name


@pytest.mark.parametrize(
    ("identifier", "prefix", "accession", "entry_name"),
    [
        ("Reverse_sp|P12345|X_HUMAN", "Reverse_sp", "P12345", "X_HUMAN"),
        ("rev_sp|P12345|X_HUMAN", "rev_sp", "P12345", "X_HUMAN"),
        ("DECOY-0-sp|P12345|X_HUMAN", "DECOY-0-sp", "P12345", "X_HUMAN"),
        ("contam_sp|P12345|X_HUMAN", "contam_sp", "P12345", "X_HUMAN"),
        ("Reverse_tr|A0A024R161", "Reverse_tr", "A0A024R161", None),
        ("sp|P12345|X_HUMAN", "sp", "P12345", "X_HUMAN"),
    ],
)
def test_decoy_and_contaminant_prefixes_keep_the_accession(identifier, prefix, accession, entry_name):
    [e] = _read_str(f">{identifier} Name OS=Homo sapiens\nA\n")
    assert (e.prefix, e.accession, e.entry_name) == (prefix, accession, entry_name)


# Leading whitespace is ignored on every line, so an indented '>' starts a header
# (as in FastaFrames) in every parser path: text handles, paths, chunked reads, FastaIndex.
_INDENTED = {
    "mid_file": ">a\nMK\n  >b\nPE\n",
    "first_header": "  >a\nMK\n>b\nPE\n",
    "tabs": "\t>a\n\tMK\t\n \t>b\n PE \n",
    "no_trailing_newline": ">a\nMK\n  >b\nPE",
    "comment": ">a\nMK\n  ;note\n\t>b\n  PE\n",
}


@pytest.mark.parametrize("text", list(_INDENTED.values()), ids=list(_INDENTED))
def test_indented_header_text_handle(text):
    assert [(e.identifier, e.sequence) for e in _read_str(text)] == [("a", "MK"), ("b", "PE")]


@pytest.mark.parametrize("text", list(_INDENTED.values()), ids=list(_INDENTED))
def test_indented_header_path_and_reader(tmp_path, text):
    path = tmp_path / "x.fasta"
    path.write_text(text)
    expected = [("a", "MK"), ("b", "PE")]
    assert [(e.identifier, e.sequence) for e in read_fasta(path)] == expected
    with FastaReader(path) as reader:
        assert [(e.identifier, e.sequence) for e in reader] == expected
    with path.open("rb") as fh:
        assert [(e.identifier, e.sequence) for e in read_fasta(fh)] == expected


@pytest.mark.parametrize("chunk", range(1, 12))
@pytest.mark.parametrize("text", list(_INDENTED.values()), ids=list(_INDENTED))
def test_indented_header_across_chunk_boundaries(tmp_path, monkeypatch, text, chunk):
    from fastatacular import _parser

    monkeypatch.setattr(_parser, "_CHUNK_CHARS", chunk)
    path = tmp_path / "x.fasta"
    path.write_text(text)
    assert [(e.identifier, e.sequence) for e in read_fasta(path)] == [("a", "MK"), ("b", "PE")]


def test_indented_header_line_number_in_errors():
    with pytest.raises(FastaParseError) as exc:
        _read_str(">a\nMK\n  >\nPE\n")
    assert exc.value.line == 3
