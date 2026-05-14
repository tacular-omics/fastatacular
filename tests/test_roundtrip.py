"""Round-trip tests: parse then re-emit."""

from __future__ import annotations

import io

from fastatacular import read_fasta, write_fasta


def test_roundtrip_preserves_sequences_and_ids():
    src = (
        ">sp|P12345|EX_HUMAN Example protein OS=Homo sapiens OX=9606 GN=EXMP PE=1 SV=2\n"
        "MKTIIALSYIFCLVFA\n"
        "ACDEFGHIKLMNPQRS\n"
        ">tr|Q99999|UNK Second entry OS=Mus musculus OX=10090\n"
        "MAGICSEQ\n"
    )
    entries = read_fasta(io.StringIO(src))
    buf = io.StringIO()
    write_fasta(entries, buf)
    rebuilt = read_fasta(io.StringIO(buf.getvalue()))

    assert len(rebuilt) == len(entries) == 2
    for original, after in zip(entries, rebuilt, strict=True):
        assert original.identifier == after.identifier
        assert original.sequence == after.sequence
        assert original.accession == after.accession
        assert original.gname == after.gname
        assert original.ncbi_tax_id == after.ncbi_tax_id
        assert original.pe == after.pe
        assert original.sv == after.sv


def test_roundtrip_via_file(tmp_path):
    src = ">id desc\nACDEFG\n>other thing\nHHHH\n"
    entries = read_fasta(io.StringIO(src))

    p = tmp_path / "out.fasta"
    write_fasta(entries, p)

    rebuilt = read_fasta(p)
    assert [e.identifier for e in rebuilt] == ["id", "other"]
    assert [e.sequence for e in rebuilt] == ["ACDEFG", "HHHH"]
