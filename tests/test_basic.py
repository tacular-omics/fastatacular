"""Smoke test: verify public API is importable."""

from fastatacular import (
    FastaParseError,
    FastaReader,
    FastaWriteError,
    SequenceEntry,
    read_fasta,
    write_fasta,
)


def test_public_api_importable():
    assert FastaParseError is not None
    assert FastaReader is not None
    assert FastaWriteError is not None
    assert SequenceEntry is not None
    assert read_fasta is not None
    assert write_fasta is not None
