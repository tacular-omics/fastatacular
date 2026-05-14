"""Frozen dataclass models for FASTA file structures."""

from dataclasses import dataclass, field


@dataclass(frozen=True, slots=True)
class SequenceEntry:
    """A single sequence entry in a FASTA file.

    Always populated:
        identifier:  the token immediately after ``>`` (before any whitespace).
        sequence:    the concatenated sequence with whitespace stripped.

    Parsed from common header conventions when available:
        prefix / accession / entry_name:
            From UniProt-style ``db|ACCESSION|ENTRY_NAME`` identifiers, or
            NCBI-style ``db|ID|...`` identifiers.
        description: free text after the identifier, before any ``KEY=value``.
        pname:       protein name (the description text minus UniProt keys).
        gname:       gene name (``GN=``).
        os_name:     organism name (``OS=``).
        ncbi_tax_id: NCBI taxonomy ID (``OX=``).
        pe:          protein existence level (``PE=``).
        sv:          sequence version (``SV=``).
        extra:       any other ``KEY=value`` pairs found in the header.
        raw_header:  the original header line text (without the leading ``>``).
    """

    identifier: str
    sequence: str
    prefix: str | None = None
    accession: str | None = None
    entry_name: str | None = None
    description: str | None = None
    pname: str | None = None
    gname: str | None = None
    os_name: str | None = None
    ncbi_tax_id: int | None = None
    pe: int | None = None
    sv: int | None = None
    extra: dict[str, str] = field(default_factory=dict)
    raw_header: str = ""
