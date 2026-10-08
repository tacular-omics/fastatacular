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
            NCBI-style ``db|ID|...`` identifiers. For legacy NCBI
            ``gi|NUMBER|db|ACCESSION|NAME`` ids the accession is ``ACCESSION``,
            the entry name ``NAME`` (when present) and the prefix ``gi``; the
            GenInfo number stays in ``identifier``. PDB chains
            (``pdb|1MBA|A``, ``gi|229552|pdb|1MBA|A``) give the accession
            ``1MBA_A`` and no entry name. NCBI ``pir||ENTRY`` and
            ``prf||NAME`` (also after ``gi|N|``) give the third field.
        description: free text after the identifier, before any ``KEY=value``.
        pname:       protein name (the description text minus UniProt keys).
        gname:       gene name (``GN=``).
        os_name:     organism name (``OS=``).
        ncbi_tax_id: NCBI taxonomy ID (``OX=``).
        pe:          protein existence level (``PE=``).
        sv:          sequence version (``SV=``).
        extra:       any other ``KEY=value`` pairs found in the header.
        raw_header:  the original header line text (without the leading ``>``).
                     The writer emits it verbatim only while it still parses to
                     the fields above; after an edit the header is rebuilt.

    ``SequenceEntry`` is frozen and compares by value, but it is not hashable:
    ``extra`` is a ``dict``. ``hash(entry)`` raises ``TypeError``, so entries
    cannot be set members or dict keys; key on ``identifier`` instead.
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

    # Explicitly unhashable (``extra`` is a dict); dataclass keeps an explicit ``None``.
    __hash__ = None  # type: ignore[assignment]

    def to_record(self) -> dict[str, str | int | None]:
        """Return this entry as a flat ``dict`` (the keys of :func:`fastatacular.to_records`)."""
        from fastatacular._records import entry_to_record

        return entry_to_record(self)
