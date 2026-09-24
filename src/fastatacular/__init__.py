"""fastatacular — A pure-Python FASTA parsing and writing library."""

from fastatacular._models import SequenceEntry
from fastatacular._parser import FastaReader, read_fasta
from fastatacular._writer import write_fasta
from fastatacular.errors import FastaError, FastaParseError, FastaWriteError

__version__ = "1.0.0"

__all__ = [
    "FastaError",
    "FastaParseError",
    "FastaReader",
    "FastaWriteError",
    "SequenceEntry",
    "read_fasta",
    "write_fasta",
]
