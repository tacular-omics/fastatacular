"""fastatacular — A pure-Python FASTA parsing and writing library."""

from fastatacular._models import SequenceEntry
from fastatacular._parser import FastaReader, read_fasta
from fastatacular._writer import write_fasta
from fastatacular.errors import FastaParseError, FastaWriteError

__version__ = "0.1.1"

__all__ = [
    "FastaParseError",
    "FastaReader",
    "FastaWriteError",
    "SequenceEntry",
    "read_fasta",
    "write_fasta",
]
