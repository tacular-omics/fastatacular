"""fastatacular — A pure-Python FASTA parsing and writing library."""

from fastatacular._index import FastaIndex
from fastatacular._models import SequenceEntry
from fastatacular._parser import FastaReader, read_fasta
from fastatacular._records import RECORD_KEYS, to_records
from fastatacular._writer import write_fasta
from fastatacular.decoys import (
    MarkovModel,
    is_decoy,
    load_markov_model,
    make_decoy_sequence,
    make_decoys,
    train_markov_model,
    write_decoy_fasta,
)
from fastatacular.errors import DecoyError, FastaError, FastaKeyError, FastaParseError, FastaWriteError

__version__ = "1.0.0"

__all__ = [
    "DecoyError",
    "FastaError",
    "FastaIndex",
    "FastaKeyError",
    "FastaParseError",
    "FastaReader",
    "FastaWriteError",
    "MarkovModel",
    "RECORD_KEYS",
    "SequenceEntry",
    "is_decoy",
    "load_markov_model",
    "make_decoy_sequence",
    "make_decoys",
    "read_fasta",
    "to_records",
    "train_markov_model",
    "write_decoy_fasta",
    "write_fasta",
]
