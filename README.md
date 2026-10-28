# fastatacular

[![PyPI](https://img.shields.io/pypi/v/fastatacular)](https://pypi.org/project/fastatacular/)
[![Python Package](https://github.com/tacular-omics/fastatacular/actions/workflows/ci.yml/badge.svg)](https://github.com/tacular-omics/fastatacular/actions/workflows/ci.yml)
[![License](https://img.shields.io/github/license/tacular-omics/fastatacular)](https://github.com/tacular-omics/fastatacular/blob/main/LICENSE)
[![Python](https://img.shields.io/pypi/pyversions/fastatacular)](https://pypi.org/project/fastatacular/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22926358.svg)](https://doi.org/10.5281/zenodo.22926358)

fastatacular is a pure-Python library with no dependencies for reading and writing
[FASTA](https://en.wikipedia.org/wiki/FASTA_format) sequence files. Each record becomes a
frozen `SequenceEntry` whose header is parsed into typed fields: the identifier, the
database prefix, accession and entry name of pipe-delimited identifiers
(`sp|P12345|EX_HUMAN`, `gi|12345|ref|NP_000001.1|`), the UniProt description keys
`OS=`, `OX=`, `GN=`, `PE=` and `SV=`, and the accession and description fields of NCBI,
Ensembl and UniRef headers (`NP_000549.1`, `ENSP00000322421.5`, `UniRef90_P69905`). The library also reads compressed files, writes
headers back unchanged, builds target-decoy databases, gives random access by identifier
or accession through a samtools-compatible `.fai` index. It is intended for proteomics and genomics pipelines
and requires Python 3.12 or later.

## Documentation

- This README: installation, quick start and the reference for reading, writing,
  header parsing and indexing.
- [docs/decoys.md](https://github.com/tacular-omics/fastatacular/blob/main/docs/decoys.md):
  decoy methods, options, the Markov models and decoy quality on the human proteome.
- [llms-full.txt](https://github.com/tacular-omics/fastatacular/blob/main/llms-full.txt):
  the complete usage guide, with every public name, signature and runnable example.
- [CHANGELOG.md](https://github.com/tacular-omics/fastatacular/blob/main/CHANGELOG.md):
  release history.
- [FASTA Toolkit](https://tacular-omics.github.io/fastatacular/): a web app built on
  fastatacular that makes decoys, checks decoy quality, cleans up, merges and converts
  FASTA files in the browser (files are not uploaded). To browse entries, use the
  [FASTA Viewer](https://pgarrett-scripps.github.io/fastaviewer/).

## Install

```bash
pip install fastatacular
```

## Quick start

Write a FASTA file, then read it back:

```python
from fastatacular import SequenceEntry, read_fasta, write_fasta

write_fasta(
    [
        SequenceEntry(
            identifier="sp|P12345|EX_HUMAN",
            sequence="MKTIIALSYIFCLVFA",
            pname="Example protein",
            os_name="Homo sapiens",
            ncbi_tax_id=9606,
            gname="EXMP",
            pe=1,
            sv=2,
        ),
        SequenceEntry(identifier="seq2", sequence="MPEPTIDEK"),
    ],
    "proteins.fasta",
)

entries = read_fasta("proteins.fasta")
for entry in entries:
    print(entry.identifier, entry.accession, len(entry.sequence))
# sp|P12345|EX_HUMAN P12345 16
# seq2 None 9
```

`read_fasta` returns a list. `FastaReader` yields one entry at a time without loading
the whole file:

```python
from fastatacular import FastaReader

with FastaReader("proteins.fasta") as reader:
    for entry in reader:
        print(entry.identifier, entry.sequence)
```

A `str` argument is always a file path. To parse FASTA text, wrap it in `io.StringIO`:

```python
import io

from fastatacular import read_fasta

entries = read_fasta(io.StringIO(">a\nMKV\n>b\n  >c\nPE\n"))
[(e.identifier, e.sequence) for e in entries]
# [('a', 'MKV'), ('b', ''), ('c', 'PE')]
```

As the example shows, a header followed directly by another header (or by the end of
the file) gives an entry with an empty sequence, and a header line with leading
whitespace is still read as a header. `write_fasta` writes an entry with an empty
sequence as its header line alone, so it reads back identically.

**Compressed files** are read transparently: gzip, bzip2 and xz, detected from the magic
bytes, not the file name. Pipes, FIFOs and `/dev/stdin` work too:

```python
entries = read_fasta("uniprot_sprot.fasta.gz")
```

Every reader and writer takes a pandas-style `compression=` keyword (type alias
`fastatacular.Compression`): `"infer"` (default), `"gzip"`, `"bz2"`, `"xz"` or `None`.
On read, `"infer"` sniffs the magic bytes, an explicit format is forced (input in another
format raises `FastaParseError`) and `None` reads plain text. A binary handle is read
as plain UTF-8 unless you pass an explicit format:

```python
entries = read_fasta("download.bin", compression="gzip")
with open("proteins.fasta.xz", "rb") as fh:
    entries = read_fasta(fh, compression="xz")
```

A PEFF file also reads as plain FASTA: its `#` file header lines are skipped and each
entry keeps its identifier, sequence and raw description. Use
[pefftacular](https://github.com/tacular-omics/pefftacular) to parse the PEFF annotations.

## Data model

Each entry is a `SequenceEntry`:

| Field | Type | Description |
|---|---|---|
| `identifier` | `str` | Header text after `>` up to the first whitespace (e.g. `sp\|P12345\|EX_HUMAN`); for an RCSB PDB entity header (`1MBA_1\|Chain A\|...`) the first pipe field |
| `sequence` | `str` | Sequence lines joined, with all whitespace removed; `""` for a header with no sequence lines |
| `prefix` | `str \| None` | Database prefix (`sp`, `tr`, `gi`, ...) when the identifier is pipe-delimited |
| `accession` | `str \| None` | Second pipe field (e.g. `P12345` in `sp\|P12345\|EX_HUMAN`); the fourth on `gi\|N\|db\|ACCESSION` identifiers; `ENTRY_CHAIN` for PDB chains (`1MBA_A` from `pdb\|1MBA\|A` or `gi\|229552\|pdb\|1MBA\|A`); the third field of NCBI `pir\|\|ENTRY` and `prf\|\|NAME` (`S71500` from `pir\|\|S71500` or `gi\|7428543\|pir\|\|S71500`); the identifier itself when it is an NCBI protein, Ensembl protein or UniRef accession (`NP_000549.1`, `ENSP00000322421.5`, `UniRef90_P69905`) |
| `entry_name` | `str \| None` | Third pipe field on UniProt identifiers (e.g. `EX_HUMAN`), or the fifth on `gi\|N\|sp\|P12345\|EX_HUMAN`; none for PDB chains |
| `description` | `str \| None` | Header text after the identifier, including any `KEY=value` pairs |
| `pname` | `str \| None` | Protein name (the description before its first `KEY=value` pair) |
| `gname` | `str \| None` | Gene name (`GN=`, or Ensembl's `gene_symbol:`) |
| `os_name` | `str \| None` | Organism name (`OS=`, UniRef's `Tax=`, or the trailing `[Homo sapiens]` of an NCBI defline) |
| `ncbi_tax_id` | `int \| None` | NCBI taxonomy ID (`OX=`, or UniRef's `TaxID=`) |
| `pe` | `int \| None` | Protein existence level (`PE=`) |
| `sv` | `int \| None` | Sequence version (`SV=`) |
| `extra` | `dict[str, str]` | Any other `KEY=value` pairs found in the header |
| `raw_header` | `str` | The header line as read, without the leading `>`; `""` on an entry you construct |

## UniProt-style headers

```python
from fastatacular import read_fasta

[entry] = read_fasta("one.fasta")
# >sp|P12345|EX_HUMAN Example protein OS=Homo sapiens OX=9606 GN=EXMP PE=1 SV=2

entry.prefix         # "sp"
entry.accession      # "P12345"
entry.entry_name     # "EX_HUMAN"
entry.pname          # "Example protein"
entry.os_name        # "Homo sapiens"
entry.ncbi_tax_id    # 9606
entry.gname          # "EXMP"
entry.pe             # 1
entry.sv             # 2
```

Other `KEY=value` pairs are kept in `entry.extra`. A header without `KEY=value` pairs
fills `description` and `pname` and leaves `extra` empty. When a key appears more than
once (`OS=Human ... OS=Mouse`), the first value is used; `description` and `raw_header`
still hold the full text.

## NCBI, Ensembl and UniRef headers

Identifiers that are a bare accession are recognized by their shape, and the fields their
database writes in the description are read into the same typed fields:

```python
[ncbi] = read_fasta("refseq.fasta")
# >NP_000549.1 hemoglobin subunit alpha [Homo sapiens]
ncbi.accession       # "NP_000549.1"
ncbi.os_name         # "Homo sapiens" (the defline's trailing brackets)

[ensembl] = read_fasta("ensembl.pep.fasta")
# >ENSP00000322421.5 pep chromosome:GRCh38:16:176680:177522:1 gene:ENSG00000206172.8 \
#  transcript:ENST00000320868.9 gene_biotype:protein_coding gene_symbol:HBA1
ensembl.accession    # "ENSP00000322421.5"
ensembl.gname        # "HBA1" (from gene_symbol:)

[uniref] = read_fasta("uniref90.fasta")
# >UniRef90_P69905 Hemoglobin subunit alpha n=12 Tax=Primates TaxID=9443 RepID=HBA_HUMAN
uniref.accession     # "UniRef90_P69905"
uniref.os_name       # "Primates" (from Tax=)
uniref.ncbi_tax_id   # 9443 (from TaxID=)
uniref.extra         # {"n": "12", "RepID": "HBA_HUMAN"}
```

Only what each format defines is parsed, and an explicit UniProt key always wins:

- NCBI: RefSeq protein (`NP_`, `XP_`, `YP_`, `WP_`, `AP_`, `ZP_`) and versioned GenPept
  (`AAA59174.1`) accessions. The organism is the defline's trailing `[...]` only, so an
  NCBI `[key=value]` qualifier is left as text, and the brackets stay in `pname` as well.
  A defline carries no taxonomy id, so `ncbi_tax_id` stays `None`.
- Ensembl: protein ids of an Ensembl `pep` file (`ENSP...`, `ENSMUSP...`). Only
  `gene_symbol:` is a gene name; `gene:` and `transcript:` are Ensembl ids, which no field
  takes, and Ensembl `pep` headers name no organism.
- UniRef: `UniRef50_`, `UniRef90_` and `UniRef100_` clusters. The cluster size (`n=`) and
  the representative member (`RepID=`) stay in `extra`.
- A decoy or contaminant tag in front of the identifier keeps the target's accession
  (`DECOY_NP_000549.1`), as it does for pipe-delimited identifiers.
- Any other first word is an identifier only: `accession` stays `None` and no description
  field is guessed from it.

## Writing

Construct entries and write them out:

```python
from fastatacular import SequenceEntry, write_fasta

entries = [
    SequenceEntry(
        identifier="sp|P12345|EX_HUMAN",
        sequence="MKTIIALSYIFCLVFA",
        pname="Example protein",
        os_name="Homo sapiens",
        ncbi_tax_id=9606,
        gname="EXMP",
        pe=1,
        sv=2,
    ),
]

write_fasta(entries, "output.fasta")
```

`dest` accepts a path string, a `pathlib.Path`, or a text or binary file object (a
binary one gets UTF-8 bytes; a handle is never closed). With the default
`compression="infer"`, a path ending in `.gz`, `.bz2` or `.xz` (any case) is written gzip,
bzip2 or xz compressed and anything else, including a handle, is plain text.
`compression="gzip"`, `"bz2"` or `"xz"` forces that format whatever the suffix (a handle
must then be binary), and `compression=None` always writes plain text. gzip output has
mtime 0 and no file name in its header, so the same entries give the same bytes.

```python
write_fasta(entries, "output.fa", compression="xz")
with open("output.fasta.gz", "wb") as fh:
    write_fasta(entries, fh, compression="gzip")
```

Sequence lines wrap at 60 characters by default. Override with `line_width=` (pass `0` to disable wrapping):

```python
write_fasta(entries, "output.fasta", line_width=80)
write_fasta(entries, "single-line.fasta", line_width=0)
```

If `raw_header` is set on an entry (as it is on every entry produced by `read_fasta`) and still matches the entry's structured fields, the writer round-trips it verbatim. If you changed a field (for example `dataclasses.replace(entry, gname="XYZ")`), or `raw_header` is empty, the header is rebuilt from the structured fields, so your edit is written.

## Decoy databases

Build target-decoy databases for FDR estimation with `reverse`, `pseudo_reverse`,
`shuffle`, `debruijn` (repeat-preserving, Moosa et al. 2020) or `markov` (order-2
chains shipped for human, mouse, yeast and E. coli, trained on UniProt, CC BY 4.0).
Every method takes `keep_residues` (e.g. `"KR"`), `keep_nterm` and `keep_cterm`, and the
same `seed` always gives the same decoys. Pure Python, no extra dependencies.

```python
from fastatacular import is_decoy, make_decoys, read_fasta, write_decoy_fasta

write_decoy_fasta("human.fasta", "human_td.fasta", method="pseudo_reverse", seed=1)

decoys = make_decoys(read_fasta("human.fasta"), method="markov", model="human", keep_residues="KR", seed=1)
next(decoys).identifier   # "DECOY_sp|..."; is_decoy(entry) checks the prefix
```

See [docs/decoys.md](docs/decoys.md) for every option, the Markov model data and
per-method quality numbers on the human proteome.

## Random access by identifier or accession

`FastaIndex(path)` reads the file once and keeps only each entry's key and byte range;
`index[key]` then reads and parses just that entry.

```python
from fastatacular import FastaIndex

index = FastaIndex("human.fasta")
entry = index["sp|P31946|1433B_HUMAN"]   # SequenceEntry, read from disk on demand
"sp|P31946|1433B_HUMAN" in index, len(index)   # no file access
index.write_fai()                        # human.fasta.fai, samtools-compatible
index = FastaIndex.from_fai("human.fasta")   # later: load the .fai instead of scanning

by_acc = FastaIndex("human.fasta", key="accession")
by_acc["P31946"]
```

- By default keys are full identifiers (the first header word, the `.fai` name, as in
  samtools). This works on target-decoy databases, where `sp|P1|X` and `DECOY_sp|P1|X`
  or `Reverse_sp|P1|X` share an accession.
- `key="accession"` keys entries by accession (`P31946` for `sp|P31946|1433B_HUMAN`, or
  the whole identifier when it has no `|`). It needs unique accessions.
- A repeated key raises `FastaError` naming the first one. Real databases do repeat
  identifiers (IP2 exports, merged databases): `FastaIndex(path, duplicates="first")`
  keeps the first entry for each key and skips the rest, as `samtools faidx` does, and
  logs one warning with the number skipped. `from_fai` takes the same option.
- A missing key raises `FastaKeyError`, which is also a `KeyError`.
- The file must be a regular, uncompressed file. A FIFO, pipe or directory raises
  `FastaError`: save the input to a file first, or read it once with `FastaReader`. gzip (including bgzip), bzip2 and xz raise `FastaError`:
  decompress first (`gunzip -k human.fasta.gz`). bgzip/`.gzi` is not supported.
  `compression="gzip"`/`"bz2"`/`"xz"` raises `FastaError`; `compression=None` skips the
  magic-byte check.
- `.fai` caveats: the name column is the first header word (mapped to the accession on
  load with `key="accession"`). Like samtools, `write_fai()` needs every sequence line of an
  entry but the last to have the same length, and no comment or blank lines inside a
  sequence; otherwise it raises `FastaError` (rewrite the file with `write_fasta` first).
  `from_fai` checks each entry's header against the `.fai` name but does not re-count
  residues, so rebuild the `.fai` whenever the FASTA changes. An entry with no sequence
  gets a row with length 0. samtools and pysam refuse a FASTA file that starts with a
  byte-order mark, so the `.fai` of such a file works only with `FastaIndex`.

## Error handling

Parse errors raise `FastaParseError`:

```python
from fastatacular import FastaParseError, read_fasta

try:
    entries = read_fasta("malformed.fasta")
except FastaParseError as e:
    print(e.line)     # offending line number
    print(e.context)  # surrounding line content
```

Input that is not UTF-8, or a corrupt compressed file, also raises `FastaParseError`
(chained to the underlying `UnicodeDecodeError` or `OSError`).

Write errors raise `FastaWriteError`, whose `index` names the bad entry. Every entry is
validated before anything is written, so a failed `write_fasta` leaves no partial file.
Both errors subclass `FastaError` (a `ValueError`), so `except FastaError` catches either.

`SequenceEntry` is frozen but not hashable (its `extra` field is a dict), so key sets and
dicts by `entry.identifier`, not by the entry. A `FastaReader` is single-pass: iterate it
once, or open a new one to read the file again.

## Development

```bash
just install      # install dependencies
just test         # run tests
just test-v       # run tests (verbose)
just cov          # run tests with coverage
just lint         # ruff lint
just format       # ruff format
just check        # lint + type check + test
just build        # build the package
just clean        # remove cache files
```

## Citation

If you use fastatacular in research, please cite the archived software release. Machine-readable citation metadata is available in [`CITATION.cff`](https://github.com/tacular-omics/fastatacular/blob/main/CITATION.cff); GitHub's **Cite this repository** menu can render it as APA or BibTeX. DOI: [10.5281/zenodo.22926358](https://doi.org/10.5281/zenodo.22926358).

## License

[MIT](https://github.com/tacular-omics/fastatacular/blob/main/LICENSE)
