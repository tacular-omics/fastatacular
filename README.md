# fastatacular

[![PyPI](https://img.shields.io/pypi/v/fastatacular)](https://pypi.org/project/fastatacular/)
[![Python Package](https://github.com/tacular-omics/fastatacular/actions/workflows/ci.yml/badge.svg)](https://github.com/tacular-omics/fastatacular/actions/workflows/ci.yml)
[![License](https://img.shields.io/github/license/tacular-omics/fastatacular)](https://github.com/tacular-omics/fastatacular/blob/main/LICENSE)
[![Python](https://img.shields.io/pypi/pyversions/fastatacular)](https://pypi.org/project/fastatacular/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22926358.svg)](https://doi.org/10.5281/zenodo.22926358)

A small, dependency-free library for reading and writing [FASTA](https://en.wikipedia.org/wiki/FASTA_format) sequence files in Python. It's built for proteomics and genomics pipelines that need fast, predictable FASTA parsing without pulling in a bioinformatics megapackage.

It understands UniProt-style description keys (`OS=`, `OX=`, `GN=`, `PE=`, `SV=`) and pipe-delimited identifiers (`sp|P12345|EX_HUMAN`, `gi|12345|ref|NP_000001.1|`) out of the box, so you get structured fields instead of a header string to parse yourself.

**Web app:** the [FASTA Toolkit](https://tacular-omics.github.io/fastatacular/) makes decoys, checks decoy quality, cleans up, merges and converts FASTA files in your browser (nothing is uploaded). To browse entries, use the [FASTA Viewer](https://pgarrett-scripps.github.io/fastaviewer/).

## Highlights

- **Zero dependencies** — pure Python, nothing else to install.
- **Two ways to read** — `read_fasta` for the whole file at once, `FastaReader` to stream entries lazily without loading everything into memory.
- **Compressed input** — `.gz`, `.bz2` and `.xz` files are read directly (detected from the file's magic bytes).
- **UniProt headers parsed for you** — accession, organism, gene name, protein existence, and sequence version come back as typed fields, not a string you have to split yourself.
- **Round-trip safe** — entries produced by `read_fasta` write back out byte-for-byte compatible headers.
- **Actionable parse errors** — `FastaParseError` reports the offending line number and surrounding context.
- **Shares its API shape with [pefftacular](https://github.com/tacular-omics/pefftacular)**, the PEFF (PSI Extended FASTA) sibling library, so switching formats doesn't mean relearning the interface.

## Install

```bash
pip install fastatacular
```

Dev install:

```bash
just install
```

## Quick start

**read_fasta** — load everything into memory at once:

```python
from fastatacular import read_fasta

entries = read_fasta("proteins.fasta")
for entry in entries:
    print(entry.identifier, len(entry.sequence))
```

**FastaReader** — iterate lazily without loading the full file:

```python
from fastatacular import FastaReader

with FastaReader("proteins.fasta") as reader:
    for entry in reader:
        process(entry)
```

**Compressed files** are read transparently: gzip, bzip2 and xz, detected from the magic
bytes, not the file name. Pipes, FIFOs and `/dev/stdin` work too:

```python
entries = read_fasta("uniprot_sprot.fasta.gz")
```

A PEFF file also reads as plain FASTA: its `#` file header lines are skipped and each
entry keeps its identifier, sequence and raw description. Use
[pefftacular](https://github.com/tacular-omics/pefftacular) to parse the PEFF annotations.

## Data model

Each entry is a `SequenceEntry`:

| Field | Type | Description |
|---|---|---|
| `identifier` | `str` | Token immediately after `>` (e.g. `sp|P12345|EX_HUMAN`) |
| `sequence` | `str` | Concatenated sequence with whitespace stripped |
| `prefix` | `str \| None` | Database prefix (`sp`, `tr`, `gi`, ...) when the id is pipe-delimited |
| `accession` | `str \| None` | Second pipe field (e.g. `P12345` in `sp\|P12345\|EX_HUMAN`) |
| `entry_name` | `str \| None` | Third pipe field on UniProt ids (e.g. `EX_HUMAN`) |
| `description` | `str \| None` | Free text after the identifier |
| `pname` | `str \| None` | Protein name (description text, minus `KEY=value` pairs) |
| `gname` | `str \| None` | Gene name (`GN=`) |
| `os_name` | `str \| None` | Organism name (`OS=`) |
| `ncbi_tax_id` | `int \| None` | NCBI taxonomy ID (`OX=`) |
| `pe` | `int \| None` | Protein existence level (`PE=`) |
| `sv` | `int \| None` | Sequence version (`SV=`) |
| `extra` | `dict[str, str]` | Any other `KEY=value` pairs found in the header |
| `raw_header` | `str` | The original header line (without leading `>`) |

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

Non-standard `KEY=value` pairs are captured in `entry.extra`. Headers with no `KEY=value` tokens leave `description` and `pname` populated and `extra` empty.

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

`dest` accepts a path string, a `pathlib.Path`, or a text-mode file object.

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
- `.fai` caveats: the name column is the first header word (mapped to the accession on
  load with `key="accession"`). Like samtools, `write_fai()` needs every sequence line of an
  entry but the last to have the same length, and no comment or blank lines inside a
  sequence; otherwise it raises `FastaError` (rewrite the file with `write_fasta` first).
  `from_fai` checks each entry's header against the `.fai` name but does not re-count
  residues, so rebuild the `.fai` whenever the FASTA changes. Entries with no sequence
  raise (samtools skips them). samtools and pysam refuse a FASTA file that starts with a
  byte-order mark, so the `.fai` of such a file works only with `FastaIndex`.

## Tables with pandas or polars

`to_records(source)` returns one plain `dict` per entry, so any data-frame library can
take the result directly. fastatacular does not ship or require pandas or polars;
install whichever you use. (The test suite runs these examples only when the library is
installed.)

```python
import pandas as pd
import polars as pl

from fastatacular import to_records

records = to_records("human.fasta")   # a path, an open text handle, or entries
df = pd.DataFrame(records)
human = pl.DataFrame(records).filter(pl.col("ncbi_tax_id") == 9606)
```

`FastaReader.to_records()` and `SequenceEntry.to_record()` give the same dicts. Every
record has these keys, in this order (`fastatacular.RECORD_KEYS`):

| key | type | value |
|---|---|---|
| `identifier` | str | text after `>` up to the first whitespace |
| `prefix`, `accession`, `entry_name` | str or None | `sp`, `P12345`, `NAME` from `sp\|P12345\|NAME` |
| `pname` | str or None | protein name (description before the first `KEY=`) |
| `gname`, `os_name` | str or None | `GN=`, `OS=` |
| `ncbi_tax_id`, `pe`, `sv` | int or None | `OX=`, `PE=`, `SV=` |
| `description` | str or None | everything after the identifier |
| `extra` | str or None | other `KEY=value` pairs, space-separated |
| `raw_header` | str | the header line without `>` |
| `length` | int | sequence length |
| `sequence` | str | the residues |

The same fields in [pefftacular](https://github.com/tacular-omics/pefftacular) records have
other names where each package follows its own model. Rename these to put both in one
frame:

| field | fastatacular | pefftacular |
|---|---|---|
| database prefix | `prefix` | `prefix` |
| accession | `accession` | `db_unique_id` |
| entry name | `entry_name` | `id` |
| protein name, gene | `pname`, `gname` | `pname`, `gname` |
| organism name | `os_name` | `tax_name` |
| taxon id, PE, SV | `ncbi_tax_id`, `pe`, `sv` | `ncbi_tax_id`, `pe`, `sv` |
| other keys | `extra`, `KEY=value` pairs | `extra`, `\Key=value` pairs |
| length, residues | `length`, `sequence` | `length`, `sequence` |

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
