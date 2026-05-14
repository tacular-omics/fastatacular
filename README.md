# fastatacular

[![License](https://img.shields.io/github/license/tacular-omics/fastatacular)](LICENSE)

Pure-Python library for reading and writing [FASTA](https://en.wikipedia.org/wiki/FASTA_format) sequence files, with optional parsing of UniProt-style description keys (`OS=`, `OX=`, `GN=`, `PE=`, `SV=`) and pipe-delimited identifiers (`sp|P12345|EX_HUMAN`, `gi|12345|ref|NP_000001.1|`).

It's the plain-FASTA companion to [pefftacular](https://github.com/tacular-omics/pefftacular) and ships with the same `read_*` / `*Reader` / `write_*` shape.

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

## Data model

Each entry is a `SequenceEntry`:

| Field | Type | Description |
|---|---|---|
| `identifier` | `str` | Token immediately after `>` (e.g. `sp|P12345|EX_HUMAN`) |
| `sequence` | `str` | Concatenated sequence with whitespace stripped |
| `prefix` | `str \| None` | Database prefix (`sp`, `tr`, `gi`, ...) when the id is pipe-delimited |
| `accession` | `str \| None` | First pipe field (e.g. `P12345`) |
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

If `raw_header` is set on an entry (as it is on every entry produced by `read_fasta`), the writer round-trips it verbatim. Otherwise the header is rebuilt from the structured fields.

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

Write errors raise `FastaWriteError`.

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

## License

[MIT](LICENSE)
