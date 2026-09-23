# fastatacular — Claude Code Guide

## Project Overview

fastatacular is a small, pure-Python (no runtime dependencies, Python >= 3.12) library
for reading and writing FASTA sequence files. It parses UniProt-style description keys
(`OS=`, `OX=`, `GN=`, `PE=`, `SV=`) and pipe-delimited identifiers
(`sp|P12345|EX_HUMAN`, `gi|12345|ref|NP_000001.1|`) into typed fields on a frozen
`SequenceEntry` dataclass, and writes entries back out with round-trip-exact headers.

Place in the tacular-omics graph: tier 0, no sibling dependencies, and no core package
currently depends on it. Its API shape deliberately mirrors the sibling
[`pefftacular`](https://github.com/tacular-omics/pefftacular) (PEFF, PSI Extended
FASTA): `read_*`, a lazy `*Reader` context manager, `write_*`, a `*ParseError` with
`line`/`context`. Keep the two consistent when changing either.

## Commands

```bash
just install          # uv sync
just test             # uv run pytest tests
just test-v           # pytest -v
just test-file FILE   # pytest one file, verbose
just cov              # pytest with term-missing coverage
just test-cov         # coverage XML + junit.xml (what CI uploads to Codecov)
just lint             # uv run ruff check src
just format           # ruff isort fix + ruff format, src only
just ty               # uv run ty check src
just check            # lint + ty + test
just build            # uv build
just clean            # remove caches / .pyc
just check-version    # python scripts/release_version.py check
```

CI (`.github/workflows/ci.yml`) is stricter than `just lint`/`just format`: it runs
`ruff check src tests`, `ruff format --check src tests`, `ty check src`,
`scripts/release_version.py check`, then pytest on 3.12-3.14 (Linux) plus macOS and
Windows, a lowest-direct-deps run, and a built-wheel import. Before pushing, run:

```bash
uv run ruff check src tests && uv run ruff format --check src tests && uv run ty check src && uv run pytest tests
```

## Architecture

```
src/fastatacular/
  __init__.py   # public re-exports + __version__ (the version source, see Releasing)
  _models.py    # SequenceEntry: frozen, slotted dataclass, one per FASTA record
  _parser.py    # header regexes, _parse_header_line, _iter_entries (the streaming core),
                # FastaReader (context manager) and read_fasta (eager list)
  _writer.py    # _build_header_line (raw_header if it still matches the fields, else rebuild), write_fasta
  errors.py     # FastaParseError(ValueError) with .line/.context; FastaWriteError(ValueError)
tests/
  test_basic.py      # version smoke test
  test_reader.py     # header parsing, comments, blank lines, error cases
  test_writer.py     # wrapping, header rebuild, write errors
  test_roundtrip.py  # read -> write -> read equality
scripts/release_version.py   # version sync/check used by just set-version / check-version and CI
```

Data flow (read): text lines -> `_iter_entries` skips blank lines and `;` comments ->
`>` line goes to `_parse_header_line` (identifier = text up to the first whitespace;
`_UNIPROT_ID` / `_PIPE_ID` regexes fill `prefix`/`accession`/`entry_name`;
`_KV_PATTERN` pulls `KEY=value` pairs) -> all whitespace is removed from sequence lines and they are joined ->
`_build_entry` makes a `SequenceEntry`, raising if the sequence is empty.

Data flow (write): `_write_entry` validates identifier/sequence -> header is
`raw_header` verbatim if non-empty and it still parses to the entry's fields, else rebuilt as
`identifier [pname|description] OS= OX= GN= PE= SV= extra...` -> sequence wrapped at
`line_width` (default 60, `<= 0` means one line).

## Public API

Exported from `fastatacular` (`__all__`):

- `read_fasta(source)` -> `list[SequenceEntry]`: read a whole file (path or text handle).
- `FastaReader(source)`: lazy iterator; must be used as `with FastaReader(p) as r: for e in r: ...`.
- `write_fasta(entries, dest, *, line_width=60)`: write entries to a path or text handle.
- `SequenceEntry`: frozen dataclass (`identifier`, `sequence`, `prefix`, `accession`,
  `entry_name`, `description`, `pname`, `gname`, `os_name`, `ncbi_tax_id`, `pe`, `sv`,
  `extra`, `raw_header`).
- `FastaParseError`: `ValueError` subclass; `.line` (1-based) and `.context`; message
  is prefixed `Line N: `.
- `FastaWriteError`: `ValueError` subclass for unwritable entries.
- `__version__`.

## Conventions

- Docstrings: short Google style; type hints carry the types. Model fields are
  documented in the `SequenceEntry` class docstring.
- Typing: full annotations, `from __future__ import annotations` in `_parser.py` and
  `_writer.py`; `ty check src` must pass. Ruff rules E, W, F, I, B, UP; line length 120.
- Errors: parsing raises `FastaParseError` with `line=` and `context=`; writing raises
  `FastaWriteError`. Non-integer `OX`/`PE`/`SV` values do not raise, they fall into
  `extra`. No logging and no warnings anywhere in the package.
- Private modules are underscore-prefixed; add new public names to `__init__.__all__`.
- Tests: plain pytest functions in `tests/test_<area>.py`, `tmp_path` for files,
  `io.StringIO` for in-memory input.

## Gotchas

- **A `str` argument is always a file path.** `read_fasta(">x\nA\n")` raises
  `FileNotFoundError`; wrap FASTA text in `io.StringIO`.
- **`raw_header` is written only while it matches.** The writer re-parses
  `raw_header`; if that gives the entry's current fields it is written verbatim
  (byte-exact round trip), otherwise the header is rebuilt, so
  `dataclasses.replace(entry, gname="X")` writes `GN=X`. A hand-built `raw_header`
  that disagrees with the fields is ignored.
- **Header rebuild from `description`.** When `raw_header` is empty and `pname` is
  `None`, the writer uses the name text of `description` (before its first `KEY=`),
  then the structured fields and `extra`; `KEY=value` pairs in `description` are only
  written if no structured field or `extra` key covers them.
- **Identifier ends at the first whitespace** (space or tab), like UniProt, BLAST,
  Biopython and samtools.
- **All whitespace is removed from sequence lines** (`A C D` -> `ACD`). Case and `*`
  are kept.
- `description` is the full text after the identifier *including* `KEY=value` pairs;
  `pname` is the part before the first key (or all of it when there are no keys).
- `accession` is the *second* pipe field (`sp|P12345|...` -> `P12345`); `entry_name`
  is only set for exactly three-field `db|ACC|NAME` ids.
- `SequenceEntry` is frozen but not hashable (`extra` is a `dict`).
- `FastaReader` errors are raised lazily, at the bad entry during iteration, not on
  open. Calling `iter()` outside `with` raises `RuntimeError`.
- An entry with no sequence lines (including the last one in the file) raises
  `FastaParseError`; `;` comment lines and blank or whitespace-only lines are skipped anywhere; a leading UTF-8 BOM is ignored.
- `just lint`/`just format` cover `src` only, but CI also checks `tests`.

## Releasing

Only the tacular-omics overseer bumps versions or publishes. See `just --list`
(`set-version`, `sync-version`, `check-version`) and the workspace CLAUDE.md release
checklist. Version source: `__version__` in `src/fastatacular/__init__.py`
(`[tool.hatch.version]`), mirrored in `CITATION.cff` by `scripts/release_version.py`.
Changelog: `CHANGELOG.md` (`## [X.Y.Z] (YYYY-MM-DD)`, newest first, `[Unreleased]` on
top). Publishing is `publish.yml` on a GitHub release (PyPI trusted publishing); Zenodo
archives each release (`.zenodo.json`).

## Workspace note

This repo is also developed inside the tacular-omics uv workspace
(`~/Repos/tacular-omics/packages/fastatacular`); there `uv run` uses the shared
`.venv` and root `uv.lock`, not this repo's `uv.lock`. See the workspace CLAUDE.md.
