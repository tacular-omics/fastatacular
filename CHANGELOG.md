# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added

- Decoy databases (`fastatacular.decoys`, re-exported at the top level): `make_decoys`,
  `make_decoy_sequence`, `write_decoy_fasta`, `is_decoy`, `DecoyError`. Methods
  `reverse`, `pseudo_reverse`, `shuffle`, `debruijn` (repeat-preserving, Moosa et al.
  2020) and `markov`. All take `keep_residues`, `keep_nterm`, `keep_cterm` and are
  reproducible with `seed`. Guide: `docs/decoys.md`.
- Markov models: `MarkovModel`, `train_markov_model`, `load_markov_model`, and order-2
  models for human, mouse, yeast and E. coli trained on UniProtKB/Swiss-Prot release
  2026_03 (CC BY 4.0, see `src/fastatacular/data/markov/NOTICE.md`).
- Compressed input: `read_fasta` and `FastaReader` read gzip, bzip2 and xz files given
  as paths, detected from the magic bytes. Paths are opened once, so pipes and FIFOs work; `bz2` and
  `lzma` are imported only when needed.
- PEFF files read as plain FASTA: `#` lines before the first `>` header (the PEFF file
  header) are skipped instead of raising "Sequence data appears before any '>' header".

### Performance

- Reading a path is about 1.3-1.5x faster (human reference proteome, 63k entries): whole
  records are split out of 1 MiB chunks instead of handling each line in Python, and
  `KEY=value` pairs are found from each `=` instead of with a lazy regex. Output is
  byte-identical on a 15-file UniProt/NCBI corpus.

### Fixed

- Input that is not valid UTF-8 raises `FastaParseError` (chained to the
  `UnicodeDecodeError`) instead of a bare `UnicodeDecodeError`. A corrupt or truncated
  compressed file raises `FastaParseError` too.

## [1.0.0] (2026-09-23)

First stable release: the public API is now stable and follows semantic versioning.

### Added

- `py.typed` marker in the wheel (`Typing :: Typed`); Python 3.14 classifier; SPDX
  `license = "MIT"` metadata.
- `FastaError(ValueError)`, exported base class of `FastaParseError` and
  `FastaWriteError`; `except FastaError` catches every fastatacular error.
- Errors carry a `hint` (also shown as a traceback note). `FastaWriteError` has an
  `index` (0-based entry position) and its message is prefixed `Entry N: `.

### Changed

- `SequenceEntry` is explicitly unhashable (`__hash__ = None`): its `extra` field is a
  dict. `hash(entry)` still raises `TypeError`, now before looking at the fields, and
  `isinstance(entry, collections.abc.Hashable)` is `False`.
- `FastaReader.__enter__` is typed to return `Self`. Iterating a reader a second time is
  documented as continuing where the first iteration stopped (empty after a full pass).
- Classifier `Development Status :: 5 - Production/Stable`.
- `KEY=value` keys are only recognised at the start of the description or after
  whitespace. `>x Protein(EC=2.7.1) OS=Homo sapiens` now gives
  `pname='Protein(EC=2.7.1)'` (was `None` with `extra={'EC': '2.7.1)'}`), and NCBI
  `[organism=...] [gene=...]` text stays in `pname` instead of becoming `extra` keys.
  A key that starts the description is still a key (`>x pH=7 sensor` gives
  `extra={'pH': '7 sensor'}`). Real UniProt, UniRef, UniParc and NCBI headers parse as before.
- A pipe identifier's `prefix` may be any text without `|` or whitespace (was
  `[A-Za-z0-9]+`), so decoy and contaminant ids such as `Reverse_sp|P1|X_HUMAN`,
  `rev_sp|...`, `DECOY-0-sp|...` and `contam_sp|...` now get `prefix`, `accession` and
  `entry_name` instead of `None`. Other ids with such a prefix (e.g.
  `NW_001494075|IGHJ1-1*03|Bos`) are now split into those fields too.

### Fixed

- `write_fasta` raises `FastaWriteError` (with `index` and `hint`) instead of writing
  a header that reads back differently: free text holding a `KEY=` token
  (`os_name="Homo sapiens GN=FAKE"`, `pname="Protein X=1 like"`), leading/trailing
  whitespace in a written field, an `extra` key that is not one
  `[A-Za-z_][A-Za-z0-9_]*` word (`extra={"bad key": "v"}`) or that a typed field owns
  (`extra={"OS": "Mouse"}`). Unedited parsed entries are written verbatim as before.
- `write_fasta` raises `FastaWriteError` instead of `TypeError` when `line_width` is
  not an int (`None`, `60.0`, `True`).
- `scripts/release_version.py sync --set X.Y.Z` also sets CITATION.cff
  `date-released` to today (adding the field if missing); CITATION.cff now has
  `date-released`.
- `write_fasta` validates every entry before writing, so an unwritable entry no longer
  leaves a partial file: a path `dest` is not created or truncated and nothing is
  written to a handle.
- Editing a parsed entry (`dataclasses.replace(entry, gname="XYZ")`) now writes the
  edit. `raw_header` is written verbatim only while it still parses to the entry's
  fields, so unedited entries still round-trip byte-exact; a `raw_header` that
  disagrees with the fields is ignored and the header is rebuilt.
- Files starting with a UTF-8 byte-order mark, and whitespace-only lines before the
  first header, no longer raise "Sequence data appears before any '>' header".
- Clearing a field on a parsed entry that has no protein name (for example
  `dataclasses.replace(entry, gname=None)` on `>x OS=Homo sapiens GN=A`) now keeps it
  cleared. The writer no longer reads the stale parsed `description` back in; an edited
  `description` is still used.
- `write_fasta` raises `FastaWriteError` instead of writing a corrupt file when a header
  field contains a line break, the identifier or sequence contains whitespace, or a
  wrapped sequence line would start with `>` or `;`.

### Tests

- Headers are checked against UniProt's FASTA header help page examples, real
  UniProtKB headers (fields taken from the entries' JSON records) and NCBI headers
  (`tests/reference/`).
- Hypothesis property tests: write/read round trip, edits surviving a write, byte-exact
  raw headers, BOM/CRLF/blank lines, and typed errors for malformed input.

## [0.1.3] (2026-09-23)

### Fixed

- Writing an entry with an empty `raw_header` and no `pname` no longer repeats the
  `KEY=value` pairs already present in `description`.
- A tab now separates the identifier from the description, as in UniProt, BLAST,
  Biopython and samtools (`>x\tdesc` gives identifier `x`).
- Internal whitespace inside sequence lines is removed, as documented.
- `just lint` and `just format` now cover `tests`, matching CI.

## [0.1.2] (2026-09-23)

### Added

- Releases are archived on Zenodo (`.zenodo.json`); no code changes.

## [0.1.1] (2026-09-23)

* Publish from GitHub Actions with PyPI trusted publishing (`publish.yml`);
  release metadata is checked against the tag.
* Keep `__version__` and `CITATION.cff` in sync with
  `scripts/release_version.py` (`just set-version X.Y.Z`).
* CI tests Python 3.12-3.14 on Linux plus macOS and Windows, the lowest
  direct dependency versions, and the built wheel.
* Standardize citation and package metadata.
* Add `CHANGELOG.md` and `CITATION.cff`; fix the README CI badge.

## [0.1.0] (2026-05-15)

### Added
- Initial release: a pure-Python library for reading and writing FASTA sequence files, with optional parsing of UniProt-style description keys and pipe-delimited identifiers.
