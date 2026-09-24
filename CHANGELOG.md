# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added

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

### Fixed

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
