# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

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
