# Contributing to fastatacular

Bug reports, fixes and new features are welcome. Use the
[issue tracker](https://github.com/tacular-omics/fastatacular/issues) for bugs and
questions, with a minimal example and the fastatacular and Python versions. Report
security problems privately as described in [SECURITY.md](SECURITY.md).

## Claiming an issue

Comment on the issue to say you are taking it before you start, and check that nobody has
already done so. One issue, one contributor: when two pull requests arrive for the same issue,
the first one that was claimed is reviewed and the other is closed, however good it is. If a
claim has had no pull request or update for 14 days, comment to take it over.

Issues labelled `good first issue` need no knowledge of the codebase. Issues labelled
`help wanted` do: each one names the files to read first and the test that shows it is done.

## Development setup

You need Python 3.12 or newer, [uv](https://docs.astral.sh/uv/) and
[just](https://just.systems/).

```bash
git clone https://github.com/tacular-omics/fastatacular.git
cd fastatacular
just install
```

`just install` runs `uv sync`. fastatacular declares no extras, so that single command
installs the whole development environment, including the dev group the tests need
(pytest, Hypothesis, pandas and polars). `uv sync --all-extras` adds nothing, because
there is no extra to add.

Check the install by reading the passed and skipped counts at the end of a test run, never
the collected total. The pandas and polars tests use `pytest.importorskip`, which collects
a test and then skips it, so a run missing those two libraries collects exactly as many
tests as a complete one. A complete environment skips nothing here.

`just --list` shows the other recipes.

Install the commit hooks once per clone with `uvx pre-commit install`. On every
commit they run `ruff check`, `ruff format --check` and a few file checks
(`uvx pre-commit run --all-files` runs them on everything). Type checks and
tests are not in the hooks, so run `just check` for those.

## Running checks

CI runs these on every pull request. Run them before you push:

```bash
uv run ruff check src tests
uv run ruff format --check src tests
uv run ty check src
uv run pytest tests            # fast; CI uses HYPOTHESIS_PROFILE=ci for more examples
python scripts/release_version.py check
```

## Pull requests

- Branch from `main` and keep each pull request to one change.
- Add tests for new behaviour and fixed bugs.
- Add a line under `## [Unreleased]` in `CHANGELOG.md` for user-visible changes.
- Do not change the version number. The maintainers set it when releasing.
- CI must pass: lint, type check, tests on Python 3.12 to 3.14, and the wheel build.

## CLAUDE.md and AI agents

`CLAUDE.md` is the maintained developer guide: commands, architecture,
conventions and known gotchas. `AGENTS.md` points other coding agents to it.
Read it before a larger change, and update it in the same pull request when your
change makes part of it wrong. AI-assisted contributions are welcome, and you are
responsible for reviewing and testing what you submit.
