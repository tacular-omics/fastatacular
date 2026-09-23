# Copilot instructions

The canonical guide for AI agents in this repo is [`CLAUDE.md`](../CLAUDE.md). Read it
first. The rules that matter most:

1. No runtime dependencies. Pure Python, `requires-python >= 3.12`.
2. Before pushing, run what CI runs:
   `uv run ruff check src tests && uv run ruff format --check src tests && uv run ty check src && uv run pytest tests`.
3. Parse failures raise `FastaParseError(msg, line=..., context=...)`; write failures
   raise `FastaWriteError`. Both subclass `ValueError`.
4. Keep the public API (`read_fasta`, `FastaReader`, `write_fasta`, `SequenceEntry`,
   errors) in step with the sibling `pefftacular` package, and list new public names in
   `src/fastatacular/__init__.py` `__all__`.
5. Never bump the version, tag or publish: only the tacular-omics overseer releases.
