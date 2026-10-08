"""The ``compression`` parameter shared by the readers and writers."""

from __future__ import annotations

from typing import Literal

from fastatacular.errors import FastaError

Compression = Literal["infer", "gzip", "bz2", "xz"] | None
"""Compression of a FASTA file, as in pandas.

``"infer"`` (the default): readers detect gzip, bzip2 or xz from the magic bytes,
writers compress by the path suffix (``.gz``, ``.bz2``, ``.xz``). An explicit
format forces it; ``None`` means plain text.
"""

_VALID: tuple[Compression, ...] = ("infer", "gzip", "bz2", "xz", None)

# ``compression`` value -> internal kind ("gz", "bz2", "xz").
_KINDS = {"gzip": "gz", "bz2": "bz2", "xz": "xz"}
_NAMES = {kind: name for name, kind in _KINDS.items()}


def _check_compression(compression: object) -> Compression:
    """Return ``compression`` if it is valid, else raise :class:`FastaError` listing the valid values."""
    if compression is None or (isinstance(compression, str) and compression in _VALID):
        return compression
    err = FastaError(f"compression must be one of 'infer', 'gzip', 'bz2', 'xz' or None, got {compression!r}")
    err.add_note("hint: pass 'infer' (the default) to use the magic bytes or suffix, or None for plain text")
    raise err


def _is_text_handle(handle: object) -> bool:
    """Whether an open file object reads/writes ``str`` (has an ``encoding``)."""
    import io

    return isinstance(handle, io.TextIOBase) or hasattr(handle, "encoding")


__all__ = ["Compression"]
