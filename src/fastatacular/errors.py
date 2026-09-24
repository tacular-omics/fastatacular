"""FASTA-specific error types.

Every exception raised by fastatacular for bad input or unwritable entries
derives from :class:`FastaError`, a :class:`ValueError` subclass, so one
``except FastaError`` catches both parse and write failures.
"""


class FastaError(ValueError):
    """Base class for all fastatacular errors."""


class FastaParseError(FastaError):
    """Raised when FASTA input cannot be parsed.

    Attributes:
        line: 1-based line number the failure was detected on, if known.
        context: The offending line text.
        hint: A short suggestion for how to fix the input, if any.
    """

    def __init__(
        self,
        message: str,
        *,
        line: int | None = None,
        context: str | None = None,
        hint: str | None = None,
    ) -> None:
        self.line = line
        self.context = context
        self.hint = hint
        super().__init__(message if line is None else f"Line {line}: {message}")
        if hint is not None:
            self.add_note(f"hint: {hint}")


class FastaWriteError(FastaError):
    """Raised when a model object cannot be serialized to FASTA.

    Attributes:
        index: 0-based position of the offending entry in the input, if known.
        hint: A short suggestion for how to fix the entry, if any.
    """

    def __init__(self, message: str, *, index: int | None = None, hint: str | None = None) -> None:
        self.index = index
        self.hint = hint
        super().__init__(message if index is None else f"Entry {index}: {message}")
        if hint is not None:
            self.add_note(f"hint: {hint}")
