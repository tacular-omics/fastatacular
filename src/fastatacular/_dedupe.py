"""Remove entries with repeated sequences and report the dropped identifiers."""

import logging
from collections.abc import Iterable, Iterator

from fastatacular._index import Duplicates, _check_duplicates
from fastatacular._models import SequenceEntry
from fastatacular.errors import FastaError

_log = logging.getLogger(__name__)


def dedupe_sequences(entries: Iterable[SequenceEntry], *, duplicates: Duplicates = "first") -> Iterator[SequenceEntry]:
    """Yield the first entry for each exact sequence, unchanged and in input order.

    Comparison is case-sensitive, includes trailing ``*`` characters and treats an
    empty sequence like any other sequence. Entries are consumed lazily; only the
    distinct sequence strings and their first identifiers are retained.

    Args:
        entries: Entries, for example from a :class:`~fastatacular.FastaReader`.
        duplicates: ``"first"`` drops later copies and logs one warning with the
            dropped count when the iterator is exhausted (none if all are unique).
            ``"error"`` raises at the first entry whose sequence has appeared before.

    Raises:
        FastaError: For an invalid option, at call time before touching ``entries``;
            or for a repeated sequence during iteration with ``duplicates="error"``.
    """
    duplicates = _check_duplicates(duplicates)

    def generate() -> Iterator[SequenceEntry]:
        first_by_sequence: dict[str, str] = {}
        dropped = 0
        for entry in entries:
            if entry.sequence in first_by_sequence:
                if duplicates == "error":
                    error = FastaError(
                        f"Entry {entry.identifier!r} repeats the sequence of {first_by_sequence[entry.sequence]!r}"
                    )
                    error.add_note('hint: duplicates="first" keeps the first entry for each sequence')
                    raise error
                dropped += 1
            else:
                first_by_sequence[entry.sequence] = entry.identifier
                yield entry
        if dropped:
            _log.warning(
                "dropped %d entr%s whose sequence repeats an earlier one (kept the first of each)",
                dropped,
                "y" if dropped == 1 else "ies",
            )

    return generate()


def duplicate_groups(entries: Iterable[SequenceEntry]) -> dict[str, list[str]]:
    """Report later identifiers under the first identifier with that exact sequence.

    Consume ``entries`` once, comparing sequences just as :func:`dedupe_sequences`
    does. Include only groups with duplicates; each list contains the dropped
    identifiers in input order, without the retained first identifier. No warning
    is logged. Like deduplication, this retains each distinct sequence string.

    The key is an identifier, not a sequence, so a database that repeats an
    identifier across two different sequences reports both their duplicates under
    that one key: entries ``("a", "MK") ("a", "PE") ("x", "MK") ("y", "PE")`` give
    ``{"a": ["x", "y"]}``, and which dropped identifier belongs to which sequence is
    no longer visible. Read the kept entries from :func:`dedupe_sequences` when
    identifiers are not unique.
    """
    first_by_sequence: dict[str, str] = {}
    groups: dict[str, list[str]] = {}
    for entry in entries:
        if entry.sequence in first_by_sequence:
            first_identifier = first_by_sequence[entry.sequence]
            groups.setdefault(first_identifier, []).append(entry.identifier)
        else:
            first_by_sequence[entry.sequence] = entry.identifier
    return groups
