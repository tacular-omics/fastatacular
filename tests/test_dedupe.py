"""Exact-sequence deduplication, duplicate reports and reader/writer composition."""

import io
import logging

import pytest

import fastatacular
from fastatacular import (
    FastaError,
    FastaReader,
    SequenceEntry,
    dedupe_sequences,
    duplicate_groups,
    read_fasta,
    write_fasta,
)


@pytest.fixture
def entries():
    return read_fasta(
        io.StringIO(
            ">sp|P69905|HBA_HUMAN  Hemoglobin OS=Homo sapiens GN=HBA1\nMVLSPADK\n"
            ">sp|Q9XYZ1|HBA_MOUSE Different header OS=Mus musculus\nMVLSPADK\n"
            ">tr|A0A1|A0A1_HUMAN Lowercase copy\nmvlspadk\n"
            ">sp|P2|B_HUMAN\tAnother protein\nPEPTIDE\n"
            ">NP_000549.1 Another identifier [Homo sapiens]\nPEPTIDE\n"
            ">sp|P3|C_HUMAN With stop marker\nPEPTIDE*\n"
            ">sp|P4|D_HUMAN Third copy\nMVLSPADK\n"
        )
    )


def test_public_exports():
    assert "dedupe_sequences" in fastatacular.__all__
    assert "duplicate_groups" in fastatacular.__all__


def test_keeps_first_entry_unchanged_in_input_order(entries):
    kept = list(dedupe_sequences(entries))
    assert [entry.identifier for entry in kept] == [
        "sp|P69905|HBA_HUMAN",
        "tr|A0A1|A0A1_HUMAN",
        "sp|P2|B_HUMAN",
        "sp|P3|C_HUMAN",
    ]
    assert all(actual is original for actual, original in zip(kept, [entries[i] for i in (0, 2, 3, 5)], strict=True))


def test_groups_report_only_dropped_identifiers_in_input_order(entries, caplog):
    with caplog.at_level(logging.WARNING, logger="fastatacular._dedupe"):
        groups = duplicate_groups(iter(entries))
    assert groups == {
        "sp|P69905|HBA_HUMAN": ["sp|Q9XYZ1|HBA_MOUSE", "sp|P4|D_HUMAN"],
        "sp|P2|B_HUMAN": ["NP_000549.1"],
    }
    assert list(groups) == ["sp|P69905|HBA_HUMAN", "sp|P2|B_HUMAN"]
    assert caplog.records == []


def test_empty_sequence_is_deduplicated():
    entries = [SequenceEntry("a", ""), SequenceEntry("b", ""), SequenceEntry("c", "MK")]
    assert list(dedupe_sequences(iter(entries))) == [entries[0], entries[2]]
    assert duplicate_groups(iter(entries)) == {"a": ["b"]}


@pytest.mark.parametrize("duplicates", ["first", "error"])
@pytest.mark.parametrize("sequences", [[], [""], ["MK", "mk", "MK*", ""]])
def test_empty_and_unique_inputs_need_no_warning(duplicates, sequences, caplog):
    entries = [SequenceEntry(str(i), sequence) for i, sequence in enumerate(sequences)]
    with caplog.at_level(logging.WARNING, logger="fastatacular._dedupe"):
        kept = list(dedupe_sequences(iter(entries), duplicates=duplicates))
        assert duplicate_groups(iter(entries)) == {}
    assert kept == entries
    assert all(actual is original for actual, original in zip(kept, entries, strict=True))
    assert caplog.records == []


def test_sequences_are_compared_even_when_identifiers_repeat():
    entries = [SequenceEntry("same", "MK"), SequenceEntry("same", "PE"), SequenceEntry("copy", "MK")]
    assert list(dedupe_sequences(entries)) == entries[:2]
    assert duplicate_groups(entries) == {"same": ["copy"]}


def test_groups_merge_under_a_repeated_first_identifier():
    """The key is an identifier, so two sequences first seen under one identifier share a list."""
    entries = [
        SequenceEntry("a", "MK"),
        SequenceEntry("a", "PE"),
        SequenceEntry("x", "MK"),
        SequenceEntry("y", "PE"),
    ]
    assert list(dedupe_sequences(entries)) == entries[:2]
    assert duplicate_groups(entries) == {"a": ["x", "y"]}


@pytest.mark.parametrize("dropped", [1, 3])
def test_first_policy_logs_one_summary_after_exhaustion(dropped, caplog):
    entries = [SequenceEntry(str(i), "MK") for i in range(dropped + 1)]
    with caplog.at_level(logging.WARNING, logger="fastatacular._dedupe"):
        result = dedupe_sequences(entries, duplicates="first")
        assert next(result) is entries[0]
        assert caplog.records == []
        assert list(result) == []
    noun = "entry" if dropped == 1 else "entries"
    assert [(record.name, record.levelname, record.getMessage()) for record in caplog.records] == [
        (
            "fastatacular._dedupe",
            "WARNING",
            f"dropped {dropped} {noun} whose sequence repeats an earlier one (kept the first of each)",
        )
    ]


def test_error_policy_stops_at_the_first_duplicate(entries, caplog):
    consumed = []

    def source():
        for entry in entries:
            consumed.append(entry)
            yield entry

    result = dedupe_sequences(source(), duplicates="error")
    assert consumed == []
    assert next(result) is entries[0]
    assert consumed == entries[:1]
    with pytest.raises(FastaError) as exc:
        next(result)
    assert str(exc.value) == "Entry 'sp|Q9XYZ1|HBA_MOUSE' repeats the sequence of 'sp|P69905|HBA_HUMAN'"
    assert exc.value.__notes__ == ['hint: duplicates="first" keeps the first entry for each sequence']
    assert consumed == entries[:2]
    assert caplog.records == []


@pytest.mark.parametrize("duplicates", ["nope", None, 1, []])
def test_invalid_policy_raises_at_call_time_without_touching_input(duplicates):
    class UnreadableEntries:
        def __iter__(self):
            pytest.fail("invalid options must be rejected before the input is accessed")

    with pytest.raises(FastaError) as exc:
        dedupe_sequences(UnreadableEntries(), duplicates=duplicates)
    assert str(exc.value) == f"duplicates must be 'error' or 'first', got {duplicates!r}"


def test_first_policy_reads_only_as_far_as_the_next_unique_entry(entries):
    consumed = []

    def source():
        for entry in entries:
            consumed.append(entry)
            yield entry

    result = dedupe_sequences(source())
    assert consumed == []
    assert next(result) is entries[0]
    assert consumed == entries[:1]
    assert next(result) is entries[2]
    assert consumed == entries[:3]
    assert list(result) == [entries[3], entries[5]]
    assert consumed == entries


def test_reader_writer_preserves_kept_headers_byte_for_byte(tmp_path):
    headers = [
        b">sp|P69905|HBA_HUMAN  Hemoglobin\tOS=Homo sapiens GN=HBA1",
        b">copy Different header",
        ">sp|P2|B_HUMAN\tProtein café OS=Mus musculus".encode(),
        b">empty  An entry without residues",
        b">empty_copy Another empty sequence",
    ]
    source = tmp_path / "merged.fasta"
    source.write_bytes(
        b"".join(
            header + b"\n" + sequence
            for header, sequence in zip(headers, [b"MK\n", b"MK\n", b"PE*\n", b"", b""], strict=True)
        )
    )
    destination = tmp_path / "unique.fasta"
    with FastaReader(source) as reader:
        write_fasta(dedupe_sequences(reader), destination)
    output_headers = [line for line in destination.read_bytes().splitlines() if line.startswith(b">")]
    assert output_headers == [headers[i] for i in (0, 2, 3)]
    assert [entry.sequence for entry in read_fasta(destination)] == ["MK", "PE*", ""]
