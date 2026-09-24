"""Decoy generation: methods, options, reproducibility, headers, Markov models, quality."""

from __future__ import annotations

import gzip
import io
import json
import re
from collections import Counter
from pathlib import Path

import pytest
from hypothesis import given
from hypothesis import strategies as st

from fastatacular import (
    DecoyError,
    FastaError,
    MarkovModel,
    SequenceEntry,
    is_decoy,
    load_markov_model,
    make_decoy_sequence,
    make_decoys,
    read_fasta,
    train_markov_model,
    write_decoy_fasta,
    write_fasta,
)
from fastatacular.decoys import AMINO_ACIDS, METHODS, MODELS

RANDOM_METHODS = ("shuffle", "debruijn", "markov")


def _proteome(n: int = 300, seed: int = 7) -> list[SequenceEntry]:
    """A synthetic proteome sampled from the human model (lengths 50-800)."""
    lengths = [50 + (i * 97) % 750 for i in range(n)]
    seqs = [
        make_decoy_sequence("A" * length, method="markov", seed=seed * 1000 + i) for i, length in enumerate(lengths)
    ]
    text = "".join(f">sp|P{i:05d}|PROT{i}_HUMAN Protein {i} OS=Homo sapiens OX=9606\n{s}\n" for i, s in enumerate(seqs))
    return read_fasta(io.StringIO(text))


@pytest.fixture(scope="module")
def proteome() -> list[SequenceEntry]:
    return _proteome()


# ---------------------------------------------------------------------------
# Methods and options
# ---------------------------------------------------------------------------


def test_reverse_and_pseudo_reverse() -> None:
    assert make_decoy_sequence("MPEPTIDEK", method="reverse") == "KEDITPEPM"
    assert make_decoy_sequence("MPEPTIDEK", method="reverse", keep_nterm=1, keep_cterm=1) == "MEDITPEPK"
    # Stretches between K/R are reversed; K and R stay put.
    assert make_decoy_sequence("ABCKDEFRGH", method="pseudo_reverse") == "CBAKFEDRHG"
    assert make_decoy_sequence("ABCKDEFRGH", method="reverse", keep_residues="KR") == "CBAKFEDRHG"
    assert make_decoy_sequence("ABCKDEFRGH", method="pseudo_reverse", keep_residues="") == "HGRFEDKCBA"


_residues = st.text(AMINO_ACIDS + "XU", min_size=0, max_size=80)


@given(
    seq=_residues,
    method=st.sampled_from(METHODS),
    keep=st.sampled_from(["", "KR", "KRP", "C"]),
    nterm=st.integers(0, 3),
    cterm=st.integers(0, 3),
    seed=st.integers(0, 10),
)
def test_kept_positions_stay_and_are_never_introduced(
    seq: str, method: str, keep: str, nterm: int, cterm: int, seed: int
) -> None:
    decoy = make_decoy_sequence(
        seq,
        method=method,
        seed=seed,
        keep_residues=keep,
        keep_nterm=nterm,
        keep_cterm=cterm,  # type: ignore[arg-type]
    )
    n = len(seq)
    assert len(decoy) == n
    lo, hi = min(nterm, n), max(min(nterm, n), n - cterm)
    assert decoy[:lo] == seq[:lo]
    assert decoy[hi:] == seq[hi:]
    for i in range(lo, hi):
        if seq[i] in keep:
            assert decoy[i] == seq[i]
        else:
            assert decoy[i] not in keep
    if method in ("reverse", "pseudo_reverse", "shuffle"):
        assert Counter(decoy) == Counter(seq)
    if method in ("markov", "debruijn"):
        # Residues outside the 20 standard amino acids stay in place.
        assert all(decoy[i] == seq[i] for i in range(n) if seq[i] not in AMINO_ACIDS)


@pytest.mark.parametrize("method", RANDOM_METHODS)
def test_same_seed_same_output_and_different_seed_differs(proteome: list[SequenceEntry], method: str) -> None:
    a = [e.sequence for e in make_decoys(proteome, method=method, seed=42)]  # type: ignore[arg-type]
    b = [e.sequence for e in make_decoys(proteome, method=method, seed=42)]  # type: ignore[arg-type]
    c = [e.sequence for e in make_decoys(proteome, method=method, seed=43)]  # type: ignore[arg-type]
    assert a == b
    assert a != c
    assert make_decoy_sequence("MKTAYIAKQRQISFVK", method=method, seed="s") == make_decoy_sequence(  # type: ignore[arg-type]
        "MKTAYIAKQRQISFVK",
        method=method,
        seed="s",  # type: ignore[arg-type]
    )


def test_unseeded_runs_differ() -> None:
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMK" * 2
    assert len({make_decoy_sequence(seq, method="shuffle") for _ in range(3)}) > 1


@pytest.mark.parametrize("method", ["reverse", "pseudo_reverse", "shuffle", "markov"])
def test_streaming_decoy_depends_only_on_its_target(proteome: list[SequenceEntry], method: str) -> None:
    full = [e.sequence for e in make_decoys(proteome, method=method, seed=3)]  # type: ignore[arg-type]
    alone = [e.sequence for e in make_decoys(proteome[5:9], method=method, seed=3)]  # type: ignore[arg-type]
    assert alone == full[5:9]


def test_make_decoys_is_lazy_for_streaming_methods() -> None:
    seen: list[int] = []

    def gen():  # type: ignore[no-untyped-def]
        for i in range(3):
            seen.append(i)
            yield SequenceEntry(identifier=f"p{i}", sequence="MKTAYIAK")

    it = make_decoys(gen(), method="shuffle", seed=1)
    assert seen == []
    next(it)
    assert seen == [0]


def test_debruijn_preserves_repeats() -> None:
    k = 2
    repeat = "WLPQHNGFEDCATSVYMIRK"
    a = "MSTNE" + repeat + "GGDDEE"
    b = "AAQQLLP" + repeat + "HHTTS"
    entries = [SequenceEntry(identifier="a", sequence=a), SequenceEntry(identifier="b", sequence=b)]
    da, db = (e.sequence for e in make_decoys(entries, method="debruijn", seed=5, k=k))
    ia, ib = a.index(repeat), b.index(repeat)
    # Positions whose (k+1)-mer lies wholly inside the repeat get the same label.
    assert da[ia + k : ia + len(repeat)] == db[ib + k : ib + len(repeat)]
    assert da != a


def test_debruijn_composition_tracks_target(proteome: list[SequenceEntry]) -> None:
    decoys = list(make_decoys(proteome, method="debruijn", seed=1))
    t = Counter("".join(e.sequence for e in proteome))
    d = Counter("".join(e.sequence for e in decoys))
    total = sum(t.values())
    assert sum(abs(t[a] - d[a]) for a in t) / total < 0.02


def test_markov_backs_off_and_keeps_nonstandard(tmp_path: Path) -> None:
    fasta = tmp_path / "tiny.fasta"
    fasta.write_text(">a\nACDACDACD\n>b\nXXACXX\n")
    model = train_markov_model(fasta, order=2)
    assert model.metadata["entries"] == 2
    decoy = make_decoy_sequence("MWYACDXW", method="markov", model=model, seed=1)
    assert len(decoy) == 8
    assert set(decoy) <= set("ACDX")
    assert decoy[6] == "X"


def test_markov_order_zero(tmp_path: Path) -> None:
    fasta = tmp_path / "t.fasta"
    fasta.write_text(">a\nAAAAAAAAAC\n")
    model = train_markov_model([fasta], order=0)
    assert model.composition["A"] == pytest.approx(0.9)
    assert set(make_decoy_sequence("MMMMMMMMMMMMMMMMMMMM", method="markov", model=model, seed=2)) <= {"A", "C"}


def test_all_kept_returns_target() -> None:
    assert make_decoy_sequence("KRKRKR", method="shuffle", keep_residues="KR", seed=1) == "KRKRKR"
    assert make_decoy_sequence("MK", method="markov", keep_nterm=1, keep_cterm=1, seed=1) == "MK"
    assert make_decoy_sequence("", method="markov", seed=1) == ""


# ---------------------------------------------------------------------------
# Headers, is_decoy, write_decoy_fasta
# ---------------------------------------------------------------------------


def test_decoy_headers(proteome: list[SequenceEntry]) -> None:
    target = proteome[0]
    [decoy] = make_decoys([target], method="reverse", prefix="rev_")
    assert decoy.identifier == "rev_" + target.identifier
    assert decoy.accession == target.accession
    assert (decoy.pname, decoy.os_name, decoy.ncbi_tax_id) == (target.pname, target.os_name, target.ncbi_tax_id)
    assert decoy.raw_header == "rev_" + target.raw_header
    assert is_decoy(decoy, prefix="rev_") and not is_decoy(target, prefix="rev_")
    assert is_decoy("DECOY_x") and not is_decoy("x")
    buf = io.StringIO()
    write_fasta([decoy], buf)
    assert buf.getvalue().startswith(">rev_sp|P00000|PROT0_HUMAN Protein 0 OS=Homo sapiens OX=9606\n")


def test_decoy_header_for_edited_entry_uses_current_fields() -> None:
    entry = SequenceEntry(identifier="x1", sequence="MKLV", gname="ABC")
    [decoy] = make_decoys([entry], method="reverse")
    assert decoy.raw_header == "DECOY_x1 GN=ABC"
    assert decoy.gname == "ABC"


@pytest.mark.parametrize("concatenate", [True, False])
def test_write_decoy_fasta(tmp_path: Path, proteome: list[SequenceEntry], concatenate: bool) -> None:
    src = tmp_path / "t.fasta"
    write_fasta(proteome, src)
    dst = tmp_path / "td.fasta"
    n = write_decoy_fasta(src, dst, method="pseudo_reverse", concatenate=concatenate)
    assert n == len(proteome)
    out = read_fasta(dst)
    decoys = [e for e in out if is_decoy(e)]
    assert len(decoys) == len(proteome)
    assert len(out) == (2 * len(proteome) if concatenate else len(proteome))
    if concatenate:
        assert out[: len(proteome)] == proteome
    buf = io.StringIO()
    write_decoy_fasta(io.StringIO(src.read_text()), buf, method="pseudo_reverse", concatenate=concatenate)
    assert buf.getvalue() == dst.read_text()


def test_write_decoy_fasta_refuses_a_decoy_database(tmp_path: Path) -> None:
    src = io.StringIO(">a\nMKV\n>DECOY_a\nVKM\n")
    with pytest.raises(DecoyError, match="already has decoy"):
        write_decoy_fasta(src, io.StringIO(), method="reverse")


# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"method": "flip"}, "Unknown decoy method"),
        ({"method": "reverse", "prefix": ""}, "prefix"),
        ({"method": "reverse", "prefix": "DE COY"}, "prefix"),
        ({"method": "reverse", "keep_nterm": -1}, "keep_nterm"),
        ({"method": "reverse", "keep_cterm": 1.5}, "keep_cterm"),
        ({"method": "reverse", "keep_residues": 3}, "keep_residues"),
        ({"method": "debruijn", "k": 0}, "k must be"),
        ({"method": "shuffle", "seed": 1.5}, "seed"),
        ({"method": "markov", "model": "martian"}, "Unknown Markov model"),
        ({"method": "markov", "keep_residues": AMINO_ACIDS}, None),
    ],
)
def test_bad_options_raise_decoy_error_eagerly(kwargs: dict, match: str | None) -> None:
    method = kwargs.pop("method")
    if match is None:  # every residue kept: valid, and the decoy is the target
        assert [e.sequence for e in make_decoys([SequenceEntry("a", "MKV")], method=method, **kwargs)] == ["MKV"]
        return
    with pytest.raises(DecoyError, match=match) as info:
        make_decoys(iter(()), method=method, **kwargs)  # raises before iteration
    assert isinstance(info.value, FastaError)


def test_model_that_excludes_every_free_residue_raises(tmp_path: Path) -> None:
    fasta = tmp_path / "k.fasta"
    fasta.write_text(">a\nKKKKRRRR\n")
    model = train_markov_model(fasta, order=1)
    with pytest.raises(DecoyError, match="zero probability"):
        make_decoy_sequence("MKV", method="markov", model=model, keep_residues="KR")


# ---------------------------------------------------------------------------
# Markov models: shipped, train, save, load, validation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", MODELS)
def test_shipped_models_load_with_attribution(name: str) -> None:
    model = load_markov_model(name)
    assert model.order == 2
    assert model.metadata["license"] == "CC BY 4.0"
    assert "UniProt" in str(model.metadata["attribution"])
    assert model.metadata["proteome"].startswith("UP")  # type: ignore[union-attr]
    assert sum(model.composition.values()) == pytest.approx(1.0)


@pytest.mark.parametrize("suffix", [".json", ".json.gz"])
def test_train_save_load_round_trip(tmp_path: Path, suffix: str) -> None:
    fasta = tmp_path / "t.fasta.gz"
    fasta.write_bytes(gzip.compress(b">a\nMKTAYIAKQR\n>b\nmktxayi\n"))
    model = train_markov_model(str(fasta), order=1, metadata={"note": "test"})
    assert model.metadata["note"] == "test"
    assert model.metadata["residues"] == 10 + 6
    idx = AMINO_ACIDS.index
    assert model.counts[1][idx("M") * 20 + idx("K")] == 2
    path = tmp_path / f"m{suffix}"
    model.save(path)
    back = load_markov_model(path)
    assert (back.order, back.counts, back.metadata, back.alphabet) == (
        model.order,
        model.counts,
        model.metadata,
        model.alphabet,
    )


def test_load_errors(tmp_path: Path) -> None:
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"format": "other"}))
    with pytest.raises(DecoyError, match="Cannot load"):
        load_markov_model(bad)
    bad.write_text("{not json")
    with pytest.raises(DecoyError, match="Cannot load"):
        load_markov_model(bad)
    with pytest.raises(DecoyError, match="Unknown Markov model"):
        load_markov_model(tmp_path / "missing.json")


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"order": 9, "counts": ()}, "order"),
        ({"order": 0, "counts": ((1,) * 20,), "alphabet": "AA"}, "alphabet"),
        ({"order": 1, "counts": ((1,) * 20,)}, "count tables"),
        ({"order": 0, "counts": ((1,) * 19,)}, "entries"),
        ({"order": 0, "counts": ((-1,) * 20,)}, "non-negative"),
        ({"order": 0, "counts": ((0,) * 20,)}, "no training data"),
    ],
)
def test_markov_model_validation(kwargs: dict, match: str) -> None:
    with pytest.raises(DecoyError, match=match):
        MarkovModel(**kwargs)


def test_train_errors(tmp_path: Path) -> None:
    with pytest.raises(DecoyError, match="order"):
        train_markov_model([], order=5)
    empty = tmp_path / "e.fasta"
    empty.write_text(">a\nXXXX\n")
    with pytest.raises(DecoyError, match="No standard"):
        train_markov_model(empty)


# ---------------------------------------------------------------------------
# Quality on a synthetic proteome (the human-proteome numbers: scripts/decoy_quality.py)
# ---------------------------------------------------------------------------

_TRYPSIN = re.compile(r"(?<=[KR])(?!P)")


def _peptides(seqs: list[str]) -> set[str]:
    return {p for s in seqs for p in _TRYPSIN.split(s) if 7 <= len(p) <= 40}


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("keep", [None, "KR"])
def test_quality(proteome: list[SequenceEntry], method: str, keep: str | None) -> None:
    targets = [e.sequence for e in proteome]
    decoys = [e.sequence for e in make_decoys(proteome, method=method, seed=11, keep_residues=keep)]  # type: ignore[arg-type]
    assert [len(d) for d in decoys] == [len(t) for t in targets]
    assert all(d != t for d, t in zip(decoys, targets, strict=True))
    t_comp, d_comp = Counter("".join(targets)), Counter("".join(decoys))
    total = sum(t_comp.values())
    assert sum(abs(t_comp[a] - d_comp[a]) for a in t_comp | d_comp) / total < 0.03
    d_peps = _peptides(decoys)
    assert len(d_peps & _peptides(targets)) / len(d_peps) < 0.01
    if keep == "KR" or method == "pseudo_reverse":
        # K/R stay put, so the decoy has the same cleavage-site positions (ignoring the proline rule).
        assert all(
            [i for i, c in enumerate(d) if c in "KR"] == [i for i, c in enumerate(t) if c in "KR"]
            for d, t in zip(decoys, targets, strict=True)
        )


def test_method_and_prefix_are_keyword_only() -> None:
    entry = SequenceEntry("a", "MKV")
    with pytest.raises(TypeError):
        make_decoys([entry], "reverse")  # type: ignore[misc]
    with pytest.raises(TypeError):
        make_decoy_sequence("MKV", "reverse")  # type: ignore[misc]
    with pytest.raises(TypeError):
        is_decoy(entry, "DECOY_")  # type: ignore[misc]
