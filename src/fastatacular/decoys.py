"""Decoy protein databases for target-decoy false discovery rate estimation.

Five methods, all streaming and reproducible with ``seed``:

``reverse``
    Reverse the sequence.
``pseudo_reverse``
    Reverse each stretch between cleavage sites, keeping the sites (``K``/``R`` by
    default) in place, so tryptic decoy peptides keep the target peptides' lengths,
    compositions and masses.
``shuffle``
    Randomly permute the residues (per stretch between kept residues).
``debruijn``
    The repeat-preserving de Bruijn method of Moosa et al. (2020): every distinct
    (k+1)-mer in the target database is given one random replacement residue, so a
    sequence repeated in the target database is repeated in the decoys.
``markov``
    Sample a new sequence from an order-k Markov chain trained on a proteome. The
    package ships order-2 models for human, mouse, yeast and E. coli.

Every method takes ``keep_residues`` (residues that stay at their positions, e.g.
``"KR"`` or ``"KRP"``), ``keep_nterm`` and ``keep_cterm`` (how many terminal residues
stay). Replacement residues are never drawn from ``keep_residues``, so with
``keep_residues="KR"`` every method keeps the target's cleavage sites exactly.

References:
    Elias JE, Gygi SP (2007). Target-decoy search strategy for increased confidence in
    large-scale protein identifications by mass spectrometry. Nat Methods 4:207-214.

    Moosa JM, Guan S, Moran MF, Ma B (2020). Repeat-preserving decoy database for false
    discovery rate estimation in peptide identification. J Proteome Res 19(3):1029-1036.
    doi:10.1021/acs.jproteome.9b00555
"""

from __future__ import annotations

import gzip
import hashlib
import itertools
import json
import os
import random
import re
from bisect import bisect
from collections import Counter
from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path
from typing import IO, Literal, get_args

from fastatacular._models import SequenceEntry
from fastatacular._parser import _build_entry, _parse_header_line, read_fasta
from fastatacular._writer import _SEQ_LINE_WIDTH, _build_header_line, write_fasta
from fastatacular.errors import DecoyError

DecoyMethod = Literal["reverse", "pseudo_reverse", "shuffle", "debruijn", "markov"]
"""The decoy methods accepted by :func:`make_decoys`."""

METHODS: tuple[str, ...] = get_args(DecoyMethod)
"""Names of the decoy methods."""

MODELS: tuple[str, ...] = ("human", "mouse", "yeast", "ecoli")
"""Names of the pretrained Markov models shipped with the package."""

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
"""The 20 standard amino acids, the alphabet of the Markov models."""

DEFAULT_PREFIX = "DECOY_"
_MAX_ORDER = 4
_MAX_TRIES = 10
_MODEL_FORMAT = "fastatacular-markov/1"


# ---------------------------------------------------------------------------
# Markov models
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True, eq=False)
class MarkovModel:
    """An order-k Markov chain over the 20 standard amino acids.

    ``counts[o]`` holds the (o+1)-mer counts for every order ``o`` from 0 to
    ``order``, in lexicographic order of ``alphabet``: ``counts[0]`` is the residue
    composition, and with indices ``a, b, c`` into ``alphabet`` ``counts[2][a*400 + b*20 + c]``
    is the count of the 3-mer ``abc``. Lower orders are used at the start of a sequence and
    for contexts that never occurred in training (back-off).

    Build one with :func:`train_markov_model`, load a shipped or saved one with
    :func:`load_markov_model`, and save one with :meth:`save`.
    """

    order: int
    counts: tuple[tuple[int, ...], ...]
    metadata: dict[str, object] = field(default_factory=dict)
    alphabet: str = AMINO_ACIDS
    _tables: dict[frozenset[str], dict[str, tuple[str, list[float]]]] = field(
        default_factory=dict, init=False, repr=False
    )

    def __post_init__(self) -> None:
        if not isinstance(self.order, int) or not 0 <= self.order <= _MAX_ORDER:
            raise DecoyError(f"Markov order must be an int from 0 to {_MAX_ORDER}, got {self.order!r}")
        size = len(self.alphabet)
        if size == 0 or len(set(self.alphabet)) != size:
            raise DecoyError(f"Markov alphabet must be non-empty with unique letters, got {self.alphabet!r}")
        if len(self.counts) != self.order + 1:
            raise DecoyError(f"An order-{self.order} model needs {self.order + 1} count tables, got {len(self.counts)}")
        for o, table in enumerate(self.counts):
            if len(table) != size ** (o + 1):
                raise DecoyError(f"Count table {o} must have {size ** (o + 1)} entries, got {len(table)}")
            if any(not isinstance(c, int) or c < 0 for c in table):
                raise DecoyError(f"Count table {o} must hold non-negative ints")
        if sum(self.counts[0]) == 0:
            raise DecoyError("The Markov model has no training data (all counts are zero)")

    @property
    def composition(self) -> dict[str, float]:
        """Residue frequencies from the order-0 counts."""
        total = sum(self.counts[0])
        return {aa: c / total for aa, c in zip(self.alphabet, self.counts[0], strict=True)}

    def save(self, path: str | Path) -> None:
        """Write the model as JSON (gzip-compressed if ``path`` ends in ``.gz``)."""
        data = {
            "format": _MODEL_FORMAT,
            "order": self.order,
            "alphabet": self.alphabet,
            "metadata": self.metadata,
            "counts": [list(t) for t in self.counts],
        }
        text = json.dumps(data, separators=(",", ":"), sort_keys=True)
        path = Path(path)
        if path.suffix == ".gz":
            path.write_bytes(gzip.compress(text.encode(), mtime=0))
        else:
            path.write_text(text, encoding="utf-8")

    def _sampling_tables(self, exclude: frozenset[str]) -> dict[str, tuple[str, list[float]]]:
        """``{context: (letters, cumulative probabilities)}`` without ``exclude``d letters."""
        cached = self._tables.get(exclude)
        if cached is not None:
            return cached
        size = len(self.alphabet)
        allowed = [i for i, aa in enumerate(self.alphabet) if aa not in exclude]
        tables: dict[str, tuple[str, list[float]]] = {}
        for o, table in enumerate(self.counts):
            for ctx_index, ctx in enumerate(itertools.product(self.alphabet, repeat=o)):
                base = ctx_index * size
                weights = [table[base + i] for i in allowed]
                total = sum(weights)
                if total == 0:
                    continue
                running = 0
                cum: list[float] = []
                for w in weights:
                    running += w
                    cum.append(running / total)
                cum[-1] = 1.0
                tables["".join(ctx)] = ("".join(self.alphabet[i] for i in allowed), cum)
        self._tables[exclude] = tables
        return tables


def _as_paths(fasta_paths: str | Path | Iterable[str | Path]) -> list[Path]:
    if isinstance(fasta_paths, (str, Path)):
        return [Path(fasta_paths)]
    return [Path(p) for p in fasta_paths]


def train_markov_model(
    fasta_paths: str | Path | Iterable[str | Path],
    *,
    order: int = 2,
    metadata: dict[str, object] | None = None,
) -> MarkovModel:
    """Count (k+1)-mers in one or more FASTA files and return a :class:`MarkovModel`.

    Only the 20 standard amino acids are counted (case-insensitive); a k-mer that
    spans any other character (``X``, ``U``, ``*``, ...) is skipped. Paths may be
    gzip, bzip2 or xz compressed.

    Args:
        fasta_paths: A FASTA path or an iterable of paths.
        order: Context length, 0 to 4. Order 2 (the shipped models) needs about
            8,000 counts; each order multiplies that by 20.
        metadata: Extra provenance to store with the model (source, date, license).

    Raises:
        DecoyError: For a bad ``order`` or input without any standard residue.
    """
    if not isinstance(order, int) or isinstance(order, bool) or not 0 <= order <= _MAX_ORDER:
        raise DecoyError(f"order must be an int from 0 to {_MAX_ORDER}, got {order!r}")
    paths = _as_paths(fasta_paths)
    split = re.compile(f"[^{AMINO_ACIDS}]+")
    counters = [Counter[str]() for _ in range(order + 1)]
    entries = 0
    for path in paths:
        for entry in read_fasta(path):
            entries += 1
            for run in split.split(entry.sequence.upper()):
                for o, counter in enumerate(counters):
                    if len(run) > o:
                        counter.update(run[i : i + o + 1] for i in range(len(run) - o))
    counts = tuple(
        tuple(counter["".join(kmer)] for kmer in itertools.product(AMINO_ACIDS, repeat=o + 1))
        for o, counter in enumerate(counters)
    )
    if sum(counts[0]) == 0:
        raise DecoyError("No standard amino acids found in the training FASTA", hint="Train on protein sequences")
    meta: dict[str, object] = {"sources": [p.name for p in paths], "entries": entries, "residues": sum(counts[0])}
    meta.update(metadata or {})
    return MarkovModel(order=order, counts=counts, metadata=meta)


def load_markov_model(name_or_path: str | Path) -> MarkovModel:
    """Load a shipped model by name (see :data:`MODELS`) or a model saved with :meth:`MarkovModel.save`.

    The shipped models are order 2, trained on the UniProtKB/Swiss-Prot entries of
    the UniProt reference proteomes (human UP000005640, mouse UP000000589, yeast
    UP000002311, E. coli K-12 UP000000625), release 2026_03. ``model.metadata`` has
    the details. UniProt data is CC BY 4.0 (https://www.uniprot.org/help/license).

    Raises:
        DecoyError: For an unknown name or a file that is not a saved model.
    """
    if isinstance(name_or_path, str) and name_or_path in MODELS:
        raw = resources.files("fastatacular").joinpath("data", "markov", f"{name_or_path}.json.gz").read_bytes()
        source = f"shipped model {name_or_path!r}"
    else:
        path = Path(name_or_path)
        if not path.is_file():
            raise DecoyError(
                f"Unknown Markov model {str(name_or_path)!r}",
                hint=f"Use one of {', '.join(MODELS)}, or the path of a model saved with MarkovModel.save()",
            )
        raw = path.read_bytes()
        source = str(path)
    try:
        if raw[:2] == b"\x1f\x8b":
            raw = gzip.decompress(raw)
        data = json.loads(raw)
        if data.get("format") != _MODEL_FORMAT:
            raise ValueError(f"format is {data.get('format')!r}, expected {_MODEL_FORMAT!r}")
        return MarkovModel(
            order=data["order"],
            counts=tuple(tuple(t) for t in data["counts"]),
            metadata=data.get("metadata", {}),
            alphabet=data["alphabet"],
        )
    except DecoyError:
        raise
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as err:
        raise DecoyError(f"Cannot load Markov model from {source}: {err}") from err


# ---------------------------------------------------------------------------
# Sequence generation
# ---------------------------------------------------------------------------


def _seed_key(seed: int | str | bytes | None) -> bytes:
    if seed is None:
        return os.urandom(16)
    if isinstance(seed, bool) or not isinstance(seed, (int, str, bytes)):
        raise DecoyError(f"seed must be an int, str, bytes or None, got {type(seed).__name__}")
    material = seed if isinstance(seed, bytes) else f"{type(seed).__name__}:{seed}".encode()
    return hashlib.blake2b(material, digest_size=16).digest()


def _check_count(name: str, value: object) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise DecoyError(f"{name} must be a non-negative int, got {value!r}")
    return value


class _Decoyer:
    """Validated options plus per-run state (the de Bruijn label map)."""

    def __init__(
        self,
        method: str,
        *,
        seed: int | str | bytes | None,
        keep_residues: str | None,
        keep_nterm: int,
        keep_cterm: int,
        k: int,
        model: str | Path | MarkovModel,
    ) -> None:
        if method not in METHODS:
            raise DecoyError(f"Unknown decoy method {method!r}", hint=f"Use one of {', '.join(METHODS)}")
        if keep_residues is None:
            keep_residues = "KR" if method == "pseudo_reverse" else ""
        if not isinstance(keep_residues, str):
            raise DecoyError(f"keep_residues must be a str, got {type(keep_residues).__name__}")
        self.method = method
        self.keep = frozenset(keep_residues)
        self.nterm = _check_count("keep_nterm", keep_nterm)
        self.cterm = _check_count("keep_cterm", keep_cterm)
        if not isinstance(k, int) or isinstance(k, bool) or k < 1:
            raise DecoyError(f"k must be a positive int, got {k!r}")
        self.k = k
        self.seed = _seed_key(seed)
        self.free_run = re.compile(f"[^{re.escape(keep_residues)}]+") if keep_residues else None
        self.model: MarkovModel | None = None
        self.free = frozenset(AMINO_ACIDS) - self.keep
        if method == "markov":
            self.model = model if isinstance(model, MarkovModel) else load_markov_model(model)
            self.free = frozenset(self.model.alphabet) - self.keep
            self.tables = self.model._sampling_tables(self.keep)
            if self.free and "" not in self.tables:
                raise DecoyError(
                    f"The model gives zero probability to every residue not in keep_residues={keep_residues!r}"
                )
        self.labels: dict[str, str] | None = None

    def rng(self, sequence: str) -> random.Random:
        digest = hashlib.blake2b(sequence.encode("utf-8", "surrogatepass"), digest_size=16, key=self.seed).digest()
        return random.Random(int.from_bytes(digest))

    def __call__(self, seq: str) -> str:
        n = len(seq)
        lo = min(self.nterm, n)
        hi = max(lo, n - self.cterm)
        method = self.method
        if method in ("reverse", "pseudo_reverse"):
            return self._rearrange(seq, lo, hi, lambda s: s[::-1])
        if method == "debruijn":
            return self._debruijn(seq, lo, hi)
        rng = self.rng(seq)
        decoy = seq
        for _ in range(_MAX_TRIES):
            if method == "shuffle":
                decoy = self._rearrange(seq, lo, hi, lambda s: "".join(rng.sample(s, len(s))))
            else:
                decoy = self._markov(seq, lo, hi, rng)
            if decoy != seq:
                break
        return decoy

    def _rearrange(self, seq: str, lo: int, hi: int, fn: Callable[[str], str]) -> str:
        core = seq[lo:hi]
        new = fn(core) if self.free_run is None else self.free_run.sub(lambda m: fn(m.group()), core)
        return seq[:lo] + new + seq[hi:]

    def _markov(self, seq: str, lo: int, hi: int, rng: random.Random) -> str:
        model = self.model
        assert model is not None  # set for method="markov"
        order = model.order
        tables = self.tables
        free = self.free
        standard = frozenset(model.alphabet)
        rand = rng.random
        out = list(seq)
        # Context: the last ``order`` decoy residues, reset by a non-standard residue.
        ctx = ""
        for ch in seq[max(0, lo - order) : lo]:
            ctx = ctx + ch if ch in standard else ""
        for i in range(lo, hi):
            ch = seq[i]
            if ch in free:
                table = tables.get(ctx)
                while table is None:
                    ctx = ctx[1:]
                    table = tables.get(ctx)
                letters, cum = table
                ch = letters[bisect(cum, rand())] if len(letters) > 1 else letters
                out[i] = ch
            if ch in standard:
                ctx = (ctx + ch)[1:] if len(ctx) == order else ctx + ch
            else:
                ctx = ""
        return "".join(out)

    def prepare_debruijn(self, sequences: Iterable[str]) -> None:
        """Relabel the de Bruijn graph of ``sequences`` (Moosa et al. 2020, section 2.1).

        Every distinct (k+1)-mer ending in a replaceable residue is an edge; ``uses``
        counts how many replaceable positions follow it. Edges are relabelled one by
        one in random order: a residue ``a`` is drawn with probability proportional
        to its remaining count ``n(a)`` in the target, then ``n(a)`` drops by the
        edge's uses, so the decoy composition tracks the target's.
        """
        k = self.k
        free = self.free
        uses: Counter[str] = Counter()
        remaining: Counter[str] = Counter()
        for seq in sequences:
            n = len(seq)
            lo = min(self.nterm, n)
            hi = max(lo, n - self.cterm)
            padded = "-" * k + seq
            for i in range(lo, hi):
                if seq[i] in free:
                    uses[padded[i : i + k + 1]] += 1
            remaining.update(ch for ch in seq[lo:hi] if ch in free)
        rng = random.Random(int.from_bytes(self.seed))
        edges = sorted(uses)
        rng.shuffle(edges)
        letters = sorted(free)
        start = [remaining[a] for a in letters]
        left = list(start)
        labels: dict[str, str] = {}
        for edge in edges:
            weights = [max(w, 0) for w in left]
            if not any(weights):
                weights = start
            j = rng.choices(range(len(letters)), weights=weights)[0]
            labels[edge] = letters[j]
            left[j] -= uses[edge]
        self.labels = labels

    def _debruijn(self, seq: str, lo: int, hi: int) -> str:
        labels = self.labels
        if labels is None:
            self.prepare_debruijn([seq])
            labels = self.labels
            assert labels is not None
        k = self.k
        padded = "-" * k + seq
        free = self.free
        out = list(seq)
        for i in range(lo, hi):
            if seq[i] in free:
                out[i] = labels[padded[i : i + k + 1]]
        return "".join(out)


def make_decoy_sequence(
    sequence: str,
    *,
    method: DecoyMethod,
    seed: int | str | bytes | None = None,
    keep_residues: str | None = None,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    k: int = 2,
    model: str | Path | MarkovModel = "human",
) -> str:
    """Return the decoy of one sequence. Options are those of :func:`make_decoys`."""
    decoyer = _Decoyer(
        method,
        seed=seed,
        keep_residues=keep_residues,
        keep_nterm=keep_nterm,
        keep_cterm=keep_cterm,
        k=k,
        model=model,
    )
    return decoyer(sequence)


def _check_prefix(prefix: str) -> None:
    if not isinstance(prefix, str) or not prefix or any(c.isspace() for c in prefix) or prefix.startswith(">"):
        raise DecoyError(
            f"Invalid decoy prefix {prefix!r}",
            hint="Use a non-empty prefix without whitespace, such as 'DECOY_' or 'rev_'",
        )


def _decoy_entry(entry: SequenceEntry, sequence: str, prefix: str) -> SequenceEntry:
    header, _ = _build_header_line(entry)
    return _build_entry(_parse_header_line(">" + prefix + header[1:], 0), [sequence], 0)


def make_decoys(
    entries: Iterable[SequenceEntry],
    *,
    method: DecoyMethod,
    prefix: str = DEFAULT_PREFIX,
    seed: int | str | bytes | None = None,
    keep_residues: str | None = None,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    k: int = 2,
    model: str | Path | MarkovModel = "human",
) -> Iterator[SequenceEntry]:
    """Yield one decoy entry per target entry.

    The decoy's header is the target's header with ``prefix`` in front of the
    identifier (``>sp|P12345|X_HUMAN ...`` becomes ``>DECOY_sp|P12345|X_HUMAN ...``),
    re-parsed, so ``identifier`` and ``prefix`` carry the decoy prefix while
    ``accession`` and the description fields are the target's.

    Args:
        entries: Target entries, e.g. from :func:`~fastatacular.read_fasta` or a
            :class:`~fastatacular.FastaReader`. Consumed lazily.
        method: One of :data:`METHODS`.
        prefix: Put in front of each decoy identifier.
        seed: Makes the output reproducible: the same seed, options and sequence
            always give the same decoy, whatever else is in the database.
            ``None`` draws a fresh random seed.
        keep_residues: Residues that stay at their positions (case-sensitive).
            Defaults to ``"KR"`` for ``pseudo_reverse`` and ``""`` otherwise.
            Replacement residues are never drawn from this set.
        keep_nterm: Number of N-terminal residues kept (e.g. 1 for the initiator Met).
        keep_cterm: Number of C-terminal residues kept.
        k: de Bruijn k-mer length (``debruijn`` only); each residue is replaced
            according to itself and the k residues before it. Moosa et al. use k=2.
        model: Markov model for ``markov``: a name from :data:`MODELS`, a
            saved-model path, or a :class:`MarkovModel`.

    ``shuffle`` and ``markov`` retry up to 10 times when the decoy equals the target;
    ``reverse``/``pseudo_reverse`` return a palindromic sequence unchanged, and any
    method returns a sequence unchanged when every position is kept.
    Residues outside the 20 standard amino acids (``X``, ``U``, ``*``, lowercase, ...)
    stay in place for ``markov`` and ``debruijn``.

    ``debruijn`` labels the de Bruijn graph of the whole input, so it reads every
    entry before yielding the first decoy, and a decoy depends on the other
    entries (that is how repeats are preserved). The other methods stream, and a
    decoy depends only on its target sequence, the seed and the options.

    Raises:
        DecoyError: For an unknown method or model, or an invalid option. Raised
            when called, before any entry is read.
    """
    _check_prefix(prefix)
    decoyer = _Decoyer(
        method,
        seed=seed,
        keep_residues=keep_residues,
        keep_nterm=keep_nterm,
        keep_cterm=keep_cterm,
        k=k,
        model=model,
    )

    def generate() -> Iterator[SequenceEntry]:
        items: Iterable[SequenceEntry] = entries
        if decoyer.method == "debruijn":
            items = list(entries)
            decoyer.prepare_debruijn(e.sequence for e in items)
        for entry in items:
            yield _decoy_entry(entry, decoyer(entry.sequence), prefix)

    return generate()


def is_decoy(entry: SequenceEntry | str, *, prefix: str = DEFAULT_PREFIX) -> bool:
    """Whether an entry (or identifier) starts with the decoy ``prefix``."""
    identifier = entry if isinstance(entry, str) else entry.identifier
    return identifier.startswith(prefix)


def write_decoy_fasta(
    src: str | Path | IO[str],
    dst: str | Path | IO[str],
    *,
    method: DecoyMethod,
    concatenate: bool = True,
    prefix: str = DEFAULT_PREFIX,
    seed: int | str | bytes | None = None,
    keep_residues: str | None = None,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    k: int = 2,
    model: str | Path | MarkovModel = "human",
    line_width: int = _SEQ_LINE_WIDTH,
) -> int:
    """Read a target FASTA and write its decoys, by default after the targets.

    With ``concatenate=True`` the output is every target entry (unchanged, headers
    verbatim) followed by every decoy; with ``False`` only the decoys. Decoy options
    are those of :func:`make_decoys`. Returns the number of decoys written.

    Raises:
        DecoyError: For invalid options, or when ``src`` already holds entries
            starting with ``prefix`` (it looks like a decoy database already).
        FastaParseError: For unreadable ``src``.
    """
    targets = read_fasta(src)
    already = next((e.identifier for e in targets if is_decoy(e, prefix=prefix)), None)
    if already is not None:
        raise DecoyError(
            f"The input already has decoy entries (e.g. {already!r} starts with {prefix!r})",
            hint="Pass the target-only FASTA, or choose another prefix",
        )
    decoys = list(
        make_decoys(
            targets,
            method=method,
            prefix=prefix,
            seed=seed,
            keep_residues=keep_residues,
            keep_nterm=keep_nterm,
            keep_cterm=keep_cterm,
            k=k,
            model=model,
        )
    )
    write_fasta(itertools.chain(targets, decoys) if concatenate else decoys, dst, line_width=line_width)
    return len(decoys)


__all__ = [
    "AMINO_ACIDS",
    "DEFAULT_PREFIX",
    "METHODS",
    "MODELS",
    "DecoyMethod",
    "MarkovModel",
    "is_decoy",
    "load_markov_model",
    "make_decoy_sequence",
    "make_decoys",
    "train_markov_model",
    "write_decoy_fasta",
]
