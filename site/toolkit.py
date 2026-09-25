"""FASTA toolkit: the Python side of the web app.

Runs inside Pyodide in a Web Worker (see worker.js), and imports in plain CPython too, so
scripts/site_smoke.py can compute the expected numbers with the same code paths. Every
public function takes and returns plain JSON-able values; state lives in ``STATE``.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import re
import statistics
import time
from collections import Counter
from pathlib import Path

import peptacular as pt

from fastatacular import (
    RECORD_KEYS,
    SequenceEntry,
    is_decoy,
    load_markov_model,
    make_decoys,
    read_fasta,
    to_records,
    train_markov_model,
    write_fasta,
)

STANDARD = "ACDEFGHIKLMNPQRSTVWY"
AMBIGUOUS = "BJOUXZ"
# Identifier prefixes that mark decoys in databases from common tools.
COMMON_DECOY_PREFIXES = ("DECOY_", "decoy_", "REV_", "rev_", "XXX_", "SHUFFLED_", "shuffled_", "RANDOM_", "random_")
WORK = Path("/tmp/fasta-toolkit")

STATE: dict = {
    "entries": [],  # working set (targets after clean-up)
    "loaded": [],  # as loaded, for "reset clean-up"
    "files": [],
    "dropped": [],  # [(identifier, reason, detail)] from the last clean-up
    "decoy": None,  # {"targets", "decoys", "output", "params"}
    "model": None,  # trained MarkovModel
    "qc_targets": None,  # (key, digested targets) reused while targets and digest settings are unchanged
}


# ---------------------------------------------------------------------------
# progress
# ---------------------------------------------------------------------------

_progress_fn = None
_last_progress = 0.0


def set_progress(fn) -> None:
    global _progress_fn
    _progress_fn = fn


def progress(stage: str, done: int, total: int, force: bool = False) -> None:
    """Report progress, at most ~10 times a second."""
    global _last_progress
    if _progress_fn is None:
        return
    now = time.monotonic()
    if force or done >= total or now - _last_progress > 0.1:
        _last_progress = now
        _progress_fn(stage, done, total)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _histogram(values: list[float], bins: int = 40, lo: float | None = None, hi: float | None = None) -> dict:
    """Equal-width bins; values above ``hi`` go into the last bin (flagged ``overflow``)."""
    if not values:
        return {"edges": [], "counts": [], "overflow": 0}
    lo = min(values) if lo is None else lo
    hi = max(values) if hi is None else hi
    if hi <= lo:
        hi = lo + 1
    width = (hi - lo) / bins
    counts = [0] * bins
    overflow = 0
    for v in values:
        if v > hi:
            overflow += 1
            counts[-1] += 1
            continue
        i = int((v - lo) / width)
        counts[min(max(i, 0), bins - 1)] += 1
    return {"edges": [lo + i * width for i in range(bins + 1)], "counts": counts, "overflow": overflow}


def _int_histogram(values: list[int], lo: int, hi: int) -> dict:
    counts = Counter(values)
    xs = list(range(lo, hi + 1))
    return {"x": xs, "counts": [counts.get(x, 0) for x in xs]}


def _percentile(sorted_vals: list[int], q: float) -> int:
    if not sorted_vals:
        return 0
    return sorted_vals[min(len(sorted_vals) - 1, int(q * (len(sorted_vals) - 1)))]


def _organism(e: SequenceEntry) -> str:
    if e.os_name:
        return f"{e.os_name} (OX={e.ncbi_tax_id})" if e.ncbi_tax_id is not None else e.os_name
    if e.ncbi_tax_id is not None:
        return f"OX={e.ncbi_tax_id}"
    return "(no OS=/OX=)"


def decoy_prefixes(entries: list[SequenceEntry], prefix: str = "DECOY_") -> dict[str, int]:
    """Count entries whose identifier starts with ``prefix`` or a common decoy prefix."""
    found: Counter[str] = Counter()
    prefixes = (prefix, *[p for p in COMMON_DECOY_PREFIXES if p != prefix])
    for e in entries:
        for p in prefixes:
            if is_decoy(e, prefix=p):
                found[p] += 1
                break
    return dict(found)


def _split_decoys(entries: list[SequenceEntry], prefixes: list[str]) -> tuple[list[SequenceEntry], list[SequenceEntry]]:
    targets, decoys = [], []
    for e in entries:
        (decoys if any(is_decoy(e, prefix=p) for p in prefixes) else targets).append(e)
    return targets, decoys


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------


def load_files(paths: list[str], names: list[str], replace: bool = True) -> str:
    """Read FASTA files (plain, .gz, .bz2 or .xz: detected from the bytes) and merge them."""
    if replace:
        STATE["entries"], STATE["files"], STATE["decoy"], STATE["dropped"] = [], [], None, []
    merged = list(STATE["loaded"]) if not replace else []
    for i, (path, name) in enumerate(zip(paths, names, strict=True)):
        progress(f"Parsing {name}", i, len(paths), force=True)
        t0 = time.perf_counter()
        entries = read_fasta(path)
        size = Path(path).stat().st_size
        STATE["files"].append(
            {"name": name, "entries": len(entries), "bytes": size, "seconds": round(time.perf_counter() - t0, 3)}
        )
        merged.extend(entries)
        Path(path).unlink(missing_ok=True)
    STATE["loaded"] = merged
    STATE["entries"] = list(merged)
    STATE["decoy"] = None
    progress("Parsed", len(paths), len(paths), force=True)
    return json.dumps(summary())


def summary() -> dict:
    entries = STATE["entries"]
    return {
        "files": STATE["files"],
        "loaded": len(STATE["loaded"]),
        "entries": len(entries),
        "residues": sum(len(e.sequence) for e in entries),
        "decoys_detected": decoy_prefixes(entries),
        "preview": [e.raw_header or e.identifier for e in entries[:5]],
    }


def clear() -> str:
    STATE.update(entries=[], loaded=[], files=[], dropped=[], decoy=None, model=None)
    return json.dumps(summary())


# ---------------------------------------------------------------------------
# clean-up
# ---------------------------------------------------------------------------


def cleanup(opts_json: str) -> str:
    """Filter and de-duplicate the loaded entries. Always starts from the loaded set."""
    o = json.loads(opts_json)
    entries = STATE["loaded"]
    dropped: list[tuple[str, str, str]] = []
    reasons: Counter[str] = Counter()

    header_re = None
    if o.get("header_regex"):
        try:
            header_re = re.compile(o["header_regex"], 0 if o.get("regex_case") else re.IGNORECASE)
        except re.error as err:
            raise ValueError(f"Bad header regex: {err}") from None
    header_keep = o.get("regex_mode", "include") == "include"

    org_text = (o.get("organism") or "").strip().lower()
    taxa = {int(t) for t in re.split(r"[,\s]+", o.get("taxa") or "") if t.strip().isdigit()}
    org_keep = o.get("organism_mode", "include") == "include"
    org_active = bool(org_text or taxa)

    min_len = int(o.get("min_len") or 0)
    max_len = int(o.get("max_len") or 0)

    alphabet = None
    if o.get("invalid"):
        alphabet = frozenset((o.get("alphabet") or STANDARD).upper() + (AMBIGUOUS if o.get("allow_ambiguous") else ""))

    kept: list[SequenceEntry] = []
    n = len(entries)
    for i, e in enumerate(entries):
        if i % 2000 == 0:
            progress("Filtering", i, n)
        header = e.raw_header or e.identifier
        if header_re is not None and bool(header_re.search(header)) != header_keep:
            reason = "header regex"
            detail = "matched" if not header_keep else "no match"
        elif org_active and _org_match(e, org_text, taxa) != org_keep:
            reason, detail = "organism", _organism(e)
        elif min_len and len(e.sequence) < min_len:
            reason, detail = "too short", str(len(e.sequence))
        elif max_len and len(e.sequence) > max_len:
            reason, detail = "too long", str(len(e.sequence))
        elif alphabet is not None and not alphabet.issuperset(e.sequence.upper()):
            bad = "".join(sorted(set(e.sequence.upper()) - alphabet))
            reason, detail = "invalid residues", bad
        else:
            kept.append(e)
            continue
        reasons[reason] += 1
        dropped.append((e.identifier, reason, detail))

    if o.get("dedup_id"):
        seen: set[str] = set()
        out = []
        for e in kept:
            if e.identifier in seen:
                reasons["duplicate identifier"] += 1
                dropped.append((e.identifier, "duplicate identifier", "later copy dropped"))
            else:
                seen.add(e.identifier)
                out.append(e)
        kept = out

    if o.get("dedup_seq"):
        il = bool(o.get("il_equivalent"))
        first: dict[str, tuple[str, str]] = {}  # key -> (identifier, sequence) of the kept copy
        out = []
        for e in kept:
            seq = e.sequence.upper()
            key = seq.replace("I", "L") if il else seq
            if key in first:
                kept_id, kept_seq = first[key]
                reason = "duplicate sequence" if seq == kept_seq else "duplicate sequence (I=L)"
                reasons[reason] += 1
                dropped.append((e.identifier, reason, f"same as {kept_id}"))
            else:
                first[key] = (e.identifier, seq)
                out.append(e)
        kept = out

    progress("Filtering", n, n, force=True)
    STATE["entries"] = kept
    STATE["dropped"] = dropped
    STATE["decoy"] = None
    return json.dumps(
        {
            "before": n,
            "after": len(kept),
            "dropped": len(dropped),
            "reasons": dict(reasons),
            "rows": [list(r) for r in dropped[:500]],
            "summary": summary(),
        }
    )


def _org_match(e: SequenceEntry, text: str, taxa: set[int]) -> bool:
    if taxa and e.ncbi_tax_id in taxa:
        return True
    return bool(text) and text in (e.os_name or "").lower()


def reset_cleanup() -> str:
    STATE["entries"] = list(STATE["loaded"])
    STATE["dropped"] = []
    STATE["decoy"] = None
    return json.dumps(summary())


# ---------------------------------------------------------------------------
# digestion (peptacular)
# ---------------------------------------------------------------------------

_RESIDUE_MASS: dict[str, float] = {}
_WATER = 0.0


def _mass_table() -> None:
    """Residue masses from peptacular: mass(XX) - mass(X), and water = mass(X) - residue."""
    global _WATER
    if _RESIDUE_MASS:
        return
    for aa in STANDARD:
        _RESIDUE_MASS[aa] = pt.mass(aa + aa) - pt.mass(aa)
    _WATER = pt.mass("G") - _RESIDUE_MASS["G"]


def peptide_mass(peptide: str) -> float | None:
    """Neutral monoisotopic mass; peptacular itself for residues outside the 20 standard."""
    _mass_table()
    try:
        return _WATER + sum(map(_RESIDUE_MASS.__getitem__, peptide))
    except KeyError:
        try:
            return float(pt.mass(peptide))
        except Exception:  # noqa: BLE001 - ambiguous residues (X, B, Z) have no mass
            return None


_FAST_OK = re.compile(r"[A-Z]+")


def _enzyme_pattern(enzyme: str) -> re.Pattern[str] | None:
    """The protease's cleavage regex from tacular; None means cleave everywhere (unspecific)."""
    from tacular import PROTEASE_LOOKUP

    info = PROTEASE_LOOKUP.get(enzyme)
    if info is None:
        raise ValueError(f"Unknown enzyme {enzyme!r}")
    pattern = info.pattern
    return None if pattern.pattern in ("()", "") else pattern


def fast_digest(seq: str, pattern: re.Pattern[str] | None, missed: int, min_len: int, max_len: int) -> list[str]:
    """Peptide strings for one sequence: the same set as ``pt.digest`` (checked for every enzyme by
    scripts/site_smoke.py), without building a ProForma annotation and spans per protein."""
    n = len(seq)
    if pattern is None:  # unspecific: every substring in the length range (missed cleavages do not apply)
        return [seq[a : a + ln] for a in range(n) for ln in range(max(min_len, 1), min(max_len, n - a) + 1)]
    found = [m.start() for m in pattern.finditer(seq)]
    sites = sorted({0, n, *found}) if found else [0, n]
    out = []
    k = len(sites)
    for i in range(k - 1):
        a = sites[i]
        for j in range(i + 1, min(i + 2 + missed, k)):
            length = sites[j] - a
            if length > max_len:
                break
            if length >= min_len:
                out.append(seq[a : sites[j]])
    return out


def digest_entries(
    entries: list[SequenceEntry], enzyme: str, missed: int, min_len: int, max_len: int, stage: str
) -> tuple[set[str], int, int]:
    """Digest every entry; return (distinct peptides, total peptides with multiplicity, failures).

    Plain upper-case sequences use ``fast_digest``; anything else goes through ``pt.digest``,
    which parses ProForma and reports sequences it cannot read as failures."""
    pattern = _enzyme_pattern(enzyme)
    distinct: set[str] = set()
    total = 0
    failed = 0
    n = len(entries)
    for i, e in enumerate(entries):
        if i % 500 == 0:
            progress(stage, i, n)
        seq = e.sequence
        if _FAST_OK.fullmatch(seq):
            peps = fast_digest(seq, pattern, missed, min_len, max_len)
        else:
            try:
                peps = [p for p, _ in pt.digest(seq, enzyme, missed_cleavages=missed, min_len=min_len, max_len=max_len)]
            except Exception:  # noqa: BLE001 - a sequence peptacular cannot parse
                failed += 1
                continue
        total += len(peps)
        distinct.update(peps)
    progress(stage, n, n, force=True)
    return distinct, total, failed


def enzymes() -> str:
    from tacular import PROTEASE_LOOKUP

    return json.dumps(sorted(PROTEASE_LOOKUP.keys()))


# ---------------------------------------------------------------------------
# stats
# ---------------------------------------------------------------------------


def stats(opts_json: str) -> str:
    o = json.loads(opts_json)
    entries = STATE["entries"]
    lengths = sorted(len(e.sequence) for e in entries)
    comp: Counter[str] = Counter()
    n = len(entries)
    for i, e in enumerate(entries):
        if i % 2000 == 0:
            progress("Counting residues", i, n)
        comp.update(e.sequence.upper())
    total = sum(comp.values()) or 1
    std = [{"aa": aa, "count": comp.get(aa, 0), "fraction": comp.get(aa, 0) / total} for aa in STANDARD]
    other = sorted(((k, v) for k, v in comp.items() if k not in STANDARD), key=lambda kv: -kv[1])
    orgs = Counter(_organism(e) for e in entries)
    ids = Counter(e.identifier for e in entries)
    seqs = Counter(e.sequence for e in entries)
    p99 = _percentile(lengths, 0.99)
    out = {
        "entries": n,
        "residues": sum(lengths),
        "length": {
            "min": lengths[0] if lengths else 0,
            "median": statistics.median(lengths) if lengths else 0,
            "mean": round(statistics.fmean(lengths), 1) if lengths else 0,
            "max": lengths[-1] if lengths else 0,
            "hist": _histogram(lengths, 40, 0, max(p99, 1)),
        },
        "composition": std,
        "other_residues": [{"aa": k, "count": v} for k, v in other],
        "entries_with_nonstandard": sum(1 for e in entries if not set(e.sequence.upper()) <= set(STANDARD)),
        "organisms": [{"name": k, "count": v} for k, v in orgs.most_common(20)],
        "organism_count": len(orgs),
        "duplicate_ids": sum(c - 1 for c in ids.values() if c > 1),
        "duplicate_sequences": sum(c - 1 for c in seqs.values() if c > 1),
        "decoys_detected": decoy_prefixes(entries),
        "peptides": None,
    }
    progress("Counting residues", n, n, force=True)
    if o.get("digest"):
        peps, total, failed = digest_entries(
            entries, o["enzyme"], int(o["missed"]), int(o["min_len"]), int(o["max_len"]), "Digesting"
        )
        out["peptides"] = {"total": total, "distinct": len(peps), "failed": failed}
    return json.dumps(out)


# ---------------------------------------------------------------------------
# decoys
# ---------------------------------------------------------------------------


def _parse_seed(seed: str):
    seed = (seed or "").strip()
    if not seed:
        return None
    return int(seed) if re.fullmatch(r"-?\d+", seed) else seed


def train_model(order: int) -> str:
    """Train a Markov model on the current targets (existing decoys excluded)."""
    targets, _ = _split_decoys(STATE["entries"], list(COMMON_DECOY_PREFIXES))
    if not targets:
        raise ValueError("Load a FASTA file first")
    WORK.mkdir(parents=True, exist_ok=True)
    path = WORK / "train.fasta"
    progress("Training Markov model", 0, 1, force=True)
    write_fasta(targets, path)
    model = train_markov_model(path, order=int(order), metadata={"trained_by": "fastatacular web toolkit"})
    path.unlink()
    STATE["model"] = model
    progress("Training Markov model", 1, 1, force=True)
    return json.dumps({"order": model.order, "metadata": model.metadata})


def load_model_file(path: str) -> str:
    model = load_markov_model(path)
    Path(path).unlink(missing_ok=True)
    STATE["model"] = model
    return json.dumps({"order": model.order, "metadata": model.metadata})


def decoys(opts_json: str) -> str:
    o = json.loads(opts_json)
    prefix = o.get("prefix") or "DECOY_"
    entries = STATE["entries"]
    if not entries:
        raise ValueError("Load a FASTA file first")
    skip = [prefix, *[p for p in COMMON_DECOY_PREFIXES if p != prefix]]
    targets, existing = _split_decoys(entries, skip)
    if not targets:
        raise ValueError("Every entry already looks like a decoy; nothing to decoy")

    method = o["method"]
    model = o.get("model") or "human"
    if method == "markov" and model in ("trained", "uploaded"):
        if STATE["model"] is None:
            raise ValueError("Train or upload a Markov model first")
        model = STATE["model"]
    keep = o.get("keep_residues")
    kwargs = {
        "method": method,
        "prefix": prefix,
        "seed": _parse_seed(o.get("seed", "")),
        "keep_residues": None if keep in (None, "") else keep,
        "keep_nterm": 1 if o.get("keep_met") else 0,
        "keep_cterm": int(o.get("keep_cterm") or 0),
        "k": int(o.get("k") or 2),
        "model": model,
    }
    t0 = time.perf_counter()
    n = len(targets)
    progress(f"Decoys ({method})", 0, n, force=True)
    out: list[SequenceEntry] = []
    for i, d in enumerate(make_decoys(targets, **kwargs)):
        out.append(d)
        if i % 500 == 0:
            progress(f"Decoys ({method})", i, n)
    progress(f"Decoys ({method})", n, n, force=True)
    seconds = time.perf_counter() - t0

    concatenate = bool(o.get("concatenate", True))
    STATE["decoy"] = {
        "targets": targets,
        "decoys": out,
        "output": targets + out if concatenate else out,
        "params": {
            **{k: v for k, v in kwargs.items() if k != "model"},
            "model": model if isinstance(model, str) else "custom",
            "concatenate": concatenate,
        },
    }
    identical = sum(1 for t, d in zip(targets, out, strict=True) if t.sequence == d.sequence)
    tc = Counter("".join(t.sequence for t in targets))
    dc = Counter("".join(d.sequence for d in out))
    tt, dt = sum(tc.values()) or 1, sum(dc.values()) or 1
    l1 = sum(abs(tc[a] / tt - dc[a] / dt) for a in set(tc) | set(dc))
    return json.dumps(
        {
            "targets": len(targets),
            "decoys": len(out),
            "output_entries": len(STATE["decoy"]["output"]),
            "existing_decoys_skipped": len(existing),
            "existing_prefixes": decoy_prefixes(existing, prefix),
            "identical_to_target": identical,
            "composition_l1": round(l1, 5),
            "seconds": round(seconds, 2),
            "examples": [
                {"target": t.sequence[:60], "decoy": d.sequence[:60], "header": d.raw_header}
                for t, d in list(zip(targets, out, strict=True))[:3]
            ],
            "params": STATE["decoy"]["params"],
        }
    )


# ---------------------------------------------------------------------------
# decoy QC
# ---------------------------------------------------------------------------


def _masses(peps: set[str], stage: str) -> tuple[list[float], int]:
    """Monoisotopic masses of ``peps``; the count of peptides without one (X, B, Z ...)."""
    _mass_table()
    residue, water = _RESIDUE_MASS.__getitem__, _WATER
    out: list[float] = []
    bad = 0
    n = len(peps)
    for i, p in enumerate(peps):
        if i % 50000 == 0:
            progress(stage, i, n)
        try:
            out.append(water + sum(map(residue, p)))
        except KeyError:
            m = peptide_mass(p)
            if m is None:
                bad += 1
            else:
                out.append(m)
    return out, bad


def _digest_side(entries: list[SequenceEntry], settings: tuple, stage: str) -> dict:
    peps, total, failed = digest_entries(entries, *settings, stage)
    masses, unmassed = _masses(peps, stage.replace("Digesting", "Masses of"))
    return {"set": peps, "total": total, "failed": failed, "masses": masses, "unmassed": unmassed}


def decoy_qc(opts_json: str) -> str:
    o = json.loads(opts_json)
    d = STATE["decoy"]
    if d is None:
        raise ValueError("Make decoys first")
    enzyme, missed, lo, hi = o["enzyme"], int(o["missed"]), int(o["min_len"]), int(o["max_len"])
    settings = (enzyme, missed, lo, hi)
    t0 = time.perf_counter()
    # The target side only depends on the targets and the digest settings, so trying another
    # decoy method re-digests only the decoys.
    key = (hash(tuple(e.sequence for e in d["targets"])), len(d["targets"]), settings)
    cached = STATE.get("qc_targets")
    target_cached = cached is not None and cached[0] == key
    if target_cached:
        t = cached[1]
    else:
        t = _digest_side(d["targets"], settings, "Digesting targets")
        t["il"] = {p.replace("I", "L") for p in t["set"]}
        STATE["qc_targets"] = (key, t)
    dd = _digest_side(d["decoys"], settings, "Digesting decoys")
    tset, dset = t["set"], dd["set"]
    shared = tset & dset
    tset_il = t["il"]
    shared_il = sum(1 for p in dset if p.replace("I", "L") in tset_il)
    tm, tbad, dm, dbad = t["masses"], t["unmassed"], dd["masses"], dd["unmassed"]
    tpep_total, tfail, dpep_total, dfail = t["total"], t["failed"], dd["total"], dd["failed"]
    mlo = min((min(x) for x in (tm, dm) if x), default=0.0)
    mhi = max((max(x) for x in (tm, dm) if x), default=1.0)
    mlo, mhi = 100 * (mlo // 100), 100 * (mhi // 100 + 1)
    bins = max(1, min(60, int((mhi - mlo) / 100)))
    return json.dumps(
        {
            "params": {"enzyme": enzyme, "missed": missed, "min_len": lo, "max_len": hi},
            "target": {
                "proteins": len(d["targets"]),
                "peptides": tpep_total,
                "distinct": len(tset),
                "failed": tfail,
                "unmassed": tbad,
            },
            "decoy": {
                "proteins": len(d["decoys"]),
                "peptides": dpep_total,
                "distinct": len(dset),
                "failed": dfail,
                "unmassed": dbad,
            },
            "shared": len(shared),
            "shared_fraction": len(shared) / len(dset) if dset else 0.0,
            "shared_il": shared_il,
            "shared_il_fraction": shared_il / len(dset) if dset else 0.0,
            "balance": len(dset) / len(tset) if tset else 0.0,
            "shared_examples": sorted(shared)[:20],
            "length": {
                "target": _int_histogram([len(p) for p in tset], lo, hi),
                "decoy": _int_histogram([len(p) for p in dset], lo, hi),
            },
            "mass": {"target": _histogram(tm, bins, mlo, mhi), "decoy": _histogram(dm, bins, mlo, mhi)},
            "target_cached": target_cached,
            "seconds": round(time.perf_counter() - t0, 2),
        }
    )


# ---------------------------------------------------------------------------
# export
# ---------------------------------------------------------------------------


def _entries_for(which: str) -> list[SequenceEntry]:
    if which == "working":
        return STATE["entries"]
    if STATE["decoy"] is None:
        raise ValueError("Make decoys first")
    return STATE["decoy"]["output"]


def _peff(entries: list[SequenceEntry], fallback_prefix: str, decoy_prefix: str | None = None) -> str:
    """PEFF text; databases whose prefix starts with ``decoy_prefix`` or a common decoy prefix get Decoy=true."""
    from pefftacular import DatabaseHeader, FileHeader, PeffError, write_peff
    from pefftacular import SequenceEntry as PeffEntry

    converted = []
    for e in entries:
        header = e.raw_header or e.identifier
        try:
            pe = PeffEntry.from_fasta(header, e.sequence)
        except PeffError:
            pe = PeffEntry.from_fasta(header, e.sequence, prefix=fallback_prefix)
        converted.append(pe)
    per_prefix = Counter(p.prefix for p in converted)
    markers = (*COMMON_DECOY_PREFIXES, decoy_prefix) if decoy_prefix else COMMON_DECOY_PREFIXES
    decoy_prefixes_ = [p for p in per_prefix if p.startswith(markers)]
    sources = tuple(f["name"] for f in STATE["files"]) or ("FASTA file",)
    header = FileHeader(
        peff_version="1.0",
        general_comments=("Converted from FASTA by the fastatacular web toolkit",),
        databases=tuple(
            DatabaseHeader(
                prefix=p,
                db_name=p,
                db_version="unknown",
                db_sources=sources,
                number_of_entries=c,
                sequence_type="AA",
                decoy=True if p in decoy_prefixes_ else None,
            )
            for p, c in per_prefix.items()
        ),
    )
    buf = io.StringIO()
    write_peff(header, converted, buf)
    return buf.getvalue()


def export(opts_json: str) -> bytes:
    o = json.loads(opts_json)
    fmt = o.get("format", "fasta")
    which = o.get("which", "working")
    progress("Exporting", 0, 1, force=True)
    if which == "dropped":
        buf = io.StringIO()
        w = csv.writer(buf)
        w.writerow(["identifier", "reason", "detail"])
        w.writerows(STATE["dropped"])
        text = buf.getvalue()
    elif which == "model":
        if STATE["model"] is None:
            raise ValueError("No trained model")
        WORK.mkdir(parents=True, exist_ok=True)
        path = WORK / "model.json.gz"
        STATE["model"].save(path)
        data = path.read_bytes()
        path.unlink()
        progress("Exporting", 1, 1, force=True)
        return data
    else:
        entries = _entries_for(which)
        if fmt == "fasta":
            buf = io.StringIO()
            write_fasta(entries, buf, line_width=int(o.get("line_width", 60)))
            text = buf.getvalue()
        elif fmt == "csv":
            buf = io.StringIO()
            w = csv.DictWriter(buf, fieldnames=list(RECORD_KEYS))
            w.writeheader()
            w.writerows(to_records(entries))
            text = buf.getvalue()
        elif fmt == "peff":
            decoy_prefix = STATE["decoy"]["params"]["prefix"] if which == "decoy" else None
            text = _peff(entries, o.get("peff_prefix") or "gen", decoy_prefix)
        else:
            raise ValueError(f"Unknown format {fmt!r}")
    data = text.encode()
    if o.get("gzip"):
        data = gzip.compress(data, mtime=0)
    progress("Exporting", 1, 1, force=True)
    return data
