"""Decoy quality report: composition, target overlap and speed per method.

    uv run python scripts/decoy_quality.py PROTEOME.fasta[.gz] [--seed 1]

For each method (with its default options and with keep_residues="KR") prints:
CPU seconds, decoys identical to their target, total length ratio, residue
composition L1 distance to the target, and the fraction of unique fully tryptic
decoy peptides (cleave after K/R, not before P; 7-40 aa, no missed cleavages)
that are also target peptides.
"""

from __future__ import annotations

import argparse
import re
import time
from collections import Counter

from fastatacular import read_fasta
from fastatacular.decoys import METHODS, make_decoys

_TRYPSIN = re.compile(r"(?<=[KR])(?!P)")


def tryptic(seqs: list[str], lo: int = 7, hi: int = 40) -> set[str]:
    return {p for s in seqs for p in _TRYPSIN.split(s) if lo <= len(p) <= hi}


def composition(seqs: list[str]) -> dict[str, float]:
    counts = Counter("".join(seqs))
    total = sum(counts.values())
    return {aa: c / total for aa, c in counts.items()}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("fasta")
    ap.add_argument("--seed", default=1, type=int)
    args = ap.parse_args()
    targets = read_fasta(args.fasta)
    t_seqs = [e.sequence for e in targets]
    t_peps = tryptic(t_seqs)
    t_comp = composition(t_seqs)
    print(f"{len(targets)} targets, {sum(map(len, t_seqs))} residues, {len(t_peps)} unique tryptic peptides (7-40 aa)")
    print("| method | keep | CPU s | identical | length ratio | composition L1 | decoy peptides in target |")
    print("|---|---|---|---|---|---|---|")
    for method in METHODS:
        for keep in (None, "KR"):
            if method == "pseudo_reverse" and keep == "KR":
                continue  # the default already keeps K/R
            t0 = time.process_time()
            decoys = list(make_decoys(targets, method, seed=args.seed, keep_residues=keep))
            cpu = time.process_time() - t0
            d_seqs = [e.sequence for e in decoys]
            same = sum(d == t for d, t in zip(d_seqs, t_seqs, strict=True))
            ratio = sum(map(len, d_seqs)) / sum(map(len, t_seqs))
            d_comp = composition(d_seqs)
            l1 = sum(abs(d_comp.get(a, 0) - t_comp.get(a, 0)) for a in set(d_comp) | set(t_comp))
            d_peps = tryptic(d_seqs)
            shared = len(d_peps & t_peps) / len(d_peps)
            label = keep if keep is not None else ("KR" if method == "pseudo_reverse" else "-")
            print(f"| {method} | {label} | {cpu:.1f} | {same} | {ratio:.3f} | {l1:.4f} | {shared:.2%} |")


if __name__ == "__main__":
    main()
