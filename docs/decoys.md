# Decoy databases

Target-decoy searching (Elias and Gygi 2007) estimates the false discovery rate of
peptide identifications by searching spectra against the real ("target") proteins plus
the same number of fake ("decoy") proteins that cannot be present in the sample.
`fastatacular.decoys` builds those decoys. It is pure Python with no dependencies, and
every method is reproducible with `seed`.

```python
from fastatacular import write_decoy_fasta

# targets followed by DECOY_ entries, the usual input for a search engine
n = write_decoy_fasta("human.fasta", "human_td.fasta", method="pseudo_reverse", seed=1)
```

```python
from fastatacular import is_decoy, make_decoys, read_fasta

targets = read_fasta("human.fasta")
decoys = list(make_decoys(targets, method="markov", model="human", keep_residues="KR", seed=1))
decoys[0].identifier   # "DECOY_sp|P31946|1433B_HUMAN"
decoys[0].accession    # "P31946": fields other than the identifier are the target's
is_decoy(decoys[0])    # True
```

## Methods

| method | what it does | reference |
|---|---|---|
| `reverse` | reverse the sequence | Elias and Gygi 2007 |
| `pseudo_reverse` | reverse each stretch between cleavage sites; `K`/`R` stay in place (the default `keep_residues="KR"`) | Elias and Gygi 2007 |
| `shuffle` | randomly permute the residues | Elias and Gygi 2007 |
| `debruijn` | give every distinct (k+1)-mer of the target database one random replacement residue, so repeated target sequence is repeated in the decoys | Moosa et al. 2020 |
| `markov` | sample a new sequence from an order-2 Markov chain trained on a proteome | - |

- **reverse / pseudo_reverse** keep the composition exactly, and pseudo-reversal also
  keeps tryptic peptide lengths and masses. They are deterministic, so `seed` has no
  effect.
- **shuffle** keeps the composition exactly. When a shuffle happens to give back the
  target it is redrawn, up to 10 times.
- **debruijn** is the repeat-preserving method of Moosa et al.:
  1. Each protein is padded with `k` gap symbols. The residue at position *i* is
     replaced by a label chosen for the (k+1)-mer ending at *i*.
  2. Labels are drawn in proportion to the target amino-acid counts still unused.
     Each draw removes that (k+1)-mer's number of occurrences from the pool, so the
     decoy composition closely tracks the target's.
  3. `k=2` is the paper's value.

  A sequence shared by two targets is shared by their decoys (all but its first `k`
  residues). Decoy peptides are therefore about as redundant as target peptides, which
  is what the method is for. Because labels depend on the whole database, `make_decoys`
  reads all entries before yielding the first decoy.
- **markov** samples residue by residue. The next residue depends on the previous 2
  decoy residues, backing off to a shorter context when a context was never seen in
  training. Composition and local dipeptide/tripeptide statistics follow the training
  proteome rather than each protein. Every decoy depends only on its own target, so the
  output streams.

## Options (every method)

- `keep_residues`: residues that stay at their positions, e.g. `"KR"`, `"KRP"` or `"C"`.
  Replacement residues are never drawn from this set, so `keep_residues="KR"` keeps
  every tryptic cleavage site of every method exactly where it was. The default is
  `"KR"` for `pseudo_reverse` and `""` otherwise.
- `keep_nterm`, `keep_cterm`: the number of terminal residues left unchanged. For
  example, `keep_nterm=1` keeps the initiator methionine.
- `seed`: an `int`, `str` or `bytes`. The same seed, method and options always give the
  same decoys, on any platform and Python version (it uses `random.Random` seeded from a
  BLAKE2 hash of the seed and the target sequence). `seed=None` draws a fresh random
  seed.
- `prefix`: prepended to the identifier (default `DECOY_`), e.g.
  `>DECOY_sp|P12345|EX_HUMAN ...`. All other header text is kept, so tools that group by
  accession still see the target's accession on the decoy entry. The prefix must be
  non-empty and contain no whitespace.
- `k` (debruijn only): the k-mer length, default 2.
- `model` (markov only): `"human"`, `"mouse"`, `"yeast"`, `"ecoli"`, a path to a saved
  model, or a `MarkovModel`.

Residues outside the 20 standard amino acids (`X`, `U`, `B`, lowercase and so on) stay at
their positions for `debruijn` and `markov`. `reverse`, `pseudo_reverse` and `shuffle`
move them with the rest of the sequence.

`write_decoy_fasta(src, dst, *, method, concatenate=True, ...)` writes the targets
followed by the decoys (or, with `concatenate=False`, the decoys only) and returns the
number of decoys. It raises `DecoyError` if `src` already contains entries with the
prefix. `DecoyError` is a `FastaError`, and every bad option raises it when
`make_decoys` is called, before any entry is read.

## Markov models

The shipped models are order 2. They were trained on the reviewed (Swiss-Prot)
reference proteomes of each organism, UniProt release 2026_03:

| name | proteome | entries |
|---|---|---|
| `human` | UP000005640, *Homo sapiens* | 20,416 |
| `mouse` | UP000000589, *Mus musculus* | 17,277 |
| `yeast` | UP000002311, *Saccharomyces cerevisiae* S288C | 6,067 |
| `ecoli` | UP000000625, *Escherichia coli* K-12 | 4,403 |

The UniProt data are licensed CC BY 4.0 (The UniProt Consortium, UniProt: the Universal
Protein Knowledgebase in 2025, Nucleic Acids Res 53:D609-D617). Each model's metadata
carries the attribution; see `src/fastatacular/data/markov/NOTICE.md`.
`scripts/train_markov_models.py` rebuilds them.

To train your own model:

```python
from fastatacular import load_markov_model, train_markov_model

model = train_markov_model(["proteome1.fasta.gz", "proteome2.fasta"], order=3)
model.save("my_model.json.gz")
model = load_markov_model("my_model.json.gz")
```

## Quality on the human proteome

These numbers come from `scripts/decoy_quality.py` on 20,416 reviewed human proteins
with `seed=1`, on one core:

1. Composition L1 is the sum of absolute differences between target and decoy
   amino-acid frequencies.
2. "Shared peptides" is the fraction of distinct fully tryptic decoy peptides of 7-40
   residues that are also target peptides (515,663 distinct target peptides).

| method | keep | CPU s | decoy == target | composition L1 | shared peptides |
|---|---|---|---|---|---|
| reverse | - | 0.5 | 0 | 0 | 0.07% |
| pseudo_reverse | KR | 1.3 | 1 | 0 | 0.08% |
| shuffle | - | 4.3 | 0 | 0 | 0.04% |
| shuffle | KR | 5.4 | 1 | 0 | 0.17% |
| debruijn | - | 8.3 | 0 | 0.0020 | 0.03% |
| debruijn | KR | 6.5 | 0 | 0.0010 | 0.03% |
| markov | - | 6.1 | 0 | 0.0038 | 0.04% |
| markov | KR | 5.2 | 0 | 0.0041 | 0.05% |

Every decoy has its target's length. The one identical decoy with `KR` kept is
RS32_HUMAN: every stretch between its K/R residues is a single residue, so no method can
change it without moving a cleavage site.

## References

- Elias JE, Gygi SP (2007). Target-decoy search strategy for increased confidence in
  large-scale protein identifications by mass spectrometry. *Nat Methods* 4:207-214.
  doi:10.1038/nmeth1019
- Moosa JM, Guan S, Moran MF, Ma B (2020). Repeat-preserving decoy database for false
  discovery rate estimation in peptide identification. *J Proteome Res*
  19(3):1029-1036. doi:10.1021/acs.jproteome.9b00555
- The UniProt Consortium (2025). UniProt: the Universal Protein Knowledgebase in 2025.
  *Nucleic Acids Res* 53:D609-D617. doi:10.1093/nar/gkae1010
