---
title: 'fastatacular: Typed FASTA reading, writing, indexing and decoy generation for proteomics in Python'
tags:
  - Python
  - Proteomics
  - Mass Spectrometry
  - FASTA
  - Bioinformatics
authors:
  - name: Patrick T. Garrett
    orcid: 0000-0002-8434-9693
    affiliation: 1
  - given-names: John R.
    surname: Yates
    suffix: III
    orcid: 0000-0001-5267-1672
    corresponding: true
    affiliation: 1
affiliations:
  - name: The Scripps Research Institute, United States
    index: 1
date: 24 September 2026
bibliography: paper.bib
---

# Summary

Protein sequence databases in FASTA format are the search space of almost every
mass-spectrometry proteomics experiment. **fastatacular** is a dependency-free Python
library that reads and writes these files as typed records. It splits UniProt-style
headers into named fields, streams plain and compressed files, gives random access
through samtools-compatible `.fai` indexes, exports records for data-frame libraries,
and builds target-decoy databases with five decoy methods. Entries can be passed
directly to the Peptacular sequence-analysis library [@peptacular].

# Statement of need

A FASTA header is free text. UniProt defines a convention for it,
`db|ACCESSION|ENTRY_NAME name OS=... OX=... GN=... PE=... SV=...` [@uniprot-fasta-help],
and proteomics code routinely re-implements the regular expressions that split it.
Edge cases that such parsers must handle include a protein name that itself contains
`=`, decoy prefixes such as `Reverse_sp|...`, tab separators, byte-order marks, and
Windows line endings. When such files are rewritten, for example after adding decoys or
filtering by taxon, edits can be lost or the header can come back in a form that no
longer parses to the same fields.

Target-decoy searching, the standard way to estimate the false discovery rate of
peptide identifications [@elias-2007], also depends on FASTA handling. Researchers need
decoys that can be reproduced from a seed, keep tryptic cleavage sites in place, and
avoid sharing peptides with the targets. fastatacular targets developers of search
pipelines and analysis scripts who need these steps done predictably in pure Python.

# State of the field

**Biopython** [@cock-2009] reads and writes FASTA through `Bio.SeqIO`, returning
`SeqRecord` objects whose header is kept as an identifier and a description string.
Random access is available through `SeqIO.index` [@biopython-seqio]. It does not split
UniProt fields. **Pyteomics** [@goloborodko-2013; @levitsky-2019] is the closest
proteomics tool. Its `fasta` module parses UniProt, UniRef, UniParc, UniMES, SPD, NCBI and
RefSeq headers, offers indexed readers, and generates reversed, shuffled and fused decoys
[@pyteomics-fasta]. It supports more header dialects than fastatacular.
**pyfaidx** [@shirley-2015] provides Pythonic random access to FASTA files through
samtools-style `.fai` indexes [@danecek-2021] and is mainly used for genomic sequences.
Stand-alone decoy generators such as DecoyPYrat [@wright-2016] build decoy databases from
the command line.

fastatacular does not replace these tools. It combines, in one small package with no
runtime dependencies: frozen, typed entry objects; a writer that re-parses its own output;
compressed and piped input; `.fai` indexing; and two decoy methods that we did not find in
the Python libraries above, namely repeat-preserving de Bruijn decoys [@moosa-2020] and
Markov-chain decoys. Its companion package, pefftacular [@pefftacular], has the same
reader and writer API for the PSI Extended FASTA Format.

# Software design

**Parsing.** `read_fasta` loads a file and `FastaReader` streams it lazily. Paths may be
plain, gzip, bzip2 or xz files. The format is detected from magic bytes rather than the
file name, so pipes, FIFOs and `/dev/stdin` also work. Each record becomes a frozen
`SequenceEntry` with the identifier, sequence, raw header, and parsed fields. The parser
recognises two identifier forms. UniProt-style `prefix|accession|entry_name` fills all
three fields. Other pipe identifiers, such as NCBI's legacy `gi|...|ref|...`, fill the
prefix and accession. Any other identifier is kept whole. In the description, `OS`, `OX`,
`GN`, `PE` and `SV` become typed fields, other `KEY=value` pairs go to `extra`, and the
remaining text is the protein name. A key is recognised only at the start of the
description or after whitespace, so `Protein(EC=2.7.1)` stays in the name. The prefix may
be any text without `|` or whitespace. Decoy and contaminant identifiers
(`Reverse_sp|P1|X`, `DECOY-0-sp|...`) therefore keep their accession. Modern NCBI
headers are read as an identifier and a free-text description. Their bracketed
`[organism=...]` text is not parsed into fields.

**Writing.** `write_fasta` validates every entry before writing anything, so a failed
write leaves no partial file. An unedited entry is written with its original header byte
for byte. An entry edited with `dataclasses.replace` has its header rebuilt from the
fields. The writer raises `FastaWriteError` instead of writing a header that would read
back differently, for example free text that contains a `KEY=` token.

**Random access and tables.** `FastaIndex` scans an uncompressed file once and keeps
each entry's byte range. Entries are then read on demand by identifier or accession.
Indexes can be saved to and loaded from samtools-compatible `.fai` files. `to_records`
returns one plain dictionary per entry for pandas or polars, without depending on
either.

**Decoys.** `make_decoys` and `write_decoy_fasta` implement five methods: reverse,
pseudo-reverse, and shuffle [@elias-2007]; de Bruijn [@moosa-2020]; and an order-2
Markov chain. Markov models for human, mouse, yeast, and *E. coli* ship with the package.
They were trained on UniProtKB/Swiss-Prot reference proteomes [@uniprot-2025], and users
can train their own. Every method accepts `keep_residues` (for example `"KR"` to keep
tryptic sites), `keep_nterm`, and `keep_cterm`. The same seed and options always give the
same decoys, because the random generator is seeded from a hash of the seed and the
target sequence. The repository's `scripts/decoy_quality.py` reports each method's
composition error and the fraction of decoy tryptic peptides that are also target
peptides. `docs/decoys.md` lists these results for the human proteome.

**Errors and testing.** All errors derive from `FastaError`, a `ValueError`. Parse errors
carry the line number, the offending text, and a repair hint. The test suite checks
parsed fields against expected values from independent sources:

- every example header on UniProt's FASTA-header help page;
- real UniProtKB headers, whose expected fields come from the entries' JSON records;
- NCBI headers, checked against NCBI esummary.

Hypothesis property tests [@maciver-2019] cover write-read round trips. Continuous
integration runs on Python 3.12 to 3.14 on Linux, and on macOS and Windows, with the
lowest supported dependency versions and the built wheel.

# Research impact statement

fastatacular is the FASTA reader of the tacular-omics packages. Peptacular's sequence
functions accept any object with a `sequence` attribute, so fastatacular entries can be
digested or weighed without conversion. Peptacular's documentation points users to
fastatacular for FASTA input [@peptacular]. pefftacular's documentation shows how to
convert a UniProt FASTA file read with fastatacular into PEFF [@pefftacular].

<!-- TODO(author): add evidence of research use (published analyses, other groups,
lab pipelines) if any exists. None is recorded in the repository, and JOSS requires
evidence that the software is used for research. -->

# Example usage

The examples below use a two-entry UniProt file (human haemoglobin alpha and beta),
written both plain and gzip-compressed. The output shown is the real output.

```python
from fastatacular import FastaReader

with FastaReader("human.fasta.gz") as reader:   # gzip detected from magic bytes
    for entry in reader:
        print(entry.accession, entry.entry_name, entry.gname,
              entry.ncbi_tax_id, entry.pe, len(entry.sequence))
# P69905 HBA_HUMAN HBA1 9606 1 142
# P68871 HBB_HUMAN HBB 9606 1 147
```

```python
from fastatacular import FastaIndex, is_decoy, make_decoys, read_fasta

targets = read_fasta("human.fasta")
decoys = list(make_decoys(targets, method="debruijn", keep_residues="KR", seed=1))
print(decoys[0].identifier, is_decoy(decoys[0]), decoys[0].sequence[:20])
# DECOY_sp|P69905|HBA_HUMAN True GPLVDYFKFAHKGLLWKTLD

index = FastaIndex("human.fasta", key="accession")
print(index["P68871"].pname)
# Hemoglobin subunit beta
```

```python
import peptacular as pt
from fastatacular import read_fasta

hbb = read_fasta("human.fasta")[1]
peptides = pt.digest(hbb, "trypsin")   # any object with .sequence is accepted
print(len(peptides), [p for p, _ in peptides[:3]])
# 15 ['MVHLTPEEK', 'SAVTALWGK', 'VNVDEVGGEALGR']
```

# Availability

fastatacular is released under the MIT license. It is distributed on PyPI
(<https://pypi.org/project/fastatacular/>), developed on GitHub
(<https://github.com/tacular-omics/fastatacular>), and archived on Zenodo
[@fastatacular]. It requires Python 3.12 or later.

# AI usage disclosure

During the preparation of this work the authors used Anthropic Claude large language models via the Claude Code interface for software-development assistance, including code, tests, and documentation, and for manuscript drafting and editing. The authors reviewed and edited all content and take full responsibility for the software and the publication.

# Acknowledgements

We thank Claire Delahunty, Ph.D., for a careful reading of the manuscript. This work was supported by the U.S. National Institutes of Health (grants R01 HL165168, R01 AG077046, R01 MH100175, and R01 AG075862 to J.R.Y.) and by the Skaggs Graduate School of Chemical and Biological Sciences at The Scripps Research Institute (P.T.G.).

# References
