# Markov decoy models: data source and license

The files in this directory are order-2 amino-acid Markov models (k-mer counts only,
no sequences) used by `fastatacular.decoys` (`method="markov"`). Rebuild them with `scripts/train_markov_models.py`.

| model | UniProt reference proteome | organism | Swiss-Prot entries | residues |
|---|---|---|---|---|
| `human` | UP000005640 | Homo sapiens (9606) | 20,416 | 11,414,601 |
| `mouse` | UP000000589 | Mus musculus (10090) | 17,277 | 9,825,191 |
| `yeast` | UP000002311 | Saccharomyces cerevisiae S288c (559292) | 6,067 | 2,938,058 |
| `ecoli` | UP000000625 | Escherichia coli K-12 (83333) | 4,403 | 1,354,431 |

Source: UniProtKB/Swiss-Prot (reviewed) entries of each reference proteome, query
`(proteome:<id>) AND (reviewed:true)`, UniProt release 2026_03 (2 September 2026),
downloaded 2026-09-24 from https://rest.uniprot.org.

Attribution: The UniProt Consortium. UniProt: the Universal Protein Knowledgebase,
https://www.uniprot.org. UniProt data is licensed under the Creative Commons
Attribution 4.0 International License (CC BY 4.0,
https://creativecommons.org/licenses/by/4.0/). The models are derived count tables.
