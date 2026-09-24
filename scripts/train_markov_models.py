"""Train the shipped Markov decoy models from UniProt reference proteomes.

Usage (downloads the four proteomes, ~17 MB, then writes src/fastatacular/data/markov/):

    uv run python scripts/train_markov_models.py [--download-dir DIR] [--order 2]

Each model is trained on the reviewed (UniProtKB/Swiss-Prot) entries of a UniProt
reference proteome. UniProt data is licensed CC BY 4.0; the attribution is stored in
each model's metadata and in src/fastatacular/data/markov/NOTICE.md.
"""

from __future__ import annotations

import argparse
import datetime as dt
import urllib.parse
import urllib.request
from pathlib import Path

from fastatacular.decoys import train_markov_model

PROTEOMES = {
    "human": ("UP000005640", "Homo sapiens", 9606),
    "mouse": ("UP000000589", "Mus musculus", 10090),
    "yeast": ("UP000002311", "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)", 559292),
    "ecoli": ("UP000000625", "Escherichia coli (strain K12)", 83333),
}
QUERY = "(proteome:{id}) AND (reviewed:true)"
URL = "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query={q}"
ATTRIBUTION = (
    "Trained on UniProtKB/Swiss-Prot data. The UniProt Consortium, UniProt: the Universal "
    "Protein Knowledgebase, https://www.uniprot.org. Licensed under CC BY 4.0 "
    "(https://creativecommons.org/licenses/by/4.0/)."
)
OUT = Path(__file__).resolve().parents[1] / "src" / "fastatacular" / "data" / "markov"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--download-dir", type=Path, default=Path("markov_training"))
    ap.add_argument("--order", type=int, default=2)
    ap.add_argument("--uniprot-release", help="release of already-downloaded files, e.g. 2026_03")
    ap.add_argument("--uniprot-release-date", help="its date, e.g. 02-September-2026")
    args = ap.parse_args()
    args.download_dir.mkdir(parents=True, exist_ok=True)
    for name, (proteome, organism, taxid) in PROTEOMES.items():
        query = QUERY.format(id=proteome)
        path = args.download_dir / f"{name}_{proteome}_reviewed.fasta.gz"
        release, release_date = args.uniprot_release, args.uniprot_release_date
        if not path.exists():
            url = URL.format(q=urllib.parse.quote(query))
            with urllib.request.urlopen(url) as resp:  # noqa: S310 - fixed https URL
                path.write_bytes(resp.read())
                release = resp.headers.get("X-UniProt-Release")
                release_date = resp.headers.get("X-UniProt-Release-Date")
        model = train_markov_model(
            path,
            order=args.order,
            metadata={
                "name": name,
                "organism": organism,
                "ncbi_tax_id": taxid,
                "proteome": proteome,
                "uniprot_query": query,
                "uniprot_release": release,
                "uniprot_release_date": release_date,
                "downloaded": dt.date.fromtimestamp(path.stat().st_mtime).isoformat(),
                "license": "CC BY 4.0",
                "attribution": ATTRIBUTION,
            },
        )
        model.save(OUT / f"{name}.json.gz")
        print(name, model.metadata["entries"], "entries", model.metadata["residues"], "residues")


if __name__ == "__main__":
    main()
