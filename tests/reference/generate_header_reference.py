"""Build ``header_reference.json``: FASTA headers with independently derived expected fields.

How it was produced (2026-09-23, Python 3.12, stdlib only, no fastatacular import):

    python tests/reference/generate_header_reference.py

Three sources, none of which uses fastatacular to decide the expected values:

1. ``spec``: every example header in UniProt's "FASTA headers" help page
   (https://www.uniprot.org/help/fasta-headers, lastModified 2026-09-21), with the
   fields transcribed by hand from the grammar on that page.
2. ``uniprotkb``: real headers from ``rest.uniprot.org/uniprotkb/<acc>.fasta``. The
   expected fields come from the same entry's JSON record (``<acc>.json``), not from
   the header text: primaryAccession, uniProtkbId, the RecName (else first SubName)
   plus " (Fragment)" when flagged, organism scientificName and taxonId, the first
   gene name (else OrderedLocusName, else ORFname), the protein existence level
   and entryAudit.sequenceVersion, following the rules on the help page.
3. ``ncbi``: real headers from NCBI efetch (``rettype=fasta``); the expected
   identifier and description are ``accessionversion`` and ``title`` from esummary.

Rerunning needs network access. The committed JSON is what the tests use.
"""

from __future__ import annotations

import json
import time
import urllib.request
from pathlib import Path

OUT = Path(__file__).with_name("header_reference.json")

UNIPROT_ACCESSIONS = [
    "P69905",  # plain Swiss-Prot
    "P04637",
    "P0DTC2",  # long organism name with digits
    "P0A7Y4",  # strain in parentheses
    "Q15118",  # name with brackets and a comma
    "Q6GZX4",  # isolate in parentheses, gene with a hyphen
    "P31946",  # name with a slash
    "P01308",
    "Q8WZ42",
    "P02769",
    "A0A024R161",  # TrEMBL
    "A0A0B4J2F0",
    "Q3SA23",  # TrEMBL fragment
    "Q8N2H2",  # TrEMBL, SubName, no gene
    "P04224",  # no gene name
    "Q8I6R7",  # fragment
    "P27748",  # strain with slashes
]
NCBI_ACCESSIONS = ["NP_000537.3", "XP_016880112.1", "AAA59174.1", "WP_000184067.1"]

_NONE = {
    "prefix": None,
    "accession": None,
    "entry_name": None,
    "pname": None,
    "os_name": None,
    "ncbi_tax_id": None,
    "gname": None,
    "pe": None,
    "sv": None,
    "extra": {},
}


def _exp(**kw: object) -> dict[str, object]:
    return {**_NONE, **kw}


SPEC = [
    (
        ">sp|Q8I6R7|ACN2_ACAGO Acanthoscurrin-2 (Fragment) OS=Acanthoscurria gomesiana OX=115339 GN=acantho2 PE=1 SV=1",
        _exp(
            identifier="sp|Q8I6R7|ACN2_ACAGO",
            prefix="sp",
            accession="Q8I6R7",
            entry_name="ACN2_ACAGO",
            pname="Acanthoscurrin-2 (Fragment)",
            os_name="Acanthoscurria gomesiana",
            ncbi_tax_id=115339,
            gname="acantho2",
            pe=1,
            sv=1,
        ),
    ),
    (
        ">sp|P27748|ACOX_CUPNH Acetoin catabolism protein X OS=Cupriavidus necator (strain ATCC 17699 / H16 / "
        "DSM 428 / Stanier 337) OX=381666 GN=acoX PE=4 SV=2",
        _exp(
            identifier="sp|P27748|ACOX_CUPNH",
            prefix="sp",
            accession="P27748",
            entry_name="ACOX_CUPNH",
            pname="Acetoin catabolism protein X",
            os_name="Cupriavidus necator (strain ATCC 17699 / H16 / DSM 428 / Stanier 337)",
            ncbi_tax_id=381666,
            gname="acoX",
            pe=4,
            sv=2,
        ),
    ),
    (
        ">sp|P04224|HA22_MOUSE H-2 class II histocompatibility antigen, E-K alpha chain OS=Mus musculus OX=10090 "
        "PE=1 SV=1",
        _exp(
            identifier="sp|P04224|HA22_MOUSE",
            prefix="sp",
            accession="P04224",
            entry_name="HA22_MOUSE",
            pname="H-2 class II histocompatibility antigen, E-K alpha chain",
            os_name="Mus musculus",
            ncbi_tax_id=10090,
            pe=1,
            sv=1,
        ),
    ),
    (
        # The double space before OX= is in the help page itself.
        ">tr|Q3SA23|Q3SA23_9HIV1 Protein Nef (Fragment) OS=Human immunodeficiency virus 1  OX=11676 GN=nef PE=3 SV=1",
        _exp(
            identifier="tr|Q3SA23|Q3SA23_9HIV1",
            prefix="tr",
            accession="Q3SA23",
            entry_name="Q3SA23_9HIV1",
            pname="Protein Nef (Fragment)",
            os_name="Human immunodeficiency virus 1",
            ncbi_tax_id=11676,
            gname="nef",
            pe=3,
            sv=1,
        ),
    ),
    (
        ">tr|Q8N2H2|Q8N2H2_HUMAN cDNA FLJ90785 fis, clone THYRO1001457, moderately similar to H.sapiens protein "
        "kinase C mu OS=Homo sapiens OX=9606 PE=2 SV=1",
        _exp(
            identifier="tr|Q8N2H2|Q8N2H2_HUMAN",
            prefix="tr",
            accession="Q8N2H2",
            entry_name="Q8N2H2_HUMAN",
            pname="cDNA FLJ90785 fis, clone THYRO1001457, moderately similar to H.sapiens protein kinase C mu",
            os_name="Homo sapiens",
            ncbi_tax_id=9606,
            pe=2,
            sv=1,
        ),
    ),
    (
        ">sp|Q4R572-2|1433B_MACFA Isoform Short of 14-3-3 protein beta/alpha OS=Macaca fascicularis OX=9541 GN=YWHAB",
        _exp(
            identifier="sp|Q4R572-2|1433B_MACFA",
            prefix="sp",
            accession="Q4R572-2",
            entry_name="1433B_MACFA",
            pname="Isoform Short of 14-3-3 protein beta/alpha",
            os_name="Macaca fascicularis",
            ncbi_tax_id=9541,
            gname="YWHAB",
        ),
    ),
    (
        ">UniRef50_Q9K794 Putative AgrB-like protein n=2 Tax=Bacillus TaxID=1386 RepID=AGRB_BACHD",
        _exp(
            identifier="UniRef50_Q9K794",
            pname="Putative AgrB-like protein",
            extra={"n": "2", "Tax": "Bacillus", "TaxID": "1386", "RepID": "AGRB_BACHD"},
        ),
    ),
    (">UPI0000000005 status=active", _exp(identifier="UPI0000000005", extra={"status": "active"})),
    (
        ">UPI000013B286 Multidrug resistance protein MdtE OS=Escherichia coli (strain K12) OX=83333 GN=mdtE "
        "AC=P37636 SS=EMBL:AAC76538 PC=UP000000625:Chromosome",
        _exp(
            identifier="UPI000013B286",
            pname="Multidrug resistance protein MdtE",
            os_name="Escherichia coli (strain K12)",
            ncbi_tax_id=83333,
            gname="mdtE",
            extra={"AC": "P37636", "SS": "EMBL:AAC76538", "PC": "UP000000625:Chromosome"},
        ),
    ),
    (
        ">UPI000000E135 Oligopeptide transporter subunit OS=Escherichia coli (strain K12) OX=83333 GN=oppF "
        "SS=EMBL:CQR80801 PC=UP000033172:Chromosome I",
        _exp(
            identifier="UPI000000E135",
            pname="Oligopeptide transporter subunit",
            os_name="Escherichia coli (strain K12)",
            ncbi_tax_id=83333,
            gname="oppF",
            extra={"SS": "EMBL:CQR80801", "PC": "UP000033172:Chromosome I"},
        ),
    ),
    (
        ">UPI00000000C1 Calmodulin-1|Calmodulin-2|Calmodulin-3 OS=Homo sapiens OX=9606 GN=CALM1|CALM2|CALM3 "
        "AC=P0DP23|P0DP24|P0DP25 SS=Ensembl:ENSP00000349467|Ensembl:ENSP00000272298|Ensembl:ENSP00000291295|"
        "Ensembl:ENSP00000472141 PC=UP000005640:Chromosome 14|Chromosome 2|Chromosome 19",
        _exp(
            identifier="UPI00000000C1",
            pname="Calmodulin-1|Calmodulin-2|Calmodulin-3",
            os_name="Homo sapiens",
            ncbi_tax_id=9606,
            gname="CALM1|CALM2|CALM3",
            extra={
                "AC": "P0DP23|P0DP24|P0DP25",
                "SS": "Ensembl:ENSP00000349467|Ensembl:ENSP00000272298|Ensembl:ENSP00000291295|Ensembl:ENSP00000472141",
                "PC": "UP000005640:Chromosome 14|Chromosome 2|Chromosome 19",
            },
        ),
    ),
    (
        ">sp|P05067 archived from Release 18.0 01-MAY-1991 SV=3",
        _exp(
            identifier="sp|P05067",
            prefix="sp",
            accession="P05067",
            pname="archived from Release 18.0 01-MAY-1991",
            sv=3,
        ),
    ),
    (
        ">tr|Q55167 archived from Release 17.0 01-JUN-2001 SV=1",
        _exp(
            identifier="tr|Q55167",
            prefix="tr",
            accession="Q55167",
            pname="archived from Release 17.0 01-JUN-2001",
            sv=1,
        ),
    ),
    (
        ">sp|P05067 archived from Release 9.2/51.2 28-NOV-2006 SV=3",
        _exp(
            identifier="sp|P05067",
            prefix="sp",
            accession="P05067",
            pname="archived from Release 9.2/51.2 28-NOV-2006",
            sv=3,
        ),
    ),
    (
        ">tr|A0RTJ8 archived from Release 11.0/36.0 29-MAY-2007 SV=1",
        _exp(
            identifier="tr|A0RTJ8",
            prefix="tr",
            accession="A0RTJ8",
            pname="archived from Release 11.0/36.0 29-MAY-2007",
            sv=1,
        ),
    ),
]

# Legacy NCBI ``db|id|...`` deflines (NCBI C++ Toolkit book, "FASTA identifiers":
# https://ncbi.github.io/cxx-toolkit/pages/ch_demo#ch_demo.T5) and generic headers.
# fastatacular keeps the first two pipe fields as prefix/accession.
OTHER = [
    (
        ">gi|4557757|ref|NP_000240.1| MutL protein homolog 1 [Homo sapiens]",
        _exp(
            identifier="gi|4557757|ref|NP_000240.1|",
            prefix="gi",
            accession="4557757",
            pname="MutL protein homolog 1 [Homo sapiens]",
        ),
    ),
    (
        ">ref|NP_000240.1| MutL protein homolog 1 [Homo sapiens]",
        _exp(
            identifier="ref|NP_000240.1|",
            prefix="ref",
            accession="NP_000240.1",
            pname="MutL protein homolog 1 [Homo sapiens]",
        ),
    ),
    (
        ">gb|AAA59174.1| insulin receptor precursor [Homo sapiens]",
        _exp(
            identifier="gb|AAA59174.1|",
            prefix="gb",
            accession="AAA59174.1",
            pname="insulin receptor precursor [Homo sapiens]",
        ),
    ),
    (">ENSP00000269305.4", _exp(identifier="ENSP00000269305.4")),
    (">contig_12\tassembled from 3 reads", _exp(identifier="contig_12", pname="assembled from 3 reads")),
]

# Spec examples whose structure SequenceEntry cannot represent. Kept as strict xfails.
KNOWN_DIFFERENCES = [
    (
        ">pp562|UPI0000126990 96% UP000000625 AAC76437 OX=83333 OS=Escherichia coli (strain K12) ; "
        "sp|P13001|BIOH_ECOLI Pimeloyl-[acyl-carrier protein] methyl ester esterase GN=bioH PE=1 SV=1",
        _exp(
            identifier="pp562|UPI0000126990",
            prefix="pp562",
            accession="UPI0000126990",
            pname="96% UP000000625 AAC76437",
            os_name="Escherichia coli (strain K12)",
            ncbi_tax_id=83333,
            gname="bioH",
            pe=1,
            sv=1,
        ),
        "Pan-proteome headers embed a second UniProtKB header after ' ; '. It has no field of its own, "
        "so the OS= value runs on to the next KEY= and takes the annotation block with it.",
    ),
]


def _get(url: str) -> str:
    time.sleep(0.5)  # NCBI allows 3 requests/s without an API key
    with urllib.request.urlopen(url, timeout=60) as r:  # noqa: S310 - fixed https URLs
        return r.read().decode("utf-8")


def _uniprot(acc: str) -> dict[str, object]:
    header = _get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta").splitlines()[0]
    d = json.loads(_get(f"https://rest.uniprot.org/uniprotkb/{acc}.json"))
    desc = d["proteinDescription"]
    name_block = desc.get("recommendedName") or desc["submissionNames"][0]
    name = name_block["fullName"]["value"]
    if desc.get("flag", "").startswith("Fragment"):
        name += " (Fragment)"
    gname = None
    for g in d.get("genes", []):
        for field, many in (("geneName", False), ("orderedLocusNames", True), ("orfNames", True)):
            if field in g:
                gname = g[field][0]["value"] if many else g[field]["value"]
                break
        if gname:
            break
    db = "sp" if "Swiss-Prot" in d["entryType"] else "tr"
    return {
        "source": f"uniprotkb:{acc}",
        "header": header,
        "expected": _exp(
            identifier=f"{db}|{d['primaryAccession']}|{d['uniProtkbId']}",
            prefix=db,
            accession=d["primaryAccession"],
            entry_name=d["uniProtkbId"],
            pname=name,
            os_name=d["organism"]["scientificName"],
            ncbi_tax_id=d["organism"]["taxonId"],
            gname=gname,
            pe=int(d["proteinExistence"].split(":")[0]),
            sv=d["entryAudit"]["sequenceVersion"],
        ),
    }


def _ncbi(acc: str) -> dict[str, object]:
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    header = _get(f"{base}/efetch.fcgi?db=protein&id={acc}&rettype=fasta").splitlines()[0]
    res = json.loads(_get(f"{base}/esummary.fcgi?db=protein&id={acc}&retmode=json"))["result"]
    doc = res[res["uids"][0]]
    return {
        "source": f"ncbi:{acc}",
        "header": header,
        "expected": _exp(identifier=doc["accessionversion"], pname=doc["title"]),
    }


def main() -> None:
    cases = [{"source": "uniprot-help", "header": h, "expected": e} for h, e in SPEC]
    cases += [{"source": "ncbi-legacy/generic", "header": h, "expected": e} for h, e in OTHER]
    cases += [_uniprot(a) for a in UNIPROT_ACCESSIONS]
    cases += [_ncbi(a) for a in NCBI_ACCESSIONS]
    known = [{"source": "uniprot-help", "header": h, "expected": e, "reason": r} for h, e, r in KNOWN_DIFFERENCES]
    OUT.write_text(json.dumps({"cases": cases, "known_differences": known}, indent=1, ensure_ascii=False) + "\n")
    print(f"wrote {len(cases)} cases + {len(known)} known differences to {OUT}")


if __name__ == "__main__":
    main()
