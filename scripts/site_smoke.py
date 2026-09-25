# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "playwright",
#     "fastatacular>=1.1,<2",
#     "peptacular>=5,<6",
#     "tacular>=2,<3",
#     "pefftacular>=1.1,<2",
# ]
# ///
"""Headless end-to-end check of the web app in ``site/``.

Serves ``site/`` with ``python -m http.server`` on a free port, drives it in headless
Chromium, exercises every feature with the example files, and compares the numbers the
page shows with the same computation in CPython (``site/toolkit.py`` imported directly,
against the PyPI releases the page installs). Fails on any console error or page error.
Then times a synthetic 20,000-entry proteome end to end.

    uv run scripts/site_smoke.py                 # PEP 723 deps: playwright + the PyPI packages
    uv run scripts/site_smoke.py --entries 2000  # smaller timing run
    # first time only: uv run --with playwright playwright install chromium

Kept out of ``tests/`` so the package's pytest run does not need playwright or a browser.
"""

from __future__ import annotations

import argparse
import bz2
import csv
import gzip
import io
import json
import lzma
import random
import shutil
import socket
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from playwright.sync_api import Page, sync_playwright

ROOT = Path(__file__).resolve().parent.parent
SITE = ROOT / "site"
EX1 = SITE / "examples" / "human_ecoli_mix.fasta"
EX2 = SITE / "examples" / "contaminants_with_decoys.fasta"
BOOT_MS = 300_000
STEP_MS = 600_000

sys.dont_write_bytecode = True  # keep site/ free of __pycache__ (it is deployed as is)
sys.path.insert(0, str(SITE))
import toolkit as tk  # noqa: E402  (site/toolkit.py, the same code the page runs)

from fastatacular import read_fasta  # noqa: E402

CHECKS: list[str] = []


def ok(cond: bool, what: str) -> None:
    if not cond:
        raise AssertionError(what)
    CHECKS.append(what)
    print(f"  ok  {what}")


def free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


# ----------------------------------------------------------------------------- page helpers


def wait_idle(page: Page, timeout: int = STEP_MS) -> None:
    page.wait_for_function("window.toolkit && window.toolkit.ready && !window.toolkit.busy", timeout=timeout)


def click_and_wait(page: Page, selector: str, cmd: str, timeout: int = STEP_MS) -> dict:
    page.evaluate(f"delete window.toolkit.results[{cmd!r}]")
    n_err = page.evaluate("window.toolkit.errors.length")
    page.click(selector)
    page.wait_for_function(
        f"window.toolkit.results[{cmd!r}] !== undefined || window.toolkit.errors.length > {n_err}", timeout=timeout
    )
    wait_idle(page, timeout)
    errs = page.evaluate("window.toolkit.errors")
    if len(errs) > n_err:
        raise AssertionError(f"{cmd} failed in the page: {errs[-1]}")
    return page.evaluate(f"window.toolkit.results[{cmd!r}]")


def download(page: Page, selector: str) -> tuple[str, bytes]:
    with page.expect_download(timeout=STEP_MS) as info:
        page.click(selector)
    d = info.value
    data = Path(d.path()).read_bytes()
    wait_idle(page)
    return d.suggested_filename, data


def tab(page: Page, name: str) -> None:
    page.click(f'.tabs button[data-tab="{name}"]')


def load_files(page: Page, paths: list[Path], append: bool = False) -> dict:
    tab(page, "input")
    page.set_checked("#appendFiles", append)
    page.evaluate("delete window.toolkit.results.load")
    page.set_input_files("#fileInput", [str(p) for p in paths])
    page.wait_for_function("window.toolkit.results.load !== undefined", timeout=STEP_MS)
    wait_idle(page)
    return page.evaluate("window.toolkit.results.load")


def tile(page: Page, tile_id: str) -> str:
    return page.inner_text(f"#{tile_id} .v")


# ----------------------------------------------------------------------------- CPython reference


def reference_load(paths: list[Path]) -> None:
    tmp = Path(tempfile.mkdtemp())
    copies = []
    for i, p in enumerate(paths):
        dst = tmp / f"{i}_{p.name}"
        shutil.copy(p, dst)
        copies.append(str(dst))
    tk.load_files(copies, [p.name for p in paths], True)
    shutil.rmtree(tmp)


def fasta_bytes_entries(data: bytes):
    return read_fasta(io.StringIO(data.decode()))


# ----------------------------------------------------------------------------- CPython-only checks


def check_digest_matches_peptacular() -> None:
    """toolkit.fast_digest must give exactly pt.digest's peptides for every enzyme the UI offers."""
    import peptacular as pt

    print("digest equivalence (CPython)")
    rng = random.Random(1)
    seqs = [e.sequence for f in (EX1, EX2) for e in read_fasta(str(f))]
    seqs += ["".join(rng.choices("ACDEFGHIKLMNPQRSTVWYKRPDE", k=rng.randint(1, 80))) for _ in range(300)]
    seqs += ["K", "KP", "PK", "KKKK", "RPK", "D", "DD", "DAD", "E", "AAAA", "KPKPRP"]
    seqs = [s for s in seqs if tk._FAST_OK.fullmatch(s)]  # others go through pt.digest itself
    names = json.loads(tk.enzymes())
    ok("trypsin" in names and "unspecific" in names and "no_enzyme" in names, f"{len(names)} enzymes offered")
    cases = 0
    for enzyme in names:
        pattern = tk._enzyme_pattern(enzyme)
        for missed in (0, 1, 2, 3):
            for lo, hi in ((1, 1000), (7, 40), (1, 5), (6, 6)):
                for seq in seqs:
                    if enzyme == "unspecific" and len(seq) > 60:
                        continue
                    ref = sorted(p for p, _ in pt.digest(seq, enzyme, missed_cleavages=missed, min_len=lo, max_len=hi))
                    got = sorted(tk.fast_digest(seq, pattern, missed, lo, hi))
                    if ref != got:
                        raise AssertionError(f"fast_digest != pt.digest: {enzyme} missed={missed} {lo}-{hi} {seq}")
                    cases += 1
    ok(
        True,
        f"fast_digest == pt.digest (same peptides, same multiplicity) for all {len(names)} enzymes, {cases:,} cases",
    )


def peff_blocks(text: str) -> dict[str, str]:
    """PEFF database header blocks by prefix."""
    return {b.split("# Prefix=")[1].split("\n")[0]: b for b in text.split("# //") if "# Prefix=" in b}


def check_peff_decoy_flag() -> None:
    """A custom decoy prefix still marks its PEFF database Decoy=true."""
    print("PEFF decoy flag (CPython)")
    reference_load([EX1])
    o = decoy_opts({"method": "reverse"}, prefix="XYZ_")
    tk.decoys(json.dumps(o))
    text = tk.export(json.dumps({"which": "decoy", "format": "peff"})).decode()
    blocks = peff_blocks(text)
    ok("Decoy=true" in blocks["XYZ_sp"] and "Decoy=true" not in blocks["sp"], "PEFF: XYZ_sp is Decoy=true, sp is not")


# ----------------------------------------------------------------------------- feature checks

DIGEST = {"enzyme": "trypsin", "missed": 1, "min_len": 7, "max_len": 40}
CLEAN = {
    "header_regex": "",
    "regex_mode": "include",
    "organism": "",
    "taxa": "",
    "organism_mode": "include",
    "min_len": 20,
    "max_len": 0,
    "invalid": True,
    "alphabet": "ACDEFGHIKLMNPQRSTVWY",
    "allow_ambiguous": False,
    "dedup_id": True,
    "dedup_seq": True,
    "il_equivalent": True,
}


def set_digest(page: Page) -> None:
    page.select_option("#dEnzyme", DIGEST["enzyme"])
    page.fill("#dMissed", str(DIGEST["missed"]))
    page.fill("#dMin", str(DIGEST["min_len"]))
    page.fill("#dMax", str(DIGEST["max_len"]))


def check_input(page: Page, work: Path) -> None:
    print("input")
    r = load_files(page, [EX1])
    ok(r["entries"] == 50 and tile(page, "loadEntries") == "50", "example 1 loads 50 entries")
    for ext, opener in ((".gz", gzip.open), (".bz2", bz2.open), (".xz", lzma.open)):
        p = work / f"human_ecoli_mix.fasta{ext}"
        with opener(p, "wb") as fh:
            fh.write(EX1.read_bytes())
        r = load_files(page, [p])
        ok(r["entries"] == 50, f"{ext} input reads 50 entries")
    page.click("#loadExample2")
    page.wait_for_function("window.toolkit.results.load && window.toolkit.results.load.entries === 13", timeout=STEP_MS)
    wait_idle(page)
    ok(True, "example 2 button loads 13 entries")
    r = load_files(page, [EX1, EX2])
    ok(r["entries"] == 63 and len(r["files"]) == 2, "merging both examples gives 63 entries")
    ok(r["decoys_detected"] == {"DECOY_": 4, "rev_": 2}, f"existing decoys detected: {r['decoys_detected']}")
    r = load_files(page, [EX2], append=True)
    ok(r["entries"] == 76 and r["loaded"] == 76, "'add to loaded' appends a third file (76 entries)")
    load_files(page, [EX1, EX2])


def check_cleanup(page: Page) -> None:
    print("clean-up")
    reference_load([EX1, EX2])
    ref = json.loads(tk.cleanup(json.dumps(CLEAN)))
    tab(page, "clean")
    page.fill("#cMinLen", "20")
    page.set_checked("#cInvalid", True)
    page.set_checked("#cIL", True)
    r = click_and_wait(page, "#cleanBtn", "cleanup")
    ok(r["after"] == ref["after"] == 57, f"clean-up keeps 57 of 63 (page {r['after']}, python {ref['after']})")
    ok(r["reasons"] == ref["reasons"], f"drop reasons match python: {r['reasons']}")
    expected = {
        "invalid residues",
        "too short",
        "duplicate identifier",
        "duplicate sequence",
        "duplicate sequence (I=L)",
    }
    ok(set(r["reasons"]) == expected, "every planted problem is reported")
    ok(page.locator("#droppedTable tbody tr").count() == 6, "dropped table lists 6 entries")
    name, data = download(page, "#dlDropped")
    rows = list(csv.reader(io.StringIO(data.decode())))
    ok(name.endswith("_dropped.csv") and len(rows) == 7, "dropped-list CSV has header + 6 rows")

    # organism + regex filters, checked against python, then restore the full clean-up
    opts = dict(CLEAN, organism="homo sapiens", taxa="9913", header_regex="PE=1", min_len=0, invalid=False)
    reference_load([EX1, EX2])
    ref2 = json.loads(tk.cleanup(json.dumps(opts)))
    page.fill("#cOrganism", "homo sapiens")
    page.fill("#cTaxa", "9913")
    page.fill("#cHeaderRegex", "PE=1")
    page.fill("#cMinLen", "")
    page.set_checked("#cInvalid", False)
    r = click_and_wait(page, "#cleanBtn", "cleanup")
    ok(
        r["after"] == ref2["after"] and r["reasons"] == ref2["reasons"],
        f"organism/OX/regex filters match python ({r['after']} kept, {r['reasons']})",
    )
    r = click_and_wait(page, "#cleanResetBtn", "resetCleanup")
    ok(r["entries"] == 63, "undo restores 63 entries")
    page.fill("#cOrganism", "")
    page.fill("#cTaxa", "")
    page.fill("#cHeaderRegex", "")
    page.fill("#cMinLen", "20")
    page.set_checked("#cInvalid", True)
    click_and_wait(page, "#cleanBtn", "cleanup")
    reference_load([EX1, EX2])
    tk.cleanup(json.dumps(CLEAN))


def check_stats(page: Page) -> None:
    print("stats")
    tab(page, "stats")
    set_digest(page)
    r = click_and_wait(page, "#statsPepBtn", "stats")
    ref = json.loads(tk.stats(json.dumps({"digest": True, **DIGEST})))
    ok(r["entries"] == ref["entries"] == 57, "stats: 57 entries")
    ok(r["peptides"] == ref["peptides"], f"theoretical peptides match python: {r['peptides']}")
    ok(
        r["composition"] == ref["composition"] and r["length"] == ref["length"],
        "composition and length histogram match",
    )
    orgs = {o["name"]: o["count"] for o in r["organisms"]}
    ok(r["organisms"] == ref["organisms"] and "Homo sapiens (OX=9606)" in orgs, f"organism breakdown matches: {orgs}")
    ok(page.locator("#orgTable tbody tr").count() == len(orgs), "organism table rendered")
    ok(page.locator("#statsResult svg").count() == 2, "length and composition charts render")
    ok(tile(page, "statPeptides") == f"{ref['peptides']['total']:,}", "peptide tile shows the count")


DECOY_CASES = [
    {"method": "reverse"},
    {"method": "pseudo_reverse"},
    {"method": "shuffle"},
    {"method": "debruijn", "k": 3},
    {"method": "markov", "model": "yeast"},
    {"method": "markov", "model": "trained"},
]


def decoy_opts(case: dict, prefix: str = "DECOY_", concatenate: bool = True) -> dict:
    return {
        "method": case["method"],
        "prefix": prefix,
        "seed": "7",
        "keep_residues": "",
        "keep_met": True,
        "keep_cterm": 0,
        "k": case.get("k", 2),
        "model": case.get("model", "human"),
        "concatenate": concatenate,
    }


def set_decoy_form(page: Page, o: dict) -> None:
    page.select_option("#mMethod", o["method"])
    page.fill("#mPrefix", o["prefix"])
    page.fill("#mSeed", o["seed"])
    page.fill("#mKeepRes", o["keep_residues"])
    page.set_checked("#mKeepMet", o["keep_met"])
    page.fill("#mKeepCterm", str(o["keep_cterm"]))
    if o["method"] == "debruijn":
        page.fill("#mK", str(o["k"]))
    if o["method"] == "markov":
        page.select_option("#mModel", o["model"])
    page.check('input[name="mOutput"][value="%s"]' % ("concat" if o["concatenate"] else "decoys"))


def check_decoys_and_qc(page: Page) -> None:
    print("decoys + QC")
    # train a Markov model in the page and in python
    tab(page, "decoys")
    page.select_option("#mMethod", "markov")
    page.fill("#mOrder", "2")
    r = click_and_wait(page, "#trainBtn", "trainModel")
    ref_model = json.loads(tk.train_model(2))
    ok(
        r["metadata"]["residues"] == ref_model["metadata"]["residues"] and r["order"] == 2,
        f"Markov model trained on the input ({r['metadata']['entries']} entries)",
    )
    name, model_bytes = download(page, "#dlModel")
    ok(name.endswith("_markov.json.gz") and model_bytes[:2] == b"\x1f\x8b", "trained model downloads as .json.gz")

    for i, case in enumerate(DECOY_CASES):
        o = decoy_opts(case)
        label = case["method"] + (f" ({case['model']})" if "model" in case else "")
        tab(page, "decoys")
        set_decoy_form(page, o)
        r = click_and_wait(page, "#decoyBtn", "decoys")
        ref = json.loads(tk.decoys(json.dumps(o)))
        ok(
            r["targets"] == ref["targets"] == 51 and r["decoys"] == 51 and r["output_entries"] == 102,
            f"{label}: 51 targets -> 51 decoys, 102 entries concatenated",
        )
        ok(r["existing_decoys_skipped"] == 6, f"{label}: 6 existing decoys left out")
        ok(
            r["examples"] == ref["examples"] and r["composition_l1"] == ref["composition_l1"],
            f"{label}: decoys identical to CPython (seeded)",
        )
        tab(page, "qc")
        set_digest(page)
        q = click_and_wait(page, "#qcBtn", "qc")
        refq = json.loads(tk.decoy_qc(json.dumps(DIGEST)))
        for key in ("shared", "shared_fraction", "shared_il", "balance", "length", "mass"):
            if q[key] != refq[key]:
                raise AssertionError(f"{label}: QC {key} differs: page {q[key]} python {refq[key]}")
        ok(
            True,
            f"{label}: QC matches python (shared {q['shared_fraction']:.4%}, balance {q['balance']:.3f}, "
            f"{q['target']['distinct']} target / {q['decoy']['distinct']} decoy peptides)",
        )
        ok(q["target_cached"] == (i > 0), f"{label}: target digest {'reused' if i else 'computed'}")
        ok(page.locator("#qcResult svg").count() == 2, f"{label}: length and mass charts render")
        ok(tile(page, "qcShared") == f"{100 * q['shared_fraction']:.2f}%", f"{label}: shared tile shows the fraction")

    # the shared-fraction definition, recomputed from scratch with peptacular
    import peptacular as pt

    d = tk.STATE["decoy"]
    tp = {
        p for e in d["targets"] for p, _ in pt.digest(e.sequence, "trypsin", missed_cleavages=1, min_len=7, max_len=40)
    }
    dp = {
        p for e in d["decoys"] for p, _ in pt.digest(e.sequence, "trypsin", missed_cleavages=1, min_len=7, max_len=40)
    }
    ok(
        abs(q["shared_fraction"] - len(tp & dp) / len(dp)) < 1e-12,
        "shared fraction = |T & D| / |D| from a fresh digest",
    )
    ok(abs(tk.peptide_mass("PEPTIDEK") - pt.mass("PEPTIDEK")) < 1e-9, "peptide mass table agrees with pt.mass")

    # uploaded model = trained model
    tab(page, "decoys")
    mpath = Path(tempfile.mkdtemp()) / "model.json.gz"
    mpath.write_bytes(model_bytes)
    page.evaluate("delete window.toolkit.results.loadModel")
    page.set_input_files("#modelFile", str(mpath))
    page.wait_for_function("window.toolkit.results.loadModel !== undefined", timeout=STEP_MS)
    wait_idle(page)
    o = decoy_opts({"method": "markov", "model": "uploaded"})
    set_decoy_form(page, o)
    r = click_and_wait(page, "#decoyBtn", "decoys")
    ref = json.loads(tk.decoys(json.dumps(decoy_opts({"method": "markov", "model": "trained"}))))
    ok(r["examples"] == ref["examples"], "uploaded model gives the same decoys as the trained one")

    # decoy-only output, custom prefix, keep_residues, downloaded and parsed
    o = decoy_opts({"method": "shuffle"}, prefix="XYZ_", concatenate=False)
    o["keep_residues"] = "KR"
    set_decoy_form(page, o)
    r = click_and_wait(page, "#decoyBtn", "decoys")
    ok(r["output_entries"] == 51, "decoy-only output has 51 entries")
    name, data = download(page, "#decoyDlBtn")
    entries = fasta_bytes_entries(data)
    ok(
        len(entries) == 51 and all(e.identifier.startswith("XYZ_") for e in entries),
        f"{name}: 51 entries, all with prefix XYZ_",
    )
    tref = {e.identifier: e.sequence for e in tk.STATE["entries"]}
    ok(
        all(
            [i for i, c in enumerate(e.sequence) if c in "KR"]
            == [i for i, c in enumerate(tref[e.identifier[4:]]) if c in "KR"]
            for e in entries
        ),
        "keep_residues=KR keeps every K/R in place",
    )
    ok(all(e.sequence[0] == tref[e.identifier[4:]][0] for e in entries), "N-terminal residue kept")
    ok(name.endswith("_decoy.fasta") and "target_decoy" not in name, f"decoy-only download is named {name}")
    tab(page, "export")
    page.check('input[name="xWhich"][value="decoy"]')
    page.select_option("#xFormat", "peff")
    page.set_checked("#xGzip", False)
    name, data = download(page, "#exportBtn")
    text = data.decode()
    ok(name.endswith("_decoy.peff") and "target_decoy" not in name, f"decoy-only export is named {name}")
    ok(
        "Decoy=true" in peff_blocks(text)["XYZ_sp"],
        "PEFF export marks the custom-prefix (XYZ_) database Decoy=true",
    )


def check_export(page: Page) -> None:
    print("export")
    tab(page, "decoys")
    set_decoy_form(page, decoy_opts({"method": "pseudo_reverse"}))
    click_and_wait(page, "#decoyBtn", "decoys")
    tab(page, "export")
    page.check('input[name="xWhich"][value="working"]')
    page.select_option("#xFormat", "fasta")
    page.set_checked("#xGzip", False)
    name, data = download(page, "#exportBtn")
    entries = fasta_bytes_entries(data)
    ok(len(entries) == 57 and name.endswith(".fasta"), f"{name}: current entries as FASTA (57)")
    page.select_option("#xFormat", "csv")
    name, data = download(page, "#exportBtn")
    rows = list(csv.DictReader(io.StringIO(data.decode())))
    ok(len(rows) == 57 and list(rows[0]) == list(tk.RECORD_KEYS), f"{name}: CSV with RECORD_KEYS columns, 57 rows")
    page.check('input[name="xWhich"][value="decoy"]')
    page.select_option("#xFormat", "fasta")
    page.set_checked("#xGzip", True)
    name, data = download(page, "#exportBtn")
    entries = fasta_bytes_entries(gzip.decompress(data))
    ok(
        name.endswith("_target_decoy.fasta.gz")
        and len(entries) == 102
        and sum(e.identifier.startswith("DECOY_") for e in entries) == 51,
        f"{name}: gz target+decoy, 51 decoys",
    )
    page.select_option("#xFormat", "peff")
    page.set_checked("#xGzip", False)
    name, data = download(page, "#exportBtn")
    import warnings

    from pefftacular import read_peff

    with warnings.catch_warnings():
        warnings.simplefilter("error")  # a PeffWarning (missing mandatory key etc.) fails the check
        _, peff_entries = read_peff(io.StringIO(data.decode()))
    n = len(peff_entries)
    ok(
        data.startswith(b"# PEFF 1.0") and n == 102,
        f"{name}: PEFF with 102 entries, read back by pefftacular without warnings",
    )


def check_errors_shown(page: Page) -> None:
    print("errors")
    tab(page, "decoys")
    page.fill("#mPrefix", "bad prefix")
    page.select_option("#mMethod", "reverse")
    n = page.evaluate("window.toolkit.errors.length")
    page.click("#decoyBtn")
    page.wait_for_function(f"window.toolkit.errors.length > {n}", timeout=STEP_MS)
    wait_idle(page)
    ok("prefix" in page.inner_text("#decoyResult .error").lower(), "an invalid prefix shows the DecoyError in the page")
    page.fill("#mPrefix", "DECOY_")


# ----------------------------------------------------------------------------- timing


def synthetic_proteome(n: int, path: Path) -> None:
    rng = random.Random(20260925)
    aas = "ACDEFGHIKLMNPQRSTVWY"
    weights = [7.0, 2.3, 4.7, 7.1, 3.7, 6.6, 2.6, 4.3, 5.7, 10.0, 2.1, 3.6, 6.3, 4.8, 5.6, 8.3, 5.4, 6.0, 1.2, 2.7]
    orgs = [("Homo sapiens", 9606)] * 8 + [("Mus musculus", 10090), ("Escherichia coli (strain K12)", 83333)]
    with path.open("w") as fh:
        for i in range(n):
            length = max(30, min(8000, int(rng.lognormvariate(6.1, 0.65))))
            seq = "M" + "".join(rng.choices(aas, weights, k=length - 1))
            os_name, ox = rng.choice(orgs)
            fh.write(f">sp|S{i:05d}|SYN{i}_SYNTH Synthetic protein {i} OS={os_name} OX={ox} GN=G{i} PE=1 SV=1\n")
            for j in range(0, len(seq), 60):
                fh.write(seq[j : j + 60] + "\n")


def timing(page: Page, work: Path, n: int, method: str = "pseudo_reverse") -> dict:
    print(f"timing: {n:,} synthetic entries")
    path = work / f"synthetic_{n}.fasta"
    synthetic_proteome(n, path)
    size_mb = path.stat().st_size / 1e6
    t: dict[str, float] = {}
    t0 = time.perf_counter()
    r = load_files(page, [path])
    t["load"] = time.perf_counter() - t0
    ok(r["entries"] == n, f"loaded {n:,} entries ({size_mb:.1f} MB)")
    tab(page, "stats")
    s = time.perf_counter()
    click_and_wait(page, "#statsBtn", "stats")
    t["stats"] = time.perf_counter() - s
    tab(page, "decoys")
    set_decoy_form(page, dict(decoy_opts({"method": method}), seed="1"))
    s = time.perf_counter()
    click_and_wait(page, "#decoyBtn", "decoys")
    t["decoys"] = time.perf_counter() - s
    tab(page, "qc")
    set_digest(page)
    s = time.perf_counter()
    q = click_and_wait(page, "#qcBtn", "qc")
    t["qc"] = time.perf_counter() - s
    tab(page, "decoys")
    set_decoy_form(page, dict(decoy_opts({"method": "shuffle"}), seed="1"))
    click_and_wait(page, "#decoyBtn", "decoys")
    tab(page, "qc")
    s = time.perf_counter()
    q2 = click_and_wait(page, "#qcBtn", "qc")
    t["qc_other_method"] = time.perf_counter() - s
    ok(q2["target_cached"], "QC after switching to shuffle reuses the target digest")
    tab(page, "export")
    page.check('input[name="xWhich"][value="decoy"]')
    page.select_option("#xFormat", "fasta")
    page.set_checked("#xGzip", True)
    s = time.perf_counter()
    _, data = download(page, "#exportBtn")
    t["export_gz"] = time.perf_counter() - s
    t["total"] = time.perf_counter() - t0
    ok(len(fasta_bytes_entries(gzip.decompress(data))) == 2 * n, f"exported {2 * n:,} target+decoy entries")
    return {
        "entries": n,
        "mb": round(size_mb, 1),
        "seconds": {k: round(v, 1) for k, v in t.items()},
        "qc": {
            "target_distinct": q["target"]["distinct"],
            "decoy_distinct": q["decoy"]["distinct"],
            "shared_fraction": q["shared_fraction"],
            "shuffle_shared_fraction": q2["shared_fraction"],
        },
    }


# ----------------------------------------------------------------------------- main


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--entries", type=int, default=20_000, help="synthetic entries for the timing run (0 = skip)")
    ap.add_argument("--port", type=int, default=0)
    ap.add_argument("--cpython-only", action="store_true", help="only the checks that need no browser")
    args = ap.parse_args()

    try:
        check_digest_matches_peptacular()
        check_peff_decoy_flag()
    except AssertionError as err:
        print(f"\nFAIL: {err}")
        return 1
    if args.cpython_only:
        print(f"\nPASS: {len(CHECKS)} checks")
        return 0

    port = args.port or free_port()
    server = subprocess.Popen(
        [sys.executable, "-m", "http.server", str(port), "--bind", "127.0.0.1"],
        cwd=SITE,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    problems: list[str] = []
    work = Path(tempfile.mkdtemp(prefix="fasta-toolkit-smoke-"))
    try:
        time.sleep(0.5)
        with sync_playwright() as p:
            browser = p.chromium.launch(headless=True)
            ctx = browser.new_context(accept_downloads=True)
            page = ctx.new_page()
            page.on("console", lambda m: problems.append(f"console.{m.type}: {m.text}") if m.type == "error" else None)
            page.on("pageerror", lambda e: problems.append(f"pageerror: {e}"))
            t0 = time.perf_counter()
            page.goto(f"http://127.0.0.1:{port}/", wait_until="load")
            wait_idle(page, BOOT_MS)
            boot = time.perf_counter() - t0
            versions = page.evaluate("window.toolkit.versions")
            print(f"boot {boot:.1f} s, {versions}")
            ok(page.locator("text=never uploaded").count() == 1, "privacy note is shown")
            ok(
                page.locator('a[href="https://pgarrett-scripps.github.io/fastaviewer/"]').count() == 1,
                "links to fastaviewer",
            )

            check_input(page, work)
            check_cleanup(page)
            check_stats(page)
            check_decoys_and_qc(page)
            check_export(page)
            check_errors_shown(page)
            result = {"boot_seconds": round(boot, 1), "versions": versions}
            if args.entries:
                result["timing"] = timing(page, work, args.entries)
            browser.close()
        # The deliberate bad-prefix step logs no console error (errors are shown in the page).
        ok(not problems, "no console or page errors" + ("" if not problems else f": {problems}"))
        print(json.dumps(result, indent=2))
        print(f"\nPASS: {len(CHECKS)} checks")
        return 0
    except Exception as err:  # noqa: BLE001
        print(f"\nFAIL: {err}")
        if problems:
            print("console/page errors:", *problems, sep="\n  ")
        return 1
    finally:
        server.terminate()
        server.wait(timeout=10)
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
