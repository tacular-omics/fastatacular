// UI thread: forms, rendering and downloads. All Python runs in worker.js.
(function () {
  const $ = (id) => document.getElementById(id);
  const { barChart, lineChart } = window.Charts;
  const EX1 = "examples/human_ecoli_mix.fasta";
  const EX2 = "examples/contaminants_with_decoys.fasta";
  const EXAMPLES = { loadExample1: [EX1], loadExample2: [EX2], loadBothExamples: [EX1, EX2] };
  // Test hook: the last result of every command (scripts/site_smoke.py reads it).
  const toolkit = (window.toolkit = { ready: false, busy: false, results: {}, errors: [] });

  // ---------------------------------------------------------------- worker RPC
  const worker = new Worker("worker.js");
  let nextId = 1;
  const pending = new Map();
  worker.onmessage = (ev) => {
    const m = ev.data;
    if (m.type === "status") return setStatus(m.text);
    if (m.type === "progress") return setProgress(m.stage, m.done, m.total);
    const p = pending.get(m.id);
    if (!p) return;
    pending.delete(m.id);
    m.ok ? p.resolve(m.result) : p.reject(Object.assign(new Error(m.error), { detail: m.detail }));
  };
  worker.onerror = (ev) => {
    setStatus("Could not start the Python worker: " + (ev.message || ev) + ". Reload the page to try again.");
    setState("error");
  };
  function call(cmd, args, transfer = []) {
    const id = nextId++;
    return new Promise((resolve, reject) => {
      pending.set(id, { resolve, reject });
      worker.postMessage({ id, cmd, args }, transfer);
    });
  }

  // ---------------------------------------------------------------- status line
  const setState = (s) => ($("runtime").dataset.state = s);
  function setStatus(text) { $("status").textContent = text; }
  function setProgress(stage, done, total) {
    $("progressBar").value = total ? Math.round((100 * done) / total) : 0;
    $("progressText").textContent = `${stage}: ${done.toLocaleString()} / ${total.toLocaleString()}`;
  }
  const state = { loaded: 0, entries: 0, decoys: false, decoyConcat: true, model: false, dropped: 0, baseName: "database" };
  // File name suffix for the decoy database, from how it was made (not the current form).
  const decoySuffix = () => (state.decoyConcat ? "_target_decoy" : "_decoy");

  function refreshButtons() {
    const idle = toolkit.ready && !toolkit.busy;
    const has = idle && state.entries > 0;
    for (const id of ["loadExample1", "loadExample2", "loadBothExamples"]) $(id).disabled = !idle;
    document.querySelectorAll(".load-example").forEach((b) => (b.disabled = !idle));
    $("clearBtn").disabled = !idle || state.loaded === 0;
    $("fileInput").disabled = !idle;
    for (const id of ["statsBtn", "statsPepBtn", "decoyBtn", "trainBtn"]) $(id).disabled = !has;
    $("cleanBtn").disabled = !idle || state.loaded === 0;
    $("cleanResetBtn").disabled = !idle || state.loaded === 0;
    $("dlDropped").disabled = !idle || state.dropped === 0;
    $("qcBtn").disabled = !idle || !state.decoys;
    $("decoyDlBtn").disabled = !idle || !state.decoys;
    $("dlModel").disabled = !idle || !state.model;
    const which = document.querySelector('input[name="xWhich"]:checked').value;
    $("exportBtn").disabled = !has || (which === "decoy" && !state.decoys);
    document.querySelectorAll(".needs-entries").forEach((e) => (e.hidden = state.loaded > 0));
    document.querySelectorAll(".needs-decoys").forEach((e) => (e.hidden = state.decoys));
  }

  async function run(label, cmd, args, transfer, target) {
    toolkit.busy = true;
    setState("loading");
    if (target) target.setAttribute("aria-busy", "true");
    refreshButtons();
    setStatus(label);
    $("progressBar").value = 0;
    $("progressText").textContent = "";
    const t0 = performance.now();
    try {
      const result = await call(cmd, args, transfer);
      toolkit.results[cmd] = result;
      const secs = ((performance.now() - t0) / 1000).toFixed(1);
      setStatus(`Ready. ${label.replace(/\.+$|…$/, "")} took ${secs} s.`);
      return result;
    } catch (err) {
      toolkit.errors.push({ cmd, message: err.message });
      setStatus("Ready. The last step failed; the message is shown below it.");
      if (target) target.innerHTML = alert("error", "Could not finish.", esc(err.message), "error");
      throw err;
    } finally {
      toolkit.busy = false;
      setState("ready");
      $("progressText").textContent = "";
      if (target) target.removeAttribute("aria-busy");
      refreshButtons();
    }
  }
  const quiet = (p) => p.catch(() => {});

  // ---------------------------------------------------------------- helpers
  const esc = (s) => String(s ?? "").replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]);
  const n = (x) => x.toLocaleString();
  // Key/value results: [label, value, optional element id].
  const stats = (items) =>
    `<dl class="stats">${items.map(([k, v, id]) => `<div${id ? ` id="${id}"` : ""}><dt>${esc(k)}</dt><dd>${esc(v)}</dd></div>`).join("")}</dl>`;
  // kind: success | warning | error; html is already escaped.
  const alert = (kind, lead, html, extraClass = "", id = "") =>
    `<p class="alert alert-${kind}${extraClass ? " " + extraClass : ""}"${id ? ` id="${id}"` : ""}${kind === "error" ? ' role="alert"' : ""}><strong>${lead}</strong> ${html}</p>`;
  const VERDICT = { good: ["success", "OK."], warn: ["warning", "Check."], bad: ["error", "Problem."] };
  const verdict = (cls, html, id) => alert(VERDICT[cls][0], VERDICT[cls][1], html, "", id);
  // headers: [text, cls] with cls 1 (= "num"), "seq" or none; rows: arrays of already-escaped cell HTML.
  function table(headers, rows, id, extra = "auto") {
    const cls = headers.map(([, c]) => (c === 1 ? ' class="num"' : c ? ` class="${c}"` : ""));
    const th = headers.map(([h], i) => `<th${cls[i] === ' class="num"' ? cls[i] : ""}>${h}</th>`).join("");
    const body = rows.map((r) => `<tr>${r.map((c, i) => `<td${cls[i]}>${c}</td>`).join("")}</tr>`).join("");
    return `<div class="table-wrap ${extra}"><table class="data"${id ? ` id="${id}"` : ""}><thead><tr>${th}</tr></thead><tbody>${body}</tbody></table></div>`;
  }
  const pct = (x, d = 2) => (100 * x).toFixed(d) + "%";
  const num = (id) => { const v = $(id).value.trim(); return v === "" ? 0 : parseInt(v, 10); };
  function download(bytes, name, type = "application/octet-stream") {
    const url = URL.createObjectURL(new Blob([bytes], { type }));
    const a = document.createElement("a");
    a.href = url;
    a.download = name;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 10000);
  }
  function digestOpts() {
    return { enzyme: $("dEnzyme").value, missed: num("dMissed"), min_len: num("dMin") || 1, max_len: num("dMax") || 1000 };
  }
  function prefixTable(found) {
    const keys = Object.keys(found || {});
    return keys.length ? keys.map((k) => `<code>${esc(k)}</code> ${found[k]}`).join(", ") : "none";
  }

  // ---------------------------------------------------------------- tabs
  const tabs = [...document.querySelectorAll('.tabs [role="tab"]')];
  function showTab(name) {
    tabs.forEach((b) => b.setAttribute("aria-selected", String(b.dataset.tab === name)));
    document.querySelectorAll('[role="tabpanel"]').forEach((p) => (p.hidden = p.dataset.panel !== name));
    const ds = $("digestSettings");
    if (name === "stats" || name === "qc") {
      const panel = document.querySelector(`[data-panel="${name}"] > .panel`);
      panel.insertBefore(ds, panel.querySelector(":scope > .actions"));
      ds.hidden = false;
    } else ds.hidden = true;
  }
  tabs.forEach((b) => b.addEventListener("click", () => showTab(b.dataset.tab)));
  document.querySelector(".tabs").addEventListener("keydown", (e) => {
    const i = tabs.indexOf(document.activeElement);
    if (i < 0 || (e.key !== "ArrowRight" && e.key !== "ArrowLeft")) return;
    const next = tabs[(i + (e.key === "ArrowRight" ? 1 : tabs.length - 1)) % tabs.length];
    next.focus();
    showTab(next.dataset.tab);
  });
  showTab("input");

  // ---------------------------------------------------------------- 1. input
  function renderSummary(s) {
    state.loaded = s.loaded;
    state.entries = s.entries;
    state.decoys = false;
    const files = s.files.map((f) => [esc(f.name), n(f.entries), (f.bytes / 1e6).toFixed(2), f.seconds]);
    $("loadSummary").innerHTML = !s.loaded ? "" :
      stats([["Entries", n(s.entries), "loadEntries"], ["Residues", n(s.residues)], ["Files", s.files.length]]) +
      (s.files.length ? table([["File"], ["Entries", 1], ["Size (MB)", 1], ["Parse time (s)", 1]], files) : "") +
      (Object.keys(s.decoys_detected).length ? verdict("warn", `Existing decoys found (${prefixTable(s.decoys_detected)}). The Decoys step leaves them out.`) : "") +
      (s.preview.length ? `<div><h3 class="label">First headers</h3><div class="seq-block">${s.preview.map(esc).join("<br>")}</div></div>` : "");
    $("cleanResult").innerHTML = "";
    $("decoyResult").innerHTML = "";
    $("qcResult").innerHTML = "";
  }

  async function loadFiles(files) {
    if (!files.length) return;
    state.baseName = files[0].name.replace(/\.(gz|bz2|xz)$/i, "").replace(/\.[^.]+$/, "") || "database";
    const payload = [];
    for (const f of files) payload.push({ name: f.name, buffer: f.buffer || (await f.arrayBuffer()) });
    const replace = !$("appendFiles").checked;
    const s = await run("Reading " + files.map((f) => f.name).join(", ") + "…", "load", { files: payload, replace }, payload.map((p) => p.buffer), $("loadSummary"));
    renderSummary(s);
    refreshButtons();
  }
  async function loadExamples(urls) {
    const files = [];
    for (const u of urls) files.push({ name: u.split("/").pop(), buffer: await (await fetch(u)).arrayBuffer() });
    await loadFiles(files);
  }
  $("fileInput").addEventListener("change", (e) => { quiet(loadFiles([...e.target.files])); e.target.value = ""; });
  const dz = $("dropZone");
  dz.addEventListener("dragover", (e) => { e.preventDefault(); dz.classList.add("over"); });
  dz.addEventListener("dragleave", () => dz.classList.remove("over"));
  dz.addEventListener("drop", (e) => { e.preventDefault(); dz.classList.remove("over"); if (toolkit.ready && !toolkit.busy) quiet(loadFiles([...e.dataTransfer.files])); });
  for (const [id, urls] of Object.entries(EXAMPLES)) {
    $(id).addEventListener("click", () => {
      if (id === "loadBothExamples") $("appendFiles").checked = false;
      quiet(loadExamples(urls));
    });
  }
  // "Load example" in the empty state of the other tabs: both examples, merged.
  document.querySelectorAll(".load-example").forEach((b) => b.addEventListener("click", () => quiet((async () => {
    $("appendFiles").checked = false;
    await loadExamples([EX1, EX2]);
    if (b.dataset.then === "decoys") await makeDecoys();
  })())));
  $("clearBtn").addEventListener("click", async () => { renderSummary(await run("Clearing…", "clear")); state.model = false; refreshButtons(); });

  // ---------------------------------------------------------------- 2. clean-up
  function cleanOpts() {
    return {
      header_regex: $("cHeaderRegex").value, regex_mode: $("cRegexMode").value, regex_case: $("cRegexCase").checked,
      organism: $("cOrganism").value, taxa: $("cTaxa").value, organism_mode: $("cOrgMode").value,
      min_len: num("cMinLen"), max_len: num("cMaxLen"),
      invalid: $("cInvalid").checked, alphabet: $("cAlphabet").value, allow_ambiguous: $("cAllowAmb").checked,
      dedup_id: $("cDedupId").checked, dedup_seq: $("cDedupSeq").checked, il_equivalent: $("cIL").checked,
    };
  }
  $("cleanBtn").addEventListener("click", () => quiet((async () => {
    const r = await run("Cleaning up…", "cleanup", cleanOpts(), [], $("cleanResult"));
    state.entries = r.after; state.dropped = r.dropped; state.decoys = false;
    const reasons = Object.entries(r.reasons).map(([k, v]) => [esc(k), n(v)]);
    const rows = r.rows.map((x) => [esc(x[0]), esc(x[1]), esc(x[2])]);
    $("cleanResult").innerHTML =
      stats([["Before", n(r.before)], ["Kept", n(r.after), "cleanKept"], ["Dropped", n(r.dropped), "cleanDropped"]]) +
      (r.dropped
        ? table([["Reason"], ["Entries", 1]], reasons) +
          `<h3 class="label">Dropped entries${r.dropped > r.rows.length ? ` (first ${r.rows.length}; download the CSV for all)` : ""}</h3>` +
          table([["Identifier"], ["Reason"], ["Detail"]], rows, "droppedTable", "")
        : `<p class="empty">Nothing dropped: every entry passed the filters.</p>`);
    refreshButtons();
  })()));
  $("cleanResetBtn").addEventListener("click", () => quiet((async () => {
    const s = await run("Undoing clean-up…", "resetCleanup");
    state.entries = s.entries; state.dropped = 0; state.decoys = false;
    $("cleanResult").innerHTML = `<p class="alert">Using all ${n(s.entries)} loaded entries.</p>`;
    refreshButtons();
  })()));
  $("dlDropped").addEventListener("click", () => quiet((async () => {
    download(await run("Exporting…", "export", { which: "dropped" }), `${state.baseName}_dropped.csv`, "text/csv");
  })()));

  // ---------------------------------------------------------------- 3. stats
  async function doStats(digest) {
    const r = await run(digest ? "Computing stats and digesting…" : "Computing stats…", "stats", { digest, ...digestOpts() }, [], $("statsResult"));
    const out = $("statsResult");
    const L = r.length;
    const items = [["Entries", n(r.entries), "statEntries"], ["Residues", n(r.residues)],
      ["Length min / median / max (residues)", `${L.min} / ${L.median} / ${L.max}`], ["Mean length (residues)", L.mean],
      ["Organisms", r.organism_count], ["Duplicate identifiers", r.duplicate_ids], ["Duplicate sequences", r.duplicate_sequences],
      ["Entries with non-standard residues", r.entries_with_nonstandard]];
    if (r.peptides) {
      const d = digestOpts();
      items.push([`Peptides (${d.enzyme}, ${d.missed} missed, ${d.min_len}–${d.max_len} residues)`, n(r.peptides.total), "statPeptides"]);
      items.push(["Distinct peptides", n(r.peptides.distinct), "statDistinct"]);
    }
    out.innerHTML = stats(items) +
      (Object.keys(r.decoys_detected).length ? verdict("warn", `Contains decoy entries: ${prefixTable(r.decoys_detected)}.`) : "") +
      (r.peptides && r.peptides.failed ? verdict("warn", `${r.peptides.failed} sequences could not be digested.`) : "");
    const charts = document.createElement("div");
    charts.className = "charts";
    out.appendChild(charts);
    const h = L.hist;
    barChart(charts, {
      title: "Protein length", xTitle: "Length (residues)", yTitle: "Entries",
      labels: h.counts.map((_, i) => (i === h.counts.length - 1 && h.overflow ? `${Math.round(h.edges[i])}+` : `${Math.round(h.edges[i])}`)),
      values: h.counts,
    });
    barChart(charts, {
      title: "Amino-acid composition", xTitle: "Residue", yTitle: "Share of residues (%)", labelEvery: 1,
      labels: r.composition.map((c) => c.aa), values: r.composition.map((c) => 100 * c.fraction), valueFmt: (v) => v.toFixed(2) + "%",
    });
    if (r.other_residues.length) out.insertAdjacentHTML("beforeend", `<p class="help">Other characters: ${r.other_residues.map((o) => `<code>${esc(o.aa)}</code> ${o.count}`).join(", ")}</p>`);
    const orgs = r.organisms.map((o) => [esc(o.name), n(o.count), pct(o.count / r.entries, 1)]);
    out.insertAdjacentHTML("beforeend", `<div><h3>Organisms${r.organism_count > 20 ? " (top 20)" : ""}</h3>` +
      table([["OS= (OX=)"], ["Entries", 1], ["Share of entries", 1]], orgs, "orgTable", "") + "</div>");
  }
  $("statsBtn").addEventListener("click", () => quiet(doStats(false)));
  $("statsPepBtn").addEventListener("click", () => quiet(doStats(true)));

  // ---------------------------------------------------------------- 4. decoys
  function syncMethod() {
    const m = $("mMethod").value;
    document.querySelectorAll(".only-markov").forEach((e) => (e.hidden = m !== "markov"));
    document.querySelectorAll(".only-debruijn").forEach((e) => (e.hidden = m !== "debruijn"));
  }
  $("mMethod").addEventListener("change", syncMethod);
  syncMethod();
  function modelInfo(r, how) {
    state.model = true;
    const md = r.metadata || {};
    $("modelInfo").textContent = `${how} order-${r.order} model` + (md.entries ? ` (${md.entries} entries, ${md.residues} residues)` : "");
    refreshButtons();
  }
  $("trainBtn").addEventListener("click", () => quiet((async () => {
    const r = await run("Training Markov model…", "trainModel", { order: num("mOrder") }, [], $("modelInfo"));
    $("mModel").value = "trained";
    modelInfo(r, "Trained");
  })()));
  $("modelFile").addEventListener("change", (e) => quiet((async () => {
    const f = e.target.files[0];
    if (!f) return;
    const buffer = await f.arrayBuffer();
    const r = await run("Loading model…", "loadModel", { file: { name: f.name, buffer } }, [buffer], $("modelInfo"));
    $("mModel").value = "uploaded";
    modelInfo(r, "Uploaded");
    e.target.value = "";
  })()));
  $("dlModel").addEventListener("click", () => quiet((async () => {
    download(await run("Exporting…", "export", { which: "model" }), `${state.baseName}_markov.json.gz`, "application/gzip");
  })()));
  function decoyOpts() {
    return {
      method: $("mMethod").value, prefix: $("mPrefix").value.trim(), seed: $("mSeed").value,
      keep_residues: $("mKeepRes").value.trim(), keep_met: $("mKeepMet").checked, keep_cterm: num("mKeepCterm"),
      k: num("mK") || 2, model: $("mModel").value,
      concatenate: document.querySelector('input[name="mOutput"]:checked').value === "concat",
    };
  }
  async function makeDecoys() {
    const o = decoyOpts();
    const r = await run(`Making ${o.method} decoys…`, "decoys", o, [], $("decoyResult"));
    state.decoys = true;
    state.decoyConcat = r.params.concatenate;
    const ex = r.examples.map((e) => [esc(e.header.split(/\s/)[0]), esc(e.target), esc(e.decoy)]);
    $("decoyResult").innerHTML =
      stats([["Targets", n(r.targets), "decoyTargets"], ["Decoys", n(r.decoys), "decoyCount"],
        ["Entries in output", n(r.output_entries), "decoyOutput"], ["Decoy = target", r.identical_to_target],
        ["Composition difference (L1, 0–2)", r.composition_l1], ["Time (s)", r.seconds]]) +
      (r.existing_decoys_skipped
        ? verdict("warn", `Left out ${r.existing_decoys_skipped} existing decoy entries (${prefixTable(r.existing_prefixes)}); they are not in the output.`, "decoySkipped")
        : verdict("good", "No existing decoys in the input.")) +
      (r.identical_to_target ? verdict("warn", `${r.identical_to_target} decoys equal their target (every stretch between kept residues is too short to change).`) : "") +
      `<div><h3 class="label">First decoys (first 60 residues)</h3>` +
      table([["Decoy identifier"], ["Target sequence", "seq"], ["Decoy sequence", "seq"]], ex) + "</div>" +
      `<p>Next: <a href="#" id="toQc">check the decoys in Decoy QC</a>.</p>`;
    $("toQc").addEventListener("click", (e) => { e.preventDefault(); showTab("qc"); $("tab-qc").focus(); });
    refreshButtons();
  }
  $("decoyBtn").addEventListener("click", () => quiet(makeDecoys()));
  $("decoyDlBtn").addEventListener("click", () => quiet((async () => {
    const bytes = await run("Exporting…", "export", { which: "decoy", format: "fasta", line_width: num("xLineWidth") });
    download(bytes, `${state.baseName}${decoySuffix()}.fasta`, "text/plain");
  })()));

  // ---------------------------------------------------------------- 5. QC
  $("qcBtn").addEventListener("click", () => quiet((async () => {
    const r = await run("Digesting targets and decoys…", "qc", digestOpts(), [], $("qcResult"));
    const out = $("qcResult");
    const sharedClass = r.shared_fraction < 0.005 ? "good" : r.shared_fraction < 0.02 ? "warn" : "bad";
    const bal = r.balance;
    const balClass = Math.abs(bal - 1) <= 0.05 ? "good" : Math.abs(bal - 1) <= 0.15 ? "warn" : "bad";
    const side = (name, s) => [name, n(s.proteins), n(s.peptides), n(s.distinct)];
    out.innerHTML =
      stats([["Shared decoy peptides", pct(r.shared_fraction), "qcShared"], ["Shared if I = L", pct(r.shared_il_fraction), "qcSharedIL"],
        ["Distinct target peptides", n(r.target.distinct), "qcTarget"], ["Distinct decoy peptides", n(r.decoy.distinct), "qcDecoy"],
        ["Decoy / target peptides", bal.toFixed(3), "qcBalance"], ["Time (s)" + (r.target_cached ? ", target digest reused" : ""), r.seconds]]) +
      verdict(sharedClass, `${n(r.shared)} of ${n(r.decoy.distinct)} distinct decoy peptides (${pct(r.shared_fraction)}) are also target peptides. Those can never be counted as decoy hits and bias the FDR estimate${sharedClass === "good" ? "; this is low" : ""}.`) +
      verdict(balClass, `The decoys give ${bal.toFixed(3)} times as many distinct peptides as the targets${balClass === "good" ? " (balanced)" : "; the FDR estimate assumes about 1"}.`) +
      (r.target.failed + r.decoy.failed ? verdict("warn", `${r.target.failed} target and ${r.decoy.failed} decoy sequences could not be digested.`) : "") +
      table([["Database"], ["Proteins", 1], ["Peptides", 1], ["Distinct peptides", 1]], [side("Target", r.target), side("Decoy", r.decoy)]) +
      (r.shared_examples.length ? `<details><summary>Shared peptides (first ${r.shared_examples.length})</summary><div class="details-body"><div class="seq-block">${r.shared_examples.map(esc).join(" ")}</div></div></details>` : "");
    const charts = document.createElement("div");
    charts.className = "charts";
    out.appendChild(charts);
    const T = { name: "Target", color: "var(--target)" }, D = { name: "Decoy", color: "var(--decoy)" };
    lineChart(charts, {
      title: "Distinct peptide length", xTitle: "Peptide length (residues)", yTitle: "Distinct peptides", x: r.length.target.x,
      series: [{ ...T, values: r.length.target.counts }, { ...D, values: r.length.decoy.counts }],
    });
    const m = r.mass.target;
    lineChart(charts, {
      title: "Distinct peptide mass", xTitle: "Monoisotopic mass (Da)", yTitle: "Distinct peptides",
      x: (m.edges.length ? m.edges.slice(0, -1) : []).map((e) => Math.round(e)),
      series: [{ ...T, values: m.counts }, { ...D, values: r.mass.decoy.counts }],
    });
  })()));

  // ---------------------------------------------------------------- 6. export
  document.querySelectorAll('input[name="xWhich"]').forEach((r) => r.addEventListener("change", refreshButtons));
  $("exportBtn").addEventListener("click", () => quiet((async () => {
    const which = document.querySelector('input[name="xWhich"]:checked').value;
    const format = $("xFormat").value;
    const gzip = $("xGzip").checked;
    const bytes = await run("Exporting…", "export", { which, format, gzip, line_width: num("xLineWidth"), peff_prefix: $("xPeffPrefix").value.trim() }, [], $("exportResult"));
    const name = `${state.baseName}${which === "decoy" ? decoySuffix() : ""}.${format}${gzip ? ".gz" : ""}`;
    download(bytes, name, gzip ? "application/gzip" : "text/plain");
    $("exportResult").innerHTML = alert("success", "Saved", `<code>${esc(name)}</code> (${(bytes.length / 1e6).toFixed(2)} MB).`, "", "exportDone");
  })()));

  // ---------------------------------------------------------------- boot
  (async () => {
    try {
      const info = await call("init");
      const sel = $("dEnzyme");
      for (const e of info.enzymes) sel.add(new Option(e, e, e === "trypsin", e === "trypsin"));
      $("versions").textContent = Object.entries(info.versions).map(([k, v]) => `${k} ${v}`).join(" · ");
      toolkit.ready = true;
      toolkit.versions = info.versions;
      setState("ready");
      setStatus(`Ready: fastatacular ${info.versions.fastatacular} running in your browser. Load a FASTA file or an example.`);
    } catch (err) {
      setState("error");
      setStatus("Could not start Python: " + err.message + ". Check the network connection and reload the page.");
      toolkit.errors.push({ cmd: "init", message: err.message });
    }
    refreshButtons();
  })();
  refreshButtons();
})();
