// UI thread: forms, rendering and downloads. All Python runs in worker.js.
(function () {
  const $ = (id) => document.getElementById(id);
  const { barChart, lineChart, fmt } = window.Charts;
  const EXAMPLES = {
    loadExample1: ["examples/human_ecoli_mix.fasta"],
    loadExample2: ["examples/contaminants_with_decoys.fasta"],
    loadBothExamples: ["examples/human_ecoli_mix.fasta", "examples/contaminants_with_decoys.fasta"],
  };
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
  worker.onerror = (ev) => setStatus("Worker failed: " + (ev.message || ev));
  function call(cmd, args, transfer = []) {
    const id = nextId++;
    return new Promise((resolve, reject) => {
      pending.set(id, { resolve, reject });
      worker.postMessage({ id, cmd, args }, transfer);
    });
  }

  // ---------------------------------------------------------------- status
  function setStatus(text) { $("status").textContent = text; }
  function setProgress(stage, done, total) {
    const pct = total ? Math.round((100 * done) / total) : 0;
    $("progressBar").style.width = pct + "%";
    $("progressText").textContent = `${stage}: ${done.toLocaleString()} / ${total.toLocaleString()}`;
  }
  const state = { loaded: 0, entries: 0, decoys: false, model: false, dropped: 0, baseName: "database" };

  function refreshButtons() {
    const idle = toolkit.ready && !toolkit.busy;
    const has = idle && state.entries > 0;
    for (const id of ["loadExample1", "loadExample2", "loadBothExamples"]) $(id).disabled = !idle;
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
  }

  async function run(label, cmd, args, transfer, target) {
    toolkit.busy = true;
    document.body.classList.add("busy");
    refreshButtons();
    setStatus(label);
    $("progressBar").style.width = "0";
    $("progressText").textContent = "";
    const t0 = performance.now();
    try {
      const result = await call(cmd, args, transfer);
      toolkit.results[cmd] = result;
      const secs = ((performance.now() - t0) / 1000).toFixed(1);
      setStatus(`Ready. ${label.replace(/\.+$/, "")} took ${secs} s.`);
      $("progressBar").style.width = "100%";
      return result;
    } catch (err) {
      toolkit.errors.push({ cmd, message: err.message });
      setStatus("Ready (last step failed).");
      if (target) {
        target.innerHTML = `<p class="error"></p>`;
        target.querySelector(".error").textContent = err.message;
      }
      throw err;
    } finally {
      toolkit.busy = false;
      document.body.classList.remove("busy");
      refreshButtons();
    }
  }
  const quiet = (p) => p.catch(() => {});

  // ---------------------------------------------------------------- helpers
  const esc = (s) => String(s ?? "").replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]);
  const tiles = (items) => `<div class="tiles">${items.map(([k, v, id]) => `<div class="tile"${id ? ` id="${id}"` : ""}><div class="v">${esc(v)}</div><div class="k">${esc(k)}</div></div>`).join("")}</div>`;
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
  function showTab(name) {
    document.querySelectorAll(".tabs button").forEach((b) => b.setAttribute("aria-selected", String(b.dataset.tab === name)));
    document.querySelectorAll(".panel").forEach((p) => (p.hidden = p.dataset.panel !== name));
    const ds = $("digestSettings");
    if (name === "stats" || name === "qc") {
      const panel = document.querySelector(`.panel[data-panel="${name}"]`);
      panel.insertBefore(ds, panel.querySelector(":scope > .row"));
      ds.hidden = false;
    } else ds.hidden = true;
  }
  document.querySelectorAll(".tabs button").forEach((b) => b.addEventListener("click", () => showTab(b.dataset.tab)));
  showTab("input");

  // ---------------------------------------------------------------- 1. input
  function renderSummary(s) {
    state.loaded = s.loaded;
    state.entries = s.entries;
    state.decoys = false;
    const files = s.files.map((f) => `<tr><td>${esc(f.name)}</td><td class="num">${f.entries.toLocaleString()}</td><td class="num">${(f.bytes / 1e6).toFixed(2)} MB</td><td class="num">${f.seconds} s</td></tr>`).join("");
    $("loadSummary").innerHTML =
      tiles([["entries", s.entries.toLocaleString(), "loadEntries"], ["residues", s.residues.toLocaleString()], ["files", s.files.length]]) +
      (s.files.length ? `<table><thead><tr><th>file</th><th class="num">entries</th><th class="num">size</th><th class="num">parse time</th></tr></thead><tbody>${files}</tbody></table>` : "") +
      (Object.keys(s.decoys_detected).length ? `<p class="verdict warn">Existing decoys found (${prefixTable(s.decoys_detected)}). The Decoys step leaves them out.</p>` : "") +
      (s.preview.length ? `<p class="muted">First headers:</p><div class="seq">${s.preview.map(esc).join("<br>")}</div>` : "");
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
    const s = await run("Reading " + files.map((f) => f.name).join(", ") + "...", "load", { files: payload, replace }, payload.map((p) => p.buffer), $("loadSummary"));
    renderSummary(s);
    refreshButtons();
  }
  $("fileInput").addEventListener("change", (e) => { quiet(loadFiles([...e.target.files])); e.target.value = ""; });
  const dz = $("dropZone");
  dz.addEventListener("dragover", (e) => { e.preventDefault(); dz.classList.add("over"); });
  dz.addEventListener("dragleave", () => dz.classList.remove("over"));
  dz.addEventListener("drop", (e) => { e.preventDefault(); dz.classList.remove("over"); if (toolkit.ready && !toolkit.busy) quiet(loadFiles([...e.dataTransfer.files])); });
  for (const [id, urls] of Object.entries(EXAMPLES)) {
    $(id).addEventListener("click", async () => {
      const files = [];
      for (const u of urls) files.push({ name: u.split("/").pop(), buffer: await (await fetch(u)).arrayBuffer() });
      if (id === "loadBothExamples") $("appendFiles").checked = false;
      quiet(loadFiles(files));
    });
  }
  $("clearBtn").addEventListener("click", async () => { renderSummary(await run("Clearing...", "clear")); $("loadSummary").innerHTML = ""; state.model = false; refreshButtons(); });

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
    const r = await run("Cleaning up...", "cleanup", cleanOpts(), [], $("cleanResult"));
    state.entries = r.after; state.dropped = r.dropped; state.decoys = false;
    const reasons = Object.entries(r.reasons).map(([k, v]) => `<tr><td>${esc(k)}</td><td class="num">${v}</td></tr>`).join("");
    const rows = r.rows.map((x) => `<tr><td>${esc(x[0])}</td><td>${esc(x[1])}</td><td>${esc(x[2])}</td></tr>`).join("");
    $("cleanResult").innerHTML =
      tiles([["before", r.before.toLocaleString()], ["kept", r.after.toLocaleString(), "cleanKept"], ["dropped", r.dropped.toLocaleString(), "cleanDropped"]]) +
      (r.dropped ? `<table><thead><tr><th>reason</th><th class="num">entries</th></tr></thead><tbody>${reasons}</tbody></table>
        <p class="muted">Dropped entries${r.dropped > r.rows.length ? ` (first ${r.rows.length}; download the CSV for all)` : ""}:</p>
        <div class="scroll"><table id="droppedTable"><thead><tr><th>identifier</th><th>reason</th><th>detail</th></tr></thead><tbody>${rows}</tbody></table></div>` : "<p>Nothing dropped.</p>");
    refreshButtons();
  })()));
  $("cleanResetBtn").addEventListener("click", () => quiet((async () => {
    const s = await run("Undoing clean-up...", "resetCleanup");
    state.entries = s.entries; state.dropped = 0; state.decoys = false;
    $("cleanResult").innerHTML = `<p>Using all ${s.entries.toLocaleString()} loaded entries.</p>`;
    refreshButtons();
  })()));
  $("dlDropped").addEventListener("click", () => quiet((async () => {
    download(await run("Exporting...", "export", { which: "dropped" }), `${state.baseName}_dropped.csv`, "text/csv");
  })()));

  // ---------------------------------------------------------------- 3. stats
  async function doStats(digest) {
    const r = await run(digest ? "Computing stats and digesting..." : "Computing stats...", "stats", { digest, ...digestOpts() }, [], $("statsResult"));
    const out = $("statsResult");
    const L = r.length;
    const items = [["entries", r.entries.toLocaleString(), "statEntries"], ["residues", r.residues.toLocaleString()],
      ["length min / median / max", `${L.min} / ${L.median} / ${L.max}`], ["mean length", L.mean],
      ["organisms", r.organism_count], ["duplicate IDs", r.duplicate_ids], ["duplicate sequences", r.duplicate_sequences],
      ["entries with non-standard residues", r.entries_with_nonstandard]];
    if (r.peptides) {
      const d = digestOpts();
      items.push([`peptides (${d.enzyme}, ${d.missed} missed, ${d.min_len}-${d.max_len} aa)`, r.peptides.total.toLocaleString(), "statPeptides"]);
      items.push(["distinct peptides", r.peptides.distinct.toLocaleString(), "statDistinct"]);
    }
    out.innerHTML = tiles(items) +
      (Object.keys(r.decoys_detected).length ? `<p class="verdict warn">Contains decoy entries: ${prefixTable(r.decoys_detected)}.</p>` : "") +
      (r.peptides && r.peptides.failed ? `<p class="verdict warn">${r.peptides.failed} sequences could not be digested.</p>` : "");
    const charts = document.createElement("div");
    charts.className = "charts";
    out.appendChild(charts);
    const h = L.hist;
    barChart(charts, {
      title: "Protein length", xTitle: "length (residues)", yTitle: "entries",
      labels: h.counts.map((_, i) => (i === h.counts.length - 1 && h.overflow ? `${Math.round(h.edges[i])}+` : `${Math.round(h.edges[i])}`)),
      values: h.counts,
    });
    barChart(charts, {
      title: "Amino-acid composition", xTitle: "residue", yTitle: "% of residues", labelEvery: 1,
      labels: r.composition.map((c) => c.aa), values: r.composition.map((c) => 100 * c.fraction), valueFmt: (v) => v.toFixed(2) + "%",
    });
    if (r.other_residues.length) out.insertAdjacentHTML("beforeend", `<p class="muted">Other characters: ${r.other_residues.map((o) => `<code>${esc(o.aa)}</code> ${o.count}`).join(", ")}</p>`);
    const orgs = r.organisms.map((o) => `<tr><td>${esc(o.name)}</td><td class="num">${o.count.toLocaleString()}</td><td class="num">${pct(o.count / r.entries, 1)}</td></tr>`).join("");
    out.insertAdjacentHTML("beforeend", `<h4>Organisms${r.organism_count > 20 ? " (top 20)" : ""}</h4><table id="orgTable"><thead><tr><th>OS= (OX=)</th><th class="num">entries</th><th class="num">share</th></tr></thead><tbody>${orgs}</tbody></table>`);
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
    const r = await run("Training Markov model...", "trainModel", { order: num("mOrder") }, [], $("modelInfo"));
    $("mModel").value = "trained";
    modelInfo(r, "Trained");
  })()));
  $("modelFile").addEventListener("change", (e) => quiet((async () => {
    const f = e.target.files[0];
    if (!f) return;
    const buffer = await f.arrayBuffer();
    const r = await run("Loading model...", "loadModel", { file: { name: f.name, buffer } }, [buffer], $("modelInfo"));
    $("mModel").value = "uploaded";
    modelInfo(r, "Uploaded");
    e.target.value = "";
  })()));
  $("dlModel").addEventListener("click", () => quiet((async () => {
    download(await run("Exporting...", "export", { which: "model" }), `${state.baseName}_markov.json.gz`, "application/gzip");
  })()));
  function decoyOpts() {
    return {
      method: $("mMethod").value, prefix: $("mPrefix").value.trim(), seed: $("mSeed").value,
      keep_residues: $("mKeepRes").value.trim(), keep_met: $("mKeepMet").checked, keep_cterm: num("mKeepCterm"),
      k: num("mK") || 2, model: $("mModel").value,
      concatenate: document.querySelector('input[name="mOutput"]:checked').value === "concat",
    };
  }
  $("decoyBtn").addEventListener("click", () => quiet((async () => {
    const o = decoyOpts();
    const r = await run(`Making ${o.method} decoys...`, "decoys", o, [], $("decoyResult"));
    state.decoys = true;
    const ex = r.examples.map((e) => `<tr><td class="seq">${esc(e.header)}<br>T: ${esc(e.target)}<br>D: ${esc(e.decoy)}</td></tr>`).join("");
    $("decoyResult").innerHTML =
      tiles([["targets", r.targets.toLocaleString(), "decoyTargets"], ["decoys", r.decoys.toLocaleString(), "decoyCount"],
        ["entries in output", r.output_entries.toLocaleString(), "decoyOutput"], ["decoy = target", r.identical_to_target],
        ["composition L1", r.composition_l1], ["time", r.seconds + " s"]]) +
      (r.existing_decoys_skipped ? `<p class="verdict warn" id="decoySkipped">Left out ${r.existing_decoys_skipped} existing decoy entries (${prefixTable(r.existing_prefixes)}); they are not in the output.</p>` : `<p class="verdict good">No existing decoys in the input.</p>`) +
      (r.identical_to_target ? `<p class="verdict warn">${r.identical_to_target} decoys equal their target (every stretch between kept residues is too short to change).</p>` : "") +
      `<p class="muted">First decoys (T = target, D = decoy, first 60 residues):</p><table><tbody>${ex}</tbody></table>
       <p>Next: <a href="#" id="toQc">check the decoys in Decoy QC</a>.</p>`;
    $("toQc").addEventListener("click", (e) => { e.preventDefault(); showTab("qc"); });
    refreshButtons();
  })()));
  $("decoyDlBtn").addEventListener("click", () => quiet((async () => {
    const bytes = await run("Exporting...", "export", { which: "decoy", format: "fasta", line_width: num("xLineWidth") });
    download(bytes, `${state.baseName}_${decoyOpts().concatenate ? "target_decoy" : "decoy"}.fasta`, "text/plain");
  })()));

  // ---------------------------------------------------------------- 5. QC
  $("qcBtn").addEventListener("click", () => quiet((async () => {
    const r = await run("Digesting targets and decoys...", "qc", digestOpts(), [], $("qcResult"));
    const out = $("qcResult");
    const sharedClass = r.shared_fraction < 0.005 ? "good" : r.shared_fraction < 0.02 ? "warn" : "bad";
    const bal = r.balance;
    const balClass = Math.abs(bal - 1) <= 0.05 ? "good" : Math.abs(bal - 1) <= 0.15 ? "warn" : "bad";
    out.innerHTML =
      tiles([["shared decoy peptides", pct(r.shared_fraction), "qcShared"], ["shared if I = L", pct(r.shared_il_fraction), "qcSharedIL"],
        ["distinct target peptides", r.target.distinct.toLocaleString(), "qcTarget"], ["distinct decoy peptides", r.decoy.distinct.toLocaleString(), "qcDecoy"],
        ["decoy / target peptides", bal.toFixed(3), "qcBalance"], ["time", r.seconds + " s"]]) +
      `<p class="verdict ${sharedClass}">${r.shared.toLocaleString()} of ${r.decoy.distinct.toLocaleString()} distinct decoy peptides (${pct(r.shared_fraction)}) are also target peptides. Those can never be counted as decoy hits and bias the FDR estimate${sharedClass === "good" ? "; this is low" : ""}.</p>` +
      `<p class="verdict ${balClass}">The decoys give ${bal.toFixed(3)} times as many distinct peptides as the targets${balClass === "good" ? " (balanced)" : "; the FDR estimate assumes about 1"}.</p>` +
      (r.target.failed + r.decoy.failed ? `<p class="verdict warn">${r.target.failed} target and ${r.decoy.failed} decoy sequences could not be digested.</p>` : "") +
      `<table><thead><tr><th></th><th class="num">proteins</th><th class="num">peptides</th><th class="num">distinct</th></tr></thead><tbody>
        <tr><td>target</td><td class="num">${r.target.proteins.toLocaleString()}</td><td class="num">${r.target.peptides.toLocaleString()}</td><td class="num">${r.target.distinct.toLocaleString()}</td></tr>
        <tr><td>decoy</td><td class="num">${r.decoy.proteins.toLocaleString()}</td><td class="num">${r.decoy.peptides.toLocaleString()}</td><td class="num">${r.decoy.distinct.toLocaleString()}</td></tr></tbody></table>` +
      (r.shared_examples.length ? `<details><summary>Shared peptides (first ${r.shared_examples.length})</summary><div class="seq">${r.shared_examples.map(esc).join(" ")}</div></details>` : "");
    const charts = document.createElement("div");
    charts.className = "charts";
    out.appendChild(charts);
    const T = { name: "target", color: "var(--target)" }, D = { name: "decoy", color: "var(--decoy)" };
    lineChart(charts, {
      title: "Distinct peptide length", xTitle: "length", yTitle: "peptides", x: r.length.target.x,
      series: [{ ...T, values: r.length.target.counts }, { ...D, values: r.length.decoy.counts }],
    });
    const m = r.mass.target;
    lineChart(charts, {
      title: "Distinct peptide mass (monoisotopic)", xTitle: "Da", yTitle: "peptides",
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
    const bytes = await run("Exporting...", "export", { which, format, gzip, line_width: num("xLineWidth"), peff_prefix: $("xPeffPrefix").value.trim() }, [], $("exportResult"));
    const name = `${state.baseName}${which === "decoy" ? "_target_decoy" : ""}.${format}${gzip ? ".gz" : ""}`;
    download(bytes, name, gzip ? "application/gzip" : "text/plain");
    $("exportResult").innerHTML = `<p id="exportDone">Saved <code>${esc(name)}</code> (${(bytes.length / 1e6).toFixed(2)} MB).</p>`;
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
      setStatus("Ready. Load a FASTA file or an example.");
    } catch (err) {
      setStatus("Could not start Python: " + err.message);
      toolkit.errors.push({ cmd: "init", message: err.message });
    }
    refreshButtons();
  })();
  refreshButtons();
})();
