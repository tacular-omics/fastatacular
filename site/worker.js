// Web Worker: runs Pyodide + fastatacular so the page never freezes on big files.
// Protocol: main -> {id, cmd, args}; worker -> {id, ok, result | error} and
// {type: "progress", stage, done, total} / {type: "status", text}.
const PYODIDE_VERSION = "0.29.5";
const PACKAGES = ["fastatacular>=1.1,<2", "peptacular>=5,<6", "tacular>=2,<3", "pefftacular>=1.1,<2"];

importScripts(`https://cdn.jsdelivr.net/pyodide/v${PYODIDE_VERSION}/full/pyodide.js`);

let pyodide = null;
let tk = null;

function status(text) {
  postMessage({ type: "status", text });
}

self.reportProgress = (stage, done, total) => postMessage({ type: "progress", stage, done, total });

async function boot() {
  status("Loading Python (Pyodide)...");
  pyodide = await loadPyodide({ indexURL: `https://cdn.jsdelivr.net/pyodide/v${PYODIDE_VERSION}/full/` });
  status("Installing fastatacular, peptacular, pefftacular...");
  await pyodide.loadPackage(["micropip", "lzma"]); // lzma: .xz input
  const micropip = pyodide.pyimport("micropip");
  await micropip.install(PACKAGES);
  const src = await (await fetch("toolkit.py", { cache: "no-cache" })).text();
  pyodide.FS.writeFile("/home/pyodide/toolkit.py", src);
  tk = pyodide.pyimport("toolkit");
  tk.set_progress(self.reportProgress);
  const versions = JSON.parse(
    pyodide.runPython(`
import json, importlib.metadata as md, sys
json.dumps({p: md.version(p) for p in ("fastatacular", "peptacular", "tacular", "pefftacular")} | {"python": sys.version.split()[0]})
`),
  );
  versions.pyodide = pyodide.version;
  return { versions, enzymes: JSON.parse(tk.enzymes()) };
}

const ready = boot();

function writeFiles(files) {
  pyodide.FS.mkdirTree("/tmp/in");
  return files.map((f, i) => {
    const path = `/tmp/in/${i}_${f.name.replace(/[^\w.-]/g, "_")}`;
    pyodide.FS.writeFile(path, new Uint8Array(f.buffer));
    return path;
  });
}

const COMMANDS = {
  init: async () => ready,
  load: ({ files, replace }) => {
    const paths = writeFiles(files);
    return JSON.parse(tk.load_files(paths, files.map((f) => f.name), replace));
  },
  clear: () => JSON.parse(tk.clear()),
  cleanup: (opts) => JSON.parse(tk.cleanup(JSON.stringify(opts))),
  resetCleanup: () => JSON.parse(tk.reset_cleanup()),
  stats: (opts) => JSON.parse(tk.stats(JSON.stringify(opts))),
  trainModel: ({ order }) => JSON.parse(tk.train_model(order)),
  loadModel: ({ file }) => JSON.parse(tk.load_model_file(writeFiles([file])[0])),
  decoys: (opts) => JSON.parse(tk.decoys(JSON.stringify(opts))),
  qc: (opts) => JSON.parse(tk.decoy_qc(JSON.stringify(opts))),
  export: (opts) => {
    const proxy = tk.export(JSON.stringify(opts));
    const bytes = proxy.toJs();
    proxy.destroy();
    return bytes;
  },
};

self.onmessage = async (ev) => {
  const { id, cmd, args } = ev.data;
  try {
    await ready;
    const result = await COMMANDS[cmd](args || {});
    const transfer = result instanceof Uint8Array ? [result.buffer] : [];
    postMessage({ id, ok: true, result }, transfer);
  } catch (err) {
    // Python errors: keep the last line (the exception) as the message.
    const text = String(err && err.message ? err.message : err).trim().split("\n");
    postMessage({ id, ok: false, error: text[text.length - 1], detail: text.slice(-6).join("\n") });
  }
};
