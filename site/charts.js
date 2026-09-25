// Small dependency-free SVG charts: single-series bars and multi-series lines.
// Colors are CSS variables so light/dark themes swap in one place.
(function () {
  const NS = "http://www.w3.org/2000/svg";
  const W = 560, H = 240, M = { l: 52, r: 16, t: 12, b: 40 };
  const tip = () => document.getElementById("tooltip");

  function el(name, attrs, parent) {
    const node = document.createElementNS(NS, name);
    for (const [k, v] of Object.entries(attrs || {})) node.setAttribute(k, v);
    if (parent) parent.appendChild(node);
    return node;
  }
  const fmt = (v) => (Math.abs(v) >= 1000 ? Math.round(v).toLocaleString() : +(+v).toPrecision(3) + "");

  function niceTicks(max, n = 4) {
    if (max <= 0) return [0];
    const raw = max / n, mag = 10 ** Math.floor(Math.log10(raw));
    const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= raw);
    const out = [];
    for (let v = 0; v <= max + 1e-9; v += step) out.push(v);
    if (out[out.length - 1] < max) out.push(out[out.length - 1] + step);
    return out;
  }

  function frame(container, title, xTitle, yTitle, ymax) {
    const box = document.createElement("div");
    box.className = "chart";
    box.innerHTML = `<h4></h4>`;
    box.querySelector("h4").textContent = title;
    container.appendChild(box);
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": title }, box);
    const ticks = niceTicks(ymax);
    const top = ticks[ticks.length - 1] || 1;
    const y = (v) => M.t + (H - M.t - M.b) * (1 - v / top);
    const ax = el("g", { class: "axis" }, svg);
    for (const t of ticks) {
      el("line", { x1: M.l, x2: W - M.r, y1: y(t), y2: y(t) }, ax);
      el("text", { x: M.l - 6, y: y(t) + 4, "text-anchor": "end" }, ax).textContent = fmt(t);
    }
    el("text", { x: (M.l + W - M.r) / 2, y: H - 4, "text-anchor": "middle" }, ax).textContent = xTitle;
    const yl = el("text", { x: 12, y: (M.t + H - M.b) / 2, "text-anchor": "middle", transform: `rotate(-90 12 ${(M.t + H - M.b) / 2})` }, ax);
    yl.textContent = yTitle;
    return { box, svg, y, ax };
  }

  function showTip(ev, html) {
    const t = tip();
    t.innerHTML = html;
    t.hidden = false;
    const x = Math.min(ev.clientX + 12, window.innerWidth - t.offsetWidth - 8);
    t.style.left = x + "px";
    t.style.top = ev.clientY + 12 + "px";
  }
  const hideTip = () => { tip().hidden = true; };

  function table(box, headers, rows) {
    const d = document.createElement("details");
    d.innerHTML = "<summary>Show table</summary>";
    const wrap = document.createElement("div");
    wrap.className = "details-body";
    const scroll = document.createElement("div");
    scroll.className = "table-wrap";
    wrap.appendChild(scroll);
    const tbl = document.createElement("table");
    tbl.className = "data";
    tbl.innerHTML = "<thead><tr>" + headers.map((h, i) => `<th class="${i ? "num" : ""}">${h}</th>`).join("") + "</tr></thead>";
    const tb = document.createElement("tbody");
    for (const r of rows) {
      const tr = document.createElement("tr");
      r.forEach((c, i) => {
        const td = document.createElement("td");
        if (i) td.className = "num";
        td.textContent = typeof c === "number" ? fmt(c) : c;
        tr.appendChild(td);
      });
      tb.appendChild(tr);
    }
    tbl.appendChild(tb);
    scroll.appendChild(tbl);
    d.appendChild(wrap);
    box.appendChild(d);
  }

  // labels: category labels; values: numbers.
  function barChart(container, { title, labels, values, xTitle, yTitle, color = "var(--target)", valueFmt = fmt, labelEvery }) {
    const ymax = Math.max(...values, 0);
    const { box, svg, y, ax } = frame(container, title, xTitle, yTitle, ymax);
    const n = values.length;
    const slot = (W - M.l - M.r) / n;
    const gap = n > 30 ? 1 : 2;
    const every = labelEvery || Math.ceil(n / 10);
    values.forEach((v, i) => {
      const x = M.l + i * slot;
      const h = y(0) - y(v);
      const g = el("g", {}, svg);
      el("rect", { x: x, y: M.t, width: slot, height: H - M.t - M.b, fill: "transparent" }, g);
      if (v > 0) {
        el("rect", { x: x + gap / 2, y: y(v), width: Math.max(1, slot - gap), height: h, fill: color }, g);
      }
      g.addEventListener("mousemove", (ev) => showTip(ev, `<b>${labels[i]}</b><br>${valueFmt(v)}`));
      g.addEventListener("mouseleave", hideTip);
      if (i % every === 0) el("text", { x: x + slot / 2, y: H - M.b + 14, "text-anchor": "middle" }, ax).textContent = labels[i];
    });
    table(box, [xTitle, yTitle], labels.map((l, i) => [String(l), values[i]]));
    return box;
  }

  // x: labels for each point; series: [{name, color, values}]
  function lineChart(container, { title, x, series, xTitle, yTitle }) {
    const ymax = Math.max(0, ...series.flatMap((s) => s.values));
    const { box, svg, y, ax } = frame(container, title, xTitle, yTitle, ymax);
    const n = x.length;
    const px = (i) => M.l + (n === 1 ? (W - M.l - M.r) / 2 : (i * (W - M.l - M.r)) / (n - 1));
    const every = Math.ceil(n / 10);
    x.forEach((lab, i) => { if (i % every === 0) el("text", { x: px(i), y: H - M.b + 14, "text-anchor": "middle" }, ax).textContent = lab; });
    for (const s of series) {
      const d = s.values.map((v, i) => `${i ? "L" : "M"}${px(i).toFixed(1)},${y(v).toFixed(1)}`).join("");
      el("path", { d, fill: "none", stroke: s.color, "stroke-width": 2, "stroke-linejoin": "round" }, svg);
    }
    const cross = el("line", { y1: M.t, y2: H - M.b, stroke: "var(--text-muted)", "stroke-width": 1, visibility: "hidden" }, svg);
    const dots = series.map((s) => el("circle", { r: 4, fill: s.color, stroke: "var(--surface)", "stroke-width": 2, visibility: "hidden" }, svg));
    const hit = el("rect", { x: M.l, y: M.t, width: W - M.l - M.r, height: H - M.t - M.b, fill: "transparent" }, svg);
    hit.addEventListener("mousemove", (ev) => {
      const r = svg.getBoundingClientRect();
      const sx = ((ev.clientX - r.left) / r.width) * W;
      const i = Math.max(0, Math.min(n - 1, Math.round(((sx - M.l) / (W - M.l - M.r)) * (n - 1))));
      cross.setAttribute("x1", px(i)); cross.setAttribute("x2", px(i)); cross.setAttribute("visibility", "visible");
      dots.forEach((d, k) => { d.setAttribute("cx", px(i)); d.setAttribute("cy", y(series[k].values[i])); d.setAttribute("visibility", "visible"); });
      showTip(ev, `<b>${xTitle} ${x[i]}</b><br>` + series.map((s) => `${s.name}: ${fmt(s.values[i])}`).join("<br>"));
    });
    hit.addEventListener("mouseleave", () => { hideTip(); cross.setAttribute("visibility", "hidden"); dots.forEach((d) => d.setAttribute("visibility", "hidden")); });
    const legend = document.createElement("div");
    legend.className = "legend";
    legend.innerHTML = series.map((s) => `<span><i style="background:${s.color}"></i>${s.name}</span>`).join("");
    box.insertBefore(legend, svg);
    table(box, [xTitle, ...series.map((s) => s.name)], x.map((l, i) => [String(l), ...series.map((s) => s.values[i])]));
    return box;
  }

  window.Charts = { barChart, lineChart, fmt };
})();
