"use strict";
/* LARMOR app — state, plot, sites, fit, quantify, processing, figures. */

const $ = id => document.getElementById(id);
const COLORS = ["#1f77b4", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
                "#e377c2", "#17becf", "#bcbd22", "#7f7f7f", "#ff7f0e"];

const state = {
  sourcePath: null,
  recipe: null,
  exp: null,              // {ppm: [], amp: []}
  models: [],             // registry from /api/models
  addModel: null,         // model name while in click-to-add mode
  lastFit: null,          // last /api/fit response
  undo: [], redo: [],
  hidden: new Set(),      // site indices hidden from the plot
  advOpen: new Set(),     // site indices with the constraints rows open
};

/* ---------------- helpers ---------------- */
async function api(endpoint, body) {
  const r = await fetch("/api/" + endpoint, {
    method: "POST", headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
  if (!r.ok) {
    let detail = r.statusText;
    try { detail = (await r.json()).detail; } catch (e) {}
    throw new Error(detail);
  }
  return r.json();
}
async function apiGet(endpoint) {
  const r = await fetch("/api/" + endpoint);
  if (!r.ok) throw new Error(r.statusText);
  return r.json();
}
function status(msg, kind) {
  const el = $("status");
  el.textContent = msg;
  el.className = kind || "";
}
const debounce = (fn, ms) => {
  let t; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); };
};
function fmt(v, digits = 6) {
  if (v === null || v === undefined || !isFinite(v)) return "";
  return Number(Number(v).toPrecision(digits));
}

/* ---------------- undo / redo ---------------- */
function snapshot() {
  if (!state.recipe) return;
  state.undo.push(JSON.stringify(state.recipe));
  if (state.undo.length > 60) state.undo.shift();
  state.redo.length = 0;
  updateUndoButtons();
  persist();
}
function undo() {
  if (!state.undo.length) return;
  state.redo.push(JSON.stringify(state.recipe));
  state.recipe = JSON.parse(state.undo.pop());
  renderSites(); simulate(); updateUndoButtons();
}
function redo() {
  if (!state.redo.length) return;
  state.undo.push(JSON.stringify(state.recipe));
  state.recipe = JSON.parse(state.redo.pop());
  renderSites(); simulate(); updateUndoButtons();
}
function updateUndoButtons() {
  $("btnUndo").disabled = !state.undo.length;
  $("btnRedo").disabled = !state.redo.length;
}

/* ---------------- session persistence ---------------- */
const persist = debounce(() => {
  if (!state.sourcePath) return;
  try {
    localStorage.setItem("larmor.session", JSON.stringify({
      sourcePath: state.sourcePath, recipe: state.recipe,
      winHi: $("winHi").value, winLo: $("winLo").value,
    }));
  } catch (e) {}
}, 400);

async function tryRestore() {
  let s;
  try { s = JSON.parse(localStorage.getItem("larmor.session")); } catch (e) { return; }
  if (!s || !s.sourcePath) return;
  status("restoring last session…", "busy");
  try {
    await loadSource(s.sourcePath);
    if (s.recipe) { state.recipe = s.recipe; renderSites(); simulate(); }
    if (s.winHi) $("winHi").value = s.winHi;
    if (s.winLo) $("winLo").value = s.winLo;
    status("session restored");
  } catch (e) { status(""); }
}

/* ---------------- loading ---------------- */
async function loadSource(path) {
  status("loading…", "busy");
  const data = await api("load", { path });
  state.sourcePath = path;
  state.recipe = data.recipe;
  state.exp = { ppm: data.ppm, amp: data.amp };
  state.undo = []; state.redo = []; state.hidden.clear();
  $("srcline").textContent = data.meta + (data.warnings.length ? "  ⚠ " + data.warnings.join(" • ") : "");
  $("srcline").title = path;
  $("winHi").value = Math.max(...data.ppm).toFixed(1);
  $("winLo").value = Math.min(...data.ppm).toFixed(1);
  ["btnSave", "btnFit", "btnQuant", "figTplLoad"].forEach(id => $(id).disabled = false);
  renderSites();
  if (state.recipe.sites.length) simulate(); else { drawPlot(null); status("no sites — pick a model below the plot, then click on the spectrum"); }
  persist();
}

/* ---------------- plot ---------------- */
function siteColor(i) { return COLORS[i % COLORS.length]; }

function drawPlot(model) {
  const traces = [{
    x: state.exp.ppm, y: state.exp.amp, name: "experiment",
    line: { color: "#16202a", width: 1.1 }, hoverinfo: "x+y",
  }];
  if (model) {
    traces.push({ x: model.x, y: model.total, name: "model",
                  line: { color: "#d62728", width: 1.5 } });
    if ($("showComp").checked) {
      model.sites.forEach((ys, i) => {
        if (state.hidden.has(i)) return;
        traces.push({ x: model.x, y: ys, name: `s${i} ${model.labels[i]}`,
                      line: { width: 1, dash: "dot", color: siteColor(i) }, opacity: 0.85 });
      });
    }
    if ($("showResid").checked && state.exp) {
      const yi = interp(state.exp.ppm, model.x, model.total);
      const resid = state.exp.amp.map((v, k) => v - yi[k]);
      const offset = -0.08 * Math.max(...state.exp.amp);
      traces.push({ x: state.exp.ppm, y: resid.map(v => v + offset),
                    name: "residual", line: { color: "#8a969e", width: 0.8 } });
    }
  }
  // draggable position markers, one dashed vertical line per visible site
  const shapes = [];
  if (state.recipe) {
    state.recipe.sites.forEach((s, i) => {
      if (state.hidden.has(i)) return;
      const p = s.params.isotropic_chemical_shift_ppm;
      if (!p || p.expr) return;               // linked positions are not draggable
      shapes.push({
        type: "line", x0: p.value, x1: p.value, y0: 0, y1: 1,
        yref: "paper", editable: true,
        line: { color: siteColor(i), width: 1.4, dash: "dash" },
      });
    });
  }
  const layout = {
    margin: { t: 12, r: 14, l: 55, b: 40 },
    xaxis: { title: "shift / ppm", autorange: "reversed", zeroline: false },
    yaxis: { title: "intensity", zeroline: false, exponentformat: "SI" },
    legend: { orientation: "h", y: 1.06, font: { size: 11 } },
    shapes, dragmode: "pan",
    uirevision: state.sourcePath,   // keep zoom across re-simulations
  };
  Plotly.react("plot", traces, layout, {
    responsive: true, displaylogo: false, scrollZoom: true,
    edits: { shapePosition: true },
    modeBarButtonsToRemove: ["lasso2d", "select2d", "autoScale2d"],
  });
}
function interp(xs, xp, yp) {
  // xp ascending; xs any order
  const out = new Array(xs.length);
  for (let k = 0; k < xs.length; k++) {
    const x = xs[k];
    if (x <= xp[0]) { out[k] = yp[0]; continue; }
    if (x >= xp[xp.length - 1]) { out[k] = yp[yp.length - 1]; continue; }
    let lo = 0, hi = xp.length - 1;
    while (hi - lo > 1) { const m = (hi + lo) >> 1; (xp[m] <= x) ? lo = m : hi = m; }
    const t = (x - xp[lo]) / (xp[hi] - xp[lo]);
    out[k] = yp[lo] + t * (yp[hi] - yp[lo]);
  }
  return out;
}

/* live simulation */
const simulate = debounce(async () => {
  if (!state.recipe || !state.recipe.sites.length) { drawPlot(null); return; }
  try {
    const hasKernel = state.recipe.sites.some(s => s.model === "czjzek");
    if (hasKernel) status("simulating…", "busy");
    const model = await api("simulate", { recipe: state.recipe, source_path: state.sourcePath });
    state.lastModel = model;
    drawPlot(model);
    status("");
  } catch (e) { status("simulate: " + e.message, "error"); }
}, 220);

/* plot events: click-to-add + marker drag */
function bindPlotEvents() {
  const plot = $("plot");
  plot.on("plotly_click", ev => {
    if (!state.addModel || !state.recipe) return;
    const pt = ev.points && ev.points[0];
    if (!pt) return;
    addSite(state.addModel, pt.x, Math.abs(pt.y));
    setAddMode(null);
  });
  plot.on("plotly_relayout", ev => {
    // shape drags arrive as {"shapes[N].x0": val, ...}
    let moved = false;
    for (const k of Object.keys(ev)) {
      const m = k.match(/^shapes\[(\d+)\]\.x0$/);
      if (!m) continue;
      const shapeIdx = parseInt(m[1]);
      const visible = state.recipe.sites
        .map((s, i) => ({ s, i }))
        .filter(o => !state.hidden.has(o.i) &&
                     o.s.params.isotropic_chemical_shift_ppm &&
                     !o.s.params.isotropic_chemical_shift_ppm.expr);
      const target = visible[shapeIdx];
      if (target) {
        snapshot();
        target.s.params.isotropic_chemical_shift_ppm.value = ev[k];
        moved = true;
      }
    }
    if (moved) { renderSites(); simulate(); }
  });
}

/* ---------------- add-site mode ---------------- */
function setAddMode(name) {
  state.addModel = name;
  document.querySelectorAll("#addbar button").forEach(b =>
    b.classList.toggle("active", b.dataset.model === name));
  $("plot").style.cursor = name ? "crosshair" : "";
  status(name ? `click on the spectrum to place a ${name} site (Esc to cancel)` : "");
}
function buildAddBar() {
  const bar = $("addbar");
  state.models.forEach(m => {
    const b = document.createElement("button");
    b.textContent = m.label;
    b.title = m.description;
    b.dataset.model = m.name;
    b.onclick = () => setAddMode(state.addModel === m.name ? null : m.name);
    bar.appendChild(b);
  });
}
function addSite(modelName, pos, amp) {
  const m = state.models.find(x => x.name === modelName);
  if (!m) return;
  snapshot();
  const params = {};
  m.params.forEach(p => {
    params[p.name] = { value: p.default, stderr: null, vary: p.vary !== false,
                       min: p.min, max: p.max, expr: null };
  });
  if (params.isotropic_chemical_shift_ppm) params.isotropic_chemical_shift_ppm.value = pos;
  if (params.amplitude) params.amplitude.value = amp || 1.0;
  const n = state.recipe.sites.length;
  state.recipe.sites.push({ model: modelName, label: `${m.label.split(" ")[0]}-${n}`, params });
  renderSites(); simulate();
}

/* ---------------- sites panel ---------------- */
const PARAM_LABELS = {
  isotropic_chemical_shift_ppm: "δiso (ppm)",
  sigma_Cq_MHz: "σ(Cq) (MHz)",
  Cq_MHz: "Cq (MHz)",
  eta: "η",
  zeta_ppm: "ζ CSA (ppm)",
  shift_fwhm_ppm: "FWHM (ppm)",
  amplitude: "amplitude",
  gl: "g/l fraction",
};

function renderSites() {
  const box = $("sites");
  box.innerHTML = "";
  if (!state.recipe) return;
  state.recipe.sites.forEach((site, i) => {
    const card = document.createElement("div");
    card.className = "site" + (state.hidden.has(i) ? " hidden-site" : "");

    const head = document.createElement("div");
    head.className = "site-head";
    const sw = document.createElement("span");
    sw.className = "swatch"; sw.style.background = siteColor(i);
    sw.title = "show/hide on plot";
    sw.onclick = () => { state.hidden.has(i) ? state.hidden.delete(i) : state.hidden.add(i); renderSites(); simulate(); };
    const name = document.createElement("input");
    name.className = "name"; name.value = site.label || `s${i}`;
    name.onchange = () => { snapshot(); site.label = name.value; };
    const tag = document.createElement("span");
    tag.className = "model-tag"; tag.textContent = `s${i} · ${site.model}`;
    const bGear = document.createElement("button");
    bGear.textContent = "⚙"; bGear.title = "constraints: link / min / max";
    bGear.onclick = () => { state.advOpen.has(i) ? state.advOpen.delete(i) : state.advOpen.add(i); renderSites(); };
    const bDup = document.createElement("button");
    bDup.textContent = "⧉"; bDup.title = "duplicate site";
    bDup.onclick = () => {
      snapshot();
      const copy = JSON.parse(JSON.stringify(site));
      copy.label = (site.label || "site") + "-copy";
      Object.values(copy.params).forEach(p => p.stderr = null);
      state.recipe.sites.push(copy);
      renderSites(); simulate();
    };
    const bDel = document.createElement("button");
    bDel.textContent = "✕"; bDel.className = "danger"; bDel.title = "remove site";
    bDel.onclick = () => { snapshot(); state.recipe.sites.splice(i, 1); state.hidden.delete(i); renderSites(); simulate(); };
    head.append(sw, name, tag, bGear, bDup, bDel);
    card.appendChild(head);

    const body = document.createElement("div");
    body.className = "params";
    for (const [pname, p] of Object.entries(site.params)) {
      const row = document.createElement("div");
      row.className = "prow";
      const lab = document.createElement("label");
      lab.innerHTML = (PARAM_LABELS[pname] || pname) +
        (p.expr ? ` <span class="lk" title="linked: ${p.expr}">⚭</span>` : "");
      lab.title = pname;
      const inp = document.createElement("input");
      inp.type = "number"; inp.step = "any";
      inp.value = fmt(p.value);
      inp.disabled = !!p.expr;
      inp.onchange = () => { snapshot(); p.value = parseFloat(inp.value); simulate(); };
      inp.addEventListener("wheel", ev => {
        ev.preventDefault();
        const step = (Math.abs(p.value) || 1) * (ev.shiftKey ? 0.1 : 0.02);
        snapshot();
        p.value = p.value + (ev.deltaY < 0 ? step : -step);
        inp.value = fmt(p.value);
        simulate();
      }, { passive: false });
      const vary = document.createElement("input");
      vary.type = "checkbox";
      vary.checked = p.vary && !p.expr;
      vary.disabled = !!p.expr;
      vary.title = p.expr ? "linked — follows its expression" : "checked = fitted, unchecked = fixed";
      vary.onchange = () => { snapshot(); p.vary = vary.checked; };
      const err = document.createElement("span");
      err.className = "err";
      err.textContent = p.stderr != null ? "± " + fmt(p.stderr, 3) : "";
      row.append(lab, inp, vary, err);
      body.appendChild(row);

      if (state.advOpen.has(i)) {
        const c = document.createElement("div");
        c.className = "constr";
        const ex = document.createElement("input");
        ex.type = "text"; ex.placeholder = `link: 0.5 * s0.${pname}`;
        ex.value = p.expr || "";
        ex.onchange = () => { snapshot(); p.expr = ex.value.trim() || null; renderSites(); simulate(); };
        const lo = document.createElement("input");
        lo.type = "number"; lo.placeholder = "min"; lo.step = "any";
        if (p.min != null) lo.value = p.min;
        lo.onchange = () => { snapshot(); p.min = lo.value === "" ? null : parseFloat(lo.value); };
        const hi = document.createElement("input");
        hi.type = "number"; hi.placeholder = "max"; hi.step = "any";
        if (p.max != null) hi.value = p.max;
        hi.onchange = () => { snapshot(); p.max = hi.value === "" ? null : parseFloat(hi.value); };
        c.append(ex, lo, hi);
        body.appendChild(c);
      }
    }
    card.appendChild(body);
    box.appendChild(card);
  });
}

/* ---------------- fit + quantify ---------------- */
async function runFit() {
  if (!state.recipe || !state.recipe.sites.length) { status("add at least one site first", "error"); return; }
  snapshot();
  status("fitting…", "busy");
  $("btnFit").disabled = true;
  try {
    const win = [parseFloat($("winHi").value), parseFloat($("winLo").value)];
    const res = await api("fit", {
      recipe: state.recipe, source_path: state.sourcePath,
      window: isFinite(win[0]) && isFinite(win[1]) ? win : null,
    });
    state.recipe = res.recipe;
    state.lastFit = res;
    renderSites();
    drawPlot(res);
    const bits = [`RMSD ${res.rmsd.toFixed(4)}`];
    if (res.frozen.length) bits.push(`frozen: ${res.frozen.join(", ")}`);
    if (res.at_bounds.length) bits.push(`⚠ at bounds: ${res.at_bounds.join(", ")}`);
    $("dock-summary").textContent = bits.join("  ·  ");
    $("report").style.display = "block";
    $("report").textContent = res.report;
    openDock();
    status(res.at_bounds.length ? "fit done — some parameters at bounds, check constraints" : "fit done");
    await runQuantify(false);
    persist();
  } catch (e) { status("fit: " + e.message, "error"); }
  finally { $("btnFit").disabled = false; }
}

let lastQuant = null;
async function runQuantify(openIt = true) {
  if (!state.recipe || !state.recipe.sites.length) return;
  try {
    const win = [parseFloat($("winHi").value), parseFloat($("winLo").value)];
    const q = await api("quantify", {
      recipe: state.recipe,
      window: isFinite(win[0]) && isFinite(win[1]) ? win : null,
    });
    lastQuant = q;
    const rows = q.rows.map((r, i) => `
      <tr>
        <td><span style="color:${siteColor(i)}">■</span> ${r.label} <small style="color:var(--faint)">(${r.model})</small></td>
        <td>${fmt(r.position_ppm, 5)}${r.position_err ? " ± " + fmt(r.position_err, 2) : ""}</td>
        <td>${r.integral.toExponential(3)}</td>
        <td class="frac">${r.fraction_pct.toFixed(1)}${r.fraction_err_pct != null ? " ± " + r.fraction_err_pct.toFixed(1) : ""} %</td>
      </tr>`).join("");
    $("qwrap").innerHTML = `
      <table class="qtable">
        <thead><tr><th>site</th><th>position (ppm)</th><th>integral</th><th>fraction</th></tr></thead>
        <tbody>${rows}</tbody>
      </table>
      <div style="font-size:11px;color:var(--faint);padding:4px 0">${q.note} · window ${q.window_ppm[0]} … ${q.window_ppm[1]} ppm</div>`;
    $("btnCsv").disabled = false;
    if (openIt) openDock();
  } catch (e) { status("quantify: " + e.message, "error"); }
}
function quantCsv() {
  if (!lastQuant) return;
  const head = "site,model,position_ppm,position_err,integral,integral_err,fraction_pct,fraction_err_pct";
  const lines = lastQuant.rows.map(r =>
    [r.label, r.model, r.position_ppm, r.position_err ?? "", r.integral,
     r.integral_err ?? "", r.fraction_pct, r.fraction_err_pct ?? ""].join(","));
  navigator.clipboard.writeText([head, ...lines].join("\n"));
  status("quantification table copied as CSV");
}
function openDock() { $("dock").classList.add("open"); $("dockToggle").textContent = "▾"; }

/* ---------------- processing ---------------- */
function procOps() {
  const raw = document.querySelector('input[name=procsrc]:checked').value === "raw";
  const ops = [];
  if (raw) {
    const lb = parseFloat($("procLb").value) || 0;
    if (lb > 0) ops.push({ op: "em", lb_hz: lb });
    const zf = parseInt($("procZf").value) || 1;
    if (zf > 1) ops.push({ op: "zf", factor: zf });
    ops.push({ op: "ft", offset_ppm: parseFloat($("procOff").value) || 0 });
  }
  const p0 = parseFloat($("p0v").value) || 0, p1 = parseFloat($("p1v").value) || 0;
  if (p0 || p1) ops.push({ op: "phase", p0, p1 });
  return { raw, ops };
}
async function applyProcessing(extraOps) {
  const { raw, ops } = procOps();
  if (extraOps) ops.push(...extraOps);
  status("processing…", "busy");
  try {
    const res = await api("process", { source_path: state.sourcePath, ops, use_raw: raw });
    state.exp = { ppm: res.ppm, amp: res.amp };
    drawPlot(state.lastModel || null); simulate();
    status("processing applied");
  } catch (e) { status("processing: " + e.message, "error"); }
}
async function resetProcessing() {
  status("reloading original…", "busy");
  try {
    const data = await api("load", { path: state.sourcePath });
    state.exp = { ppm: data.ppm, amp: data.amp };
    drawPlot(state.lastModel || null); simulate();
    status("original data restored");
  } catch (e) { status(e.message, "error"); }
}

/* ---------------- figure studio ---------------- */
const figStatus = (m, err) => { const el = $("figStatus"); el.textContent = m; el.style.color = err ? "var(--bad)" : "var(--muted)"; };
async function figTemplates() {
  try {
    const res = await api("figure/template", { source_path: state.sourcePath, recipe: state.recipe });
    const box = $("figTplBtns"); box.innerHTML = "";
    for (const [name, spec] of Object.entries(res.templates)) {
      const b = document.createElement("button");
      b.textContent = name;
      b.onclick = () => { $("figSpec").value = JSON.stringify(spec, null, 2); figStatus("template loaded"); };
      box.appendChild(b);
    }
  } catch (e) { figStatus(e.message, true); }
}
function parseFigSpec() {
  try { return JSON.parse($("figSpec").value); }
  catch (e) { figStatus("invalid JSON: " + e.message, true); return null; }
}
async function figPreview() {
  const spec = parseFigSpec(); if (!spec) return;
  try {
    figStatus("rendering…");
    const res = await api("figure/preview", { spec });
    $("figImg").src = "data:image/png;base64," + res.png_base64;
    $("figImg").style.display = "block";
    figStatus("");
  } catch (e) { figStatus(e.message, true); }
}
async function figExport() {
  const spec = parseFigSpec(); if (!spec) return;
  const formats = Array.from(document.querySelectorAll(".figfmt:checked")).map(c => c.value);
  if (!$("figPath").value.trim()) { figStatus("give an export base path", true); return; }
  try {
    figStatus("exporting…");
    const res = await api("figure/export", { spec, path: $("figPath").value.trim(), formats });
    figStatus("saved: " + res.saved.join("  "));
  } catch (e) { figStatus(e.message, true); }
}

/* ---------------- file browser ---------------- */
async function browse(path) {
  try {
    const res = await api("browse", { path: path || "" });
    $("bwPath").value = res.path;
    const ul = $("bwList"); ul.innerHTML = "";
    const add = (name, tag, onClick) => {
      const li = document.createElement("li");
      const t = document.createElement("span");
      t.className = "tag" + (tag === "dir" ? " dir" : "");
      t.textContent = tag;
      li.append(t, document.createTextNode(" " + name));
      li.onclick = onClick;
      ul.appendChild(li);
    };
    res.dirs.forEach(d => add(d, "dir", () => browse(res.path ? res.path + "\\" + d : d)));
    res.expnos.forEach(d => add(d, "EXPNO", () => { closeBrowser(); loadAndReport(res.path + "\\" + d); }));
    res.files.forEach(f => add(f, "file", () => { closeBrowser(); loadAndReport(res.path + "\\" + f); }));
    $("bwUp").onclick = res.parent ? () => browse(res.parent) : null;
    $("bwUp").disabled = !res.parent;
  } catch (e) { status("browse: " + e.message, "error"); }
}
async function loadAndReport(path) {
  try {
    await loadSource(path);   // backend handles .fxmla, EXPNO dirs, and .json recipes
  } catch (e) { status("load: " + e.message, "error"); }
}
function openBrowser() {
  $("modal-overlay").classList.add("open");
  browse(state.sourcePath ? state.sourcePath.replace(/[\\/][^\\/]*$/, "") : "");
}
function closeBrowser() { $("modal-overlay").classList.remove("open"); }

/* ---------------- save recipe ---------------- */
async function saveRecipe() {
  const def = (state.recipe && state.recipe.sample ? state.recipe.sample.replace(/[^\w-]+/g, "_") : "fit") + ".recipe.json";
  const path = prompt("Save recipe as (full path):",
    (localStorage.getItem("larmor.lastSaveDir") || "C:\\") + def);
  if (!path) return;
  try {
    const res = await api("save", { recipe: state.recipe, path });
    localStorage.setItem("larmor.lastSaveDir", path.replace(/[^\\/]*$/, ""));
    status("recipe saved: " + res.saved);
  } catch (e) { status("save: " + e.message, "error"); }
}

/* ---------------- wiring ---------------- */
function bindUI() {
  $("btnOpen").onclick = openBrowser;
  $("bwCancel").onclick = closeBrowser;
  $("bwGo").onclick = () => browse($("bwPath").value.trim());
  $("bwPath").addEventListener("keydown", e => { if (e.key === "Enter") browse($("bwPath").value.trim()); });
  $("modal-overlay").addEventListener("click", e => { if (e.target === $("modal-overlay")) closeBrowser(); });

  $("btnSave").onclick = saveRecipe;
  $("btnFit").onclick = runFit;
  $("btnQuant").onclick = () => runQuantify(true);
  $("btnUndo").onclick = undo;
  $("btnRedo").onclick = redo;
  $("btnCsv").onclick = quantCsv;
  $("dockToggle").onclick = () => {
    const d = $("dock");
    d.classList.toggle("open");
    $("dockToggle").textContent = d.classList.contains("open") ? "▾" : "▴";
  };
  $("showResid").onchange = () => drawPlot(state.lastModel || state.lastFit);
  $("showComp").onchange = () => drawPlot(state.lastModel || state.lastFit);
  $("btnZoomWin").onclick = () => {
    const gd = $("plot");
    const xr = gd.layout && gd.layout.xaxis && gd.layout.xaxis.range;
    if (xr) { $("winHi").value = Math.max(...xr).toFixed(1); $("winLo").value = Math.min(...xr).toFixed(1); persist(); }
  };

  // tabs
  document.querySelectorAll("#tabs button").forEach(b => {
    b.onclick = () => {
      document.querySelectorAll("#tabs button").forEach(x => x.classList.remove("on"));
      document.querySelectorAll(".tabpane").forEach(x => x.classList.remove("on"));
      b.classList.add("on");
      $("tab-" + b.dataset.tab).classList.add("on");
    };
  });

  // processing
  document.querySelectorAll('input[name=procsrc]').forEach(r =>
    r.onchange = () => { $("rawopts").style.display = r.value === "raw" && r.checked ? "flex" : "none"; });
  const syncPair = (range, num) => {
    $(range).oninput = () => { $(num).value = $(range).value; };
    $(num).onchange = () => { $(range).value = $(num).value; };
  };
  syncPair("p0", "p0v"); syncPair("p1", "p1v");
  $("btnApplyProc").onclick = () => applyProcessing();
  $("btnAutophase").onclick = () => applyProcessing([{ op: "autophase" }]);
  $("btnBaseline").onclick = () => applyProcessing([{ op: "baseline", order: parseInt($("blOrder").value) || 3 }]);
  $("btnResetProc").onclick = resetProcessing;

  // figure studio
  $("figTplLoad").onclick = figTemplates;
  $("figPreview").onclick = figPreview;
  $("figExport").onclick = figExport;

  // keyboard
  document.addEventListener("keydown", e => {
    if (e.target.tagName === "INPUT" || e.target.tagName === "TEXTAREA") return;
    if (e.key === "Escape") setAddMode(null);
    else if (e.key === "f" || e.key === "F") runFit();
    else if (e.key === "q" || e.key === "Q") runQuantify(true);
    else if ((e.ctrlKey || e.metaKey) && e.key === "z") { e.preventDefault(); undo(); }
    else if ((e.ctrlKey || e.metaKey) && e.key === "y") { e.preventDefault(); redo(); }
  });
}

/* ---------------- boot ---------------- */
(async function boot() {
  bindUI();
  Plotly.newPlot("plot", [], { margin: { t: 12, r: 14, l: 55, b: 40 } },
                 { responsive: true, displaylogo: false, scrollZoom: true });
  bindPlotEvents();
  try {
    const res = await apiGet("models");
    state.models = res.models;
    buildAddBar();
  } catch (e) { status("could not load model registry: " + e.message, "error"); }
  tryRestore();
})();
