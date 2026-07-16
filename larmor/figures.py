"""Figure studio: publication-quality NMR figures from declarative JSON specs.

Inspired by NMRVEW (T. Charpentier-era CEMHTI notebooks): composable 1D graphs
with per-trace scale/offset/window-normalization, 2D contour maps with top and
right projections, overlaid external 1D traces, F1-band sub-projections and
slope lines, and series plots (saturation recovery, REDOR). A figure is a
plain JSON dict ("spec"), so it can be saved next to the data and re-rendered
exactly -- same philosophy as the fit recipe.

Spec kinds:
  {"kind": "1d",     ...}  stacked/overlaid 1D traces (spectra, fits, components)
  {"kind": "2d",     ...}  contour map (MQMAS, HMQC, SQ-DQ, ...) + projections
  {"kind": "series", ...}  points-vs-time (satrec T1, REDOR dephasing)

Every renderer returns a matplotlib Figure; export() writes png/svg/pdf.
"""
from __future__ import annotations

import io
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

# ---------------------------------------------------------------------------
# style presets (rcParams bundles); sizes in inches follow journal conventions
STYLES: dict[str, dict] = {
    "article": {  # single column (~3.5 in wide)
        "figsize": (3.5, 2.8),
        "rc": {"font.size": 8, "axes.linewidth": 0.8, "lines.linewidth": 1.0,
               "xtick.direction": "out", "ytick.direction": "out",
               "font.family": "sans-serif", "legend.fontsize": 6.5,
               "legend.frameon": False},
    },
    "article-wide": {  # double column (~7 in wide)
        "figsize": (7.0, 3.2),
        "rc": {"font.size": 9, "axes.linewidth": 0.8, "lines.linewidth": 1.2,
               "xtick.direction": "out", "ytick.direction": "out",
               "font.family": "sans-serif", "legend.fontsize": 7.5,
               "legend.frameon": False},
    },
    "presentation": {
        "figsize": (9.0, 5.5),
        "rc": {"font.size": 16, "axes.linewidth": 1.6, "lines.linewidth": 2.2,
               "xtick.direction": "out", "ytick.direction": "out",
               "font.family": "sans-serif", "legend.fontsize": 13,
               "legend.frameon": False},
    },
    "thesis": {
        "figsize": (5.8, 4.0),
        "rc": {"font.size": 11, "axes.linewidth": 1.0, "lines.linewidth": 1.4,
               "xtick.direction": "out", "ytick.direction": "out",
               "font.family": "serif", "legend.fontsize": 9,
               "legend.frameon": False},
    },
}

#: superscripted isotope label, e.g. "27Al" -> "$^{27}$Al NMR shift (ppm)"
def nucleus_xlabel(nucleus: str) -> str:
    digits = "".join(c for c in nucleus if c.isdigit())
    symbol = "".join(c for c in nucleus if c.isalpha())
    return rf"$^{{{digits}}}${symbol} NMR shift (ppm)"


# ---------------------------------------------------------------------------
# trace sources

def _norm_window(x: np.ndarray, y: np.ndarray, window) -> np.ndarray:
    """Normalize y to max 1 within a ppm window (NMRVEW's norm_0_to_1)."""
    if window:
        sel = (x >= min(window)) & (x <= max(window))
        peak = np.abs(y[sel]).max() if sel.any() else np.abs(y).max()
    else:
        peak = np.abs(y).max()
    return y / peak if peak else y


def load_trace(t: dict) -> tuple[np.ndarray, np.ndarray, dict]:
    """Resolve one trace spec to (x_ppm, y, meta).

    Trace kinds:
      {"data": {"x": [...], "y": [...]}}                inline arrays
      {"path": "...fxmla"}                              embedded dmfit spectrum
      {"path": "<EXPNO dir>", "procno": 1}              Bruker processed 1D
      {"recipe": "...json", "part": "total"}            simulated fit total
      {"recipe": "...json", "part": "site", "site": 0}  one fit component
      {"recipe": "...json", "part": "residual"}         experiment - total
    """
    meta: dict = {}
    if "data" in t:
        return np.asarray(t["data"]["x"], float), np.asarray(t["data"]["y"], float), meta

    if "recipe" in t:
        from larmor import engine
        from larmor.recipe import Recipe

        recipe = Recipe.load(t["recipe"])
        meta["nucleus"] = recipe.nucleus
        exp_ppm = None
        if not engine.needs_kernel(recipe):
            exp_ppm, _, _ = load_trace({"path": recipe.source_path})
        x, total, per_site = engine.simulate(recipe, exp_ppm=exp_ppm)
        part = t.get("part", "total")
        if part == "total":
            return x, total, meta
        if part == "site":
            i = int(t["site"])
            meta["label"] = recipe.sites[i].label
            return x, per_site[i], meta
        if part == "residual":
            ex, ey, _ = load_trace({"path": recipe.source_path})
            yi = np.interp(ex, x, total)
            return ex, ey - yi, meta
        raise ValueError(f"unknown recipe part {part!r}")

    path = Path(t["path"])
    if path.suffix.lower() in (".fxmla", ".fxml"):
        from larmor.io import fxmla

        dm = fxmla.read(path)
        meta["nucleus"] = dm.dimensions[0].nucleus
        x, y = dm.spectrum.ppm, dm.spectrum.amplitude
    else:
        from larmor.io import bruker

        exp = bruker.read_expno(path, procno=int(t.get("procno", 1)))
        meta["nucleus"] = exp.nucleus
        x, y = exp.processed_ppm, exp.processed.astype(float)
    order = np.argsort(x)
    return np.asarray(x)[order], np.asarray(y)[order], meta


# ---------------------------------------------------------------------------
# 1D figures

def render_1d(spec: dict) -> Figure:
    style = STYLES[spec.get("style", "article")]
    with plt.rc_context(style["rc"]):
        fig, ax = plt.subplots(figsize=spec.get("figsize", style["figsize"]))
        nucleus = None
        for t in spec.get("traces", []):
            x, y, meta = load_trace(t)
            nucleus = nucleus or meta.get("nucleus")
            if t.get("normalize") is not None:
                y = _norm_window(x, y, t["normalize"] if t["normalize"] is not True else None)
            y = y * float(t.get("scale", 1.0)) + float(t.get("offset", 0.0))
            (line,) = ax.plot(x, y,
                              lw=t.get("linewidth", None),
                              ls=t.get("linestyle", "-"),
                              alpha=t.get("alpha", 1.0),
                              label=t.get("label", meta.get("label")))
            if t.get("color"):
                line.set_color(t["color"])
        xlabel = spec.get("xlabel") or (nucleus_xlabel(nucleus) if nucleus else "shift (ppm)")
        ax.set_xlabel(xlabel)
        if spec.get("xlim"):
            hi, lo = spec["xlim"]
            ax.set_xlim(max(hi, lo), min(hi, lo))  # ppm axes run high -> low
        else:
            ax.invert_xaxis()
        if spec.get("ylim"):
            ax.set_ylim(*spec["ylim"])
        ax.set_yticks([]) if spec.get("hide_yaxis", True) else None
        if spec.get("hide_yaxis", True):
            ax.spines["left"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["top"].set_visible(False)
        for a in spec.get("annotations", []):
            ax.text(a["x"], a["y"], a["text"], fontsize=a.get("fontsize"),
                    ha=a.get("ha", "left"))
        if any(t.get("label") for t in spec.get("traces", [])):
            ax.legend(loc=spec.get("legend_loc", "best"),
                      ncol=spec.get("legend_ncol", 1))
        fig.tight_layout()
        return fig


# ---------------------------------------------------------------------------
# 2D figures

def load_2d(path: str | Path, procno: int = 1):
    """Read a processed Bruker 2D (2rr) read-only. Returns (x_F2, y_F1, Z)."""
    import nmrglue as ng

    pdata = Path(path) / "pdata" / str(procno)
    dic, Z = ng.bruker.read_pdata(str(pdata))
    axes = []
    for key, npts in (("procs", Z.shape[1]), ("proc2s", Z.shape[0])):
        p = dic[key]
        si, sf = int(p["SI"]), float(p["SF"])
        offset, sw = float(p["OFFSET"]), float(p["SW_p"])
        axes.append(offset - np.arange(npts) * (sw / sf / npts))
    x_f2, y_f1 = axes
    return x_f2, y_f1, Z.astype(float)


def render_2d(spec: dict) -> Figure:
    style = STYLES[spec.get("style", "article")]
    x, y, Z = load_2d(spec["path"], int(spec.get("procno", 1)))
    Z = Z / np.abs(Z).max()

    lev = spec.get("levels", {})
    n = int(lev.get("n", 12))
    if "min_frac" in lev:
        min_frac = float(lev["min_frac"])
    else:
        # default the contour floor to the measured noise: sample the outer
        # 5% frame of the matrix, put the lowest contour at ~8 sigma
        edge = max(1, min(Z.shape) // 20)
        frame = np.concatenate([Z[:edge].ravel(), Z[-edge:].ravel(),
                                Z[:, :edge].ravel(), Z[:, -edge:].ravel()])
        min_frac = float(np.clip(8.0 * frame.std(), 0.02, 0.5))
    if lev.get("mode", "log") == "log":
        levels = np.logspace(np.log10(min_frac), 0, n)
    else:
        levels = np.linspace(min_frac, 1.0, n)

    with plt.rc_context(style["rc"]):
        base_w, base_h = spec.get("figsize", (5.2, 5.6))
        show_top = spec.get("proj_top", True)
        show_right = spec.get("proj_right", True)
        fig = plt.figure(figsize=(base_w, base_h))
        gs = fig.add_gridspec(2, 2, width_ratios=[5, 1.15], height_ratios=[1.15, 5],
                              hspace=0.04, wspace=0.04)
        ax = fig.add_subplot(gs[1, 0])
        ax_top = fig.add_subplot(gs[0, 0], sharex=ax) if show_top else None
        ax_right = fig.add_subplot(gs[1, 1], sharey=ax) if show_right else None

        cmap = spec.get("cmap", "viridis")
        ax.contour(x, y, Z, levels=levels, cmap=cmap, linewidths=0.7)
        if spec.get("negative"):
            ax.contour(x, y, -Z, levels=levels, colors="crimson",
                       linewidths=0.7, linestyles="dashed")

        xlim = spec.get("xlim") or (float(x.max()), float(x.min()))
        ylim = spec.get("ylim") or (float(y.max()), float(y.min()))
        ax.set_xlim(max(xlim), min(xlim))
        ax.set_ylim(max(ylim), min(ylim))
        ax.set_xlabel(spec.get("xlabel", "F2 shift (ppm)"))
        ax.set_ylabel(spec.get("ylabel", "F1 (ppm)"))
        ax.grid(spec.get("grid", True), ls="-.", lw=0.5, color="0.75")

        xsel = (x >= min(xlim)) & (x <= max(xlim))
        ysel = (y >= min(ylim)) & (y <= max(ylim))
        mode = spec.get("projection", "skyline")
        reducer = (lambda a, axis: a.max(axis=axis)) if mode == "skyline" else \
                  (lambda a, axis: a.sum(axis=axis))

        if ax_top is not None:
            proj = reducer(Z[np.ix_(ysel, xsel)], 0)
            proj = proj / np.abs(proj).max()
            ax_top.plot(x[xsel], proj, lw=1.0, label=spec.get("proj_label", "projection"))
            for t in spec.get("overlay_top", []):
                tx, ty, _ = load_trace(t)
                ty = _norm_window(tx, ty, t.get("normalize"))
                ty = ty * float(t.get("scale", 1.0)) + float(t.get("offset", 0.0))
                ax_top.plot(tx, ty, lw=1.0, label=t.get("label"),
                            color=t.get("color"), alpha=t.get("alpha", 0.9))
            for sp_ in spec.get("subproj", []):
                f1a, f1b = sp_["f1"]
                band = (y >= min(f1a, f1b)) & (y <= max(f1a, f1b))
                sub = Z[np.ix_(band, xsel)].sum(axis=0)
                sub = sub / np.abs(sub).max() * float(sp_.get("scale", 1.0))
                ax_top.plot(x[xsel], sub, lw=0.9, ls="--", label=sp_.get("label"))
            ax_top.set_yticks([])
            plt.setp(ax_top.get_xticklabels(), visible=False)
            for s in ("left", "right", "top"):
                ax_top.spines[s].set_visible(False)
            if any(l.get_label() and not l.get_label().startswith("_")
                   for l in ax_top.lines):
                ax_top.legend(loc=spec.get("legend_top_loc", "upper right"))

        if ax_right is not None:
            proj = reducer(Z[np.ix_(ysel, xsel)], 1)
            proj = proj / np.abs(proj).max()
            ax_right.plot(proj, y[ysel], lw=1.0)
            ax_right.set_xticks([])
            plt.setp(ax_right.get_yticklabels(), visible=False)
            for s in ("right", "top", "bottom"):
                ax_right.spines[s].set_visible(False)

        for sl in spec.get("slopes", []):
            xs = np.array([min(xlim), max(xlim)])
            ax.plot(xs, sl["slope"] * xs + sl.get("intercept", 0.0),
                    color=sl.get("color", "k"), lw=sl.get("linewidth", 0.9),
                    ls=sl.get("linestyle", "-"))

        if spec.get("annotation"):
            ax.text(0.04, 0.94, spec["annotation"], transform=ax.transAxes,
                    fontsize=style["rc"]["font.size"] + 2, va="top")
        return fig


# ---------------------------------------------------------------------------
# series figures (saturation recovery, REDOR)

def load_series(spec: dict) -> dict:
    """Extract a series from a Bruker EXPNO's TopSpin analysis files.

    mode "satrec": pdata/<procno>/t1ints.txt (delay / integral blocks)
    mode "redor":  pdata/<procno>/redor.txt  (S0 / S pairs + spinning speed)
    Inline data:   {"data": {"x": [...], "y": [...], ["yerr": ...]}}
    """
    if "data" in spec:
        d = spec["data"]
        return {"x": np.asarray(d["x"], float), "y": np.asarray(d["y"], float),
                "yerr": np.asarray(d["yerr"], float) if "yerr" in d else None}

    pdata = Path(spec["path"]) / "pdata" / str(spec.get("procno", 1))
    mode = spec.get("mode", "satrec")
    if mode == "satrec":
        vals = [float(v) for v in (pdata / "t1ints.txt").read_text().split()]
        # header = total line count; then blocks of 3 lines per point:
        # (delay,0,0) (0,0,1) (npts,integral,0); terminated by a (-1,0,0) line
        delays, integrals = [], []
        i = 1
        while i + 9 <= len(vals) + 1:
            if vals[i] < 0:  # -1 sentinel
                break
            block = vals[i:i + 9]
            if len(block) < 9:
                break
            delays.append(block[0])
            integrals.append(block[7])
            i += 9
        x, y = np.array(delays), np.array(integrals)
        y = y / np.abs(y).max()
        return {"x": x, "y": np.abs(y), "yerr": None}
    if mode == "redor":
        text = (pdata / "redor.txt").read_text().splitlines()
        masr = None
        rows = []
        for ln in text:
            parts = ln.split()
            if ln.strip().startswith("Spinning speed"):
                masr = float(parts[-1])
            if len(parts) == 5 and parts[0].isdigit():
                rows.append([float(v) for v in parts])
        arr = np.array(rows)
        n, s0, s = arr[:, 0], arr[:, 1], arr[:, 2]
        with np.errstate(divide="ignore", invalid="ignore"):
            ds = 1.0 - np.where(s0 != 0, s / s0, np.nan)
        x = n / masr * 1000.0 if masr else n  # ms of recoupling (n rotor periods)
        return {"x": x, "y": ds, "yerr": None,
                "xlabel": "recoupling time (ms)" if masr else "rotor cycles",
                "ylabel": r"$\Delta S/S_0$"}
    raise ValueError(f"unknown series mode {mode!r}")


def _fit_satrec(x: np.ndarray, y: np.ndarray, stretched: bool):
    import lmfit

    def model(p):
        beta = p["beta"].value if stretched else 1.0
        return p["a"] * (1.0 - np.exp(-((x / p["t1"]) ** beta)))

    params = lmfit.Parameters()
    params.add("a", value=float(y.max()), min=0)
    params.add("t1", value=float(x[np.searchsorted(y, 0.63 * y.max())]
                                 if y.max() > 0 else 1.0), min=1e-6)
    if stretched:
        params.add("beta", value=1.0, min=0.2, max=2.0)
    out = lmfit.minimize(lambda p: model(p) - y, params, method="leastsq")
    return out, (lambda xx: out.params["a"].value *
                 (1 - np.exp(-((xx / out.params["t1"].value) **
                               (out.params["beta"].value if stretched else 1.0)))))


def render_series(spec: dict) -> Figure:
    style = STYLES[spec.get("style", "article")]
    data = load_series(spec)
    mode = spec.get("mode", "satrec")
    with plt.rc_context(style["rc"]):
        fig, ax = plt.subplots(figsize=spec.get("figsize", style["figsize"]))
        ax.plot(data["x"], data["y"], "o", ms=4, mfc="white",
                label=spec.get("label", {"satrec": "integrals",
                                         "redor": r"$\Delta S/S_0$"}.get(mode)))
        if data.get("yerr") is not None:
            ax.errorbar(data["x"], data["y"], yerr=data["yerr"], fmt="none",
                        ecolor="0.4", capsize=2)
        note = None
        if mode == "satrec" and spec.get("fit", True):
            sel = data["x"] > 0
            out, curve = _fit_satrec(data["x"][sel], data["y"][sel],
                                     stretched=spec.get("stretched", False))
            xx = np.logspace(np.log10(max(data["x"][sel].min(), 1e-4)),
                             np.log10(data["x"].max()), 200)
            ax.plot(xx, curve(xx), "-", lw=1.2, label="fit")
            t1 = out.params["t1"]
            err = f" ± {t1.stderr:.2g}" if t1.stderr else ""
            note = f"$T_1$ = {t1.value:.3g}{err} s"
            if spec.get("stretched"):
                b = out.params["beta"]
                note += f", β = {b.value:.2f}"
            ax.set_xscale("log")
        ax.set_xlabel(spec.get("xlabel", data.get("xlabel",
                      "recovery delay (s)" if mode == "satrec" else "time")))
        ax.set_ylabel(spec.get("ylabel", data.get("ylabel",
                      "normalized integral" if mode == "satrec" else "")))
        if note:
            ax.text(0.05, 0.9, note, transform=ax.transAxes)
        if spec.get("annotation"):
            ax.text(0.05, 0.78, spec["annotation"], transform=ax.transAxes)
        ax.legend(loc=spec.get("legend_loc", "lower right"))
        fig.tight_layout()
        return fig


# ---------------------------------------------------------------------------
RENDERERS = {"1d": render_1d, "2d": render_2d, "series": render_series}


def render(spec: dict) -> Figure:
    kind = spec.get("kind")
    if kind not in RENDERERS:
        raise ValueError(f"unknown figure kind {kind!r} (valid: {list(RENDERERS)})")
    return RENDERERS[kind](spec)


def render_png_bytes(spec: dict, dpi: int = 130) -> bytes:
    fig = render(spec)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


def export(spec: dict, out_base: str | Path,
           formats: tuple[str, ...] = ("png", "svg", "pdf"),
           dpi: int = 600) -> list[str]:
    """Write the figure in each format next to out_base (no extension)."""
    fig = render(spec)
    out_base = Path(out_base)
    out_base.parent.mkdir(parents=True, exist_ok=True)
    saved = []
    for fmt in formats:
        target = out_base.with_suffix("." + fmt)
        fig.savefig(target, format=fmt, dpi=dpi, bbox_inches="tight")
        saved.append(str(target))
    plt.close(fig)
    return saved
