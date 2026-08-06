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
# A shared, literature-quality base: clean sans-serif, thin dark-grey axes, ticks
# in, generous label spacing, no chartjunk. Per-style dicts only tweak size/family.
_BASE_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Helvetica Neue", "DejaVu Sans"],
    "mathtext.fontset": "dejavusans",
    "axes.edgecolor": "#2b2b2b",
    "axes.labelcolor": "#1a1a1a",
    "axes.labelweight": "normal",
    "axes.titleweight": "bold",
    "axes.labelpad": 3.5,
    "axes.linewidth": 0.9,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.major.size": 4.0, "ytick.major.size": 4.0,
    "xtick.minor.size": 2.2, "ytick.minor.size": 2.2,
    "xtick.major.width": 0.9, "ytick.major.width": 0.9,
    "xtick.color": "#2b2b2b", "ytick.color": "#2b2b2b",
    "lines.linewidth": 1.3, "lines.solid_capstyle": "round",
    "legend.frameon": False, "legend.handlelength": 1.5,
    "legend.labelspacing": 0.3, "legend.borderaxespad": 0.4,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "savefig.dpi": 300,
    "savefig.bbox": "tight",
}


def _style(figsize, **over):
    return {"figsize": figsize, "rc": {**_BASE_RC, **over}}


STYLES: dict[str, dict] = {
    "article": _style((3.5, 2.8), **{"font.size": 8.5, "legend.fontsize": 7,
                                     "lines.linewidth": 1.2}),
    "article-wide": _style((7.0, 3.2), **{"font.size": 9, "legend.fontsize": 7.5,
                                          "lines.linewidth": 1.3}),
    "presentation": _style((9.0, 5.5), **{"font.size": 16, "axes.linewidth": 1.4,
                                          "lines.linewidth": 2.4,
                                          "legend.fontsize": 13,
                                          "xtick.major.size": 6,
                                          "ytick.major.size": 6}),
    "thesis": _style((5.8, 4.0), **{"font.size": 11, "font.family": "serif",
                                    "mathtext.fontset": "dejavuserif",
                                    "lines.linewidth": 1.5, "legend.fontsize": 9}),
    # journal presets (single-column widths + each house's typographic minimums)
    "nature": _style((3.50, 2.7), **{"font.size": 7, "axes.linewidth": 0.6,
                                     "lines.linewidth": 1.0, "legend.fontsize": 6,
                                     "axes.labelpad": 2.5}),
    "acs": _style((3.25, 2.6), **{"font.size": 8, "lines.linewidth": 1.0,
                                  "legend.fontsize": 7}),
    "rsc": _style((3.27, 2.6), **{"font.size": 8, "lines.linewidth": 1.1,
                                  "legend.fontsize": 7}),
}

#: Named, generic STRUCTURAL presets — what kind of NMR figure this is, not
#: what journal it's for (that's STYLES above; the two combine freely: pick a
#: template for the layout, a style for the typography). Each maps to a
#: partial spec the Plotting studio pre-fills; nucleus-agnostic by design —
#: name by what the figure IS, never by which nucleus happened to be in the
#: reference example. Values are the fixed parameters of that figure type;
#: user-specific parts (traces/panels/categories) are filled in afterward.
TEMPLATES: dict[str, dict] = {
    "Stacked series": {
        "description": "One 1D spectrum per sample, offset vertically, each "
                       "labelled at its own trace end instead of a legend box "
                       "— the standard way to show a composition/time series "
                       "without a fit (e.g. a raw-spectra overview panel).",
        "kind": "1d",
        "spec": {"x_is_ppm": True, "hide_yaxis": True, "legend_loc": "none"},
        "trace_defaults": {"end_label": True},
    },
    "Deconvolution grid": {
        "description": "One panel per fitted spectrum: experiment + total "
                       "fit + every component, filled and labelled by "
                       "position — the standard full-deconvolution figure "
                       "(dmfit/ssNake style) for a series of related fits.",
        "kind": "batch_grid",
        "spec": {"component_mode": "fill", "peak_labels": "position",
                 "show_total": True, "show_experiment": True, "legend": True},
    },
    "Composition series (shaded component)": {
        "description": "Same grid, but only ONE component per panel is "
                       "filled/highlighted (the rest are invisible, only the "
                       "total fit line shows them) with its position and "
                       "population % labelled — for tracking one species "
                       "across a composition series.",
        "kind": "batch_grid",
        "spec": {"component_mode": "fill", "shade_only": [0],
                 "peak_labels": "position+pct", "show_total": True,
                 "show_experiment": True, "legend": False, "cols": 1},
    },
    "Composition trend": {
        "description": "Population (%) — or any fitted quantity — vs. a "
                       "composition variable, one line+marker series per "
                       "site/species, with a legend. Not a spectrum: a "
                       "generic x–y plot (set 'x_is_ppm' off).",
        "kind": "1d",
        "spec": {"x_is_ppm": False, "hide_yaxis": False, "grid": False,
                 "legend_loc": "best"},
        "trace_defaults": {"marker": "o"},
    },
    "Species distribution": {
        "description": "A 100%-stacked bar of species population vs. "
                       "composition — the categorical alternative to a "
                       "composition-trend line plot when you want every "
                       "sample's full speciation at a glance.",
        "kind": "species_bar",
        "spec": {"normalize": True, "value_labels": True},
    },
    "2D correlation": {
        "description": "A 2D contour map (MQMAS/DQ-SQ/HMQC/…) with top and "
                       "side projections and optional diagonal/connectivity "
                       "reference lines.",
        "kind": "2d",
        "spec": {"contour_mode": "contour", "proj_top": True,
                 "proj_right": True},
    },
}


#: the app's one categorical identity palette (Okabe-Ito derived, colorblind-
#: safe, fixed order — never cycled past, never re-derived per feature). A
#: static-figure copy of the live desktop's site colors (larmor.desktop.theme.
#: LIGHT_SERIES): the core layer can't import the desktop layer, but an
#: exported figure should still show site A in the same colour the fit table
#: does. Keep the two lists' VALUES in sync if either changes.
SITE_COLORS = ["#0072b2", "#d55e00", "#009e73", "#b0568c", "#8f6e00",
              "#0e7c86", "#7f3fbf", "#117733", "#882255", "#5a5a5a"]


def site_color(i: int) -> str:
    return SITE_COLORS[i % len(SITE_COLORS)]


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
        d = t["data"]
        if d.get("yerr") is not None:
            meta["yerr"] = np.asarray(d["yerr"], float)
        return np.asarray(d["x"], float), np.asarray(d["y"], float), meta

    if "recipe" in t:
        from larmor import engine
        from larmor.recipe import Recipe

        recipe = Recipe.load(t["recipe"])
        meta["nucleus"] = recipe.nucleus
        exp_ppm = None
        if not engine.needs_kernel(recipe) and recipe.source_path:
            try:
                exp_ppm, _, _ = load_trace({"path": recipe.source_path})
            except Exception:
                exp_ppm = None       # source unreachable -> fall back below
        if exp_ppm is None and not engine.needs_kernel(recipe):
            # data-less (or unreachable-source) non-kernel fit: a grid framed
            # around the sites, same as a data-less dmfit .fxmla just below
            exp_ppm, _ = _simulate_model_curve(recipe)
        x, total, per_site = engine.simulate(recipe, exp_ppm=exp_ppm)
        part = t.get("part", "total")
        if part == "total":
            return x, total, meta
        if part == "site":
            i = int(t["site"])
            meta["label"] = recipe.sites[i].label
            return x, per_site[i], meta
        if part == "residual":
            if not recipe.source_path:
                raise ValueError(
                    "residual needs the fit's source spectrum, but this "
                    "recipe carries no source_path (a data-less fit)")
            ex, ey, _ = load_trace({"path": recipe.source_path})
            yi = np.interp(ex, x, total)
            return ex, ey - yi, meta
        raise ValueError(f"unknown recipe part {part!r}")

    path = Path(t["path"])
    if path.suffix.lower() in (".fxmla", ".fxml"):
        from larmor.io import fxmla

        dm = fxmla.read(path)
        meta["nucleus"] = dm.dimensions[0].nucleus
        if dm.spectrum is not None:
            x, y = dm.spectrum.ppm, dm.spectrum.amplitude
        else:                                    # a data-less dmfit fit → its model
            recipe, _ = fxmla.to_recipe(dm)
            x, y = _simulate_model_curve(recipe)
    else:
        # any Bruker path (1r/2rr file, pdata folder, EXPNO), CSV or Varian
        from larmor.loader import load_any

        x, y, rec, _meta, _warns = load_any(str(path))
        meta["nucleus"] = rec.get("nucleus", "")
    order = np.argsort(x)
    return np.asarray(x)[order], np.asarray(y)[order], meta


def _simulate_model_curve(recipe, n: int = 3000):
    """A model curve for a fit that carries no spectrum: simulate over a grid
    spanning the sites (± a margin from their widths, or the fit window)."""
    from larmor import engine
    win = getattr(recipe, "fit_window_ppm", None)
    if win and len(win) == 2:
        lo, hi = min(win), max(win)
    else:
        centers = [float(c.value) for s in recipe.sites
                  for c in [s.params.get("isotropic_chemical_shift_ppm")]
                  if c is not None]
        if centers:
            m = engine.site_width_margin(recipe.sites)
            lo, hi = min(centers) - m, max(centers) + m
        else:
            lo, hi = -100.0, 100.0
    x = np.linspace(lo, hi, n)
    _, total, _ = engine.simulate(recipe, exp_ppm=x)
    return x, total


# ---------------------------------------------------------------------------
# 1D figures

def _norm_factor(x, y, mode: str) -> float:
    """Normalisation divisor for a trace: by peak, by area, or by edge noise."""
    y = np.asarray(y, float)
    if mode == "max":
        return float(np.max(np.abs(y))) or 1.0
    if mode == "area":
        trap = getattr(np, "trapezoid", None) or np.trapz
        return float(trap(np.abs(y), np.asarray(x, float))) or 1.0
    if mode == "noise":
        n = max(3, y.size // 20)
        return float(np.std(np.concatenate([y[:n], y[-n:]]))) or 1.0
    return 1.0


def render_1d(spec: dict) -> Figure:
    norm_mode = spec.get("norm")                    # None|"max"|"area"|"noise"
    with plt.rc_context(_rc(spec)):
        fig, ax = plt.subplots(figsize=_figsize(spec))
        nucleus = None
        # load every trace first (so we can normalise consistently and difference
        # against a reference — for series comparison)
        loaded = []
        for t in spec.get("traces", []):
            x, y, meta = load_trace(t)
            nucleus = nucleus or meta.get("nucleus")
            factor = 1.0                             # multiplicative norm factor
            if t.get("normalize") is not None:
                y2 = _norm_window(x, y, t["normalize"] if t["normalize"] is not True else None)
                factor = (y2 / y)[np.isfinite(y) & (y != 0)][:1]
                factor = float(factor[0]) if len(factor) else 1.0
                y = y2
            elif norm_mode:
                nf = _norm_factor(x, y, norm_mode)
                factor = 1.0 / nf if nf else 1.0
                y = np.asarray(y, float) * factor
            if meta.get("yerr") is not None:         # keep error bars consistent
                meta["yerr"] = np.asarray(meta["yerr"], float) * factor
            loaded.append([np.asarray(x, float), np.asarray(y, float), meta, t])
        if spec.get("difference") and loaded:        # subtract the first trace
            rx, ry = loaded[0][0], loaded[0][1]
            for row in loaded[1:]:
                row[1] = row[1] - np.interp(row[0], rx, ry)
            loaded = loaded[1:] if spec.get("difference") == "drop_ref" else loaded
        for x, y, meta, t in loaded:
            scale = float(t.get("scale", 1.0))
            y = y * scale + float(t.get("offset", 0.0))
            (line,) = ax.plot(x, y,
                              lw=t.get("linewidth", None),
                              ls=t.get("linestyle", "-"),
                              marker=t.get("marker"),
                              alpha=t.get("alpha", 1.0),
                              label=t.get("label", meta.get("label")))
            if t.get("color"):
                line.set_color(t["color"])
            if t.get("end_label"):
                # a label at the trace's own displayed end -- the literature
                # convention for a stacked series (e.g. a composition label
                # beside each spectrum) instead of a legend box per trace
                ppm = spec.get("x_is_ppm", True)
                at_min = bool(ppm)                # ppm axis is inverted: its
                edge = int(np.argmin(x) if at_min else np.argmax(x))  # right edge is the MIN x
                txt = t.get("label") or meta.get("label") or ""
                ax.annotate(txt, (x[edge], y[edge]), textcoords="offset points",
                           xytext=(-4 if at_min else 4, 2),
                           ha="right" if at_min else "left",
                           va="bottom", fontsize=plt.rcParams["font.size"],
                           color=line.get_color())
            ye = meta.get("yerr")
            if (ye is not None and np.isfinite(ye).any()
                    and t.get("err_visible", True)):
                ew = float(t.get("err_width", 1.2))
                ax.errorbar(x, y, yerr=np.abs(np.nan_to_num(ye)) * abs(scale),
                            fmt="none",
                            ecolor=t.get("err_color") or line.get_color(),
                            capsize=t.get("err_capsize", t.get("capsize", 3.5)),
                            elinewidth=ew, capthick=ew,
                            alpha=t.get("err_alpha", t.get("alpha", 1.0)),
                            zorder=line.get_zorder() + 1)
        # ppm spectra run high→low with a hidden intensity axis; a generic x-y
        # plot (e.g. a parameter-vs-sample series) keeps a normal, upright axis
        x_is_ppm = spec.get("x_is_ppm", True)
        hide_y = spec.get("hide_yaxis", x_is_ppm)
        xlabel = spec.get("xlabel") or (
            nucleus_xlabel(nucleus) if (nucleus and x_is_ppm)
            else ("shift (ppm)" if x_is_ppm else "x"))
        ax.set_xlabel(xlabel)
        if spec.get("ylabel"):
            ax.set_ylabel(spec["ylabel"])
        if spec.get("xlim"):
            hi, lo = spec["xlim"]
            ax.set_xlim((max(hi, lo), min(hi, lo)) if x_is_ppm else (lo, hi))
        elif x_is_ppm:
            ax.invert_xaxis()
        if spec.get("ylim"):
            lo, hi = spec["ylim"]
            ax.set_ylim(min(lo, hi), max(lo, hi))
        _apply_axis_extras(ax, spec)                 # ticks, grid, spines
        for a in spec.get("annotations", []):
            ax.text(a["x"], a["y"], a["text"], fontsize=a.get("fontsize"),
                    ha=a.get("ha", "left"))
        if hide_y:
            ax.set_yticks([])
            for s in ("left", "right", "top"):
                ax.spines[s].set_visible(False)
        elif not spec.get("box"):                    # clean upright x-y plot
            for s in ("right", "top"):
                ax.spines[s].set_visible(False)
        _apply_legend(ax, spec, has_labels=any(
            t.get("label") for t in spec.get("traces", [])))
        if spec.get("title"):
            ax.set_title(spec["title"])
        fig.tight_layout()
        return fig


# ---------------------------------------------------------------------------
# shared axis customisation (used by the 1D renderer and the plotting studio)

def _rc(spec: dict) -> dict:
    """rcParams for this figure: the chosen style, with optional font-size /
    line-width / legend-font-size overrides from the spec."""
    rc = dict(STYLES.get(spec.get("style", "article"), STYLES["article"])["rc"])
    if spec.get("font_size"):
        rc["font.size"] = float(spec["font_size"])
    if spec.get("legend_fontsize"):
        rc["legend.fontsize"] = float(spec["legend_fontsize"])
    if spec.get("line_width"):
        rc["lines.linewidth"] = float(spec["line_width"])
    if spec.get("tick_direction"):
        rc["xtick.direction"] = rc["ytick.direction"] = spec["tick_direction"]
    return rc


def _figsize(spec: dict):
    return spec.get("figsize",
                    STYLES.get(spec.get("style", "article"), STYLES["article"])["figsize"])


def _apply_axis_extras(ax, spec: dict) -> None:
    """Ticks, grid and tick styling common to the customisation panel."""
    from matplotlib.ticker import MultipleLocator

    # explicit tick positions+labels (e.g. sample names) take precedence
    if spec.get("xticks"):
        ax.set_xticks([float(p) for p, _ in spec["xticks"]])
        ax.set_xticklabels([str(lab) for _, lab in spec["xticks"]],
                           rotation=spec.get("xtick_rotation", 0),
                           ha="right" if spec.get("xtick_rotation") else "center")
    elif spec.get("xtick_step"):
        ax.xaxis.set_major_locator(MultipleLocator(float(spec["xtick_step"])))
    if spec.get("yticks"):
        ax.set_yticks([float(p) for p, _ in spec["yticks"]])
        ax.set_yticklabels([str(lab) for _, lab in spec["yticks"]])
    elif spec.get("ytick_step"):
        ax.yaxis.set_major_locator(MultipleLocator(float(spec["ytick_step"])))
    if spec.get("minor_ticks"):
        ax.minorticks_on()
    if spec.get("tick_direction"):
        ax.tick_params(which="both", direction=spec["tick_direction"])
    if spec.get("tick_labelsize"):
        ax.tick_params(labelsize=float(spec["tick_labelsize"]))
    if spec.get("grid"):
        ax.grid(True, which="major", ls=spec.get("grid_style", ":"),
                lw=0.6, color="0.7", alpha=0.8)


def _apply_legend(ax, spec: dict, has_labels: bool, default_loc: str = "best") -> None:
    """Legend at the requested location; 'none'/'off' hides it entirely."""
    loc = spec.get("legend_loc", default_loc)
    if not has_labels or loc in ("none", "off", None):
        return
    ax.legend(loc=loc, ncol=spec.get("legend_ncol", 1),
              fontsize=spec.get("legend_fontsize"),
              title=spec.get("legend_title"),
              framealpha=spec.get("legend_frame_alpha", 0.0),
              frameon=bool(spec.get("legend_frame", False)))


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
        sf = float(p["SF"])
        offset, sw = float(p["OFFSET"]), float(p["SW_p"])
        axes.append(offset - np.arange(npts) * (sw / sf / npts))
    x_f2, y_f1 = axes
    return x_f2, y_f1, Z.astype(float)


def render_2d(spec: dict) -> Figure:
    style = STYLES[spec.get("style", "article")]
    rc = _rc(spec)
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

    with plt.rc_context(rc):
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
        mode_c = spec.get("contour_mode", "contour")   # contour|density|filled|both
        cs = None
        if mode_c in ("density", "filled", "both"):
            ax.contourf(x, y, Z, levels=levels, cmap=cmap)
        if mode_c in ("contour", "both"):
            colors = None if spec.get("contour_colored", True) else "0.2"
            cs = ax.contour(x, y, Z, levels=levels,
                            cmap=cmap if colors is None else None,
                            colors=colors, linewidths=0.7)
            if spec.get("contour_values"):             # print the level values
                ax.clabel(cs, fontsize=6, fmt="%.2f")
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
            if (spec.get("legend_top_loc", "upper right") not in ("none", "off")
                    and any(l.get_label() and not l.get_label().startswith("_")
                            for l in ax_top.lines)):
                ax_top.legend(loc=spec.get("legend_top_loc", "upper right"))

        if ax_right is not None:
            proj = reducer(Z[np.ix_(ysel, xsel)], 1)
            proj = proj / np.abs(proj).max()
            ax_right.plot(proj, y[ysel], lw=1.0)
            ax_right.set_xticks([])
            plt.setp(ax_right.get_yticklabels(), visible=False)
            for s in ("right", "top", "bottom"):
                ax_right.spines[s].set_visible(False)

        # reference lines: CS axis, quadrupolar-induced-shift axis, iso guides
        for sl in list(spec.get("slopes", [])) + list(spec.get("iso_lines", [])):
            xs = np.array([min(xlim), max(xlim)])
            ax.plot(xs, sl["slope"] * xs + sl.get("intercept", 0.0),
                    color=sl.get("color", "k"), lw=sl.get("linewidth", 0.9),
                    ls=sl.get("linestyle", "-"), label=sl.get("label"))
        if (spec.get("legend_loc", "lower left") not in ("none", "off")
                and any(l.get_label() and not l.get_label().startswith("_")
                        for l in ax.lines)):
            ax.legend(loc=spec.get("legend_loc", "lower left"),
                      fontsize=spec.get("legend_fontsize", 7))

        if spec.get("annotation"):
            ax.text(0.04, 0.94, spec["annotation"], transform=ax.transAxes,
                    fontsize=style["rc"]["font.size"] + 2, va="top")
        if spec.get("title"):
            (ax_top or ax).set_title(spec["title"])
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
    data = load_series(spec)
    mode = spec.get("mode", "satrec")
    with plt.rc_context(_rc(spec)):
        fig, ax = plt.subplots(figsize=_figsize(spec))
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
        if spec.get("title"):
            ax.set_title(spec["title"])
        if spec.get("ylim"):
            lo, hi = spec["ylim"]; ax.set_ylim(min(lo, hi), max(lo, hi))
        if spec.get("grid"):
            ax.grid(True, ls=spec.get("grid_style", ":"), lw=0.6, color="0.7")
        _apply_legend(ax, spec, has_labels=True, default_loc="lower right")
        fig.tight_layout()
        return fig


# ---------------------------------------------------------------------------
# batch grid: small-multiples publication figure from already-fitted spectra
# (a deconvolution grid, or a composition series highlighting one component)

def render_batch_grid(spec: dict) -> Figure:
    """One panel per already-fitted spectrum: experiment + total fit +
    (optionally) its components, laid out as a small-multiples grid.

    ``panels``: ``[{"recipe": "path/to/fit.recipe.json", "title": "…"}, …]``.
    Each recipe is resolved through ``load_trace`` exactly like a 1D trace, so
    it follows the fit's own ``source_path`` for the experimental data
    (kernel or not, data-less or not) — nothing here re-implements that.

    A panel may instead omit ``"recipe"`` and give ``"data_path"`` (a raw
    spectrum with no matching fit yet — series_grid.py produces these when a
    batch CSV names a sample it couldn't pair with a saved .recipe.json): it
    draws as experiment-only, no total/components. A panel with neither key
    draws as an empty "no data located" placeholder rather than raising, so
    one unresolved sample never blocks rendering the rest of the grid.

    ``component_mode``: ``"fill"`` (filled + outlined, the classic
    deconvolution look), ``"dashed"`` (outline only), or ``"hidden"`` (total +
    experiment only). ``shade_only``: a list of site indices to draw and skip
    the rest — combined with ``component_mode="fill"``, this is the
    "composition series highlighting one component" style (draw every other
    component invisible, shade just the one of interest). ``peak_labels``:
    ``None`` | ``"position"`` | ``"label"`` | ``"position+pct"`` (needs
    ``larmor.quantify`` — adds each shown component's integrated population %).
    """
    from pathlib import Path as _Path
    from larmor.recipe import Recipe

    panels_spec = spec.get("panels", [])
    n = len(panels_spec)
    if n == 0:
        raise ValueError("batch_grid needs at least one panel (a fitted recipe)")
    cols = int(spec.get("cols") or 0) or max(1, int(np.ceil(np.sqrt(n))))
    rows = int(np.ceil(n / cols))
    comp_mode = spec.get("component_mode", "fill")     # fill|dashed|hidden
    shade_only = set(spec.get("shade_only") or [])
    comp_alpha = float(spec.get("component_alpha", 0.35))
    show_total = spec.get("show_total", True)
    show_exp = spec.get("show_experiment", True)
    peak_labels = spec.get("peak_labels")              # None|position|label|position+pct
    label_fmt = spec.get("peak_label_fmt", "{pos:.1f}")
    x_is_ppm = spec.get("x_is_ppm", True)
    pw, ph = spec.get("panel_size", (3.0, 2.4))

    with plt.rc_context(_rc(spec)):
        fsz = plt.rcParams["font.size"]
        fig, axes = plt.subplots(rows, cols, figsize=(pw * cols, ph * rows),
                                 squeeze=False)
        legend_handles: dict[int, tuple] = {}   # site index -> (handle, label)
        for idx, p in enumerate(panels_spec):
            ax = axes[idx // cols][idx % cols]
            recipe = Recipe.load(p["recipe"]) if p.get("recipe") else None
            x_total = y_total = None
            if recipe is not None:
                x_total, y_total, _ = load_trace({"recipe": p["recipe"], "part": "total"})
            title = (p.get("title") or (recipe.sample if recipe else "")
                    or (_Path(p["recipe"]).stem if p.get("recipe") else p.get("sample", "?")))
            nucleus = recipe.nucleus if recipe is not None else p.get("nucleus", "")

            pops = None
            if recipe is not None and peak_labels == "position+pct":
                try:
                    from larmor.quantify import quantify
                    q = quantify(recipe, getattr(recipe, "fit_window_ppm", None))
                    pops = {int(r["site"][1:]): r["fraction_pct"] for r in q["rows"]}
                except Exception:
                    pops = None

            ex = ey = None
            src = recipe.source_path if recipe is not None else p.get("data_path", "")
            if show_exp and src:
                try:
                    ex, ey, _ = load_trace({"path": src})
                    ax.plot(ex, ey, color="#222", lw=0.9, label="_experiment")
                except Exception:
                    ex = ey = None

            if recipe is not None and comp_mode != "hidden":
                for i, site in enumerate(recipe.sites):
                    if shade_only and i not in shade_only:
                        continue
                    xi, yi, _ = load_trace({"recipe": p["recipe"], "part": "site",
                                            "site": i})
                    col = site_color(i)
                    if comp_mode == "fill":
                        ax.fill_between(xi, yi, color=col, alpha=comp_alpha, lw=0)
                        (line,) = ax.plot(xi, yi, color=col, lw=0.8)
                    else:
                        (line,) = ax.plot(xi, yi, color=col, lw=1.0, ls="--")
                    legend_handles.setdefault(i, (line, site.label or f"s{i}"))
                    if peak_labels:
                        pos = site.params.get("isotropic_chemical_shift_ppm")
                        if pos is None:
                            continue
                        txt = label_fmt.format(pos=float(pos.value))
                        if peak_labels == "label":
                            txt = site.label or txt
                        elif peak_labels == "position+pct" and pops and i in pops:
                            txt = f"{txt}\n{pops[i]:.0f}%"
                        ax.annotate(txt, (float(pos.value), float(np.max(yi))),
                                    textcoords="offset points", xytext=(0, 3),
                                    ha="center", fontsize=max(fsz - 2, 5))

            if recipe is not None and show_total:
                ax.plot(x_total, y_total, color="#c0392b", lw=1.2,
                        ls=(0, (4, 2)), label="_fit")

            if recipe is not None and peak_labels:  # headroom so labels clear the title
                lo, hi = ax.get_ylim()
                ax.set_ylim(lo, hi * 1.14)

            x_ref = x_total if x_total is not None else ex
            if x_is_ppm:
                if spec.get("xlim"):
                    hi, lo = spec["xlim"]
                    ax.set_xlim(max(hi, lo), min(hi, lo))
                elif x_ref is not None:
                    ax.set_xlim(float(np.max(x_ref)), float(np.min(x_ref)))
            ax.set_yticks([])
            for s in ("left", "right", "top"):
                ax.spines[s].set_visible(False)
            if x_ref is None:            # neither a fit nor data was found yet
                ax.text(0.5, 0.5, "no data located", ha="center", va="center",
                        transform=ax.transAxes, color="#999", fontsize=max(fsz - 1, 5))
                ax.set_xticks([])
            ax.set_title(title, fontsize=fsz, pad=2)
            if x_ref is not None and (spec.get("xlabel") or (x_is_ppm and nucleus)):
                ax.set_xlabel(spec.get("xlabel") or nucleus_xlabel(nucleus),
                             fontsize=fsz)

        for j in range(n, rows * cols):
            axes[j // cols][j % cols].set_visible(False)

        if spec.get("shared_scale"):
            live = [a for a in fig.axes if a.get_visible()]
            los = [a.get_xlim()[1] for a in live]; his = [a.get_xlim()[0] for a in live]
            for a in live:
                a.set_xlim(max(his), min(los))

        if spec.get("legend", True) and legend_handles:
            handles = [h for h, _ in legend_handles.values()]
            labels = [lbl for _, lbl in legend_handles.values()]
            fig.legend(handles, labels, loc="lower center",
                      ncol=min(len(handles), 8), bbox_to_anchor=(0.5, -0.02),
                      frameon=False, fontsize=max(fsz - 1, 5))
        if spec.get("suptitle"):
            fig.suptitle(spec["suptitle"])
        fig.tight_layout()
        return fig


# ---------------------------------------------------------------------------
# species distribution: a 100%-stacked bar of population vs. composition

def render_species_bar(spec: dict) -> Figure:
    """A 100%-stacked bar chart: one bar per composition point, segments are
    named species/sites in a fixed category order and colour (never re-sorted
    by value — identity stays with the category, not its rank).

    ``categories``: ``[label, …]`` — the x-axis (composition) points.
    ``series``: ``[{"label": "Q2", "values": [...]}, …]`` — one entry per
    species, ``values`` aligned with ``categories``, in **stacking order**
    (first at the bottom). Values need not already sum to 100 — each bar is
    normalized to its own total unless ``normalize`` is False.
    """
    categories = [str(c) for c in spec.get("categories", [])]
    series = spec.get("series", [])
    if not categories or not series:
        raise ValueError("species_bar needs `categories` and `series`")
    normalize = spec.get("normalize", True)
    vals = np.array([[float(v) for v in s["values"]] for s in series], dtype=float)
    if normalize:
        totals = vals.sum(axis=0)
        totals[totals == 0] = 1.0
        vals = vals / totals * 100.0

    with plt.rc_context(_rc(spec)):
        fig, ax = plt.subplots(figsize=spec.get("figsize", (5.0, 3.2)))
        x = np.arange(len(categories))
        bottom = np.zeros(len(categories))
        width = float(spec.get("bar_width", 0.7))
        for i, s in enumerate(series):
            col = s.get("color") or site_color(i)
            ax.bar(x, vals[i], width, bottom=bottom, color=col,
                  label=s.get("label", f"series {i}"),
                  edgecolor=spec.get("edge_color", "white"), linewidth=0.6)
            if spec.get("value_labels"):
                for xi, (v, b) in enumerate(zip(vals[i], bottom)):
                    if v >= float(spec.get("value_label_min_pct", 6)):
                        ax.text(xi, b + v / 2, f"{v:.0f}", ha="center",
                               va="center", fontsize=plt.rcParams["font.size"] - 2,
                               color="white")
            bottom += vals[i]
        ax.set_xticks(x)
        ax.set_xticklabels(categories, rotation=spec.get("xtick_rotation", 0),
                           ha="right" if spec.get("xtick_rotation") else "center")
        ax.set_ylim(0, 100 if normalize else spec.get("ymax") or bottom.max() * 1.05)
        ax.set_ylabel(spec.get("ylabel", "population (%)" if normalize else "value"))
        if spec.get("xlabel"):
            ax.set_xlabel(spec["xlabel"])
        for s in ("right", "top"):
            ax.spines[s].set_visible(False)
        if "legend_loc" not in spec:      # default: a row of swatches above the bars
            ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.18),
                      ncol=min(len(series), 6), frameon=False,
                      fontsize=plt.rcParams["legend.fontsize"])
        else:
            _apply_legend(ax, spec, has_labels=True)
        if spec.get("title"):
            ax.set_title(spec["title"])
        fig.tight_layout()
        return fig


# ---------------------------------------------------------------------------
RENDERERS = {"1d": render_1d, "2d": render_2d, "series": render_series,
            "batch_grid": render_batch_grid, "species_bar": render_species_bar}


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
