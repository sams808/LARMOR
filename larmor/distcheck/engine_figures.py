"""Engine stage ``engine-figures`` of the distribution check: every
renderer of larmor.figures on a minimal valid spec, the drawn geometry
checked against the numbers it was drawn from, every kind exported to PNG /
SVG / PDF and each file checked from its first to its last bytes; the
shipped two-panel figure spec rendered once its relative paths are resolved
against the examples folder. Qt-free (matplotlib Agg), every file under
``ctx["out"]``.
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path

import numpy as np

from larmor.distcheck.engines import _check, _examples, _expno, _out, _small_kernels

__all__ = ["figures_stage"]


#: what a complete file of each format starts and ends with (a truncated
#: write -- a backend dying half-way -- keeps its magic but loses its tail)
_FIG_MAGIC = {"png": (lambda b: b.startswith(b"\x89PNG\r\n\x1a\n") and b"IEND" in b[-16:],
                      "\\x89PNG ... IEND"),
              "pdf": (lambda b: b.startswith(b"%PDF-") and b"%%EOF" in b[-64:], "%PDF- ... %%EOF"),
              "svg": (lambda b: b"<svg" in b and b"</svg>" in b[-64:], "<svg ... </svg>")}


def _fig_jcamp(path: Path, values: dict) -> None:
    """A minimal TopSpin parameter file (procs / proc2s)."""
    lines = ["##TITLE= distcheck synthetic", "##JCAMPDX= 5.0"]
    lines += [f"##${k}= {v}" for k, v in values.items()]
    path.write_text("\n".join(lines + ["##END="]) + "\n", encoding="ascii")


def _fig_synthetic_2d(folder: Path, f2=(60.0, 6.0), f1=(40.0, 4.0)) -> Path:
    """A processed Bruker 2D EXPNO (pdata/1/2rr + procs + proc2s, one
    submatrix, little-endian int32) holding one Gaussian peak at (F2, F1)
    ppm: what figures.load_2d reads when the bundled 3QMAS is absent."""
    pdata = folder / "pdata" / "1"
    pdata.mkdir(parents=True, exist_ok=True)
    sf = 130.3
    dims = {"procs": (128, 150.0, 200.0), "proc2s": (64, 120.0, 160.0)}  # SI, OFFSET, span (ppm)
    axes = {}
    for name, (si, offset, span) in dims.items():
        _fig_jcamp(pdata / name, {"BYTORDP": 0, "DTYPP": 0, "NC_proc": 0, "OFFSET": offset,
                                  "SF": sf, "SI": si, "SW_p": span * sf, "XDIM": si})
        axes[name] = offset - np.arange(si) * (span / si)
    x, y = axes["procs"], axes["proc2s"]
    z = np.exp(-0.5 * (((x[None, :] - f2[0]) / f2[1]) ** 2 + ((y[:, None] - f1[0]) / f1[1]) ** 2))
    np.round(z * 2 ** 28).astype("<i4").tofile(pdata / "2rr")
    return folder


def _fig_absolute_spec(ex: Path, out: Path) -> Path:
    """The shipped two-panel figure spec, its recipes copied into ``out``
    with ABSOLUTE paths: the spec names its recipes, and the recipes their
    data, relative to the folder that holds examples/ -- not where a user's
    shell stands (the resolution gui_tools applies, Qt-free)."""
    spec = json.loads((ex / "pCABS2-4_fits.figure.json").read_text(encoding="utf-8"))
    for panel in spec.get("panels", []):
        name = Path(panel["recipe"]).name
        rec = json.loads((ex / name).read_text(encoding="utf-8"))
        sp = Path(str(rec.get("source_path") or ""))
        if not sp.is_absolute() and sp.parts and sp.parts[0] == "examples":
            sp = ex.joinpath(*sp.parts[1:])
        rec["source_path"] = str(sp.resolve())
        copy = out / name
        copy.write_text(json.dumps(rec, indent=1), encoding="utf-8")
        panel["recipe"] = str(copy)
    target = out / "pCABS2-4_fits.figure.json"
    target.write_text(json.dumps(spec, indent=1), encoding="utf-8")
    return target


def _fig_saved_fit(out: Path, sample: str, pos, amps, seed: int) -> str:
    """A saved non-kernel 11B fit with its CSV spectrum (tests/test_figures)."""
    from larmor import engine
    from larmor.recipe import Param, Recipe, SiteModel

    rec = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sample=sample, sites=[
        SiteModel(model="gauss_lor", label="AB"[i], params={
            "isotropic_chemical_shift_ppm": Param(p), "shift_fwhm_ppm": Param(5.0),
            "amplitude": Param(a), "gl": Param(1.0, vary=False)})
        for i, (p, a) in enumerate(zip(pos, amps))])
    x = np.linspace(-40.0, 40.0, 400)
    _, y, _ = engine.simulate(rec, exp_ppm=x)
    y = np.asarray(y, float) + np.random.default_rng(seed).normal(0.0, 1.0, x.size)
    csv = out / f"{sample}_raw.csv"
    csv.write_text("# nucleus = 11B\n# larmor_MHz = 160\n"
                   + "\n".join(f"{xi:.4f} {yi:.4f}" for xi, yi in zip(x, y)), encoding="utf-8")
    rec.source_path = str(csv)
    p = out / f"{sample}.recipe.json"
    rec.save(p)
    return str(p)


def figures_stage(say, ctx):
    """larmor.figures: every renderer in figures.RENDERERS (1d, 2d, series,
    batch_grid, species_bar, infinite_field) on a minimal valid spec, the
    drawn geometry checked against the numbers it was built from, then
    exported to PNG / SVG / PDF (each file's head and tail checked); the 2D
    map through figures.load_2d on the bundled 3QMAS (a synthetic 2rr
    otherwise), and the shipped pCABS2-4_fits.figure.json rendered once its
    relative paths are resolved against the examples folder."""
    # the shipped grid simulates Czjzek / amorphous fits: reduced kernels
    # in quick mode, the user's kernel disk cache untouched
    with _small_kernels(ctx, say):
        _figures(say, ctx)


def _figures(say, ctx):
    import matplotlib.pyplot as plt

    from larmor import figures
    from larmor.io import bruker
    from larmor.recipe import Recipe

    out = _out(ctx, "figures")
    ex = _examples(ctx)
    built: dict = {}

    # ---- 1d: normalise / scale / offset arithmetic on inline traces, and an
    # EXPNO trace drawn exactly as the reader returns it
    x = np.linspace(-50.0, 50.0, 501)
    spec = {"kind": "1d", "xlim": [40, -40], "traces": [
        {"data": {"x": x.tolist(), "y": (3.0 * np.exp(-(x / 8.0) ** 2)).tolist()}, "label": "a",
         "normalize": True},
        {"data": {"x": x.tolist(), "y": (2.0 * np.exp(-((x - 10.0) / 4.0) ** 2)).tolist()}, "label": "b",
         "offset": 1.2, "scale": 0.5, "normalize": True}]}
    e1 = _expno(ctx, 3616)
    if e1 is not None:
        spec["traces"].append({"path": str(e1), "label": "3616", "scale": 1e-6})
    fig = figures.render(spec)
    lines = list(fig.axes[0].lines)
    ya, yb = lines[0].get_ydata(), lines[1].get_ydata()
    _check(len(lines) == len(spec["traces"]) and abs(ya.max() - 1.0) < 1e-9 and abs(yb.max() - 1.7) < 1e-9
           and abs(lines[1].get_xdata()[np.argmax(yb)] - 10.0) < 0.2 + 1e-9,
           f"1d: {len(lines)} lines, normalised max {ya.max():.6f} (1), scaled+offset max {yb.max():.6f} (1.7)")
    _check(fig.axes[0].get_xlim() == (40.0, -40.0), f"1d: xlim {fig.axes[0].get_xlim()} (a ppm axis runs 40 -> -40)")
    note = "inline traces"
    if e1 is not None:
        ref = np.sort(np.asarray(bruker.read(e1).data, float)) * 1e-6
        got = np.sort(np.asarray(lines[2].get_ydata(), float))
        _check(got.size == ref.size and np.allclose(got, ref, rtol=1e-9, atol=1e-12),
               f"1d: the 3616 trace drawn ({got.size} points, max {got.max():.4g}) is not the 1r "
               f"bruker.read returns ({ref.size} points, max {ref.max():.4g})")
        note += f" + the 3616 EXPNO ({got.size} points = its 1r)"
    else:
        note += " (no bundled examples: no EXPNO trace)"
    plt.close(fig)
    built["1d"] = spec
    say(f"1d: {note}; normalise -> 1.000, x0.5 + 1.2 -> {yb.max():.3f}, ppm axis reversed")

    # ---- 2d: load_2d's axes against the parameter files, the projections
    # against the matrix
    e3 = _expno(ctx, 3620)
    if e3 is None:
        e3 = _fig_synthetic_2d(out / "syn2d", f2=(60.0, 6.0), f1=(40.0, 4.0))
        src2d = "a synthetic 2rr (no bundled 3QMAS): one peak at F2 60 / F1 40 ppm"
    else:
        src2d = "the bundled 27Al 3QMAS (3620)"
    x2, y2, Z = figures.load_2d(e3)
    for name, ax_ in (("procs", x2), ("proc2s", y2)):
        txt = (e3 / "pdata" / "1" / name).read_text(encoding="latin-1")
        prm = {k: float(re.search(rf"##\${k}=\s*([-+0-9.eE]+)", txt).group(1)) for k in ("OFFSET", "SW_p", "SF", "SI")}
        step = prm["SW_p"] / prm["SF"] / prm["SI"]
        _check(ax_.size == int(prm["SI"]) and abs(ax_[0] - prm["OFFSET"]) < 1e-9
               and np.allclose(np.diff(ax_), -step, rtol=1e-9),
               f"load_2d {name} axis: {ax_.size} points from {ax_[0]:.4f} ppm, step {np.diff(ax_)[0]:.5f} "
               f"(file: SI {prm['SI']:.0f}, OFFSET {prm['OFFSET']:.4f}, step {-step:.5f})")
    _check(Z.shape == (y2.size, x2.size) and np.isfinite(Z).all() and Z.max() > 0, f"load_2d matrix {Z.shape}")
    i1, i2 = np.unravel_index(int(np.argmax(Z)), Z.shape)
    if e3.name == "syn2d":
        _check(abs(x2[i2] - 60.0) <= abs(x2[1] - x2[0]) and abs(y2[i1] - 40.0) <= abs(y2[1] - y2[0]),
               f"synthetic 2rr read back with its peak at F2 {x2[i2]:.2f} / F1 {y2[i1]:.2f} ppm (60 / 40)")
    spec = {"kind": "2d", "path": str(e3), "nucleus": "27Al", "levels": {"n": 8},
            "slopes": [{"slope": 1.0, "intercept": 0.0, "label": "CS axis"}]}
    fig = figures.render(spec)
    ax, top, right = fig.axes[0], fig.axes[1], fig.axes[2]
    tx, ty = top.lines[0].get_xdata(), top.lines[0].get_ydata()
    rx, ry = right.lines[0].get_xdata(), right.lines[0].get_ydata()
    # skyline projections: their maxima sit on the column / row of the matrix maximum
    _check(abs(ty.max() - 1.0) < 1e-9 and tx[np.argmax(ty)] == x2[i2] and ry[np.argmax(rx)] == y2[i1],
           f"2d projections peak at F2 {tx[np.argmax(ty)]:.2f} / F1 {ry[np.argmax(rx)]:.2f} ppm; "
           f"the matrix maximum is at {x2[i2]:.2f} / {y2[i1]:.2f}")
    _check(ax.get_xlim() == (float(x2.max()), float(x2.min())) and "27" in ax.get_xlabel()
           and r"\delta_1" in ax.get_ylabel() and len(ax.collections) >= 1,
           f"2d axes: xlim {ax.get_xlim()}, labels {ax.get_xlabel()!r} / {ax.get_ylabel()!r}")
    plt.close(fig)
    built["2d"] = spec
    say(f"2d: {src2d}; {Z.shape} matrix, axes = OFFSET - i*SW_p/SF/SI of procs / proc2s, "
        f"projections peak at the matrix maximum ({x2[i2]:.1f}, {y2[i1]:.1f}) ppm")

    # ---- series: a saturation recovery with T1 = 2 s; the T1 the figure prints is its fit's
    tau = np.array([0.0, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0])
    yrec = (1.0 - np.exp(-tau / 2.0)) + np.random.default_rng(5).normal(0.0, 0.003, tau.size)
    spec = {"kind": "series", "mode": "satrec", "data": {"x": tau.tolist(), "y": yrec.tolist()}}
    fig = figures.render(spec)
    ax = fig.axes[0]
    pts = ax.lines[0]
    m = re.search(r"\$T_1\$ = ([0-9.eE+-]+)", " ".join(t.get_text() for t in ax.texts))
    t1 = float(m.group(1)) if m else float("nan")
    _check(np.allclose(pts.get_xdata(), tau) and np.allclose(pts.get_ydata(), yrec),
           "series: the points drawn are not the data given")
    # 0.3 % noise on 11 delays: the fitted T1 lands within ~1 % (3 % allowed)
    _check(abs(t1 - 2.0) < 0.06 and ax.get_xscale() == "log" and len(ax.lines) == 2,
           f"series: printed T1 {t1} s for a 2.000 s recovery ({len(ax.lines)} lines, {ax.get_xscale()} axis)")
    plt.close(fig)
    built["series"] = spec
    say(f"series: satrec through 11 delays (0.3 % noise) -> the figure prints T1 = {t1:g} s (truth 2 s)")

    # ---- batch_grid: the shipped two-panel spec (paths resolved), or two synthetic fits
    t0 = time.perf_counter()
    if ex is not None:
        spec_path = _fig_absolute_spec(ex, out)
        spec = json.loads(spec_path.read_text(encoding="utf-8"))
        what = f"the shipped {spec_path.name}"
    else:
        paths = [_fig_saved_fit(out, f"g{k}", [10.0, -5.0], [80.0, 40.0 + 20.0 * k], k) for k in range(2)]
        spec = {"kind": "batch_grid", "panels": [{"recipe": p, "title": f"g{k}", "xlim": [30.0, -30.0]}
                                                for k, p in enumerate(paths)],
                "cols": 2, "legend": False, "peak_labels": "label"}
        what = "two synthetic saved fits (no bundled examples)"
    fig = figures.render(spec)
    axes = [a for a in fig.axes if a.get_visible()]
    _check(len(axes) == len(spec["panels"]), f"batch_grid: {len(axes)} panels for {len(spec['panels'])} recipes")
    for ax, panel in zip(axes, spec["panels"]):
        rec = Recipe.load(panel["recipe"])
        lab = {ln.get_label(): ln for ln in ax.lines}
        _check(ax.get_title() == panel["title"] and ax.get_xlim() == (max(panel["xlim"]), min(panel["xlim"])),
               f"batch_grid panel {panel['title']!r}: title {ax.get_title()!r}, xlim {ax.get_xlim()}")
        _check("_experiment" in lab and "_fit" in lab, f"batch_grid panel {ax.get_title()!r}: lines {sorted(lab)}")
        # the experiment drawn is the recipe's own source spectrum
        if rec.source_kind == "bruker":
            raw = np.asarray(bruker.read(rec.source_path).data, float)
        else:
            raw = np.loadtxt(rec.source_path, comments="#")[:, 1]
        got = np.asarray(lab["_experiment"].get_ydata(), float)
        _check(got.size == raw.size and np.allclose(np.sort(got), np.sort(raw), rtol=1e-6, atol=1e-9 * np.abs(raw).max()),
               f"batch_grid panel {ax.get_title()!r}: experiment {got.size} points, max {got.max():.4g}; "
               f"its source has {raw.size}, max {raw.max():.4g}")
        # the components drawn add up to the total drawn
        comps = [ln for ln in ax.lines if ln.get_label() not in ("_experiment", "_fit")]
        tot = np.asarray(lab["_fit"].get_ydata(), float)
        ssum = np.sum([np.interp(lab["_fit"].get_xdata(), c.get_xdata(), c.get_ydata()) for c in comps], axis=0)
        _check(len(comps) == len(rec.sites) and np.allclose(ssum, tot, atol=1e-6 * np.abs(tot).max()),
               f"batch_grid panel {ax.get_title()!r}: {len(comps)} component lines for {len(rec.sites)} sites, "
               f"sum - total max {np.abs(ssum - tot).max():.3g}")
        texts = sorted(t.get_text() for t in ax.texts)
        _check(texts == sorted(s.label for s in rec.sites),
               f"batch_grid panel {ax.get_title()!r}: peak labels {texts}, sites {[s.label for s in rec.sites]}")
    plt.close(fig)
    built["batch_grid"] = spec
    say(f"batch_grid: {what} in {time.perf_counter() - t0:.1f} s; per panel the title and x range of the "
        f"spec, the experiment = its source data, components summing to the total, labels = the sites")

    # ---- species_bar: every bar is its values normalised to 100 %
    vals = [[10.0, 20.0, 5.0], [30.0, 20.0, 15.0], [60.0, 60.0, 30.0]]
    spec = {"kind": "species_bar", "categories": ["P-5", "P-10", "P-20"],
            "series": [{"label": f"Q{i}", "values": v} for i, v in enumerate(vals)]}
    fig = figures.render(spec)
    ax = fig.axes[0]
    heights = np.array([[p.get_height() for p in c.patches] for c in ax.containers])
    bottoms = np.array([[p.get_y() for p in c.patches] for c in ax.containers])
    want = np.array(vals) / np.array(vals).sum(axis=0) * 100.0
    _check(heights.shape == want.shape and np.allclose(heights, want) and np.allclose(bottoms[1:], np.cumsum(want, 0)[:-1])
           and np.allclose(heights.sum(0), 100.0) and ax.get_ylim() == (0.0, 100.0)
           and [t.get_text() for t in ax.get_xticklabels()] == spec["categories"],
           f"species_bar: segment heights {np.round(heights, 2).tolist()}, expected {np.round(want, 2).tolist()}")
    plt.close(fig)
    built["species_bar"] = spec
    say(f"species_bar: 3 x 3 stacked segments = values / column total x 100 (e.g. Q0 at P-5 {want[0, 0]:.2f} %)")

    # ---- infinite_field: exact lines dcg = d_iso + b / nu0^2 -> the starred intercepts
    truth = {"glass A": (-12.0, -2.0e4), "glass B": (3.5, -6.0e4)}
    nu = [96.3, 128.4, 160.5, 192.6]
    spec = {"kind": "infinite_field", "nucleus": "11B", "spin": 1.5, "eta": 0.7, "samples": [
        {"label": k, "points": [[n, d + b / n ** 2, 0.1] for n in nu]} for k, (d, b) in truth.items()]}
    fig = figures.render(spec)
    stars = [ln for ln in fig.axes[0].lines if ln.get_marker() == "*"]
    got = [float(s.get_ydata()[0]) for s in stars]
    _check(len(stars) == 2 and np.allclose(got, [d for d, _ in truth.values()], atol=1e-6)
           and all(float(s.get_xdata()[0]) == 0.0 for s in stars),
           f"infinite_field: intercepts drawn at {got} ppm, truth {[d for d, _ in truth.values()]}")
    plt.close(fig)
    built["infinite_field"] = spec
    say(f"infinite_field: two exact 4-field lines -> starred intercepts {got[0]:.3f} / {got[1]:.3f} ppm (truth -12 / 3.5)")

    # ---- every kind of the registry covered, and exported
    missing = set(figures.RENDERERS) - set(built)
    _check(not missing, f"figure kinds with no check here: {sorted(missing)} (add a spec to this stage)")
    t0 = time.perf_counter()
    sizes = []
    for kind, spec in built.items():
        fmts = ("png", "svg", "pdf") + (("json",) if kind == "infinite_field" else ())
        for p in figures.export(spec, out / f"fig_{kind}", formats=fmts, dpi=90):
            fmt = Path(p).suffix[1:]
            data = Path(p).read_bytes()
            if fmt == "json":
                _check(json.loads(data.decode("utf-8")) == json.loads(json.dumps(spec)),
                       f"{kind}: the JSON sidecar is not the spec")
                continue
            ok, magic = _FIG_MAGIC[fmt]
            _check(ok(data) and len(data) > 1000,
                   f"{kind}.{fmt}: {len(data)} bytes, starts {data[:12]!r}, ends {data[-12:]!r} (want {magic})")
            sizes.append(len(data))
    say(f"export: {len(built)} kinds x png/svg/pdf ({len(sizes)} files, {min(sizes)}-{max(sizes)} bytes, head and "
        f"tail checked) + the JSON sidecar in {time.perf_counter() - t0:.1f} s")
