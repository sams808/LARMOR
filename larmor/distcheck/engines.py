"""Engine stages of the distribution check: every Qt-free capability on
synthetic data (readers, processing, fitting, 2D, batch, series, error
tools, exports, provenance, help rendering ...). Each stage is
``(name, fn(say, ctx))``; see ``larmor.distcheck`` for the contract.

Every stage is deterministic (seeded), touches no network, writes only under
``ctx["out"]`` and uses the bundled example data (``ctx["examples"]``) when
it is present, synthetic data otherwise. In ``--quick`` mode the Czjzek /
amorphous kernels are built on a reduced (Cq, eta) grid (``_small_kernels``)
so the whole set runs in a minute or two; the first stage that needs
mrsimulator pays its import (~20 s on a cold machine) and says so.
"""
from __future__ import annotations

import contextlib
import os
import time
from pathlib import Path

import numpy as np

__all__ = ["stages"]

NUCLEUS, LARMOR_MHZ, MAS_HZ = "27Al", 130.3, 20000.0


# ------------------------------------------------------------ helpers
def _out(ctx, name: str) -> Path:
    """A fresh sub-folder of the scratch folder for one stage."""
    p = Path(ctx["out"]) / name
    p.mkdir(parents=True, exist_ok=True)
    return p


def _examples(ctx) -> Path | None:
    ex = ctx.get("examples")
    return Path(ex) if ex and (Path(ex) / "pCABS2-4").is_dir() else None


def _expno(ctx, n: int) -> Path | None:
    ex = _examples(ctx)
    p = ex / "pCABS2-4" / str(n) if ex else None
    return p if p is not None and p.is_dir() else None


@contextlib.contextmanager
def _small_kernels(ctx, say=None):
    """Czjzek-family and amorphous kernels on a reduced (Cq, eta) grid in
    quick mode, the kernel disk cache left untouched in both modes (the
    check must not write outside its scratch folder), and the in-memory
    kernel cache cleared on both sides so no stage inherits a reduced
    kernel."""
    from larmor import engine
    from larmor.models import quadrupolar

    saved = dict(engine.KERNEL_SETTINGS)
    amorph = (quadrupolar.AMORPH_N_CQ, quadrupolar.AMORPH_N_ETA)
    had_env = "LARMOR_NO_SESSION" in os.environ
    os.environ.setdefault("LARMOR_NO_SESSION", "1")
    if ctx.get("quick"):
        engine.KERNEL_SETTINGS.update({"npts": 512, "n_cq": 12, "n_eta": 3})
        quadrupolar.AMORPH_N_CQ, quadrupolar.AMORPH_N_ETA = 12, 3
        if say is not None:
            say("quick mode: Czjzek kernels on a 12 x 3 (Cq, eta) grid, 512 points; "
                "amorphous on 12 x 3")
    engine.clear_kernel_cache()
    try:
        yield
    finally:
        engine.KERNEL_SETTINGS.clear()
        engine.KERNEL_SETTINGS.update(saved)
        quadrupolar.AMORPH_N_CQ, quadrupolar.AMORPH_N_ETA = amorph
        engine.clear_kernel_cache()
        if not had_env:
            os.environ.pop("LARMOR_NO_SESSION", None)


def _import_mrsimulator(say) -> None:
    """mrsimulator is imported lazily by the quadrupolar models; its import
    is the one slow step of these stages, so it is timed and named."""
    import sys
    if "mrsimulator" in sys.modules:
        return
    t0 = time.perf_counter()
    import mrsimulator  # noqa: F401
    say(f"import mrsimulator {mrsimulator.__version__}: {time.perf_counter() - t0:.1f} s")


def _ppm_axis(n: int = 400, hi: float = 150.0, lo: float = -100.0) -> np.ndarray:
    """A Bruker-like descending ppm axis."""
    return np.linspace(hi, lo, n)


def _two_line_recipe(amp_a: float = 100.0, amp_b: float = 40.0):
    """The synthetic truth: two Gauss/Lorentz 27Al lines (no kernel needed)."""
    from larmor.recipe import Param, Recipe, SiteModel

    return Recipe(nucleus=NUCLEUS, larmor_frequency_MHz=LARMOR_MHZ, spin_rate_Hz=MAS_HZ,
                  sample="distcheck", sites=[
                      SiteModel(model="gauss_lor", label="Al4", params={
                          "isotropic_chemical_shift_ppm": Param(60.0),
                          "shift_fwhm_ppm": Param(14.0), "amplitude": Param(amp_a),
                          "gl": Param(0.6, vary=False)}),
                      SiteModel(model="gauss_lor", label="Al6", params={
                          "isotropic_chemical_shift_ppm": Param(5.0),
                          "shift_fwhm_ppm": Param(10.0), "amplitude": Param(amp_b),
                          "gl": Param(0.6, vary=False)})])


def _synthetic_1d(n: int = 400, noise: float = 0.5, seed: int = 1,
                  amp_a: float = 100.0, amp_b: float = 40.0):
    """(x_ppm descending, y, truth recipe) for a two-line 27Al spectrum."""
    from larmor import engine

    x = _ppm_axis(n)
    truth = _two_line_recipe(amp_a, amp_b)
    # simulate() returns its own (ascending) axis: resample onto x
    xs, ys, _ = engine.simulate(truth, exp_ppm=x)
    y = np.interp(x, xs, np.asarray(ys, float))
    y = y + np.random.default_rng(seed).normal(0.0, noise, x.size)
    return x, y, truth


def _start_recipe(x=None):
    """A deliberately displaced two-line starting point for the fits."""
    from larmor.recipe import Param

    rec = _two_line_recipe(70.0, 30.0)
    rec.sites[0].params["isotropic_chemical_shift_ppm"] = Param(55.0, min=30.0, max=90.0)
    rec.sites[0].params["shift_fwhm_ppm"] = Param(18.0, min=2.0, max=60.0)
    rec.sites[1].params["isotropic_chemical_shift_ppm"] = Param(8.0, min=-20.0, max=30.0)
    rec.sites[1].params["shift_fwhm_ppm"] = Param(8.0, min=2.0, max=60.0)
    return rec


def _check(cond: bool, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


# ------------------------------------------------------------ readers
def _write_csv_spectrum(path: Path, x, y, nucleus: str = "11B", larmor: float = 160.0) -> Path:
    """The two-line-header CSV the loaders accept (tests/test_series_mode_ui)."""
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"# nucleus = {nucleus}\n# larmor_MHz = {larmor:g}\n")
        for xi, yi in zip(x, y):
            f.write(f"{xi:.5f} {yi:.5f}\n")
    return path


def _synthetic_varian(folder: Path):
    """A one-FID Varian ``.fid`` directory written with nmrglue (binary fid +
    procpar), with the parameters larmor.io.varian reads appended."""
    import nmrglue as ng

    udic = ng.fileiobase.create_blank_udic(1)
    udic[0].update({"size": 256, "complex": True, "encoding": "direct", "sw": 20000.0,
                    "obs": 104.26, "car": 0.0, "label": "Al27", "time": True, "freq": False})
    dic = ng.varian.create_dic(udic)
    t = np.arange(256) / 20000.0
    fid = (np.exp(2j * np.pi * 1500.0 * t) * np.exp(-t / 0.01)).astype(np.complex64)
    ng.varian.write(str(folder), dic, fid, overwrite=True)
    extra = [("tn", '"Al27"', True), ("sfrq", "104.26", False), ("reffrq", "104.26", False),
             ("sw", "20000.0", False), ("rfl", "10000.0", False), ("rfp", "0.0", False),
             ("seqfil", '"s2pul"', True), ("comment", '"distcheck"', True)]
    lines = []
    for name, value, is_str in extra:
        lines.append(f"{name} {2 if is_str else 1} {2 if is_str else 1} 0 0 0 2 1 0 1 64")
        lines.append(f"1 {value}")
        lines.append("0")
    with open(folder / "procpar", "a", encoding="ascii") as f:
        f.write("\n".join(lines) + "\n")
    return fid


def readers(say, ctx):
    """larmor.io: the Bruker reader on the bundled EXPNOs (processed, raw fid
    through the Fourier transform, imaginary channel, acqus metadata, the
    2D ser), load_any on an EXPNO and on a CSV, the dmfit import, the
    text/CSV/JSON/fxmla exports and their round trip, the folder scanner,
    a synthetic Varian dataset."""
    # the dmfit import of amorphous lines simulates them (an mrsimulator
    # kernel): reduced in quick mode, never written to the user's disk cache
    with _small_kernels(ctx, say):
        _readers(say, ctx)


def _readers(say, ctx):
    from larmor import fourier, loader, processing
    from larmor.io import bruker, export, fxmla, scan, spectra, varian
    from larmor.recipe import Recipe

    out = _out(ctx, "readers")
    ex = _examples(ctx)
    if ex is None:
        say("no bundled example data: the Bruker / dmfit readers run on nothing here")
    else:
        e1, e2, e3 = ex / "pCABS2-4" / "3616", ex / "pCABS2-4" / "1118", ex / "pCABS2-4" / "3620"
        d = bruker.read(e1)
        _check(d.ndim == 1 and d.domain == "freq" and d.data.size == 2048,
               f"3616 processed: ndim {d.ndim} {d.domain} {d.data.shape}")
        ppm = d.axes[0].values
        _check(ppm is not None and ppm.size == d.data.size and np.ptp(ppm) > 500.0,
               "3616 ppm axis missing or implausible")
        say(f"3616 1r: {d.data.size} points, {ppm[0]:.1f}..{ppm[-1]:.1f} ppm, "
            f"nucleus {d.meta.get('nucleus')} at {d.meta.get('larmor_MHz'):.3f} MHz, "
            f"MAS {d.meta.get('spin_rate_Hz')} Hz ({d.meta.get('mas_source')})")
        imag = bruker.read_imag(e1)
        _check(imag is not None and imag.size == d.data.size, "read_imag did not return the 1i")
        meta = bruker.read_acqus_meta(e1)
        _check(meta["nucleus"] == "27Al" and abs(meta["larmor_MHz"] - 130.323) < 0.01
               and meta["spin_rate_Hz"] == 26000.0, f"acqus metadata wrong: {meta}")
        raw = bruker.read(e1 / "fid")
        _check(raw.domain == "time" and np.iscomplexobj(raw.data), "the fid did not read as a complex FID")
        say(f"3616 fid: {raw.data.size} complex points after the digital filter; warnings {raw.warnings}")
        fx, fy = fourier.ft1d(raw.data, raw.axes[0].sw_Hz, raw.axes[0].obs_MHz,
                              ops=[{"op": "em", "lb_hz": 100.0}])
        mag = np.abs(fy)
        _check(fy.size >= raw.data.size and np.all(np.isfinite(mag))
               and mag.max() > 4 * np.median(mag), "ft1d of the fid shows no line")
        s = processing.from_bruker_fid(str(e1))
        s = processing.apply(s, [{"op": "em", "lb_hz": 100.0}, {"op": "zf"}, {"op": "ft"},
                                 {"op": "autophase"}])
        _check(s.domain == "freq" and s.y.real.max() > 0.9 * np.abs(s.y).max(),
               "from_bruker_fid -> em/zf/ft/autophase did not give an absorption line")
        say(f"fid -> spectrum: {s.y.size} points, peak at {s.x_ppm[np.argmax(s.y.real)]:.1f} ppm "
            f"(axis relative to the carrier)")
        xp = bruker.read_expno(e1)
        _check(xp.fid is not None and xp.processed is not None and xp.nucleus == "27Al",
               "read_expno lost the fid or the processed data")
        two = bruker.read(e3)
        _check(two.ndim == 2 and two.data.shape == (64, 1024) and two.hyper,
               f"3620 2rr: {two.data.shape}, hyper {two.hyper is not None}")
        ser = bruker.read(e3 / "ser")
        d2 = fourier.ft2d_from_nmrdata(ser)
        _check(ser.domain == "time" and d2.z.ndim == 2 and np.isfinite(d2.z).all()
               and d2.z.max() > 0, "ft2d_from_nmrdata on the 3620 ser failed")
        say(f"3620: 2rr {two.data.shape} with quadrants {sorted(two.hyper)}; ser {ser.data.shape} "
            f"-> ft2d {d2.z.shape}")
        for p in (e1, e1 / "pdata" / "1" / "1r"):
            ppm_l, amp_l, rec, meta_l, warns = loader.load_any(str(p))
            _check(rec["nucleus"] == "27Al" and rec["spin_rate_Hz"] == 26000.0
                   and rec["source_kind"] == "bruker" and len(ppm_l) == 2048 and rec["acquisition"],
                   f"load_any({p.name}) recipe wrong: {rec.get('nucleus')} {rec.get('spin_rate_Hz')}")
        say(f"load_any(EXPNO): sample {rec['sample']!r}, sha {rec['source_sha256'][:12]}, "
            f"{len(rec['acquisition'])} acquisition keys, warnings {warns}")
        for f, nuc, models in ((e1 / "pdata" / "1" / "1r.fxml", "27Al", {"czjzek"}),
                               (e2 / "pdata" / "1" / "25_fit2021.fxml", "11B", {"amorphous", "gauss_lor"})):
            dm = fxmla.read(f)
            rec_dm, notes = fxmla.to_recipe(dm)
            got = {s_.model for s_ in rec_dm.sites}
            _check(rec_dm.nucleus == nuc and got == models and len(rec_dm.sites) == 3,
                   f"dmfit import of {f.name}: {rec_dm.nucleus} {got}")
            say(f"dmfit {f.name}: version {dm.version}, {len(rec_dm.sites)} lines {sorted(got)}, "
                f"{len(notes)} import note(s)")
        entries = scan.list_dir(ex)
        _check(any(t.is_sample and t.name == "pCABS2-4" for t in entries), "list_dir missed the sample folder")
        infos = scan.scan_sample(ex / "pCABS2-4")
        kinds = {i.expno: (i.nucleus, i.ndim, i.kind) for i in infos}
        _check(kinds == {"1118": ("11B", 1, "Single pulse"), "3616": ("27Al", 1, "Single pulse"),
                         "3620": ("27Al", 2, "MQMAS")}, f"scan_sample: {kinds}")
        labels = scan.disambiguate(["a", "a"], [str(e1), str(e2)])
        _check(len(set(labels)) == 2, f"disambiguate left duplicates: {labels}")
        say(f"scan: {len(infos)} EXPNOs {sorted(kinds)}, sample_label {scan.sample_label(str(e1), {'sample': 'pCABS2-4'})!r}, "
            f"disambiguated {labels}")

    # ---- synthetic: CSV, exports round trip, Varian
    x, y, truth = _synthetic_1d(n=300, noise=0.0)
    p_csv = spectra.write_csv(out / "syn.csv", x, y, {"nucleus": "11B", "larmor_MHz": 160.0,
                                                     "sample": "distcheck"})
    ppm_c, amp_c, meta_c = spectra.read_csv(p_csv)          # comes back ascending
    order = np.argsort(x)
    _check(np.allclose(ppm_c, x[order], atol=1e-5) and np.allclose(amp_c, y[order], atol=1e-6)
           and meta_c.get("nucleus") == "11B" and meta_c.get("larmor_MHz") == 160.0,
           f"spectra.write_csv/read_csv round trip failed ({meta_c})")
    p2 = _write_csv_spectrum(out / "hdr.csv", x, y)
    for p in (Path(p_csv), p2):
        ppm_l, amp_l, rec, meta_l, _w = loader.load_any(str(p))
        _check(rec["nucleus"] == "11B" and rec["larmor_frequency_MHz"] == 160.0 and len(ppm_l) == x.size,
               f"load_any({p.name}): {rec.get('nucleus')} {rec.get('larmor_frequency_MHz')}")
    say(f"CSV: write_csv/read_csv round trip and load_any on both header styles ({x.size} points)")
    rec = truth
    rec.sample = "distcheck"
    # the exporters interpolate the experiment onto the model axis with
    # np.interp, which needs an ASCENDING axis -- every loader delivers one;
    # a descending axis is only reported here, never asserted on
    labels, cols = export.curve_table(rec, x, y)
    i_exp, i_mod = labels.index("experiment"), labels.index("model")
    desc_ok = np.allclose(cols[i_exp], cols[i_mod], atol=1e-6 * y.max())
    say("curve_table on a descending axis keeps the experiment column: "
        + ("yes" if desc_ok else "NO (np.interp needs ascending ppm; all loaders sort ascending)"))
    x, y = x[::-1].copy(), y[::-1].copy()
    export.export_text(rec, x, y, out / "fit.txt")
    export.export_csv_params(rec, out / "fit.csv")
    export.export_fxmla(rec, x, y, out / "fit.fxmla")
    export.export_json(rec, out / "fit.json")
    export.export_curves_csv(rec, x, y, out / "fit_curves.csv")
    for name in ("fit.txt", "fit.csv", "fit.fxmla", "fit.json", "fit_curves.csv"):
        _check((out / name).stat().st_size > 100, f"export {name} is empty")
    labels, cols = export.curve_table(rec, x, y)
    _check(np.allclose(cols[labels.index("experiment")], cols[labels.index("model")], atol=1e-6 * y.max()),
           "curve_table: experiment and model columns differ on a noise-free spectrum")
    back = Recipe.load(out / "fit.json")
    _check([s_.model for s_ in back.sites] == ["gauss_lor", "gauss_lor"] and back.nucleus == NUCLEUS,
           "export_json did not write a loadable recipe")
    dm = fxmla.read(out / "fit.fxmla")
    rec_dm, _notes = fxmla.to_recipe(dm)
    pos_in = [s_.params["isotropic_chemical_shift_ppm"].value for s_ in rec.sites]
    pos_out = [s_.params["isotropic_chemical_shift_ppm"].value for s_ in rec_dm.sites]
    _check(len(rec_dm.sites) == 2 and np.allclose(pos_in, pos_out, atol=0.05),
           f"fxmla round trip moved the lines: {pos_in} -> {pos_out}")
    ppm_f, amp_f, rec_f, _m, _w = loader.load_any(str(out / "fit.fxmla"))
    _check(len(ppm_f) == x.size and abs(amp_f.max() - y.max()) < 1e-3 * y.max(),
           "load_any on the exported fxmla did not return the spectrum")
    say(f"exports: txt {(out / 'fit.txt').stat().st_size} B, params csv, json, curves csv, fxmla "
        f"(re-imported: {len(rec_dm.sites)} lines at {pos_out[0]:.2f} / {pos_out[1]:.2f} ppm, "
        f"spectrum {len(ppm_f)} points)")
    vdir = out / "syn.fid"
    fid_in = _synthetic_varian(vdir)
    _check(varian.is_varian(vdir), "the synthetic Varian folder is not recognised")
    nd = varian.read(vdir)
    _check(nd.ndim == 1 and nd.data.size == fid_in.size and nd.meta["nucleus"] == "27Al"
           and abs(nd.meta["larmor_MHz"] - 104.26) < 1e-6, f"varian.read: {nd.meta}")
    ppm_v, amp_v, meta_v = varian.read_spectrum(vdir)
    peak = float(ppm_v[np.argmax(amp_v)])
    # rfl = sw/2 at rfp = 0 puts the carrier at 0 ppm: the 1500 Hz line sits 14.4 ppm from it
    _check(np.isfinite(amp_v).all() and abs(abs(peak) - 1500.0 / 104.26) < 1.0,
           f"Varian spectrum peak at {peak:.2f} ppm (expected +-{1500.0 / 104.26:.2f})")
    ppm_l, amp_l, rec_v, meta_l, warns = loader.load_any(str(vdir))
    _check(rec_v["source_kind"] == "varian" and rec_v["nucleus"] == "27Al", f"load_any(varian): {rec_v.get('nucleus')}")
    say(f"Varian: {nd.data.size}-point synthetic FID read back as {nd.meta['nucleus']} at "
        f"{nd.meta['larmor_MHz']} MHz, spectrum peak {peak:.2f} ppm, load_any ok")


# ------------------------------------------------------------ processing
def _fid(n=512, sw=10000.0, freqs=(1000.0,), t2=0.02):
    from larmor import processing

    t = np.arange(n) / sw
    y = np.zeros(n, complex)
    for f in freqs:
        y += np.exp(2j * np.pi * f * t) * np.exp(-t / t2)
    return processing.Spectrum1D(x_ppm=None, y=y, sfo1_MHz=100.0, sw_Hz=sw, domain="time")


def _dephased(p0_deg=37.0, p1_deg=-20.0, n=512):
    from larmor import processing

    x = np.linspace(-50.0, 50.0, n)
    y = 1.0 / (1.0 + ((x - 5.0) / 3.0) ** 2) + 0.0j
    y = y + 0.6 / (1.0 + ((x + 20.0) / 4.0) ** 2)
    s = processing.from_processed(x, y, 100.0)
    s.y = s.y * np.exp(1j * np.deg2rad(p0_deg + p1_deg * (np.arange(n) / (n - 1) - 0.5)))
    return s


def processing_stage(say, ctx):
    """larmor.processing: every op of the pipeline on a synthetic FID /
    spectrum, the autophase resolution, the iterative and pybaselines
    baselines, the recipe replay (loader.apply_processing), the Fourier
    module (ft1d, States recombination, ft2d), the MAS-rate resolver, the
    phase-drag arithmetic and the display geometry."""
    import dataclasses

    from larmor import baseline, display, fourier, loader, masrate, phasedrag, processing, pybaseline
    from larmor.recipe import Recipe

    out = _out(ctx, "processing")
    peak = lambda s: float(s.x_ppm[np.argmax(np.abs(s.y))])  # noqa: E731
    # time-domain ops, each followed by zf/ft: the line must stay at +10 ppm
    for ops in ([{"op": "em", "lb_hz": 20.0}], [{"op": "gm", "lb_hz": -20.0, "gb": 0.3}],
                [{"op": "sine", "ssb": 2.0}], [{"op": "sine", "ssb": 3.0, "power": 2}],
                [{"op": "traf", "lb_hz": 15.0}], [{"op": "fcor", "factor": 0.5}],
                [{"op": "shift_fid", "points": 2}], [{"op": "tdeff", "points": 256}],
                [{"op": "lp", "n_predict": 64, "n_coeff": 8}]):
        s = processing.apply(_fid(), ops + [{"op": "zf"}, {"op": "ft"}])
        _check(s.domain == "freq" and abs(peak(s) - 10.0) < 0.4, f"{ops} moved the line to {peak(s):.2f} ppm")
    s = processing.apply(_fid(), [{"op": "tdeff", "points": 256}])
    _check(s.y.size == 256, "tdeff did not truncate")
    s = processing.apply(_fid(), [{"op": "zf", "factor": 4}, {"op": "ft", "offset_ppm": 5.0}])
    _check(s.y.size == 2048 and abs(peak(s) - 15.0) < 0.4, f"zf x4 + ft offset: {s.y.size} points, peak {peak(s):.2f}")
    say("time domain: em gm sine qsine traf fcor shift_fid tdeff lp zf ft keep the line at +10 ppm")
    # frequency-domain ops
    s = processing.apply(_fid(), [{"op": "ft"}])
    x0 = peak(s)
    processing.op_sr(s, sr_hz=200.0)
    _check(abs(peak(s) - x0 - 2.0) < 1e-6, "sr 200 Hz did not shift the axis by 2 ppm")
    two = dict(freqs=(1000.0, -1500.0))      # two lines: p0 and p1 are both determined
    s = processing.apply(_fid(**two), [{"op": "zf"}, {"op": "ft"}, {"op": "phase", "p0": 90.0}])
    _check(s.y.real.max() < 0.8 * np.abs(s.y).max(), "phase p0=90 left the absorption line unchanged")
    processing.apply(s, [{"op": "autophase"}])
    _check(s.y.real.max() > 0.9 * np.abs(s.y).max(), "autophase did not recover the absorption line")
    s = processing.apply(_fid(), [{"op": "ft"}, {"op": "phase", "p0": 90.0}, {"op": "magnitude"}])
    _check((s.y.real >= 0).all(), "magnitude left negative points")
    ref = processing.apply(_fid(**two), [{"op": "fcor", "factor": 0.5}, {"op": "zf"}, {"op": "ft"}])
    s = processing.apply(_fid(**two), [{"op": "fcor", "factor": 0.5}, {"op": "zf"}, {"op": "ft"},
                                       {"op": "phase", "p0": 70.0}])
    s.y = s.y.real + 0j
    processing.apply(s, [{"op": "hilbert"}, {"op": "autophase"}])
    _check(s.y.real.max() > 0.85 * ref.y.real.max(), "hilbert + autophase on a real spectrum failed")
    s = processing.apply(_fid(), [{"op": "ft"}])
    n0 = s.y.size
    s = processing.apply(s, [{"op": "ift"}, {"op": "em", "lb_hz": 30.0}, {"op": "ft"}])
    _check(s.y.size == n0 and abs(peak(s) - 10.0) < 0.4, "ift / re-apodise / ft round trip failed")
    s = processing.apply(_fid(), [{"op": "ft"}, {"op": "extract", "hi_ppm": 20.0, "lo_ppm": 0.0}])
    _check(s.x_ppm.min() >= -0.1 and s.x_ppm.max() <= 20.1 and abs(peak(s) - 10.0) < 0.4, "extract window wrong")
    s = processing.apply(_fid(), [{"op": "ft"}, {"op": "normalize"}, {"op": "scale", "factor": 2.0},
                                  {"op": "offset", "value": 1.0}, {"op": "real"}, {"op": "conj"}])
    _check(abs(float(np.abs(s.y).max()) - 3.0) < 1e-6, f"normalize/scale/offset chain: max {np.abs(s.y).max()}")
    _check(processing.chain_start_domain([{"op": "em", "lb_hz": 10}, {"op": "ft"}]) == "time"
           and processing.chain_start_domain([{"op": "baseline", "order": 2}]) == "freq", "chain_start_domain")
    cv = processing.channel_view(np.array([3 + 4j]), "magnitude")
    _check(abs(float(cv[0]) - 5.0) < 1e-12, "channel_view magnitude")
    say("frequency domain: sr phase autophase magnitude hilbert ift/ft extract normalize scale offset real conj")
    # autophase angles and the recorded-chain resolution
    s = _dephased()
    copy = lambda s_: dataclasses.replace(s_, y=np.array(s_.y, copy=True), x_ppm=np.array(s_.x_ppm, copy=True))  # noqa: E731
    p0, p1 = processing.autophase_angles(copy(s))
    _check(abs(p0 + 37.0) < 3.0 and abs(p1 - 20.0) < 8.0, f"autophase_angles found ({p0:.1f}, {p1:.1f}), expected (-37, 20)")
    try:
        p0a, p1a = processing.autophase_angles(copy(s), method="acme")
        say(f"autophase: scan ({p0:.1f}, {p1:.1f}) deg, acme ({p0a:.1f}, {p1a:.1f}) deg")
    except ImportError as exc:
        say(f"autophase: scan ({p0:.1f}, {p1:.1f}) deg; acme unavailable: {exc}")
    ops, folds = processing.resolve_autophase(copy(s), [{"op": "phase", "p0": 10.0, "p1": 4.0, "pivot_frac": 0.3},
                                                        {"op": "autophase"}], pivot_frac=0.3)
    _check(len(folds) == 1 and ops[0]["op"] == "phase" and len(ops) == 1, f"resolve_autophase: {ops}")
    _check(processing.resolve_autophase(copy(s), [{"op": "baseline", "order": 2}], pivot_frac=0.3)
           == ([{"op": "baseline", "order": 2}], []), "resolve_autophase changed a chain without autophase")
    # baselines: polynomial, iterative, pybaselines, flat, two-point, subtract_avg
    x = np.linspace(-50, 50, 600)
    from larmor.models.analytic import gauss_lor
    signal = gauss_lor(x, 0.0, 4.0, 10.0, 1.0)
    drift = 0.5 + 0.02 * x + 0.001 * x ** 2
    s = processing.apply(processing.from_processed(x, signal + drift, 100.0), [{"op": "baseline", "order": 2}])
    _check(np.abs(np.concatenate([s.y.real[:50], s.y.real[-50:]])).max() < 0.2
           and abs(s.y.real.max() - 10.0) < 0.6, "polynomial baseline did not remove the drift")
    rng = np.random.default_rng(0)
    xi = np.arange(512)
    g = lambda c, w, a: a * np.exp(-0.5 * ((xi - c) / w) ** 2)  # noqa: E731
    peaks = g(180, 8, 100) + g(260, 5, 60) + g(330, 12, 45)
    roll = 30 * np.sin(2 * np.pi * xi / 800) + 12
    r = baseline.iterative_baseline(peaks + roll + rng.normal(0, 1.5, xi.size))
    _check(r.converged and np.sqrt(np.mean((r.baseline - roll) ** 2)) / np.ptp(roll) < 0.1,
           f"iterative_baseline: converged {r.converged}, rms error {np.sqrt(np.mean((r.baseline - roll) ** 2)) / np.ptp(roll):.3f}")
    s2 = processing.apply(processing.from_processed(np.linspace(100, -100, 512), peaks + roll, 100.0),
                          [{"op": "iterbaseline", "dead_time_pts": 0, "smoothness": 1.0, "threshold_factor": 1.0}])
    _check(abs(s2.y.real.max() - 100.0) < 15.0, "iterbaseline op left the peak heights wrong")
    say(f"baselines: polynomial, iterative ({r.n_iter} iterations, converged {r.converged})")
    t0 = time.perf_counter()
    avail = pybaseline.available()
    _check(avail, "pybaselines is not importable in this installation")
    xb = np.linspace(200.0, -100.0, 512)
    base = 8.0 + 0.05 * (xb - 50) + 0.0008 * (xb - 50) ** 2
    yb = 100 * np.exp(-0.5 * ((xb - 60) / 2.0) ** 2) + 60 * np.exp(-0.5 * ((xb - 10) / 4.0) ** 2) + base \
        + np.random.RandomState(0).normal(0, 0.3, xb.size)
    rms = lambda a: float(np.sqrt(np.mean(np.square(a))))  # noqa: E731
    methods = ("arpls", "snip", "modpoly") if ctx["quick"] else tuple(pybaseline.METHODS)
    for key in methods:
        prm = pybaseline.default_params(key, xb.size)
        est = pybaseline.compute(xb, yb, key, **prm)
        _check(est.shape == yb.shape and np.isfinite(est).all() and rms(est - base) < 0.3 * rms(yb - base),
               f"pybaselines {key}: baseline error {rms(est - base):.2f} vs signal {rms(yb - base):.2f}")
        s3 = processing.apply(processing.from_processed(xb, yb, 100.0), [{"op": "pybaseline", "method": key, **prm}])
        _check(np.allclose(s3.y.real, yb - est), f"pybaseline op {key} differs from compute()")
    say(f"pybaselines {pybaseline.METHODS[methods[0]].label} + {len(methods) - 1} more of {len(pybaseline.METHODS)} "
        f"methods, {time.perf_counter() - t0:.1f} s (includes the import)")
    s = processing.apply(processing.from_processed(x, signal + 3.0, 100.0), [{"op": "flat_baseline"}])
    _check(abs(float(s.y.real[:20].mean())) < 0.3, "flat_baseline left the offset")
    s = processing.apply(processing.from_processed(x, signal + 0.1 * x, 100.0),
                         [{"op": "twopoint_bg", "x1": -45.0, "y1": -4.5, "x2": 45.0, "y2": 4.5}])
    _check(abs(float(s.y.real[0])) < 0.5 and abs(float(s.y.real[-1])) < 0.5, "twopoint_bg did not remove the slope")
    s = processing.apply(processing.from_processed(x, signal + 2.0, 100.0), [{"op": "subtract_avg", "hi_ppm": -30.0, "lo_ppm": -50.0}])
    _check(abs(float(s.y.real[:20].mean())) < 0.3, "subtract_avg left the offset")
    # align / combine / peaks
    a = processing.from_processed(x, gauss_lor(x, 0.0, 4.0, 1.0, 1.0), 100.0)
    b = processing.from_processed(x, gauss_lor(x, 3.0, 4.0, 1.0, 1.0), 100.0)
    shift = processing.align(a, b)
    _check(abs(shift + 3.0) < 0.2, f"align found {shift:.2f} ppm, expected -3")
    c = processing.combine(a, a, op="subtract")
    _check(abs(c.y).max() < 1e-9, "combine subtract of a spectrum from itself is not zero")
    # (centres off the sample grid: two exactly equal neighbouring samples
    # are reported as two peaks by pick_peaks when min_sep_ppm is 0)
    pk = processing.pick_peaks(x, gauss_lor(x, 0.05, 4.0, 1.0, 1.0) + gauss_lor(x, 20.05, 4.0, 0.5, 1.0),
                               threshold_frac=0.1, min_sep_ppm=1.0)
    _check(len(pk) == 2 and abs(pk[0]["ppm"] - 0.05) < 0.2 and abs(pk[1]["ppm"] - 20.05) < 0.2, f"pick_peaks: {pk}")
    say(f"align {shift:.2f} ppm, combine, pick_peaks {[(round(p['ppm'], 1), round(p['fwhm_ppm'], 1)) for p in pk]}")
    # replay of a recorded chain through the loader
    xr = np.linspace(-50, 50, 400)
    rec = Recipe(nucleus=NUCLEUS, larmor_frequency_MHz=100.0,
                 processing=[{"op": "hilbert"}, {"op": "ift"}, {"op": "em", "lb_hz": 50}, {"op": "ft"},
                             {"op": "baseline", "order": 1}])
    ppm_r, amp_r, notes = loader.apply_processing(rec, xr, gauss_lor(xr, 12.3, 1.0, 1.0, 0.0) + 3.0)
    _check(np.allclose(ppm_r, xr, atol=1e-9) and abs(ppm_r[int(np.argmax(amp_r))] - 12.3) < 0.5
           and any("replayed" in n for n in notes), f"apply_processing replay: {notes}")
    say(f"loader.apply_processing: {notes[0]}")
    # fourier: ft1d, states_recombine, ft2d on a synthetic States set
    sw, n, freq = 10000.0, 512, 1500.0
    t = np.arange(n) / sw
    ppm1, spec1 = fourier.ft1d(np.exp(2j * np.pi * freq * t) * np.exp(-t / 0.02), sw, 100.0,
                               ops=[{"op": "em", "lb_hz": 20}])
    _check(abs(ppm1[np.argmax(np.abs(spec1))] - 15.0) < 0.3, "ft1d peak position")
    ser = np.arange(8 * 4).reshape(8, 4).astype(complex)
    _check(fourier.states_recombine(ser, "States").shape == (4, 4)
           and fourier.states_recombine(ser, "QF").shape == (8, 4), "states_recombine shapes")
    n1, n2, sw1, sw2 = 16, 128, 4000.0, 20000.0
    t2 = np.arange(n2) / sw2
    t1 = np.arange(n1) / sw1
    cos = np.exp(-t1[:, None] / 0.05) * np.cos(2 * np.pi * 500.0 * t1)[:, None]
    sin = np.exp(-t1[:, None] / 0.05) * np.sin(2 * np.pi * 500.0 * t1)[:, None]
    direct = (np.exp(2j * np.pi * 2000.0 * t2) * np.exp(-t2 / 0.02))[None, :]
    ser = np.empty((2 * n1, n2), complex)
    ser[0::2] = cos * direct
    ser[1::2] = sin * direct
    prm = fourier.FT2DParams(mode="States", f2_ops=[{"op": "em", "lb_hz": 20}], f1_ops=[{"op": "em", "lb_hz": 20}])
    f2, f1, z = fourier.ft2d(ser, sw2, sw1, 100.0, 100.0, params=prm)
    i1, i2 = np.unravel_index(int(np.argmax(np.abs(z))), z.shape)
    _check(abs(f2[i2] - 20.0) < 0.6 and abs(f1[i1] - 5.0) < 0.6, f"ft2d peak at F2 {f2[i2]:.1f}, F1 {f1[i1]:.1f}")
    say(f"fourier: ft1d, states_recombine, ft2d of a {ser.shape} States set -> {z.shape}, peak F2 {f2[i2]:.1f} / F1 {f1[i1]:.1f} ppm")
    # MAS rate resolution (the store redirected into the scratch folder)
    tr = masrate.parse_title_rate("Rotor SR31648\nMASR 35.714 kHz")
    _check(tr.hz == 35714.0, f"parse_title_rate: {tr}")
    ev = masrate.MasEvidence(acqus_Hz=4200.0, title_Hz=20000.0, booking_Hz=22000.0)
    res = masrate.resolve(ev)
    _check(res.uncertain and res.rate_Hz == 22000.0, f"resolve: {res}")
    _check(masrate.resolve(masrate.MasEvidence(None, None, None), masr_zero=True).source == "static", "static resolve")
    key = ("C:/data/2026-05", "RS2427418", "27Al")
    store = out / "mas_confirmations.jsonl"
    masrate.remember(key, ev, 22000.0, path=store)
    got = masrate.lookup(key, ev, path=store)
    masrate.forget(key, ev, path=store)
    _check(got and got["rate_Hz"] == 22000.0 and masrate.lookup(key, ev, path=store) is None,
           "masrate remember/lookup/forget")
    say(f"masrate: title {masrate.format_hz(tr.hz)} Hz, majority rule -> {res.source}, store round trip in {store.name}")
    # phase drag arithmetic
    _check(phasedrag.wrap_p0(190) == -170 and phasedrag.wrap_p0(540) == 180, "wrap_p0")
    rng = np.random.default_rng(7)
    yy = rng.standard_normal(256) + 1j * rng.standard_normal(256)
    xx = np.linspace(-100.0, 100.0, 256)
    s1 = processing.op_phase(processing.from_processed(xx, yy.copy(), 100.0), 37.0, 210.0, pivot_frac=0.30)
    s2 = processing.op_phase(processing.from_processed(xx, yy.copy(), 100.0),
                             phasedrag.compensate_p0(37.0, 210.0, 0.30, 0.62), 210.0, pivot_frac=0.62)
    _check(np.allclose(s1.y, s2.y, atol=1e-9), "compensate_p0 did not keep the spectrum unchanged")
    gsd = phasedrag.DragGesture()
    _check(gsd.move(40, 2) == (10.0, 0.0) and gsd.lock == "p0", "DragGesture lock")
    # display geometry
    xd = np.linspace(200.0, -100.0, 801)
    yd = 5000.0 * np.exp(-((xd - 60.0) / 15.0) ** 2) + 800.0 * np.exp(-((xd + 20.0) / 8.0) ** 2)
    area = display.area(xd, yd)
    _check(abs(area - (5000.0 * np.sqrt(np.pi) * 15.0 + 800.0 * np.sqrt(np.pi) * 8.0)) < 1e-3 * area, "display.area")
    _check(abs(display.y_factor("max", xd, yd) - 1.0 / yd.max()) < 1e-12, "display.y_factor")
    (xlo, xhi), (ylo, yhi) = display.full_extents(xd, yd)
    _check(xlo < -100.0 < 200.0 < xhi and ylo < 0 < yd.max() < yhi, "display.full_extents")
    xo, yo = display.overlay_display(xd, yd, scale=2.5, shift=1.2, yoff=0.3, stack=0.4, span=100.0)
    _check(np.allclose(xo, xd + 1.2) and np.allclose(yo, yd * 2.5 + 70.0), "display.overlay_display")
    _check(display.snap_back((250.0, 400.0), (-100.0, 200.0)) and display.view_limits((-100.0, 200.0), (0.0, 1000.0))
           == {"xMin": -400.0, "xMax": 500.0, "yMin": -2000.0, "yMax": 3000.0}, "display.snap_back / view_limits")
    say(f"phasedrag + display: area {area:.0f}, badge {display.overlay_badge(2.5, 1.2, 0.3)!r}")


# ------------------------------------------------------------ models
def models_stage(say, ctx):
    """Every lineshape model in the registry, simulated on a small ppm axis
    through engine.simulate (27Al; a 11B and a 31P case for the isotope
    path), the kernel-axis route of make_context, and a small Czjzek
    kernel built directly (the fit path's weights @ K multiply)."""
    from larmor import engine, models
    from larmor.models.base import SimContext
    from larmor.recipe import Param, Recipe, SiteModel

    with _small_kernels(ctx, say):
        _import_mrsimulator(say)
        x = np.linspace(-150.0, 150.0, 300)
        names = list(models.REGISTRY)
        expected = {"gauss_lor", "gl_norm", "jmultiplet", "sidebands", "voigt", "exchange2", "czjzek",
                    "ext_czjzek", "czjzek_d", "czjzek_corr", "amorphous", "quad_ct", "quad_first",
                    "quad_csa", "csa_czjzek", "csa_mas", "function", "spectrum"}
        _check(expected <= set(names), f"models missing from the registry: {expected - set(names)}")
        total = 0.0
        for name in names:
            m = models.get(name)
            params = m.defaults()
            params["amplitude"].value = 100.0
            extra = {}
            if name == "function":
                extra["func"] = "exp(-((x - b) / c) ** 2)"
                params["b"].value = 10.0
                params["c"].value = 6.0
            if name == "spectrum":
                extra["ref"] = {"ppm": x.tolist(), "amp": np.exp(-((x + 30.0) / 25.0) ** 2).tolist()}
            if name == "ext_czjzek":
                np.random.seed(0)          # mrsimulator's extended-Czjzek pdf is Monte-Carlo
            site = SiteModel(model=name, label=name, params=params, **extra)
            rec = Recipe(nucleus=NUCLEUS, larmor_frequency_MHz=LARMOR_MHZ, spin_rate_Hz=MAS_HZ, sites=[site])
            t0 = time.perf_counter()
            gx, y, per = engine.simulate(rec, kernel=engine.Axis(x_ppm=x))
            dt = time.perf_counter() - t0
            total += dt
            _check(gx.shape == x.shape and len(per) == 1 and np.all(np.isfinite(y)),
                   f"{name}: non-finite output or wrong shape {np.shape(y)}")
            _check(0.0 < float(np.max(y)) <= 100.0 * 1.05,
                   f"{name}: max {float(np.max(y)):.3g} for amplitude 100 (all zeros = the model did not render)")
            say(f"{name:12s} {m.label[:36]:36s} max {float(np.max(y)):7.2f}  {dt:5.2f} s"
                + ("  (mrsimulator)" if m.needs_quadrupolar or name in ("csa_mas", "csa_czjzek") else ""))
        say(f"{len(names)} models rendered in {total:.1f} s")
        # the make_context route: a Czjzek site renders on the kernel's own axis
        rec = Recipe(nucleus=NUCLEUS, larmor_frequency_MHz=LARMOR_MHZ, spin_rate_Hz=MAS_HZ, sites=[
            SiteModel(model="czjzek", label="A", params={
                "isotropic_chemical_shift_ppm": Param(62.0), "sigma_Cq_MHz": Param(1.6),
                "shift_fwhm_ppm": Param(6.0), "line_fwhm_ppm": Param(0.5, vary=False), "amplitude": Param(100.0)})])
        _check(engine.needs_kernel(rec) and not engine.needs_kernel(_two_line_recipe()), "needs_kernel")
        kx, ky, _ = engine.simulate(rec, exp_ppm=_ppm_axis(400))
        _check(kx[0] < -300 and kx[-1] > 300 and np.all(np.diff(kx) > 0), f"kernel axis {kx[0]:.0f}..{kx[-1]:.0f} ppm")
        on_x = np.interp(_ppm_axis(400), kx, ky)
        _check(abs(float(_ppm_axis(400)[np.argmax(on_x)]) - 62.0) < 15.0, "Czjzek line is not near its isotropic shift")
        info = engine.kernel_cache_info()
        say(f"make_context: Czjzek on the kernel axis {kx[0]:.0f}..{kx[-1]:.0f} ppm ({kx.size} points); "
            f"{info['entries']} kernel(s) cached, {info['mb']:.1f} MB")
        # other isotopes: 11B (spin 3/2) and a spin-1/2 CSA manifold on 31P
        ctx_b = SimContext(nucleus="11B", larmor_MHz=160.46, spin_rate_Hz=20000.0, x_ppm=np.linspace(-40.0, 60.0, 300))
        for name, values in (("quad_ct", {"isotropic_chemical_shift_ppm": 18.0, "Cq_MHz": 2.6, "eta": 0.15,
                                          "shift_fwhm_ppm": 1.0, "amplitude": 1.0}),
                             ("czjzek", {"isotropic_chemical_shift_ppm": 0.0, "sigma_Cq_MHz": 0.3,
                                         "shift_fwhm_ppm": 3.0, "line_fwhm_ppm": 0.0, "amplitude": 1.0})):
            yb = models.get(name).render(values, ctx_b)
            _check(np.isfinite(yb).all() and 0.3 < float(yb.max()) <= 1.05, f"11B {name}: max {float(yb.max()):.3g}")
        ctx_p = SimContext(nucleus="31P", larmor_MHz=162.0, spin_rate_Hz=5000.0, x_ppm=np.linspace(-200.0, 200.0, 800))
        yp = models.get("csa_mas").render({"isotropic_chemical_shift_ppm": -12.0, "zeta_ppm": 60.0, "eta": 0.4,
                                           "shift_fwhm_ppm": 1.5, "amplitude": 1.0}, ctx_p)
        n_lines = int(np.sum((yp[1:-1] > yp[:-2]) & (yp[1:-1] > yp[2:]) & (yp[1:-1] > 0.02)))
        _check(n_lines >= 3, f"31P csa_mas at 5 kHz shows {n_lines} sidebands")
        say(f"11B quad_ct / czjzek and a 31P CSA manifold ({n_lines} sidebands at 5 kHz) rendered")
        # a kernel built directly: the fit path's weights @ K multiply
        k = engine.build_kernel("27Al", 130.32, 12500.0, sw_Hz=150000.0, npts=128, ref_offset_ppm=30.0,
                                cq_max_MHz=6.0, n_cq=3, n_eta=3)
        _check(k.K.dtype == np.float32 and k.K.shape == (9, 128), f"kernel {k.K.shape} {k.K.dtype}")
        w = k.weights(2.0)
        yk = w @ k.K
        _check(w.shape == (9,) and np.isfinite(yk).all() and yk.max() > 0, "kernel weights / multiply")
        ax = engine.kernel_axis_ppm("27Al", 130.32, sw_Hz=150000.0, npts=128, ref_offset_ppm=30.0)
        _check(np.allclose(ax, k.x_ppm), "kernel_axis_ppm disagrees with the built kernel's axis")
        _check(engine.kernel_cq_max(20.0) == 25.0 and engine.kernel_cq_max(30.0) == 50.0, "kernel_cq_max ladder")
        sw, ref = engine.kernel_window(_ppm_axis(400), LARMOR_MHZ)
        _check(sw >= engine.KERNEL_MIN_SW_HZ and -100 < ref < 150, f"kernel_window {sw} {ref}")
        say(f"build_kernel: {k.K.shape} float32 basis, weights(2.0) sum {w.sum():.3f}, axis {ax[0]:.0f}..{ax[-1]:.0f} ppm, "
            f"margin {engine.site_width_margin(rec.sites):.0f} ppm")


# ------------------------------------------------------------ physics helpers
def physics_stage(say, ctx):
    """The physics helpers around the models: sideband detection on a
    synthetic MAS manifold, Herzfeld-Berger intensities and the manifold
    measurement (+ its inversion outside quick mode), the static CT feature
    map and its inversion, the Czjzek distribution, unit conversions, the
    isotope table, literature shift ranges, the command-palette matcher,
    start-value estimation and the measurement tools."""
    from larmor import (convert, czjzek_dist, engine, estimate, herzfeld_berger as hb, measure, nuclei,
                        palette, refranges, sidebands, staticct)
    from larmor.recipe import Param, Recipe, SiteModel

    with _small_kernels(ctx, say):
        X = np.linspace(-50.0, 50.0, 2000)                  # 1H at 100 MHz: 1 kHz = 10 ppm
        r = Recipe(nucleus="1H", larmor_frequency_MHz=100.0, spin_rate_Hz=1000.0, sites=[
            SiteModel(model="sidebands", label="S", params={
                "isotropic_chemical_shift_ppm": Param(0.0), "shift_fwhm_ppm": Param(1.0),
                "amplitude": Param(1.0), "ssb_ratio": Param(0.4), "n_ssb": Param(3.0), "gl": Param(1.0)})])
        xc, y, _ = engine.simulate(r, exp_ppm=X)
        det = sidebands.detect(xc, y, 100.0, nu_rot_Hz=1000.0)
        _check(det.ok and abs(det.nu_rot_Hz - 1000.0) < 20 and det.n_orders >= 3 and det.confidence > 0.8,
               f"sideband detection: {det.message} (ok {det.ok}, confidence {det.confidence:.2f})")
        det2 = sidebands.detect(xc, y, 100.0, nu_rot_Hz=0.0, scan=True)
        _check(det2.ok and det2.scanned and abs(det2.nu_rot_Hz - 1000.0) < 50,
               f"sideband scan found {det2.nu_rot_Hz} Hz")
        ratio = sidebands.seed_ratio(det)
        _check(abs(ratio - 0.4) < 0.1, f"seed_ratio {ratio:.2f}")
        say(f"sidebands: {sidebands.describe(det)}; scan {det2.nu_rot_Hz:.0f} Hz; ratio {ratio:.2f}")
        NU0, NUR = 162.0, 5000.0
        inten = hb.sideband_intensities(60.0, 0.4, NU0, NUR, n_max=6)
        _check(abs(sum(inten.values()) - 1.0) < 1e-3 and inten[0] > inten[3] > 0, "Herzfeld-Berger intensities")
        pos = hb.sideband_positions_ppm(-12.0, NU0, NUR, inten)
        _check(abs(pos[1] - pos[0] - NUR / NU0) < 1e-9, "sideband_positions_ppm spacing")
        xh = np.linspace(-200.0, 200.0, 3001)
        sig = 1.5 / (2 * np.sqrt(2 * np.log(2)))
        yh = np.full_like(xh, 0.02)
        for n, i_n in inten.items():
            yh += i_n / (sig * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((xh - pos[n]) / sig) ** 2)
        raw = hb.measure_manifold(xh, yh, -12.0, NU0, NUR, 4)
        _check({-1, 0, 1} <= set(raw) and raw[0] > raw[2] > 0, f"measure_manifold {raw}")
        line = (f"Herzfeld-Berger: {len(inten)} orders carry intensity ({hb.n_sidebands_needed(60.0, NU0, NUR)} needed), "
                f"manifold measured over {len(raw)} orders")
        if not ctx["quick"]:
            fitted = hb.fit_sidebands(raw, NU0, NUR)
            _check(fitted.ok and abs(fitted.zeta_ppm - 60.0) < 3.0 and abs(fitted.eta - 0.4) < 0.1,
                   f"fit_sidebands: zeta {fitted.zeta_ppm:.1f} eta {fitted.eta:.2f}")
            line += f", inverted to zeta {fitted.zeta_ppm:.1f} ppm, eta {fitted.eta:.2f}"
        else:
            line += " (the 10-20 s inversion runs outside quick mode)"
        say(line)
        feats = staticct.ct_static_features(4.0, 0.3, 2.5, LARMOR_MHZ)
        dists = [min(abs(e - h) for h in feats.horns) for e in feats.edges]
        edge = feats.edges[int(np.argmax(dists))]
        read = staticct.read_cq_eta(feats.horns[0] + 10.0, feats.horns[1] + 10.0, edge + 10.0, 2.5, LARMOR_MHZ)
        _check(read.ok and abs(read.Cq_MHz - 4.0) < 0.1 and abs(read.eta - 0.3) < 0.05
               and abs(read.delta_iso_ppm - 10.0) < 1.0, f"staticct read {read}")
        say(f"staticct: horns {feats.horns[0]:.1f}/{feats.horns[1]:.1f} ppm -> Cq {read.Cq_MHz:.2f} MHz, eta {read.eta:.2f}, "
            f"d_iso {read.delta_iso_ppm:.1f} ppm")
        cq = czjzek_dist.suggested_cq_axis(1.6, 400)
        marg = czjzek_dist.marginal_cq(1.6, cq)
        _check(abs(marg.sum() - 1.0) < 1e-6 and abs(cq[np.argmax(marg)] - 3.73 * 1.6) < 0.1, "czjzek marginal")
        pdf = czjzek_dist.czjzek_pdf(1.6, cq, np.linspace(0, 1, 11))
        wts = czjzek_dist.czjzek_weights(1.6, 5.0, np.linspace(0.05, 25, 12), np.linspace(0, 1, 3))
        _check(pdf.shape == (11, 400) and wts.shape == (36,) and np.isfinite(wts).all(), "czjzek_pdf / weights")
        say(f"czjzek_dist: mode |Cq| {cq[np.argmax(marg)]:.2f} MHz for sigma 1.6, rms Pq {czjzek_dist.rms_pq(1.6):.2f}, "
            f"mode Pq {czjzek_dist.mode_pq(1.6):.2f}")
        pq = convert.pq_from_cq_eta(4.0, 0.3)
        _check(abs(convert.cq_from_pq_eta(pq, 0.3) - 4.0) < 1e-9 and convert.nu_q(4.0, 2.5) > 0, "convert pq/nu_q")
        _check(abs(convert.Hz_to_ppm(convert.ppm_to_Hz(12.5, 130.3), 130.3) - 12.5) < 1e-9, "convert ppm<->Hz")
        _check(convert.ct_second_order_shift_ppm(5.0, 2.5, 156.0) < 0, "second-order shift sign")
        d11, d22, d33 = convert.csa_principal_from_haeberlen(10.0, 60.0, 0.4)
        back = convert.csa_haeberlen_from_principal(d11, d22, d33)
        _check(np.allclose(back, (10.0, 60.0, 0.4), atol=1e-6), f"CSA Haeberlen round trip {back}")
        span, skew = convert.csa_span_skew(d11, d22, d33)
        _check(np.allclose(convert.csa_principal_from_span_skew(10.0, span, skew), sorted((d11, d22, d33), reverse=True),
                           atol=1e-6), "CSA span/skew round trip")
        d_hz = convert.dipolar_Hz(10.7084, -4.3173, 1.5)
        _check(abs(convert.distance_from_dipolar(10.7084, -4.3173, d_hz) - 1.5) < 1e-9, "dipolar round trip")
        say(f"convert: Pq {pq:.3f} MHz, 13C-15N at 1.5 A = {d_hz:.0f} Hz, CSA span {span:.1f} skew {skew:.2f}")
        al = next(i for i in nuclei.all_isotopes() if i.symbol == "27Al")
        _check(al.spin == 2.5 and abs(al.larmor_MHz(9.4) - 104.37) < 0.05 and abs(al.larmor_from_1H(500.0) - 130.39) < 0.05,
               f"nuclei 27Al: spin {al.spin}, {al.larmor_MHz(9.4):.2f} MHz at 9.4 T")
        b11 = nuclei.primary_isotope("B")
        _check(b11 is not None and b11.symbol == "11B" and nuclei.feasibility("Al") == "favorable", "nuclei table")
        say(f"nuclei: {len(nuclei.all_isotopes())} isotopes; 27Al {al.larmor_from_1H(500.13):.3f} MHz on a 500 MHz magnet "
            f"(B0 {nuclei.b0_from_1H(500.13):.2f} T)")
        hit = refranges.assign("27Al", 65.0)
        _check(hit is not None and hit["label"] == "Al[4]" and refranges.family_for("11B", 15.0) == "BO3", "refranges")
        _check(refranges.ranges_for("29Si") and refranges.positions_for("19F") and refranges.citation_for(hit),
               "refranges tables")
        say(f"refranges: 65 ppm 27Al -> {hit['label']} ({refranges.citation_for(hit)[:40]}...)")
        _check(palette.fuzzy_score("fit", "Fit") and palette.rank("fit", ["Profit", "Fit", "Refit"]) == [1, 2, 0]
               and palette.fold("χ² map") == "chi2 map" and palette.strip_mnemonic("&File") == "File", "palette matcher")
        xs = np.linspace(-150.0, 150.0, 400)
        q = Recipe(nucleus=NUCLEUS, larmor_frequency_MHz=LARMOR_MHZ, spin_rate_Hz=0.0, sites=[
            SiteModel(model="quad_ct", label="Q", params={
                "isotropic_chemical_shift_ppm": Param(20.0), "Cq_MHz": Param(6.0), "eta": Param(0.3),
                "shift_fwhm_ppm": Param(2.0), "amplitude": Param(1.0)})])
        _, ys, _ = engine.simulate(q, exp_ppm=xs)
        got = estimate.start_values("quad_ct", xs, ys, NUCLEUS, LARMOR_MHZ)
        _check(2.0 < got.get("Cq_MHz", 0.0) < 12.0, f"estimate.start_values quad_ct: {got}")
        centre, width = estimate.band_width_ppm(xs, ys)
        _check(width > 5.0 and -60 < centre < 60, f"band_width_ppm {centre:.1f} {width:.1f}")
        say(f"estimate: static 27Al Cq 6 MHz read back as {got['Cq_MHz']:.2f} MHz, band {width:.0f} ppm wide at {centre:.0f} ppm; "
            f"cq_for_width(20 ppm) {estimate.cq_for_width(20.0, NUCLEUS, LARMOR_MHZ):.2f} MHz")
        xm = np.linspace(-50, 50, 2000)
        gm = lambda c, w, a: a * np.exp(-4 * np.log(2) * ((xm - c) / w) ** 2)  # noqa: E731
        ym = gm(0, 6, 1.0) + gm(20, 6, 0.5)
        _check(abs(measure.fwhm(xm, ym, (10, -10)) - 6.0) < 0.05, "measure.fwhm")
        rows = measure.integrate_regions(xm, ym, [(10, -10), (30, 10)])
        _check(abs(rows[0]["percent"] - 66.7) < 0.5 and abs(rows[1]["percent"] - 33.3) < 0.5, f"integrate_regions {rows}")
        _check(abs(measure.centre_of_mass(xm, ym, (30, 10)) - 20.0) < 0.1
               and measure.snr(xm, ym + np.random.default_rng(0).normal(0, 1e-3, xm.size), (10, -10)) > 100, "measure")
        say(f"measure: fwhm {measure.fwhm(xm, ym, (10, -10)):.2f} ppm, regions {rows[0]['percent']:.1f} / {rows[1]['percent']:.1f} %")


# ------------------------------------------------------------ fitting
class _SelftestSay:
    """selftest.core_fits logs through ``say(...)`` and ``say.fail(stage, exc)``."""

    def __init__(self, note):
        self.note = note
        self.fails: list = []

    def __call__(self, *parts):
        self.note("selftest:", *parts)

    def fail(self, stage, exc=None):
        self.fails.append(f"{stage}: {exc!r}" if exc is not None else stage)
        self.note("selftest FAIL:", stage, repr(exc) if exc is not None else "")


def _czjzek_recipe(pos=62.0, sigma=1.6, dcs=6.0, amp=100.0, bounds=False):
    from larmor.recipe import Param, Recipe, SiteModel

    if bounds:
        params = {"isotropic_chemical_shift_ppm": Param(pos, min=30.0, max=90.0),
                  "sigma_Cq_MHz": Param(sigma, min=0.3, max=2.4), "shift_fwhm_ppm": Param(dcs, min=1.0, max=30.0),
                  "line_fwhm_ppm": Param(0.5, vary=False), "amplitude": Param(amp, min=0.0)}
    else:
        params = {"isotropic_chemical_shift_ppm": Param(pos), "sigma_Cq_MHz": Param(sigma),
                  "shift_fwhm_ppm": Param(dcs), "line_fwhm_ppm": Param(0.5, vary=False), "amplitude": Param(amp)}
    return Recipe(nucleus=NUCLEUS, larmor_frequency_MHz=LARMOR_MHZ, spin_rate_Hz=MAS_HZ, sample="distcheck",
                  sites=[SiteModel(model="czjzek", label="Al4", params=params)])


def fitting_stage(say, ctx):
    """larmor.fit on a synthetic two-line Gauss/Lorentz spectrum (free, then
    with a linked width), on a Czjzek line (the kernel path, reduced in
    quick mode), and selftest.core_fits -- the first-time user's two fits --
    reused as is."""
    from larmor import engine, fit as fitmod, selftest

    with _small_kernels(ctx, say):
        _import_mrsimulator(say)
        x, y, truth = _synthetic_1d(n=400, noise=0.5, seed=1)
        start = _start_recipe()
        calls: list = []
        res = fitmod.fit(start, x, y, window_ppm=(120.0, -60.0), iter_cb=lambda *a, **k: calls.append(1))
        sites = res.recipe.sites
        pos = [s.params["isotropic_chemical_shift_ppm"].value for s in sites]
        amps = [s.params["amplitude"].value for s in sites]
        errs = [s.params["amplitude"].stderr for s in sites]
        _check(abs(pos[0] - 60.0) < 0.5 and abs(pos[1] - 5.0) < 0.5 and abs(amps[0] - 100.0) < 3.0
               and abs(amps[1] - 40.0) < 3.0, f"gauss_lor fit landed at {pos} ppm, amplitudes {amps}")
        _check(res.rmsd < 0.02 and res.lmfit_result.nfev >= 8 and not getattr(res.lmfit_result, "aborted", False),
               f"gauss_lor fit: rmsd {res.rmsd:.4f}, nfev {res.lmfit_result.nfev}, message {res.lmfit_result.message}")
        _check(all(e is not None and np.isfinite(e) and e > 0 for e in errs), f"no error bars: {errs}")
        _check(res.recipe.fit_rmsd == res.rmsd and res.recipe.software.get("larmor") and res.recipe.fit_window_ppm,
               "the fitted recipe lacks its rmsd / software stamp / window")
        _check(res.x_ppm.size == x.size and res.y_fit.shape == res.x_ppm.shape and len(res.per_site) == 2
               and res.lmfit_result.residual.size <= x.size, "FitResult arrays")
        say(f"gauss_lor x2: {len(calls)} evaluations, rmsd {res.rmsd:.4f}; Al4 {pos[0]:.2f} ppm, {amps[0]:.1f} ± {errs[0]:.1f}; "
            f"Al6 {pos[1]:.2f} ppm, {amps[1]:.1f} ± {errs[1]:.1f}; software {res.recipe.software.get('larmor')}")
        linked = _start_recipe()
        linked.sites[1].params["shift_fwhm_ppm"].expr = "s0.shift_fwhm_ppm"
        linked.sites[1].params["amplitude"].expr = "0.4 * s0.amplitude"
        res2 = fitmod.fit(linked, x, y, window_ppm=(120.0, -60.0))
        w = [s.params["shift_fwhm_ppm"].value for s in res2.recipe.sites]
        a2 = [s.params["amplitude"].value for s in res2.recipe.sites]
        _check(abs(w[0] - w[1]) < 1e-9 and abs(a2[1] - 0.4 * a2[0]) < 1e-9 and res2.rmsd < 0.1,
               f"linked fit: widths {w}, amplitudes {a2}, rmsd {res2.rmsd:.3f}")
        _check(fitmod.translate_expr("0.4 * s0.amplitude", linked) == "0.4 * s0_amp", "translate_expr")
        say(f"linked: B width = A width ({w[1]:.2f} ppm), B amplitude = 0.4 A ({a2[1]:.1f}), rmsd {res2.rmsd:.4f}, "
            f"stderr of the derived amplitude {res2.recipe.sites[1].params['amplitude'].stderr}")
        kx, ky, _ = engine.simulate(_czjzek_recipe(), exp_ppm=x)
        ycz = np.interp(x, kx, np.asarray(ky, float)) + np.random.default_rng(4).normal(0.0, 0.4, x.size)
        t0 = time.perf_counter()
        res3 = fitmod.fit(_czjzek_recipe(pos=58.0, sigma=1.0, dcs=8.0, amp=60.0, bounds=True), x, ycz,
                          window_ppm=(120.0, -60.0))
        dt = time.perf_counter() - t0
        s0 = res3.recipe.sites[0].params
        _check(abs(s0["isotropic_chemical_shift_ppm"].value - 62.0) < 3.0 and abs(s0["sigma_Cq_MHz"].value - 1.6) < 0.4
               and res3.rmsd < 0.06, f"czjzek fit: pos {s0['isotropic_chemical_shift_ppm'].value:.1f}, "
               f"sigma {s0['sigma_Cq_MHz'].value:.2f}, rmsd {res3.rmsd:.4f}")
        say(f"czjzek: pos {s0['isotropic_chemical_shift_ppm'].value:.1f} ppm, sigma {s0['sigma_Cq_MHz'].value:.2f} MHz, "
            f"rmsd {res3.rmsd:.4f}, {res3.lmfit_result.nfev} evaluations in {dt:.1f} s; fit axis {res3.x_ppm.size} points "
            f"(the kernel's)")
        wrapped = _SelftestSay(say)
        t0 = time.perf_counter()
        ok = selftest.core_fits(wrapped)
        _check(ok and not wrapped.fails, f"selftest.core_fits failed: {wrapped.fails}")
        say(f"selftest.core_fits: both first-time-user fits pass in {time.perf_counter() - t0:.1f} s"
            + (" (reduced kernel)" if ctx["quick"] else ""))


# ------------------------------------------------------------ error tools
def error_tools_stage(say, ctx):
    """The error estimators and verdicts on a fitted synthetic spectrum: the
    multi-start (auto_fit), the chi-square profile, Monte-Carlo refits
    (sequential), the chi2 map, identifiability (correlations), the physical
    sanity checks, the residual diagnostics and the fit-health verdict."""
    from larmor import (autofit, chi2map, diagnostics, fit as fitmod, fithealth, identifiability, quantify,
                        sanity)
    from larmor.recipe import Param, Recipe

    quick = ctx["quick"]
    window = (120.0, -60.0)
    x, y, truth = _synthetic_1d(n=400, noise=0.5, seed=2)
    rec = _start_recipe()
    res = fitmod.fit(rec, x, y, window_ppm=window)
    t0 = time.perf_counter()
    auto = autofit.auto_fit(_start_recipe(), x, y, window_ppm=window, n_starts=2 if quick else 6, seed=1)
    _check(auto.best_rmsd <= res.rmsd + 1e-6 and len(auto.trials) >= 2 and abs(auto.trials[0] - auto.best_rmsd) < 1e-12
           and auto.result is not None, f"auto_fit: best {auto.best_rmsd:.4f} vs plain {res.rmsd:.4f}, trials {auto.trials}")
    say(f"auto_fit: {len(auto.trials)} trials in {time.perf_counter() - t0:.1f} s, best rmsd {auto.best_rmsd:.4f} "
        f"({auto.n_improved} improved on the start)")
    t0 = time.perf_counter()
    prof = autofit.error_profile(rec, x, y, site=0, param="isotropic_chemical_shift_ppm", window_ppm=window,
                                 n_points=5 if quick else 11, span=3.0, parallel=False)
    _check(prof.chi2.min() == prof.chi2_min and prof.ci68[0] is not None and prof.ci68[1] is not None
           and prof.ci68[0] < prof.best_value < prof.ci68[1] and not prof.ran_parallel
           and 100 < prof.dof < x.size, f"error_profile: ci68 {prof.ci68}, best {prof.best_value}, dof {prof.dof}")
    say(f"error_profile: {prof.values.size} points in {time.perf_counter() - t0:.1f} s, Al4 position "
        f"{prof.best_value:.2f} ppm, 68 % [{prof.ci68[0]:.2f}, {prof.ci68[1]:.2f}], noise {np.sqrt(prof.noise_var):.2f}")
    t0 = time.perf_counter()
    mc = autofit.monte_carlo_errors(Recipe.from_dict(rec.to_dict()), x, y, window_ppm=window, n_trials=4, seed=3,
                                    parallel=False)
    _check(mc.n_ok == 4 and mc.params and all(np.isfinite(p.std) for p in mc.params) and not mc.ran_parallel
           and abs(mc.noise - 0.5) < 0.25 and mc.site_integrals is not None and mc.site_integrals.shape == (4, 2),
           f"monte_carlo_errors: n_ok {mc.n_ok}, noise {mc.noise}, integrals {None if mc.site_integrals is None else mc.site_integrals.shape}")
    report = mc.report()
    _check(isinstance(report, str) and len(report) > 50, "MonteCarloResult.report")
    say(f"monte_carlo_errors: {mc.trials} trials in {time.perf_counter() - t0:.1f} s, noise {mc.noise:.2f}, "
        f"{len(mc.params)} parameters, Al4 position std {next(p.std for p in mc.params if p.param == 'isotropic_chemical_shift_ppm'):.3f} ppm")
    A, B, Z, (a0, b0) = chi2map.chi2_surface(rec.to_dict(), x, y, window, axis_a=(0, "isotropic_chemical_shift_ppm"),
                                             axis_b=(0, "shift_fwhm_ppm"), n=5 if quick else 11)
    ib, ia = np.unravel_index(int(np.argmin(Z)), Z.shape)
    _check(Z.shape == (A.size, B.size) and abs(A[ia] - a0) <= (A[1] - A[0]) + 1e-9
           and abs(B[ib] - b0) <= (B[1] - B[0]) + 1e-9, "chi2_surface minimum is not at the fitted values")
    varying = chi2map.varying_params(rec.to_dict())
    _check(any(v[0] == 0 and v[1] == "isotropic_chemical_shift_ppm" for v in varying)
           and not any(v[1] == "gl" for v in varying), f"varying_params {varying}")
    say(f"chi2_surface: {Z.shape} grid, minimum at ({A[ia]:.2f} ppm, {B[ib]:.2f} ppm) = fitted ({a0:.2f}, {b0:.2f}); "
        f"{len(varying)} varying parameters")
    names, corr = identifiability.corr_matrix(res.lmfit_result)
    pairs = identifiability.unidentifiable_pairs(res.lmfit_result)
    _check(corr is not None and corr.shape == (len(names), len(names)) == (6, 6) and pairs == [],
           f"identifiability: {len(names)} names, pairs {pairs}")
    warns = sanity.check_recipe(res.recipe, window=window)
    bad = Recipe.from_dict(res.recipe.to_dict())
    bad.sites[0].params["gl"] = Param(1.4, vary=False)
    bad.sites[1].params["shift_fwhm_ppm"] = Param(-2.0)
    bad_warns = sanity.check_recipe(bad, window=window)
    _check(warns == [] and len(bad_warns) >= 2 and sanity.summarize(bad_warns), f"sanity: {warns} / {bad_warns}")
    say(f"identifiability: {len(names)} parameters, strongest |r| {np.abs(corr - np.eye(len(names))).max():.2f}, "
        f"no degenerate pair; sanity: clean fit {warns}, broken recipe -> {sanity.summarize(bad_warns)!r}")
    resid = np.asarray(res.lmfit_result.residual, float)
    white = diagnostics.residual_structure(resid)
    wavy = diagnostics.residual_structure(np.sin(np.linspace(0, 20, 400)))
    _check(not white["structured"] and wavy["structured"] and abs(diagnostics.lag1_autocorr(resid)) < 0.3
           and diagnostics.runs_test(resid)["n"] == resid.size, f"diagnostics: white {white}, wavy {wavy}")
    quant_rows = quantify.quantify(res.recipe)["rows"]
    h = fithealth.assess(res.recipe, res.y_exp, res.y_fit, lmfit_result=res.lmfit_result, at_bounds=res.at_bounds or [],
                         frozen=res.frozen_sites or [], window=res.recipe.fit_window_ppm, rmsd=res.rmsd,
                         quant_rows=quant_rows, fitted=True, ppm=x, x_fit=res.x_ppm)
    _check(h.fitted and h.level in ("ok", "check") and h.pill_text() and h.chi_text().startswith("RMSD")
           and h.redchi is not None, f"fithealth on a clean fit: {h.level} {h.flags}")
    hb = fithealth.assess(bad, res.y_exp, res.y_fit, lmfit_result=res.lmfit_result, window=res.recipe.fit_window_ppm,
                          rmsd=res.rmsd, ppm=x, x_fit=res.x_ppm)
    _check(hb.level == "bad" and any(f.kind == "physical" for f in hb.flags), f"fithealth on an impossible value: {hb.level}")
    say(f"diagnostics: residual runs z {white['runs_z']:.2f}, lag1 {white['lag1']:.2f} -> {white['message']}; "
        f"fithealth: {h.pill_text()} · {h.chi_text()} · {h.summary()[:60]}; impossible gl -> {hb.pill_text()}")


# ------------------------------------------------------------ constraints
ISO = "isotropic_chemical_shift_ppm"


def _dict_sites(exprs) -> list:
    out = []
    for i, e in enumerate(exprs):
        out.append({"model": "gauss_lor", "label": "ABCDE"[i], "params": {
            ISO: {"value": 10.0 * i + 5.0, "stderr": None, "vary": True, "min": None, "max": None, "expr": e},
            "shift_fwhm_ppm": {"value": 6.0, "stderr": None, "vary": True, "min": 0.1, "max": None, "expr": None},
            "amplitude": {"value": 100.0, "stderr": None, "vary": True, "min": 0.0, "max": None, "expr": None},
            "gl": {"value": 1.0, "stderr": None, "vary": False, "min": 0.0, "max": 1.0, "expr": None}}})
    return out


def constraints_stage(say, ctx):
    """The constraint helpers: sanitising / remapping links on recipe-dict
    sites, the glass-protocol bounds, the constraint library (capture on
    one recipe, apply on another), the lines table's cell grammar and the
    parameter-status read-out (fixed / linked / at a bound)."""
    from larmor import cellparse, constraint_library, constraints_util, paramstatus
    from larmor.recipe import Param

    iso = lambda s: s["params"][ISO]["expr"]  # noqa: E731
    s = _dict_sites([None, f"s1.{ISO} + 5.3"])
    _check(constraints_util.sanitize_constraints(s) == [f"s1.{ISO}"] and iso(s[1]) is None, "self-reference not dropped")
    s = _dict_sites([f"s9.{ISO}", None])
    _check(constraints_util.sanitize_constraints(s) == [f"s0.{ISO}"], "out-of-range reference not dropped")
    s = _dict_sites([f"s1.{ISO}", f"s0.{ISO}"])
    _check(len(constraints_util.sanitize_constraints(s)) == 1 and sum(iso(t) is not None for t in s) == 1, "cycle not broken")
    s = _dict_sites([None, f"s0.{ISO} + 5.3"])
    _check(constraints_util.sanitize_constraints(s) == [], "a valid link was dropped")
    s = _dict_sites([None, None, f"s1.{ISO} + 2"])
    s.pop(0)
    dropped = constraints_util.remap_exprs_after_delete(s, 0)
    _check(iso(s[1]) == f"s0.{ISO} + 2" and dropped == [], f"remap after delete: {iso(s[1])} {dropped}")
    s = _dict_sites([f"s2.{ISO}", None, None])
    s = [s[2], s[1], s[0]]
    constraints_util.remap_exprs_after_move(s, {0: 2, 1: 1, 2: 0})
    _check(iso(s[2]) == f"s0.{ISO}" and iso(s[0]) is None, "remap after move")
    _check(constraints_util.param_refs("0.5 * s0.amplitude + s1.shift_fwhm_ppm") == [(0, "amplitude"), (1, "shift_fwhm_ppm")]
           and constraints_util.site_refs("0.5 * s0.amplitude + s1.shift_fwhm_ppm") == {0, 1}
           and constraints_util.references_self(f"s1.{ISO}", 1, ISO), "param_refs / site_refs / references_self")
    s = _dict_sites([None, None])
    notes = constraints_util.restrict_glass_protocol(s)
    p0 = s[0]["params"]
    _check(p0[ISO]["min"] == 2.0 and p0[ISO]["max"] == 8.0 and p0["shift_fwhm_ppm"]["min"] == 4.0 and notes,
           f"restrict_glass_protocol: {p0[ISO]} {p0['shift_fwhm_ppm']} {notes}")
    say(f"constraints_util: sanitize (self / out-of-range / cycle), remap after delete and move, "
        f"glass protocol -> {notes[0] if notes else ''}")
    rec = _start_recipe()
    rec.sites[1].params["shift_fwhm_ppm"].expr = "s0.shift_fwhm_ppm"
    cset = constraint_library.capture(rec)
    fresh = _two_line_recipe()
    applied = constraint_library.apply(fresh, cset)
    _check(fresh.sites[1].params["shift_fwhm_ppm"].expr == "s0.shift_fwhm_ppm" and fresh.sites[0].params[ISO].min == 30.0
           and fresh.sites[0].params["gl"].vary is False and "s1.shift_fwhm_ppm" in applied
           and "link" in constraint_library.describe(cset), f"constraint library: applied {applied}")
    say(f"constraint_library: {constraint_library.describe(cset)} -> {len(applied)} constraints re-applied")
    base = {"param_name": ISO, "param_unit": "ppm", "this_index": 2, "n_sites": 5, "larmor_MHz": 130.323}
    P = lambda text, **kw: cellparse.parse_cell(text, **{**base, **kw})  # noqa: E731
    r = P("62.6266")
    _check(r.set_value and abs(r.value - 62.6266) < 1e-9 and r.error is None, f"cell '62.6266': {r}")
    _check(P("A").expr == f"s0.{ISO}" and P("A+20").expr == f"s0.{ISO} + 20", "cell links")
    _check(P("0.5B", param_name="amplitude", param_unit="").expr == "0.5*s1.amplitude", "cell '0.5B'")
    _check(P("A + 500 kHz", param_name="Cq_MHz", param_unit="MHz").expr == "s0.Cq_MHz + 0.5", "cell Hz conversion")
    r = P("62.6 [0, 100]")
    _check(r.set_value and abs(r.value - 62.6) < 1e-9 and r.set_min and r.min == 0.0 and r.set_max and r.max == 100.0,
           f"cell with bounds: {r}")
    r = P("C")
    _check(r.error, "a link to the line itself was accepted")
    _check(cellparse.format_link(f"s0.{ISO} + 20", ISO) == "A+20" and cellparse.format_link("0.5*s1.amplitude", "amplitude") == "0.5B"
           and cellparse.index_to_letter(27) == "AB" and cellparse.letter_to_index("AB") == 27, "format_link / letters")
    say(f"cellparse: values, links (A, A+20, 0.5B, A + 500 kHz -> MHz), bounds, self-link refused ({r.error!r})")
    rec.sites[0].params[ISO] = Param(30.0, min=30.0, max=90.0)
    st = paramstatus.recipe_statuses(rec)
    _check(st[(0, "gl")].kind == "fixed" and st[(1, "shift_fwhm_ppm")].kind == "linked"
           and st[(0, ISO)].kind == "at_bound" and st[(0, ISO)].side == "min" and st[(1, ISO)].kind == "free",
           f"paramstatus: {st}")
    summ = paramstatus.summary(rec)
    foot = paramstatus.footnote([("A d_iso", st[(0, ISO)]), ("B FWHM", st[(1, "shift_fwhm_ppm")]), ("A gl", st[(0, "gl")])], "text")
    _check("fixed" in summ and "linked" in summ and "bound" in summ and "‡" in foot and "§" in foot and "†" in foot,
           f"paramstatus summary {summ!r} / footnote {foot!r}")
    cells = paramstatus.csv_fields(rec.sites[1].params["shift_fwhm_ppm"], st[(1, "shift_fwhm_ppm")])
    _check(len(cells) == 5 and cells[3] == "s0.shift_fwhm_ppm", f"csv_fields {cells}")
    _check(paramstatus.effective_bounds("gauss_lor", "shift_fwhm_ppm", Param(5.0)) == (0.1, None), "effective_bounds")
    say(f"paramstatus: {summ}; footnote {foot[:70]}...")


# ------------------------------------------------------------ several spectra
def _b_recipe(sample: str, larmor: float = 160.0):
    from larmor.recipe import Param, Recipe, SiteModel

    return Recipe(nucleus="11B", larmor_frequency_MHz=larmor, spin_rate_Hz=0.0, sample=sample, sites=[
        SiteModel(model="gauss_lor", label="A", params={
            ISO: Param(14.0, min=0.0, max=30.0), "shift_fwhm_ppm": Param(5.0, min=0.1),
            "amplitude": Param(80.0, min=0.0), "gl": Param(1.0, vary=False)}),
        SiteModel(model="gauss_lor", label="B", params={
            ISO: Param(1.0, min=-10.0, max=15.0), "shift_fwhm_ppm": Param(3.5, min=0.1),
            "amplitude": Param(50.0, min=0.0), "gl": Param(1.0, vary=False)})])


def _b_spectrum(x, pos_shift: float, amps, seed: int):
    from larmor import engine
    from larmor.recipe import Param, Recipe, SiteModel

    tr = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, spin_rate_Hz=0.0, sites=[
        SiteModel(model="gauss_lor", label="A", params={
            ISO: Param(15.0 + pos_shift), "shift_fwhm_ppm": Param(6.0), "amplitude": Param(amps[0]), "gl": Param(1.0)}),
        SiteModel(model="gauss_lor", label="B", params={
            ISO: Param(1.0), "shift_fwhm_ppm": Param(3.5), "amplitude": Param(amps[1]), "gl": Param(1.0)})])
    xs, m, _ = engine.simulate(tr, exp_ppm=x)
    return np.interp(x, xs, m) + np.random.default_rng(seed).normal(0.0, 1.0, x.size)


def multi_stage(say, ctx):
    """Several spectra at once: the shared-model fit over two fields
    (multifit), the amplitude-only batch fit with released positions, its
    covariance error analysis, the shared / error CSVs and the publication
    bundle; the sequential series fit through the series-mode plumbing
    (entries_for_sweep, apply_sweep_result, carry_into, components); and
    the publication report (larmor.batch) on saved fits."""
    import copy

    from larmor import batch, batchfit, components, multifit, seqfit, seriesmode
    from larmor.io import bundle

    out = _out(ctx, "multi")
    quick = ctx["quick"]
    x = np.linspace(-20.0, 60.0, 400)
    window = (40.0, -10.0)
    shifts, amps = [0.0, 0.3, -0.3], [[100, 60], [80, 90], [120, 40]]
    entries = [(_b_recipe(f"g{k}"), x, _b_spectrum(x, sh, am, k), window)
               for k, (sh, am) in enumerate(zip(shifts, amps))]
    # the two-field co-fit (positions and widths shared, amplitudes free)
    t0 = time.perf_counter()
    mres = multifit.fit_multi([(_b_recipe("f1", 104.26), x, _b_spectrum(x, 0.0, [100, 60], 11)),
                               (_b_recipe("f2", 195.48), x, _b_spectrum(x, 0.0, [80, 90], 12))],
                              share=(ISO, "shift_fwhm_ppm"), windows=[window, window])
    pA = [r.sites[0].params[ISO].value for r in mres.recipes]
    _check(len(mres.recipes) == 2 and abs(pA[0] - pA[1]) < 1e-9 and abs(pA[0] - 15.0) < 0.3
           and all(r < 0.1 for r in mres.rmsd) and [d["kind"] for d in mres.per_dataset] == ["1d", "1d"],
           f"fit_multi: positions {pA}, rmsd {mres.rmsd}")
    say(f"multifit over 104.26 / 195.48 MHz: shared A at {pA[0]:.2f} ppm, rmsd {mres.rmsd[0]:.3f} / {mres.rmsd[1]:.3f} "
        f"in {time.perf_counter() - t0:.1f} s")
    # batch fit: amplitudes free, positions released by 20 %, widths held
    t0 = time.perf_counter()
    bres = batchfit.batch_fit(entries, release=(ISO,), release_frac=0.2)
    posA = [r.sites[0].params[ISO].value for r in bres.recipes]
    fwhm = [r.sites[0].params["shift_fwhm_ppm"].value for r in bres.recipes]
    ampA = [r.sites[0].params["amplitude"].value for r in bres.recipes]
    # the shared width stays at the model's 5 ppm while the synthetic truth
    # is 6 ppm, so every absolute amplitude absorbs the same mismatch (about
    # +9 %); what a batch fit is FOR -- the amplitudes relative to each
    # other across the series -- must come out right (truth 100 : 80 : 120)
    ratios = [a / ampA[0] for a in ampA]
    _check(bres.released == (ISO,) and len(bres.recipes) == 3 and all(abs(w - 5.0) < 1e-9 for w in fwhm)
           and posA[1] > posA[0] > posA[2]
           and all(abs(r - t) < 0.03 for r, t in zip(ratios, (1.0, 0.8, 1.2))),
           f"batch_fit: positions {posA}, widths {fwhm}, amplitudes {ampA} (ratios {ratios})")
    data = [(p, a, w) for _r, p, a, w in entries]
    batchfit.batch_error_analysis(bres, data, method="covariance")
    rows = batchfit.error_table(bres, method="covariance")
    amp_rows = [r for r in rows if r["param"] == "amplitude"]
    _check(bres.error_method == "covariance" and amp_rows and all(r["stderr"] is not None for r in amp_rows),
           "batch covariance errors missing")
    line = f"batch_fit x3 in {time.perf_counter() - t0:.1f} s: A at {[round(p, 2) for p in posA]} ppm (released), covariance errors"
    if not quick:
        batchfit.batch_error_analysis(bres, data, method="montecarlo", n_trials=8, seed=1, parallel=False)
        _check(bres.error_method == "montecarlo" and "covariance" in bres.error_detail, "batch Monte-Carlo errors")
        line += " + Monte-Carlo (8 trials)"
    say(line)
    batchfit.write_shared_csv(bres, out / "batch_table.csv")
    batchfit.write_error_csv(bres, out / "batch_errors.csv", "covariance")
    for name in ("batch_table.csv", "batch_errors.csv"):
        head = (out / name).read_text(encoding="utf-8").splitlines()[0]
        _check(head.startswith("scope,") and (out / name).stat().st_size > 200, f"{name}: {head}")
    shared = batchfit.shared_table(bres)
    cols, piv = batchfit.pivot_by_spectrum(shared, 3)
    back = batchfit.result_from_dict(batchfit.result_to_dict(bres))
    _check(any(r["scope"] == "shared" for r in shared) and cols and back.labels == bres.labels and len(back.recipes) == 3,
           "shared_table / pivot / result dict round trip")
    bout = bundle.write_bundle(bres, data, out / "bundle", kind="batch")
    names = {p.name for p in (out / "bundle").iterdir()}
    _check({"README.txt", "manifest.csv", "batch_table.csv"} <= names and len(bout.manifest) == 3 and bout.warnings == []
           and any(n.endswith("_curves.csv") for n in names) and any(n.endswith(".recipe.json") for n in names),
           f"write_bundle: {sorted(names)} warnings {bout.warnings}")
    say(f"batch CSVs + publication bundle: {len(names)} files ({bout.summary[:60]})")
    # the sequential series fit through the series-mode plumbing
    model = _b_recipe("series").to_dict()
    members = [{"recipe": copy.deepcopy(model), "ppm": x, "amp": _b_spectrum(x, sh, am, 20 + k), "name": f"s{k}"}
               for k, (sh, am) in enumerate(zip(shifts, amps))]
    sweep = seriesmode.entries_for_sweep(members, default_window=window)
    _check(len(sweep) == 3 and sweep[0][3] == window, "entries_for_sweep")
    t0 = time.perf_counter()
    sres = seqfit.run_sequential(sweep, passes=1, start="first", propagate=(ISO, "shift_fwhm_ppm"))
    applied = seriesmode.apply_sweep_result(members, sres)
    pos = [a["recipe"]["sites"][0]["params"][ISO]["value"] for a in applied]
    _check(all(a is not None for a in applied) and sres.passes == 1 and len(sres.history) == 1
           and all(abs(p - (15.0 + sh)) < 0.5 for p, sh in zip(pos, shifts))
           and all(a["recipe"]["fit_rmsd"] is not None and len(a["y_fit"]) == len(a["x"]) for a in applied),
           f"run_sequential / apply_sweep_result: positions {pos}")
    say(f"seqfit (1 pass, warm-started) in {time.perf_counter() - t0:.1f} s: A at {[round(p, 2) for p in pos]} ppm, "
        f"mean rmsd {sres.history[0]['mean']:.3f}")
    new, note = seriesmode.carry_into(members[2]["recipe"], applied[0]["recipe"], (ISO,),
                                      seriesmode.amp_max(members[2]["amp"]), seriesmode.amp_max(members[0]["amp"]),
                                      src_name="s0")
    _check(abs(new["sites"][0]["params"][ISO]["value"] - pos[0]) < 1e-9 and note, f"carry_into: {note}")
    cands = seriesmode.carry_candidates([m["recipe"] for m in members])
    _check("shift_fwhm_ppm" in cands and "amplitude" not in cands
           and seriesmode.member_status(applied[0]["recipe"]) == "fitted"
           and seriesmode.member_status(model) == "unfitted", "carry_candidates / member_status")
    _check(components.pair_sites([{"label": "A"}, {"label": "B"}, {"label": "C"}], [{"label": "C"}, {"label": "A"}])
           == [1, None, 0], "components.pair_sites")
    comps = components.component_map([m["recipe"] for m in members])
    _check(components.display_names(comps) == ["A", "B"], f"component_map: {components.display_names(comps)}")
    say(f"seriesmode: carry_into -> {note!r}; carry candidates {cands}; components {components.display_names(comps)}")
    # the publication report on saved fits (CSV spectra + .recipe.json)
    paths = []
    for k in range(2):
        yk = _b_spectrum(x, shifts[k], amps[k], 30 + k)
        csv = _write_csv_spectrum(out / f"glass{k}.csv", x, yk)
        rec = _b_recipe(f"glass{k}")
        rec.source_path, rec.source_kind, rec.fit_window_ppm = str(csv), "csv", window
        rec.save(out / f"glass{k}.recipe.json")
        paths.append(str(out / f"glass{k}.recipe.json"))
    t0 = time.perf_counter()
    rep = batch.run_batch(paths, out / "report", error_method="covariance", make_plots=not quick,
                          formats=("csv", "latex", "markdown"))
    got = {Path(f).name for f in rep.files}
    _check(rep.n_fits == 2 and rep.n_sites == 4 and {"table.csv", "table.tex", "report.md"} <= got, f"run_batch: {rep}")
    loaded, warns = batch.load_entries(paths)
    table = batch.build_table(loaded, error_method="covariance")
    _check(len(loaded) == 2 and len(table.rows) == 4 and table.headers and table.error_method == "covariance",
           f"load_entries / build_table: {len(loaded)} entries, {len(table.rows)} rows")
    say(f"batch report in {time.perf_counter() - t0:.1f} s: {sorted(got)}; table {len(table.rows)} rows x "
        f"{len(table.headers)} columns" + ("" if quick else " (+ overlay figures)"))


# ------------------------------------------------------------ the list
def stages() -> list:
    return [("engine-readers", readers), ("engine-processing", processing_stage),
            ("engine-models", models_stage), ("engine-physics", physics_stage),
            ("engine-fitting", fitting_stage), ("engine-error-tools", error_tools_stage),
            ("engine-constraints", constraints_stage), ("engine-multi", multi_stage)]
