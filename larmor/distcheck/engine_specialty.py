"""Engine stage ``engine-specialty`` of the distribution check: the
specialist tools on synthetic data with a known truth -- relaxation series
(satrec / series), REDOR, variable temperature, VOCS stitching, QCPMG and the
multi-field centre-of-gravity extrapolation, the shielding calibration, DFT
(.magres) import, the SIMPSON bridge, the 2D MQMAS engine, the correlation
decomposition, and what engine-physics / engine-error-tools leave unchecked in
the physics helpers. The truths are physical laws and constants (the exact
isolated-pair REDOR curve, Arrhenius, e Q Vzz / h, the second-order
quadrupolar shifts, Maricq-Waugh moments, IUPAC receptivities), closed forms
and numpy. Qt-free, deterministic, every file under ``ctx["out"]``.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from larmor.distcheck.engines import _check, _expno, _import_mrsimulator, _out

__all__ = ["specialty_stage"]


def specialty_stage(say, ctx):
    """The relaxation / dipolar / series tools, then the structure / 2D /
    correlation tools and the physics helpers' remainder."""
    _specialty_series(say, ctx)
    _specialty_struct(say, ctx)


#: IUPAC 2001 gyromagnetic ratios (Harris et al., Pure Appl. Chem. 73, 1795),
#: rad s^-1 T^-1 -- typed here, independently of larmor.redor's MHz/T table
_SSER_GAMMA = {"13C": 6.728284e7, "15N": -2.712618e7}
#: mu0/4pi (T m A^-1), hbar (J s) and R (J mol^-1 K^-1), CODATA 2018
_SSER_MU0_4PI, _SSER_HBAR, _SSER_R = 1.0e-7, 1.054571817e-34, 8.314462618


def _specialty_series(say, ctx):
    """The relaxation / dipolar / series tools on synthetic data with a known
    truth: the saturation-recovery T1 path (a Bruker pseudo-2D ser + vdlist
    written with nmrglue -> process_slices -> integrate_zones -> fit_buildup,
    satrec / series analyze, per-site T1 of two overlapping lines, the keep
    mask, stretched / inversion / CPMG kinds, TopSpin's ct1t2.txt), REDOR
    (the dipolar constant against the physical formula, the pair and
    short-time fits, the redor.txt route), the Arrhenius / VFT fits, the VOCS
    stitch, QCPMG (spikelet comb at 1/tau, the period, T2, the sum-echo
    absorption), the multi-field CT centre-of-gravity extrapolation and the
    shielding-to-shift calibration."""
    out = _out(ctx, "specialty_series")
    _sser_relaxation(say, out)
    _sser_redor(say, out)
    _sser_vt(say)
    _sser_vocs(say)
    _sser_qcpmg(say)
    _sser_fields(say)
    _sser_shiftcal(say, out)


# ------------------------------------------------------------ relaxation
def _sser_write_ser(folder: Path, delays_tok, ser, sfo1: float, sw: float) -> None:
    """A Bruker pseudo-2D EXPNO (acqus/acqu2s/ser/vdlist), float64 data, no
    digital filter (GRPDLY 0, DSPFVS 20), written with nmrglue."""
    import nmrglue as ng

    folder.mkdir(parents=True, exist_ok=True)
    udic = ng.fileiobase.create_blank_udic(2)
    udic[0].update({"size": ser.shape[0], "complex": False, "encoding": "states", "sw": 1.0,
                    "obs": sfo1, "car": 0.0, "label": "F1", "time": True, "freq": False})
    udic[1].update({"size": ser.shape[1], "complex": True, "encoding": "direct", "sw": sw,
                    "obs": sfo1, "car": 0.0, "label": "1H", "time": True, "freq": False})
    dic = ng.bruker.create_dic(udic)
    dic["acqus"].update({"NUC1": "<1H>", "SFO1": sfo1, "BF1": sfo1, "SW_h": sw, "O1": 0.0,
                         "PULPROG": "<satrec_zg>", "DECIM": 1, "DSPFVS": 20, "GRPDLY": 0,
                         "DTYPA": 2, "BYTORDA": 0, "NS": 8})
    dic["acqu2s"].update({"TD": ser.shape[0]})
    ng.bruker.write(str(folder), dic, ser, overwrite=True)
    (folder / "vdlist").write_text("\n".join(delays_tok) + "\n", encoding="ascii")


def _sser_relaxation(say, out: Path) -> None:
    from larmor import satrec, series
    from larmor.recipe import Param, Recipe, SiteModel

    sfo1, sw, npts, lb, scale, phi = 100.0, 10000.0, 1024, 20.0, 1.0e6, 30.0
    # (label, offset Hz, T2 s, M0, T1 s): A alone; B and C 0.9 ppm apart
    # (FWHM 0.6 ppm each after LB), so one window cannot separate them
    lines = (("A", 1500.0, 0.010, 1.0, 0.40), ("B", -2000.0, 0.008, 0.8, 0.25),
             ("C", -2090.0, 0.008, 0.6, 4.00))
    tok = ("5m", "20m", "50m", "100m", "250m", "0.5", "1s", "2.5", "6", "15s")
    delays = np.array([0.005, 0.02, 0.05, 0.1, 0.25, 0.5, 1.0, 2.5, 6.0, 15.0])
    t = np.arange(npts) / sw
    rng = np.random.default_rng(3)
    ser = np.zeros((delays.size, npts), complex)
    for k, tau in enumerate(delays):
        for _lab, f, t2, m0, t1 in lines:
            ser[k] += m0 * (1.0 - np.exp(-tau / t1)) * np.exp(2j * np.pi * f * t - t / t2)
    noise = 2.0e-3 * (rng.standard_normal(ser.shape) + 1j * rng.standard_normal(ser.shape))
    ser = (ser * np.exp(1j * np.deg2rad(phi)) + noise) * scale     # + a receiver phase
    exp = out / "satrec" / "11"
    _sser_write_ser(exp, tok, ser, sfo1, sw)

    vals, src = series.read_delays(exp)
    _check(src == "vdlist" and vals.shape == delays.shape and np.allclose(vals, delays, rtol=1e-12, atol=0.0),
           f"read_delays: {src} {vals} (truth {delays})")
    # the series kind from acqus PULPROG (what kind=None callers rely on): a
    # CPMG program beside it, so a fallback to 'satrec' cannot pass for detection
    cp = out / "satrec" / "12"
    cp.mkdir(parents=True, exist_ok=True)
    (cp / "acqus").write_text("##TITLE= distcheck\n##JCAMPDX= 5.0\n##$PULPROG= <cpmg1d>\n##$TD= 64\n##END=\n",
                              encoding="ascii")
    kinds = (series.detect_kind(exp), series.detect_kind(cp))
    _check(kinds == ("satrec", "cpmg"), f"detect_kind from PULPROG: satrec_zg -> {kinds[0]!r}, cpmg1d -> {kinds[1]!r}")
    x, slices = satrec.process_slices(exp, lb_hz=lb)
    _check(slices.shape == (delays.size, x.size) and np.all(np.diff(x) > 0),
           f"process_slices: slices {slices.shape}, axis {x.size} points")
    step = float(x[1] - x[0])
    near = np.abs(x - 15.0) < 3.0
    peak_a = float(x[near][np.argmax(slices[-1][near])])
    _check(abs(peak_a - 15.0) <= step, f"line A at {peak_a:.3f} ppm, truth 15.000 (bin {step:.3f} ppm)")
    # DFT sum rule: the points of a spectrum add up to N x the first FID point,
    # which process_slices halves (FCOR) -- so a correctly phased relaxed slice
    # sums to N/2 x the total magnetisation it started from
    m_last = scale * sum(m0 * (1.0 - np.exp(-delays[-1] / t1)) for *_r, m0, t1 in lines)
    sum_ratio = float(slices[-1].sum() / (0.5 * slices.shape[1] * m_last))
    _check(abs(sum_ratio - 1.0) < 0.01, f"relaxed-slice sum = {sum_ratio:.4f} x N/2 x M (phase or scale lost)")
    say(f"satrec EXPNO: nmrglue ser {ser.shape} + vdlist in m/s units read back exactly; process_slices "
        f"{slices.shape}, line A at {peak_a:.3f} ppm (truth 15.000), relaxed slice sums to "
        f"{sum_ratio:.4f} x N/2 x M0 (DFT sum rule)")

    # window integration: zone A alone, one zone over the B/C pair
    zone_a, zone_bc = (21.0, 9.0), (-14.45, -26.45)
    ig = series.integrate_zones(x, slices, [zone_a, zone_bc])
    _check(ig.shape == (2, delays.size), f"integrate_zones shape {ig.shape}")
    ra = series.fit_buildup(delays, ig[0], kind="satrec")
    m0_fit = float(ra["curve"](1.0e6) * ra["norm"])
    # the zone holds the Lorentzian (FWHM 1/(pi T2) + LB) out to +-600 Hz:
    # its fraction of the area is (2/pi) atan(600 / HWHM)
    hwhm = 0.5 * (1.0 / (np.pi * lines[0][2]) + lb)
    m0_true = 0.5 * slices.shape[1] * scale * lines[0][3] * (2.0 / np.pi) * np.arctan(600.0 / hwhm)
    out_t1, _curve = satrec.fit_t1(delays, ig[0])
    t1_fit, i0_fit = float(out_t1.params["t1"].value), float(out_t1.params["i0"].value)
    an = satrec.analyze(exp, window_ppm=zone_a, lb_hz=lb)
    an_mag = satrec.analyze(exp, window_ppm=zone_a, lb_hz=lb, mode="magnitude")
    sa = series.analyze(exp, kind="satrec", window_ppm=zone_a, lb_hz=lb)
    got = {"fit_buildup": ra["tau"], "satrec.fit_t1": t1_fit, "satrec.analyze": an.t1_s,
           "satrec.analyze(magnitude)": an_mag.t1_s, "series.analyze": sa.tau}
    for name, v in got.items():   # FID noise at 0.2 % of M0 moves T1 by ~0.3 %: 3 % is a wide margin
        _check(abs(v / 0.40 - 1.0) < 0.03, f"{name}: T1 {v:.4f} s, truth 0.400 s")
    _check(abs(m0_fit / m0_true - 1.0) < 0.015 and abs(i0_fit - 1.0) < 0.02,
           f"zone A M0 {m0_fit:.4g} vs analytic {m0_true:.4g}; fit_t1 i0 {i0_fit:.4f} (truth 1)")
    _check(an.summary.startswith("T1 = ") and sa.kind == "satrec" and sa.tau_name == "T1",
           f"result labels: {an.summary!r} / {sa.kind} {sa.tau_name}")
    say("T1 zone A (truth 0.400 s): " + ", ".join(f"{k} {v:.4f}" for k, v in got.items())
        + f"; M0 {m0_fit / m0_true:.4f} x the analytic zone area, fit_t1 i0 {i0_fit:.4f}")

    # per-site: NNLS on the fitted lineshapes separates B and C
    fwhm = {lab: (1.0 / (np.pi * t2) + lb) / sfo1 for lab, _f, t2, _m, _t in lines}
    rec = Recipe(nucleus="1H", larmor_frequency_MHz=sfo1, spin_rate_Hz=0.0, sample="distcheck", sites=[
        SiteModel(model="gauss_lor", label=lab, params={
            "isotropic_chemical_shift_ppm": Param(f / sfo1), "shift_fwhm_ppm": Param(fwhm[lab]),
            "amplitude": Param(1.0), "gl": Param(0.0, vary=False)})
        for lab, f, _t2, _m, _t in lines])
    per = series.analyze_per_site(exp, rec, kind="satrec", lb_hz=lb)
    _check([r.label for r in per] == ["A", "B", "C"], f"analyze_per_site labels {[r.label for r in per]}")
    for r, (lab, _f, _t2, _m, t1) in zip(per, lines):
        _check(abs(r.tau / t1 - 1.0) < 0.03, f"per-site {lab}: T1 {r.tau:.4f} s, truth {t1:.3f} s")
    mixed = series.fit_buildup(delays, ig[1], kind="satrec")["tau"]
    say(f"analyze_per_site (B/C overlap, 0.9 ppm apart): A {per[0].tau:.4f}, B {per[1].tau:.4f}, "
        f"C {per[2].tau:.3f} s (truth 0.40 / 0.25 / 4.00); one window over B+C gives {mixed:.3f} s")

    # the keep mask (the dialog's click-to-exclude): a corrupted point at 0.5 s
    bad = ig[0].copy()
    bad[5] *= 0.6
    keep = np.ones(delays.size, bool)
    keep[5] = False
    t_bad = series.fit_buildup(delays, bad, kind="satrec")["tau"]
    t_keep = series.fit_buildup(delays, bad, keep=keep, kind="satrec")["tau"]
    _check(abs(t_bad / 0.40 - 1.0) > 0.10 and abs(t_keep / 0.40 - 1.0) < 0.03,
           f"keep mask: corrupted fit {t_bad:.4f} s, with the point excluded {t_keep:.4f} s (truth 0.400)")
    say(f"fit_buildup keep mask: a 40 %-low point at 0.5 s pulls T1 to {t_bad:.4f} s; "
        f"excluded -> {t_keep:.4f} s (truth 0.400)")

    # the other series kinds on seeded build-up / decay curves
    ts = np.logspace(-2, 1.5, 14)
    ys = 1.0 - np.exp(-((ts / 3.0) ** 0.6)) + rng.normal(0.0, 2.0e-3, ts.size)
    st = series.fit_buildup(ts, ys, kind="satrec", stretched=True)
    st2, _c2 = satrec.fit_t1(ts, ys, stretched=True)
    ti = np.array([0.01, 0.1, 0.3, 0.6, 1.0, 1.5, 2.5, 4.0, 8.0, 20.0])
    inv = series.fit_buildup(ti, 1.0 - 2.0 * np.exp(-ti / 1.6) + rng.normal(0.0, 2.0e-3, ti.size), kind="invrec")
    null = float(inv["curve"](1.6 * np.log(2.0)))          # inversion recovery crosses zero at T1 ln 2
    tc = np.linspace(0.001, 0.06, 12)
    t2r = series.fit_buildup(tc, np.exp(-tc / 0.012) + rng.normal(0.0, 2.0e-3, tc.size), kind="cpmg")
    _check(abs(st["tau"] / 3.0 - 1.0) < 0.03 and abs(st["beta"] - 0.6) < 0.03
           and abs(st2.params["t1"].value / 3.0 - 1.0) < 0.03 and abs(st2.params["beta"].value - 0.6) < 0.03,
           f"stretched: series tau {st['tau']:.3f} beta {st['beta']:.3f}, satrec t1 "
           f"{st2.params['t1'].value:.3f} beta {st2.params['beta'].value:.3f} (truth 3.0 / 0.6)")
    _check(abs(inv["tau"] / 1.6 - 1.0) < 0.02 and abs(null) < 0.01 and abs(t2r["tau"] / 0.012 - 1.0) < 0.02,
           f"invrec T1 {inv['tau']:.4f} s (1.6), null {null:.4f}; cpmg T2 {t2r['tau'] * 1e3:.3f} ms (12)")
    say(f"kinds: stretched tau {st['tau']:.3f} s, beta {st['beta']:.3f} (truth 3.0 / 0.6; fit_t1 "
        f"{st2.params['t1'].value:.3f} / {st2.params['beta'].value:.3f}); invrec T1 {inv['tau']:.4f} s (1.6), "
        f"curve at T1 ln2 {null:+.4f}; cpmg T2 {t2r['tau'] * 1e3:.3f} ms (12.0)")

    # TopSpin's ct1t2.txt (the quantitativity chip's T1 source), both layouts
    legacy = out / "satrec" / "ct1t2_legacy.txt"
    legacy.write_text(
        "Dataset :\n distcheck/11/pdata/1\nAREA fit :\n I[t]=I[0]+P*exp(-t/T1)\n\n"
        "10 points for Integral 1,  Integral Region from 21.000 to 9.000 ppm\nResults     Comp. 1\n\n"
        f"I[0]  =    1.000e+00\nP     =   -1.000e+00\nT1    =     {ra['tau'] * 1e3:.3f}m\n"
        "SD    =    1.000e-03\n\n", encoding="ascii")
    simfit = out / "satrec" / "ct1t2_simfit.txt"
    simfit.write_text(
        "SIMFIT RESULTS\n==============\n\nAREA fit : Saturation-Recovery(T1) : I[t]=I[0](1-exp(-t/T1))\n\n"
        "10 points for Integral 2,  Integral Region from -14.450 to -26.450 ppm\n\n"
        "Results     Comp. 1       Comp. 2\n\nI[0]  =    8.000e-01     6.000e-01\n"
        f"T1    =      {per[1].tau * 1e3:.3f}m        {per[2].tau:.3f}s\n\nRSS   =    1.0e-04\n"
        "SD    =    1.0e-03\n", encoding="ascii")
    reg_l, reg_s = satrec.read_ct1t2(legacy), satrec.read_ct1t2(simfit)
    _check(len(reg_l) == 1 and reg_l[0].index == 1 and (reg_l[0].hi_ppm, reg_l[0].lo_ppm) == (21.0, 9.0)
           and abs(reg_l[0].t1_s - ra["tau"]) < 1e-6, f"read_ct1t2 legacy: {reg_l}")
    _check(len(reg_s) == 1 and len(reg_s[0].components_s) == 2 and abs(reg_s[0].t1_s - per[2].tau) < 1e-3
           and abs(min(reg_s[0].components_s) - per[1].tau) < 1e-6, f"read_ct1t2 SIMFIT: {reg_s}")
    say(f"read_ct1t2: legacy T1 {reg_l[0].t1_s:.4f} s over {reg_l[0].hi_ppm:g}..{reg_l[0].lo_ppm:g} ppm; "
        f"SIMFIT components {[round(c, 4) for c in reg_s[0].components_s]} s -> governing {reg_s[0].t1_s:.3f} s")


# ------------------------------------------------------------ REDOR
def _sser_redor(say, out: Path) -> None:
    from larmor import convert, redor

    r0 = 2.50
    d = redor.dipolar_constant_hz("13C", "15N", r0)
    g1, g2 = _SSER_GAMMA["13C"], _SSER_GAMMA["15N"]
    d_phys = _SSER_MU0_4PI * abs(g1 * g2) * _SSER_HBAR / (2.0 * np.pi * (r0 * 1e-10) ** 3)
    d_conv = abs(convert.dipolar_Hz(g1 / (2e6 * np.pi), g2 / (2e6 * np.pi), r0))   # signed b/2pi there
    ratio8 = redor.dipolar_constant_hz("13C", "15N", r0 / 2.0) / d
    r_back = redor.distance_angstrom("13C", "15N", d)
    # the module's table carries 5-6 significant digits of gamma: 1e-4 relative
    _check(abs(d / d_phys - 1.0) < 1e-4 and abs(d_conv / d_phys - 1.0) < 1e-9,
           f"dipolar constant 13C-15N at {r0} A: redor {d:.3f} Hz, convert {d_conv:.3f} Hz, "
           f"mu0/4pi gI gS hbar / (2 pi r^3) = {d_phys:.3f} Hz")
    _check(abs(ratio8 - 8.0) < 1e-9 and abs(r_back - r0) < 1e-9,
           f"1/r^3 law: D(r/2)/D(r) = {ratio8:.6f}; distance round trip {r_back:.6f} A")
    say(f"REDOR constants: 13C-15N at {r0:.2f} A, D = {d:.2f} Hz = (mu0/4pi) gI gS hbar/(2 pi r^3) "
        f"{d_phys:.2f} Hz (IUPAC gammas) = convert.dipolar_Hz; D(r/2)/D(r) = {ratio8:.4f}; r round trip exact")

    # the exact isolated-pair curve, from outside the module (Mueller et al.,
    # J. Magn. Reson. A 113, 81, 1995): S/S0 = (sqrt2 pi/4) J_1/4(sqrt2 lam)
    # J_-1/4(sqrt2 lam), lam = D N Tr with D the dipolar constant in Hz
    from scipy.special import jv

    def mueller(lam):
        x = np.sqrt(2.0) * np.asarray(lam, float)
        with np.errstate(invalid="ignore"):
            s = (np.sqrt(2.0) * np.pi / 4.0) * jv(0.25, x) * jv(-0.25, x)
        return np.where(x > 0, 1.0 - s, 0.0)

    lam = np.linspace(0.0, 3.0, 61)
    dev = float(np.abs(redor.redor_pair_curve(1.0, lam) - mueller(lam)).max())
    # the module's 200 x 100 powder quadrature reproduces it to ~5e-4; a curve
    # twice too fast in lambda (the I-spin line at +-D instead of +-D/2)
    # deviates by ~0.6
    _check(dev < 2e-3, f"REDOR powder curve vs the exact isolated-pair form: max deviation {dev:.2e}")
    # a dephasing made from the exact form for the 2.50 A pair, + noise
    rng = np.random.default_rng(13)
    masr = 10000.0
    n_cyc = np.arange(2, 42, 2)
    ntr = n_cyc / masr
    ds = mueller(d_phys * ntr) + rng.normal(0.0, 0.003, ntr.size)
    pair = redor.analyze(ntr, ds, pair=("13C", "15N"), regime="pair")
    m2_vv = 0.75 * (4.0 / 15.0) * (2.0 * np.pi * d_phys) ** 2     # van Vleck, one spin-1/2 partner
    _check(pair.regime == "isolated pair" and pair.n_used == ntr.size and abs(pair.d_hz / d_phys - 1.0) < 0.015
           and abs(pair.distance_A - r0) < 0.01 and abs(pair.m2 / m2_vv - 1.0) < 0.03,
           f"REDOR pair fit: D {pair.d_hz:.2f} Hz (truth {d_phys:.2f}), r {pair.distance_A:.4f} A (truth {r0}), "
           f"M2 {pair.m2:.4g} (van Vleck {m2_vv:.4g})")
    low = ds < 0.2
    short = redor.analyze(ntr, ds, pair=("13C", "15N"), regime="short")
    # the parabola is the second-order expansion, ~6 % above the exact curve
    # at dS/S0 = 0.2, so D comes out up to ~3 % low; + noise
    _check(short.regime == "short-time parabola" and short.n_used == int(low.sum()) >= 3
           and abs(short.d_hz / d_phys - 1.0) < 0.05, f"REDOR short-time fit: D {short.d_hz:.2f} Hz on "
           f"{short.n_used} points (truth {d_phys:.2f})")
    auto = redor.analyze(ntr, ds, regime="auto")
    up = ~low
    auto_up = redor.analyze(ntr[up], ds[up], regime="auto")
    _check(auto.regime == "short-time parabola" and any("auto-selected" in n for n in auto.notes)
           and auto_up.regime == "isolated pair" and abs(auto_up.d_hz / d_phys - 1.0) < 0.015,
           f"REDOR auto regime: {auto.regime} / upper branch {auto_up.regime} D {auto_up.d_hz:.2f} Hz")
    lam = np.array([0.005, 0.2])
    para = redor.short_time_curve(1.0, lam) / redor.redor_pair_curve(1.0, lam)
    _check(abs(para[0] - 1.0) < 1e-3 and para[1] > para[0],
           f"the parabola is not the small-lambda limit of the powder curve: ratio {para}")
    say(f"REDOR: powder curve = the exact isolated-pair form (Mueller 1995) to {dev:.1e}; dephasing of the "
        f"{r0:.2f} A pair from that form (D {d_phys:.2f} Hz, noise 0.003): pair fit {pair.d_hz:.2f} Hz -> r "
        f"{pair.distance_A:.4f} A, M2 {pair.m2:.4g} rad^2/s^2 (van Vleck {m2_vv:.4g}); short-time {short.d_hz:.2f} Hz "
        f"on {short.n_used} points; auto -> {auto.regime}, upper branch only -> {auto_up.regime} "
        f"{auto_up.d_hz:.2f} Hz; parabola/curve {para[0]:.5f} at lambda 0.005")

    # the TopSpin redor.txt route (the dialog and `larmor redor`)
    pdata = out / "redor" / "9" / "pdata" / "1"
    pdata.mkdir(parents=True, exist_ok=True)
    s0 = 2.0e10
    rows = "\n".join(f"{n:10d} {s0:14.6e} {s0 * (1.0 - v):14.6e} {0:13d} {0:13d}" for n, v in zip(n_cyc, ds))
    (pdata / "redor.txt").write_text(
        f"Dataset :\ndistcheck/9\n\nSpinning speed : {masr:.6f}\n\nPeak 1 (xy 172.1 ppm)\n\n"
        "    number integral(S0) integral(S*) intensity(S0) intensity(S*)\n\n" + rows + "\n", encoding="ascii")
    n_b, ds_b, masr_b = redor.read_redor_txt(pdata / "redor.txt")
    _check(masr_b == masr and np.array_equal(n_b, n_cyc) and np.allclose(ds_b, ds, atol=2e-6),
           f"read_redor_txt: MAS {masr_b}, {n_b.size} rows, max |dS/S0 error| {np.max(np.abs(ds_b - ds)):.2g}")
    fe = redor.analyze_expno(out / "redor" / "9", pair=("13C", "15N"), regime="pair")
    _check(abs(fe.d_hz / d - 1.0) < 0.015 and abs(fe.distance_A - r0) < 0.01,
           f"analyze_expno: D {fe.d_hz:.2f} Hz, r {fe.distance_A:.4f} A (truth {d:.2f} Hz, {r0} A)")
    say(f"REDOR redor.txt ({n_b.size} rows, MAS {masr_b:.0f} Hz) -> analyze_expno D {fe.d_hz:.2f} Hz, "
        f"r {fe.distance_A:.4f} A; {fe.summary}")


# ------------------------------------------------------------ VT
def _sser_vt(say) -> None:
    from larmor import vt

    ea, tau0 = 35.0e3, 2.0e-13                      # J/mol, s
    temps = np.linspace(220.0, 360.0, 9)
    tau = tau0 * np.exp(ea / (_SSER_R * temps))     # tau_c(T); the rate is 1/tau
    arr = vt.fit_arrhenius(temps, 1.0 / tau)
    t_mid = 287.3
    k_mid = float(arr["curve"](t_mid)) * tau0 * np.exp(ea / (_SSER_R * t_mid))
    _check(abs(arr["Ea_kJmol"] / 35.0 - 1.0) < 1e-6 and abs(arr["A"] * tau0 - 1.0) < 1e-6
           and abs(k_mid - 1.0) < 1e-9, f"Arrhenius: Ea {arr['Ea_kJmol']:.6f} kJ/mol (35), A {arr['A']:.6e} "
           f"(truth {1.0 / tau0:.6e}), curve at {t_mid} K off by {k_mid - 1.0:.2g}")
    k = 1.0e12 * np.exp(-900.0 / (temps - 150.0))
    vft = vt.fit_vft(temps, k)
    _check(abs(vft["B_K"] / 900.0 - 1.0) < 0.005 and abs(vft["T0_K"] - 150.0) < 0.5
           and abs(vft["A"] / 1.0e12 - 1.0) < 0.02, f"VFT: B {vft['B_K']:.2f} K (900), T0 {vft['T0_K']:.2f} K "
           f"(150), A {vft['A']:.4g} (1e12)")
    say(f"vt: Arrhenius on 1/tau_c(T) at {temps.size} temperatures {temps[0]:.0f}-{temps[-1]:.0f} K: Ea "
        f"{arr['Ea_kJmol']:.4f} kJ/mol (truth 35), A {arr['A']:.4e} s^-1 (truth {1.0 / tau0:.4e}); VFT B "
        f"{vft['B_K']:.2f} K (900), T0 {vft['T0_K']:.2f} K (150)")


# ------------------------------------------------------------ VOCS
def _sser_vocs(say) -> None:
    from larmor import vocs

    def truth(x):
        return np.exp(-0.5 * ((x + 300.0) / 110.0) ** 2) + 1.3 * np.exp(-0.5 * ((x - 160.0) / 70.0) ** 2)

    def piece(centre, half, n, roll=0.10):
        """The truth inside one offset window, with a sin^2 excitation
        roll-off over the outer `roll` of each side."""
        x = np.linspace(centre - half, centre + half, n)
        y = truth(x)
        edge = roll * half
        for sgn in (-1.0, 1.0):
            dist = np.abs(x - (centre + sgn * half))
            m = dist < edge
            y[m] *= np.sin(0.5 * np.pi * dist[m] / edge) ** 2
        return x, y

    # offsets 330 ppm apart, 560 ppm wide, the third piece sampled twice as
    # finely; trimming 15 % per edge removes the 10 % roll-off and still
    # leaves ~60 ppm of clean overlap at each junction
    pieces = [piece(-480.0, 280.0, 601), piece(-150.0, 280.0, 601), piece(180.0, 280.0, 1201)]
    finest = 560.0 / 1200
    errs = {}
    for method in ("skyline", "average"):
        res = vocs.stitch(pieces, method=method, trim_frac=0.15)
        cov = res.coverage > 0
        tr = truth(res.ppm)
        err = float(np.max(np.abs(res.amp[cov] - tr[cov])) / tr.max())
        area = float(np.trapezoid(res.amp, res.ppm) / np.trapezoid(tr, res.ppm))
        dx = float(np.median(np.diff(res.ppm)))
        _check(cov.all() and res.coverage.max() == 2 and err < 0.01 and abs(area - 1.0) < 0.005
               and dx <= 1.05 * finest and len(res.notes) == 1,
               f"vocs {method}: max error {err:.3g} of the peak, area {area:.4f} x truth, grid {dx:.3f} ppm "
               f"(finest {finest:.3f}), coverage max {res.coverage.max()}, notes {res.notes}")
        errs[method] = (err, area)
    far = vocs.stitch([piece(-500.0, 150.0, 301), piece(500.0, 150.0, 301)], trim_frac=0.1)
    gap = (far.ppm > -300.0) & (far.ppm < 300.0)
    _check(any("gap" in n.lower() for n in far.notes) and np.all(far.amp[gap] == 0.0)
           and np.all(far.coverage[gap] == 0), f"vocs gap: notes {far.notes}")
    say(f"vocs: 3 offset pieces (10 % roll-off, trim 15 %): skyline max error {errs['skyline'][0]:.1e} of "
        f"the peak, area {errs['skyline'][1]:.5f} x truth; average {errs['average'][0]:.1e}, "
        f"{errs['average'][1]:.5f}; grid at the finest piece's spacing; non-tiling offsets -> gap warned, "
        "zero filled")


# ------------------------------------------------------------ QCPMG
def _sser_qcpmg(say) -> None:
    from scipy.signal import find_peaks

    from larmor import qcpmg

    sw, sfo, period, n_echo = 128000.0, 100.0, 128, 32
    tau = period / sw                                # echo period 1 ms: spikelets every 1 kHz
    f0, fwhm, t2, phi = 2350.0, 4000.0, 0.015, 25.0
    sig_t = 2.0 * np.sqrt(np.log(2.0)) / (np.pi * fwhm)     # Gaussian echo -> Gaussian line of this FWHM (Hz)
    tt = (np.arange(period) - period // 2) / sw              # echo top at the block centre
    echo = np.exp(-((tt / sig_t) ** 2) + 2j * np.pi * f0 * tt)
    rng = np.random.default_rng(5)
    fid = np.concatenate([echo * np.exp(-k * tau / t2) for k in range(n_echo)]) * np.exp(1j * np.deg2rad(phi))
    fid = fid + 0.002 * (rng.standard_normal(fid.size) + 1j * rng.standard_normal(fid.size))

    ppm_s, spec_s = qcpmg.spikelet_spectrum(fid, sw, sfo, lb_Hz=20.0, zf=2)
    mag = np.abs(spec_s)
    bin_hz = float(ppm_s[1] - ppm_s[0]) * sfo
    pk, _ = find_peaks(mag, height=0.1 * mag.max(), distance=int(0.5 / tau / bin_hz))
    f_pk = ppm_s[pk] * sfo
    off_comb = np.abs(f_pk - np.round(f_pk * tau) / tau)
    _check(pk.size >= 5 and np.all(np.abs(np.diff(f_pk) - 1.0 / tau) <= 0.5 * bin_hz) and off_comb.max() <= 0.5 * bin_hz,
           f"spikelets at {np.round(f_pk, 1)} Hz: not a comb at multiples of 1/tau = {1.0 / tau:.1f} Hz")
    _check(abs(qcpmg.spikelet_spacing_ppm(period, sw, sfo) - 1.0 / tau / sfo) < 1e-12,
           f"spikelet_spacing_ppm {qcpmg.spikelet_spacing_ppm(period, sw, sfo)} vs 1/tau {1.0 / tau / sfo} ppm")
    found, score = qcpmg.find_period_by_correlation(fid)
    auto = qcpmg.detect_period(fid)
    cnst = [1.0] * 16
    cnst[7] = 1.0 / tau                              # Bruker CNST7 = the spikelet spacing (Hz)
    per_meta = qcpmg.echo_period_from_meta({"sw_Hz": sw, "cnst": cnst}, n_points=fid.size)
    _check(found == period and score > 0.9 and auto == period and per_meta == (float(period), "CNST7"),
           f"echo period: correlation {found} ({score:.3f}), autocorrelation {auto}, from CNST7 "
           f"{per_meta}; truth {period}")
    say(f"qcpmg spikelets: {pk.size} at {np.round(f_pk).astype(int).tolist()} Hz, spacing "
        f"{np.median(np.diff(f_pk)):.1f} Hz = 1/tau {1.0 / tau:.1f} Hz; period {found} points from the data "
        f"(score {score:.4f}) and from CNST7")

    echoes = qcpmg.split_echoes(fid, period)
    top = qcpmg.echo_top_point(echoes)
    # ssNake's C + B exp(-t/T2) on 32 noisy echo tops: the free offset C lets
    # T2 scatter by ~1 %, hence 4 %
    f = qcpmg.fit_t2(tau, qcpmg.echo_decay(echoes, top), period=period)
    _check(echoes.shape == (n_echo, period) and top == period // 2 and f.ok and abs(f.T2_s / t2 - 1.0) < 0.04
           and abs(f.lb_Hz - 1.0 / (np.pi * f.T2_s)) < 1e-9, f"echo train: split {echoes.shape}, top {top} "
           f"(truth {period // 2}), T2 {f.T2_s * 1e3:.3f} ms (truth {t2 * 1e3:.1f}), ok {f.ok}")
    ppm_e, spec_e = qcpmg.sum_echo_spectrum(fid, period, sw, sfo, t2_weight_s=f.T2_s)
    p0 = qcpmg.autophase0(spec_e)
    ph = qcpmg.phase_spectrum(spec_e, p0)
    re = ph.real
    purity = float(np.abs(ph.imag).max() / re.max())
    win = ((f0 - 2.0 * fwhm) / sfo, (f0 + 2.0 * fwhm) / sfo)
    meas = qcpmg.measure_cg(ppm_e, re, win)
    width = qcpmg.fwhm_hz(ppm_e, re, sfo, win)
    bin_e = float(ppm_e[1] - ppm_e[0]) * sfo
    _check(abs(p0 - phi) < 1.0 and purity < 0.02, f"sum echo: autophase0 {p0:.2f} deg (receiver {phi}), "
           f"|imag|/real {purity:.3f}")
    _check(abs(meas.cg_ppm - f0 / sfo) < 0.05 and not meas.flags and abs(width - fwhm) <= 2.0 * bin_e,
           f"sum echo line: CG {meas.cg_ppm:.4f} ppm (truth {f0 / sfo:.4f}), flags {meas.flags}, FWHM "
           f"{width:.1f} Hz (truth {fwhm:.0f}, bin {bin_e:.1f})")
    ppm_c, env_c = qcpmg.coadd_spectrum(fid, period, sw, sfo, lb_Hz=0.0)
    cg_c, _sig = qcpmg.centre_of_gravity(ppm_c, env_c, win)
    # magnitude rectifies the noise: a positive floor across the +-80 ppm
    # window moves the CG by a fraction of a ppm
    _check(abs(cg_c - f0 / sfo) < 0.5, f"coadd_spectrum CG {cg_c:.3f} ppm (truth {f0 / sfo:.3f})")
    # every echo has the same shape, so the spikelet heights sample the one-echo
    # spectrum -- the sum-echo absorption -- at the comb frequencies
    env_at = np.interp(ppm_s[pk], ppm_e, re / re.max())
    h = mag[pk] / mag[pk].max()
    trace = float(np.max(np.abs(h - env_at / env_at.max())))
    _check(trace < 0.02, f"spikelet heights do not trace the sum-echo envelope: max deviation {trace:.3f}")
    say(f"qcpmg sum echo: T2 {f.T2_s * 1e3:.2f} ms (truth 15.00), matched LB {f.lb_Hz:.1f} Hz; autophase0 "
        f"{p0:.2f} deg (receiver 25), |imag|/real {purity:.4f}; CG {meas.cg_ppm:.3f} ppm (truth "
        f"{f0 / sfo:.3f}), FWHM {width:.0f} Hz (truth 4000, bin {bin_e:.1f}); coadd CG {cg_c:.2f} ppm; "
        f"spikelet heights trace the envelope within {trace:.1e}")


# ------------------------------------------------------------ multi-field CG
def _sser_fields(say) -> None:
    from larmor import convert, qcpmg, qcpmg_fields as qf

    spin, diso, cq, eta = 1.5, -70.0, 3.2, 0.7          # a 35Cl-like site
    pq = cq * np.sqrt(1.0 + eta ** 2 / 3.0)
    spin_f = (spin * (spin + 1.0) - 0.75) / (spin ** 2 * (2.0 * spin - 1.0) ** 2)

    def qis(nu):
        """CT second-order quadrupolar-induced shift of the centre of gravity
        (ppm): -(3/40) (PQ/nu0)^2 [I(I+1) - 3/4] / [I^2 (2I-1)^2] 1e6."""
        return -(3.0 / 40.0) * (pq / nu) ** 2 * spin_f * 1.0e6

    rng = np.random.default_rng(11)
    sig_q_ref, nu_ref, sig_csd = 5.786, 58.79, 6.0      # quadrupolar width ~ 1/nu0^2, shift spread constant
    points, moments, dcgs = [], [], []
    for nu in (58.79, 78.354, 107.811):                  # 35Cl at 14.1 / 18.8 / 25.8 T
        dcg = diso + qis(nu)
        alt = (qf.dcg_at_field(diso, cq, nu, spin, eta),
               diso + convert.ct_second_order_shift_ppm(convert.pq_from_cq_eta(cq, eta), spin, nu))
        _check(all(abs(a - dcg) < 1e-9 for a in alt), f"dcg at {nu} MHz: Eq. 1 {dcg:.6f}, module {alt}")
        sig = float(np.hypot(sig_q_ref * (nu_ref / nu) ** 2, sig_csd))
        x = np.linspace(dcg - 400.0, dcg + 400.0, 4001)
        y = np.exp(-0.5 * ((x - dcg) / sig) ** 2) + rng.normal(0.0, 1e-3, x.size)
        win = (dcg - 6.0 * sig, dcg + 6.0 * sig)
        meas = qcpmg.measure_cg(x, y, win)
        points.append(qf.FieldPoint.from_measurement(nu, meas, source=f"synthetic {nu} MHz"))
        moments.append((nu, qcpmg.second_moment_ppm(x, y, win)))
        dcgs.append((dcg, meas.cg_ppm))
    _check(all(abs(m - t) < 0.05 for t, m in dcgs), f"measure_cg: {dcgs} (truth, measured)")
    res = qf.infinite_field_diso(points, spin=spin, eta=eta)
    _check(abs(res.delta_iso_ppm - diso) < 0.5 and abs(res.cq_MHz / cq - 1.0) < 0.01
           and abs(res.pq_MHz / pq - 1.0) < 0.01 and res.dof == 1 and not res.note,
           f"infinite_field_diso: diso {res.delta_iso_ppm:.3f} ppm (truth {diso}), CQ {res.cq_MHz:.4f} MHz "
           f"(truth {cq}), PQ {res.pq_MHz:.4f} ({pq:.4f}), dof {res.dof}, note {res.note!r}")
    ws = qf.second_moment_split(moments)
    wq_true = qf.GAUSS_FWHM_PER_SIGMA * sig_q_ref
    wcsd_true = qf.GAUSS_FWHM_PER_SIGMA * sig_csd
    _check(ws.ok and abs(ws.wq_lo_ppm / wq_true - 1.0) < 0.02 and abs(ws.wcsd_ppm / wcsd_true - 1.0) < 0.02,
           f"second_moment_split: W_q {ws.wq_lo_ppm:.3f} (truth {wq_true:.3f}), W_csd {ws.wcsd_ppm:.3f} "
           f"(truth {wcsd_true:.3f}) ppm, ok {ws.ok}")
    xs = np.array([1.0, 2.5, 3.0, 4.2, 6.0])
    es = np.array([0.05, 0.08, 0.05, 0.10, 0.06])
    ys = 0.7 * xs - 3.0 + rng.normal(0.0, 1.0, xs.size) * es
    b, a, cov = qf.weighted_line(xs, ys, es)
    ref, ref_cov = np.polyfit(xs, ys, 1, w=1.0 / es, cov="unscaled")
    _check(np.allclose([b, a], ref, rtol=1e-9, atol=0.0) and np.allclose(cov, ref_cov, rtol=1e-9, atol=0.0),
           f"weighted_line ({b}, {a}) / {cov.tolist()} vs numpy polyfit {ref} / {ref_cov.tolist()}")
    say(f"qcpmg_fields: CG {', '.join(f'{m:.3f}' for _t, m in dcgs)} ppm at 58.8/78.4/107.8 MHz (Eq. 1: "
        f"{', '.join(f'{t:.3f}' for t, _m in dcgs)}) -> diso {res.delta_iso_ppm:.3f} ppm (truth -70), CQ "
        f"{res.cq_MHz:.4f} MHz (3.2), PQ {res.pq_MHz:.4f} ({pq:.4f}), chi2/dof {res.chi2_red:.2g}; second-moment "
        f"split W_q {ws.wq_lo_ppm:.2f} / W_csd {ws.wcsd_ppm:.2f} ppm (truth {wq_true:.2f} / {wcsd_true:.2f}); "
        "weighted_line = numpy polyfit")


# ------------------------------------------------------------ shift calibration
def _sser_shiftcal(say, out: Path) -> None:
    from larmor import shiftcal
    from larmor.shiftcal import CalibrationPoint, ShiftCalibration

    sig = np.array([210.0, 262.5, 331.0, 395.2])
    sigma_ref = 560.0                                # delta = sigma_ref - sigma: slope -1
    exact = shiftcal.fit_calibration([CalibrationPoint(f"ref{i}", s, sigma_ref - s, 0.2)
                                      for i, s in enumerate(sig)], "19F")
    _check(abs(exact.slope + 1.0) < 1e-9 and abs(exact.intercept - sigma_ref) < 1e-6 and exact.rmse_ppm < 1e-9
           and exact.dof == 2 and exact.flags == [], f"exact sigma_ref line: slope {exact.slope:.12f}, intercept "
           f"{exact.intercept:.9f}, rmse {exact.rmse_ppm:.2g}, flags {exact.flags}")
    sr = shiftcal.from_sigma_ref(sigma_ref, "19F")
    _check(sr.delta(300.0) == 260.0 and sr.kind == "sigma_ref" and sr.zeta_scale == 1.0,
           f"from_sigma_ref: delta(300) = {sr.delta(300.0)}")
    # a general line delta = a sigma + b, weighted by the reference errors
    a, b = -0.85, 140.0
    err = np.array([0.2, 0.3, 0.2, 0.4])
    d = a * sig + b + np.random.default_rng(9).normal(0.0, 1.0, sig.size) * err
    cal = shiftcal.fit_calibration([CalibrationPoint(f"w{i}", s, v, e)
                                    for i, (s, v, e) in enumerate(zip(sig, d, err))], "19F")
    ref, ref_cov = np.polyfit(sig, d, 1, w=1.0 / err, cov="unscaled")
    scale = max(1.0, cal.chi2_red)                  # inflated by chi2_red > 1 only
    _check(np.allclose([cal.slope, cal.intercept], ref, rtol=1e-9, atol=0.0)
           and np.allclose(cal.cov, ref_cov * scale, rtol=1e-9, atol=0.0),
           f"weighted calibration ({cal.slope}, {cal.intercept}) cov {cal.cov} vs numpy {ref} / {ref_cov.tolist()}")
    _check(abs(cal.slope - a) < 4.0 * cal.slope_err and abs(cal.intercept - b) < 4.0 * cal.intercept_err,
           f"calibration a {cal.slope:.4f} +- {cal.slope_err:.4f} (truth {a}), b {cal.intercept:.2f} +- "
           f"{cal.intercept_err:.2f} (truth {b})")
    centre = float(np.average(sig, weights=1.0 / err ** 2))
    _check(cal.delta_err(centre) < cal.delta_err(sig.min() - 50.0) and cal.delta_err(centre) < cal.delta_err(sig.max() + 50.0),
           f"delta_err does not grow away from the weighted centroid: {cal.delta_err(centre)}")
    path = out / "19F.shiftcal.json"
    cal.save(path)
    back = ShiftCalibration.load(path)
    _check(back.slope == cal.slope and back.intercept == cal.intercept and back.cov == cal.cov
           and back.describe() == cal.describe() and [p.name for p in back.points] == ["w0", "w1", "w2", "w3"],
           f"shiftcal save/load round trip: {back.describe()}")
    two = shiftcal.fit_calibration([CalibrationPoint("NaF", 402.07, -225.0, 0.3),
                                    CalibrationPoint("CaF2", 232.92, -108.82, 0.3)], "19F")
    _check(two.dof == 0 and shiftcal.FLAG_ZERO_DOF in two.flags and abs(two.delta(402.07) + 225.0) < 1e-9,
           f"two-point calibration: dof {two.dof}, flags {two.flags}")
    say(f"shiftcal: exact delta = {sigma_ref:g} - sigma -> slope {exact.slope:.6f}, intercept {exact.intercept:.4f}; "
        f"weighted a {cal.slope:.4f} +- {cal.slope_err:.4f} (truth -0.85), b {cal.intercept:.2f} +- "
        f"{cal.intercept_err:.2f} (140), covariance = numpy polyfit; .shiftcal.json round trip; 2 points -> zero-dof flag")


# ------------------------------------------------------------ specialty: structure, 2D, physics remainder
def _specialty_struct(say, ctx):
    """DFT tensors, the SIMPSON bridge, the 2D MQMAS engine, the correlation
    decomposition and what engine-physics / engine-error-tools leave
    unchecked in the physics helpers -- each against an independent truth:
    tensors built here and CODATA constants for the .magres import (sigma_iso
    = trace/3, Haeberlen zeta/eta, C_Q = eQVzz/h, eta_Q from the
    eigenvalues), the analytic second-order shifts of a 27Al 3QMAS map (F2 =
    d_iso + d2, F1 = d_iso - 10/17 d2) and 2D fits that recover a position and
    an F1 referencing offset, the Herzfeld-Berger moments (Maricq-Waugh M2
    and the sign-carrying M3), the static-CT powder average and its eta = 0
    singularities, mrsimulator's Czjzek Monte-Carlo, IUPAC receptivities,
    the closed-form correlation of two co-centred Gaussians, and injected
    defects for the sanity checks."""
    _import_mrsimulator(say)
    out = _out(ctx, "specialty_struct")
    _sstr_dft(say, out)
    _sstr_simpson(say, out)
    _sstr_twod(say, ctx)
    _sstr_correlate(say)
    _sstr_convert_nuclei(say)
    _sstr_staticct(say)
    _sstr_herzfeld_berger(say)
    _sstr_czjzek(say)
    _sstr_estimate(say)
    _sstr_measure(say)
    _sstr_refranges(say)
    _sstr_sanity(say)
    _sstr_diagnostics(say)
    _sstr_identifiability(say)


#: CODATA 2018 (e and h exact): elementary charge, Planck constant, Hartree
#: energy, Bohr radius -- the first-principles EFG and cutoff conversions
_SSTR_E, _SSTR_H, _SSTR_EH, _SSTR_A0 = 1.602176634e-19, 6.62607015e-34, 4.3597447222071e-18, 5.29177210903e-11


def _sstr_rot(a: float, b: float, c: float) -> np.ndarray:
    """A ZYZ rotation (radians): puts a tensor in a known, non-trivial orientation."""
    def rz(t):
        return np.array([[np.cos(t), -np.sin(t), 0.0], [np.sin(t), np.cos(t), 0.0], [0.0, 0.0, 1.0]])
    ry = np.array([[np.cos(b), 0.0, np.sin(b)], [0.0, 1.0, 0.0], [-np.sin(b), 0.0, np.cos(b)]])
    return rz(a) @ ry @ rz(c)


def _sstr_dft(say, out: Path) -> None:
    """larmor.dft on a synthetic CASTEP-style .magres: two symmetry-equivalent
    27Al atoms (same principal values, two orientations, one shielding tensor
    with an antisymmetric part) and one 17O."""
    import hashlib
    import json

    from larmor import dft, engine, nuclei, sanity, shiftcal
    from larmor.recipe import Recipe

    v_au = _SSTR_EH / (_SSTR_E * _SSTR_A0 ** 2)          # V m^-2 per atomic unit of EFG
    ry_ev = _SSTR_EH / (2.0 * _SSTR_E)                   # 1 Ry in eV
    ms_pv, efg_pv = (540.0, 555.0, 615.0), (-0.08, -0.12, 0.20)      # sigma (ppm), V (a.u.): 27Al, C_Q ~7 MHz
    o_ms, o_efg = (210.0, 230.0, 290.0), (-0.30, -0.50, 0.80)        # 17O, C_Q ~ -4.8 MHz
    skew = np.array([[0.0, 3.0, -2.0], [-3.0, 0.0, 1.0], [2.0, -1.0, 0.0]])  # antisymmetric: no NMR effect
    r1, r2, r3 = _sstr_rot(0.3, 0.9, -0.4), _sstr_rot(1.7, 2.2, 0.6), _sstr_rot(-0.8, 0.5, 2.9)
    tens = {("Al", 1): (r1 @ np.diag(ms_pv) @ r1.T + skew, r1 @ np.diag(efg_pv) @ r1.T),
            ("Al", 2): (r2 @ np.diag(ms_pv) @ r2.T, r2 @ np.diag(efg_pv) @ r2.T),
            ("O", 1): (r3 @ np.diag(o_ms) @ r3.T, r3 @ np.diag(o_efg) @ r3.T)}
    lines = ["#$magres-abinitio-v1.0", "[calculation]", "calc_code CASTEP", "calc_code_version 23.1",
             "calc_xcfunctional PBE", "calc_cutoffenergy 816.34 eV", "[/calculation]", "[atoms]",
             "units lattice Angstrom", "units atom Angstrom"]
    lines += [f"atom {el} {el} {i} {0.9 * k:.3f} 0.0 0.0" for k, (el, i) in enumerate(tens)]
    lines += ["[/atoms]", "[magres]", "units ms ppm", "units efg au"]
    for (el, i), (ms, efg) in tens.items():
        lines.append(f"ms {el} {i} " + " ".join(f"{v:.12f}" for v in ms.ravel()))
        lines.append(f"efg {el} {i} " + " ".join(f"{v:.12f}" for v in efg.ravel()))
    path = out / "synthetic.magres"
    path.write_text("\n".join(lines + ["[/magres]", ""]), encoding="ascii")

    mf = dft.read_magres_file(path)
    labels = [s.label for s in mf.sites]
    _check(labels == ["Al1", "Al2", "O1"] and mf.sha256 == hashlib.sha256(path.read_bytes()).hexdigest(),
           f"read_magres_file: sites {labels}, sha256 {mf.sha256[:12]}")
    _check(mf.calc.code == "CASTEP" and mf.calc.xc == "PBE" and mf.calc.cutoff_wfc_Ry is not None
           and abs(mf.calc.cutoff_wfc_Ry - 816.34 / ry_ev) < 1e-6,
           f"[calculation]: {mf.calc.to_dict()} (816.34 eV = {816.34 / ry_ev:.6f} Ry)")
    warns = dft.assign_isotopes(mf.sites)
    _check(warns == [] and [s.isotope for s in mf.sites] == ["27Al", "27Al", "17O"],
           f"assign_isotopes: {[s.isotope for s in mf.sites]}, {warns}")
    al1, _al2, o1 = mf.sites
    # shielding: sigma_iso = trace / 3 (the antisymmetric part has none); Haeberlen
    # |zz - iso| >= |xx - iso| >= |yy - iso|, zeta = zz - iso, eta = (yy - xx) / zeta
    sh = al1.shielding()
    iso_true = float(np.mean(ms_pv))
    zz, xx, yy = sorted(np.array(ms_pv) - iso_true, key=abs, reverse=True)
    eta_cs = (yy - xx) / zz
    _check(abs(sh["iso_ppm"] - float(np.trace(tens[("Al", 1)][0])) / 3.0) < 1e-6
           and abs(sh["iso_ppm"] - iso_true) < 1e-6 and abs(sh["zeta_ppm"] - zz) < 1e-6 and abs(sh["eta"] - eta_cs) < 1e-6,
           f"Al1 shielding {sh}; built sigma_iso {iso_true}, zeta {zz}, eta {eta_cs:.6f}")
    # EFG: C_Q = e Q Vzz / h (signed), eta = (Vxx - Vyy) / Vzz with |Vxx| <= |Vyy| <= |Vzz|
    qtab = {i.symbol: i.quad_moment_barn for i in nuclei.all_isotopes() if i.symbol in ("27Al", "17O")}
    cq_al = qtab["27Al"] * 1e-28 * efg_pv[2] * v_au * _SSTR_E / _SSTR_H / 1e6
    cq_o = qtab["17O"] * 1e-28 * o_efg[2] * v_au * _SSTR_E / _SSTR_H / 1e6
    vxx, vyy, vzz = sorted(efg_pv, key=abs)
    eta_q = (vxx - vyy) / vzz
    qa, qo = al1.quadrupolar(), o1.quadrupolar()
    _check(qa is not None and abs(qa["Cq_MHz"] / cq_al - 1.0) < 1e-5 and abs(qa["eta"] - eta_q) < 1e-6,
           f"27Al EFG -> {qa}; e Q Vzz / h = {cq_al:.6f} MHz (Q {qtab['27Al']} b), eta {eta_q:.6f}")
    _check(qo is not None and abs(qo["Cq_MHz"] / cq_o - 1.0) < 1e-5 and qo["Cq_MHz"] < 0.0,
           f"17O EFG -> {qo}; e Q Vzz / h = {cq_o:.6f} MHz (Q {qtab['17O']} b < 0: the signed value is negative)")
    # the shipped Q against Pyykko 2018 (Mol. Phys. 116, 1328): compilations
    # differ by a few % (27Al: 0.140-0.150 b in common tables), so 5 %
    lit = {"27Al": 0.1466, "17O": -0.02558}
    dev = {k: qtab[k] / v - 1.0 for k, v in lit.items()}
    _check(all(abs(d) < 0.05 for d in dev.values()), f"quadrupole moments {qtab} vs Pyykko 2018 {lit}")
    say(f"dft: synthetic.magres (CASTEP, 816.34 eV = {mf.calc.cutoff_wfc_Ry:.4f} Ry, sha256 ok) -> {labels}; "
        f"Al1 sigma_iso {sh['iso_ppm']:.4f} ppm = trace/3, zeta {sh['zeta_ppm']:.4f}, eta {sh['eta']:.4f} "
        f"(built 570 / {zz:.0f} / {eta_cs:.4f})")
    say(f"dft: EFG -> C_Q = eQVzz/h: 27Al {qa['Cq_MHz']:.5f} MHz (CODATA {cq_al:.5f}), eta {qa['eta']:.4f} "
        f"(built {eta_q:.4f}); 17O {qo['Cq_MHz']:.5f} MHz (CODATA {cq_o:.5f}, Q < 0); table Q vs Pyykko 2018 "
        f"{dev['27Al']:+.1%} / {dev['17O']:+.1%}")
    groups = dft.group_equivalent(dft.sites_for_isotope(mf.sites, "27Al"))
    g = groups[0] if groups else None
    _check(len(groups) == 1 and g.multiplicity == 2 and g.members == ["Al1", "Al2"]
           and abs(g.shielding()["zeta_ppm"] - zz) < 1e-6 and abs(g.quadrupolar()["Cq_MHz"] - cq_al) < 1e-4,
           f"group_equivalent: {[(x.label, x.multiplicity, x.members) for x in groups]}")
    # seeding: the lineshape depends on |C_Q| only and the models bound Cq_MHz
    # at >= 0.01, so a negative computed C_Q (here the 17O site: Q < 0) must
    # seed its magnitude -- a signed seed would be clamped to 0.01 MHz
    o_seed = o1.to_site_dict("quad_ct", reference_ppm=290.0)["params"]["Cq_MHz"]["value"]
    _check(abs(o_seed - abs(qo["Cq_MHz"])) < 1e-9 * abs(qo["Cq_MHz"]),
           f"to_site_dict seeds the 17O site (computed C_Q {qo['Cq_MHz']:.5f} MHz) at {o_seed} MHz, not |C_Q|")
    sigma_ref = 630.0
    p = g.to_site_dict("quad_ct", reference_ppm=sigma_ref)["params"]
    _check(abs(p["isotropic_chemical_shift_ppm"]["value"] - (sigma_ref - iso_true)) < 1e-6
           and abs(p["Cq_MHz"]["value"] - cq_al) < 1e-4 and abs(p["eta"]["value"] - eta_q) < 1e-6,
           f"to_site_dict(quad_ct, sigma_ref {sigma_ref}): "
           f"{ {k: v['value'] for k, v in p.items()} }")
    # the path the import dialog and the CLI share: the site dict must load,
    # pass the physical checks and simulate where its tensor says -- the CT
    # centre of gravity of a powder (static here) is d_iso + d2(P_Q)
    dicts, notes = dft.sites_to_recipe_dicts(groups, "quad_ct", reference_ppm=sigma_ref)
    rec = Recipe.from_dict({"nucleus": "27Al", "larmor_frequency_MHz": 130.3, "spin_rate_Hz": 0.0, "sites": dicts})
    flags = sanity.check_recipe(rec)
    gx, _y, per = engine.simulate(rec, exp_ppm=np.linspace(250.0, -150.0, 4001))
    ct = np.clip(np.asarray(per[0], float), 0.0, None)
    centroid = float((np.asarray(gx) * ct).sum() / ct.sum())
    pq = cq_al * np.sqrt(1.0 + eta_q ** 2 / 3.0)
    d2 = -(3.0 / 40.0) * (2.5 * 3.5 - 0.75) / (2.5 ** 2 * 4.0 ** 2) * (pq / 130.3) ** 2 * 1e6
    expect = sigma_ref - iso_true + d2
    _check(len(dicts) == 1 and set(dicts[0]) == {"model", "label", "params"} and any("×2 (Al1, Al2)" in n for n in notes)
           and flags == [] and abs(centroid - expect) < 0.5,
           f"sites_to_recipe_dicts: {len(dicts)} site(s), notes {notes}, sanity {flags}; static CT centroid "
           f"{centroid:.2f} ppm vs d_iso + d2 = {expect:.2f}")
    prov = dft.import_provenance(mf, "27Al", "quad_ct", groups, shiftcal.from_sigma_ref(sigma_ref, "27Al"))
    ps = prov["sites"][0]
    _check(prov["sha256"] == mf.sha256 and ps["multiplicity"] == 2 and ps["members"] == ["Al1", "Al2"]
           and abs(ps["delta_pred_ppm"] - (sigma_ref - iso_true)) < 1e-6 and abs(ps["Cq_MHz"] - cq_al) < 1e-4
           and prov["calculation"]["code"] == "CASTEP" and json.loads(json.dumps(prov)) == prov,
           f"import_provenance: {prov}")
    say(f"dft: Al1/Al2 (two orientations) grouped x{g.multiplicity}, zeta kept {g.shielding()['zeta_ppm']:.3f} ppm; "
        f"quad_ct seed d_iso {p['isotropic_chemical_shift_ppm']['value']:.3f} ppm (= {sigma_ref:g} - {iso_true:g}), "
        f"C_Q {p['Cq_MHz']['value']:.4f} MHz, eta {p['eta']['value']:.4f}; the imported recipe is physical and its "
        f"static CT centroid {centroid:.2f} ppm = d_iso + d2 {expect:.2f}; provenance record JSON-clean")


def _sstr_simpson(say, out: Path) -> None:
    """larmor.simpson: the spin-system / REDOR input text, the SIMP reader on
    a known FID, and the detection (absent here: run() must refuse cleanly)."""
    import os

    from larmor import simpson
    from larmor.recipe import Param, Recipe, SiteModel

    rec = Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, spin_rate_Hz=20000.0, sites=[
        SiteModel(model="quad_ct", label="A", params={
            "isotropic_chemical_shift_ppm": Param(60.0), "Cq_MHz": Param(4.2), "eta": Param(0.35),
            "shift_fwhm_ppm": Param(2.0), "amplitude": Param(1.0)})])
    block = simpson.spinsys_block(rec)
    rows = {ln.split()[0]: ln.split() for ln in block.splitlines()[1:-1] if ln.strip()}
    quad, shift = rows.get("quadrupole"), rows.get("shift")
    _check(quad is not None and abs(float(quad[3]) - 4.2e6) < 1e-3 and abs(float(quad[4]) - 0.35) < 1e-12
           and shift is not None and shift[2] == "60.0p" and rows.get("channels") == ["channels", "27Al"],
           f"spinsys_block: C_Q must be written in Hz (4.2 MHz -> 4.2e6):\n{block}")
    text = simpson.redor_input(rec, "19F", -1234.0, spin_rate_hz=12500.0)
    tr = next((float(ln.split()[2]) for ln in text.splitlines() if ln.strip().startswith("variable tr")), None)
    _check(tr is not None and abs(tr - 1.0 / 12500.0) < 1e-15 and "dipole 1 2 -1234" in text
           and "nuclei 27Al 19F" in text, f"redor_input: rotor period {tr} (expected {1.0 / 12500.0})")
    n, sw, f0 = 64, 2000.0, 250.0
    t = np.arange(n) / sw
    fid = np.exp(2j * np.pi * f0 * t) * np.exp(-t / 0.01)
    fpath = out / "sim.fid"
    fpath.write_text(f"SIMP\nNP={n}\nSW={sw:g}\nTYPE=FID\nDATA\n"
                     + "".join(f"{v.real:.12f} {v.imag:.12f}\n" for v in fid) + "END\n", encoding="ascii")
    res = simpson.parse_fid(fpath)
    peak = float(np.fft.fftfreq(n, d=float(res.x[1] - res.x[0]))[int(np.argmax(np.abs(np.fft.fft(res.y))))])
    _check(res.y.size == n and np.allclose(res.y, fid, atol=1e-9) and abs(res.x[1] - 1.0 / sw) < 1e-15
           and abs(peak - f0) < 1e-6, f"parse_fid: {res.y.size} points, dwell {res.x[1]}, line at {peak} Hz (wrote {f0})")
    exe = simpson.simpson_available()
    if exe is None:
        try:
            simpson.run("spinsys {\n}\n")
        except RuntimeError as exc:
            _check("not on PATH" in str(exc), f"simpson.run without SIMPSON: unclear message {exc}")
        else:
            raise AssertionError("simpson.run returned a result with no SIMPSON on PATH")
        found = "not on PATH -> run() refuses with a clear RuntimeError, nothing written"
    else:
        _check(Path(exe).is_file() and os.access(exe, os.X_OK), f"simpson_available() -> {exe!r} is not an executable")
        found = f"found at {exe} (not run here: the recoupling block is lab-specific)"
    say(f"simpson: spinsys C_Q {float(quad[3]):.4g} Hz, REDOR rotor period {tr * 1e6:.1f} us, SIMP reader "
        f"{res.y.size} points with the {peak:.0f} Hz line in place; SIMPSON {found}")


def _sstr_mq_recipe(nu0: float, diso, cq: float, eta: float, fwhm: float, amp: float = 1.0, fixed: bool = False):
    """One 27Al quad_ct site for the 2D checks: C_Q, eta and the width held (the
    2D basis is discrete in C_Q), d_iso free in [0, 60] unless ``fixed``."""
    from larmor.recipe import Param, Recipe, SiteModel

    pos = Param(diso, min=0.0, max=60.0) if not fixed else Param(diso, vary=False)
    return Recipe(nucleus="27Al", larmor_frequency_MHz=nu0, spin_rate_Hz=0.0, sites=[
        SiteModel(model="quad_ct", label="A", params={
            "isotropic_chemical_shift_ppm": pos, "Cq_MHz": Param(cq, vary=False), "eta": Param(eta, vary=False),
            "shift_fwhm_ppm": Param(fwhm, vary=False), "amplitude": Param(amp, min=0.0)})])


def _sstr_twod(say, ctx) -> None:
    """larmor.twod: a tiny 27Al 3QMAS kernel, the site positions against the
    analytic second-order relations, two few-parameter fits, the 2D
    processing replay, and the bundled 3620 3QMAS read through twod."""
    import time

    from larmor import convert, czjzek_dist, twod
    from larmor.recipe import Param, Recipe, SiteModel

    nu0, spin = 130.3, 2.5
    try:
        t0 = time.perf_counter()
        k = twod.build_mqmas_kernel("27Al", nu0, f2_window=(110.0, -70.0), f1_window=(110.0, -30.0),
                                    n2=64, n1=48, n_cq=3, n_eta=2, cq_max_MHz=8.0)
        t_k = time.perf_counter() - t0
        # sheared 3QMAS, I = 5/2: echo ratio k = 19/12, so a pure shift lands in
        # the native F1 at (k - 3) / (1 + k) = -17/31 of its value
        ratio = 19.0 / 12.0
        c_true = (ratio - 3.0) / (1.0 + ratio)
        _check(k.K.shape == (6, 48, 64) and abs(k.f1_cs_scale - c_true) < 5e-3
               and np.all(np.diff(k.f1_ppm) > 0) and np.all(np.diff(k.f2_ppm) > 0),
               f"MQMAS kernel {k.K.shape}, F1 shift scale {k.f1_cs_scale:.5f} (3QMAS I=5/2: {c_true:.5f})")
        # second-order isotropic shifts: F2 centroid = d_iso + d2 (the CT MAS
        # centre of gravity), F1 = d_iso - (10/17) d2 in the d1-isotropic
        # convention, d2 = -(3/40) [I(I+1) - 3/4] / [I^2 (2I-1)^2] (P_Q / nu0)^2;
        # the measured agreement is ~0.05 ppm on this 2.8 ppm grid, 0.5 ppm
        # leaves margin while a wrong shear, slope or P_Q factor moves >= 2 ppm
        pos = []
        for cq, eta in ((8.0, 0.0), (8.0, 1.0)):            # both on the kernel's (C_Q, eta) grid
            z, _per = twod.simulate_2d(_sstr_mq_recipe(nu0, 30.0, cq, eta, 3.0), k)
            zc = np.clip(z, 0.0, None)
            p2, p1 = zc.sum(axis=0), zc.sum(axis=1)
            f2c, f1c = float((p2 * k.f2_ppm).sum() / p2.sum()), float((p1 * k.f1_ppm).sum() / p1.sum())
            pq = cq * np.sqrt(1.0 + eta ** 2 / 3.0)
            d2 = -(3.0 / 40.0) * (spin * (spin + 1.0) - 0.75) / (spin ** 2 * (2.0 * spin - 1.0) ** 2) * (pq / nu0) ** 2 * 1e6
            got = convert.ct_second_order_shift_ppm(convert.pq_from_cq_eta(cq, eta), spin, nu0)
            _check(abs(got - d2) < 1e-9 * abs(d2), f"convert.ct_second_order_shift_ppm {got} vs {d2}")
            _check(abs(f2c - (30.0 + d2)) < 0.5 and abs(f1c - (30.0 - 10.0 / 17.0 * d2)) < 0.5,
                   f"3QMAS C_Q {cq} eta {eta}: centroid F2 {f2c:.2f} / F1 {f1c:.2f} ppm, analytic "
                   f"{30.0 + d2:.2f} / {30.0 - 10.0 / 17.0 * d2:.2f}")
            pos.append(f"C_Q {cq:g} eta {eta:g}: F2 {f2c:.2f} / F1 {f1c:.2f} (analytic {30.0 + d2:.2f} / "
                       f"{30.0 - 10.0 / 17.0 * d2:.2f})")
        slope = twod.qis_slope("27Al", nu0)
        _check(abs(slope + 10.0 / 17.0) < 0.02, f"qis_slope {slope:.4f}, expected -10/17 = {-10.0 / 17.0:.4f}")
        say(f"twod: 3QMAS kernel {k.K.shape} in {t_k:.2f} s, F1 shift scale {k.f1_cs_scale:.5f} (-17/31); "
            f"d_iso 30 ppm sites, {'; '.join(pos)}; QIS slope {slope:.3f} (-10/17)")
        # a Czjzek site is the kernel's weights times its (C_Q, eta) basis: its
        # centroids must be the WEIGHTED second-order shifts of the grid points
        # in the kernel's documented eta-major order (a weight/basis mismatch
        # moves them), and the weights must be czjzek_dist's d = 5 density
        sigma = 1.6
        w = k.weights(sigma)
        cq_g, eta_g = np.meshgrid(k.cq_grid_MHz, k.eta_grid, indexing="xy")
        coef = -(3.0 / 40.0) * (spin * (spin + 1.0) - 0.75) / (spin ** 2 * (2.0 * spin - 1.0) ** 2) * 1e6 / nu0 ** 2
        d2w = float((w * coef * (cq_g ** 2 * (1.0 + eta_g ** 2 / 3.0)).ravel()).sum() / w.sum())
        cz = Recipe(nucleus="27Al", larmor_frequency_MHz=nu0, sites=[SiteModel(model="czjzek", label="G", params={
            "isotropic_chemical_shift_ppm": Param(30.0), "sigma_Cq_MHz": Param(sigma), "shift_fwhm_ppm": Param(3.0),
            "line_fwhm_ppm": Param(0.0), "amplitude": Param(1.0)})])
        zc = np.clip(twod.simulate_2d(cz, k)[0], 0.0, None)
        p2, p1 = zc.sum(axis=0), zc.sum(axis=1)
        f2c, f1c = float((p2 * k.f2_ppm).sum() / p2.sum()), float((p1 * k.f1_ppm).sum() / p1.sum())
        wdev = float(np.abs(w / w.sum() - czjzek_dist.czjzek_weights(sigma, 5.0, k.cq_grid_MHz, k.eta_grid)).max())
        _check(wdev < 1e-9 and abs(f2c - (30.0 + d2w)) < 0.5 and abs(f1c - (30.0 - 10.0 / 17.0 * d2w)) < 0.5,
               f"2D Czjzek sigma {sigma}: centroid F2 {f2c:.2f} / F1 {f1c:.2f}, weighted analytic {30.0 + d2w:.2f} / "
               f"{30.0 - 10.0 / 17.0 * d2w:.2f}; weights vs czjzek_dist {wdev:.1e}")
        say(f"twod Czjzek (sigma {sigma} MHz): weights = czjzek_dist d=5 to {wdev:.0e}; centroid F2 {f2c:.2f} / "
            f"F1 {f1c:.2f} ppm = the weighted second-order shifts {30.0 + d2w:.2f} / {30.0 - 10.0 / 17.0 * d2w:.2f}")
        # fit 1: the experiment referenced like the kernel (beta held at 0),
        # d_iso started 5 ppm off; fit 2: the experiment's F1 axis 6 ppm off,
        # the site held, the F1 referencing offset beta fitted
        truth = _sstr_mq_recipe(nu0, 30.0, 8.0, 0.0, 3.0)
        z, _ = twod.simulate_2d(truth, k)
        noisy = z + np.random.default_rng(12).normal(0.0, 0.005 * float(z.max()), z.shape)
        t0 = time.perf_counter()
        start = _sstr_mq_recipe(nu0, 25.0, 8.0, 0.0, 3.0, amp=0.5)
        start.mqmas_f1_ref_ppm, start.mqmas_f1_ref_vary = 0.0, False
        res = twod.fit_2d(start, twod.Data2D(f2_ppm=k.f2_ppm, f1_ppm=k.f1_ppm, z=noisy, nucleus="27Al",
                                             larmor_MHz=nu0), kernel=k)
        p_iso = start.sites[0].params["isotropic_chemical_shift_ppm"]
        _check(abs(p_iso.value - 30.0) < 0.3 and res.rmsd < 0.01 and start.mqmas_f1_ref_ppm == 0.0
               and p_iso.stderr is not None and res.z_fit.shape == k.shape,
               f"fit_2d from d_iso 25: {p_iso.value:.3f} ppm (truth 30), rmsd {res.rmsd:.4f}, beta {start.mqmas_f1_ref_ppm}")
        held = _sstr_mq_recipe(nu0, 30.0, 8.0, 0.0, 3.0, amp=0.5, fixed=True)
        res2 = twod.fit_2d(held, twod.Data2D(f2_ppm=k.f2_ppm, f1_ppm=k.f1_ppm - 6.0, z=noisy, nucleus="27Al",
                                             larmor_MHz=nu0), kernel=k)
        _check(abs(held.mqmas_f1_ref_ppm + 6.0) < 0.5 and res2.rmsd < 0.01,
               f"fit_2d F1 referencing: beta {held.mqmas_f1_ref_ppm:.3f} ppm (truth -6), rmsd {res2.rmsd:.4f}")
        say(f"twod.fit_2d: d_iso 25 -> {p_iso.value:.3f} +- {p_iso.stderr:.3f} ppm (truth 30, beta held), rmsd "
            f"{res.rmsd:.4f}; F1 referenced 6 ppm off -> beta {held.mqmas_f1_ref_ppm:.3f} ppm, rmsd {res2.rmsd:.4f}; "
            f"{time.perf_counter() - t0:.2f} s")
        # 2D processing replayed from a recorded op list: exact hypercomplex
        # phasing (a 50 deg rotation of absorption/dispersion undone), the
        # calibration shift, the shear F1' = F1 + factor (F2 - ref), and the
        # Hilbert path (real quadrant only): the live row preview must equal the
        # applied 2D phasing
        f2 = np.linspace(-100.0, 100.0, 201)                # 1 ppm steps
        f1 = np.linspace(-40.0, 40.0, 41)                   # 2 ppm steps
        xl = (f2 - 10.0) / 2.0
        absn = np.tile(1.0 / (1.0 + xl ** 2), (f1.size, 1))
        disp = np.tile(xl / (1.0 + xl ** 2), (f1.size, 1))
        th = np.deg2rad(50.0)
        d = twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=absn * np.cos(th) - disp * np.sin(th),
                        ri=absn * np.sin(th) + disp * np.cos(th), ir=np.zeros_like(absn), ii=np.zeros_like(absn))
        done = twod.replay_ops(d, [{"op": "phase", "axis": "f2", "p0": 50.0, "p1": 0.0, "pivot": 10.0},
                                   {"op": "shift", "d2": 1.5, "d1": -2.0}])
        dot = np.zeros((f1.size, f2.size))
        j26, i0 = int(np.argmin(np.abs(f2 - 26.0))), int(np.argmin(np.abs(f1)))
        dot[i0, j26] = 1.0                                  # a point at F2 = 26, F1 = 0
        sh = twod.replay_ops(twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=dot), [{"op": "shear", "factor": 1.0, "ref_ppm": 0.0}])
        moved = float(f1[int(np.argmax(sh.z[:, j26]))])
        try:
            twod.replay_ops(d, [{"op": "warp"}])
            refused = False
        except ValueError as exc:
            refused = "warp" in str(exc)
        real = twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=d.z)
        applied = real.phased("f2", 40.0, 25.0, pivot_ppm=5.0)
        preview = real.phase_line("f2", 7, 40.0, 25.0, pivot_ppm=5.0)
        perr = float(np.abs(done.z - absn).max())
        _check(perr < 1e-9 and np.allclose(done.f2_ppm, f2 + 1.5) and np.allclose(done.f1_ppm, f1 - 2.0)
               and abs(moved - 26.0) < 1e-9 and refused and np.abs(applied.z[7] - preview).max() < 1e-12
               and not np.allclose(applied.z, real.z),
               f"replay_ops: phase error {perr:.2e}, shear moved the point to F1 {moved:.2f} (expected 26), unknown op "
               f"refused {refused}, Hilbert preview vs applied {np.abs(applied.z[7] - preview).max():.1e}")
        say(f"twod 2D processing: hypercomplex p0 50 deg recovers the absorption to {perr:.1e}, shift (+1.5, -2.0) "
            f"relabels the axes, shear x1 moves an F2 = 26 ppm point to F1 {moved:.1f}, an unknown op is refused; "
            f"Hilbert-path row preview = applied phasing")
        e3 = _expno(ctx, 3620)
        if e3 is None:
            say("no bundled example data: twod.read_bruker_2d on the 3620 3QMAS is skipped "
                "(the synthetic maps above cover the 2D engine)")
            return
        t0 = time.perf_counter()
        m = twod.read_bruker_2d(e3)
        t_r = time.perf_counter() - t0
        _check(m.z.shape == (64, 1024) and m.has_hyper and m.nucleus == "27Al" and abs(m.larmor_MHz - 130.323) < 0.01
               and m.spin_rate_Hz == 26000.0 and np.all(np.diff(m.f2_ppm) > 0) and np.all(np.diff(m.f1_ppm) > 0),
               f"read_bruker_2d(3620): {m.z.shape}, hyper {m.has_hyper}, {m.nucleus} {m.larmor_MHz} MHz, "
               f"MAS {m.spin_rate_Hz}")
        back = m.phased("f2", 30.0, 12.0, pivot_ppm=60.0).phased("f2", -30.0, -12.0, pivot_ppm=60.0)
        rot = m.phased("f1", 47.0, 0.0)
        norm0 = sum(float(np.sum(np.square(a))) for a in (m.z, m.ri, m.ir, m.ii))
        norm1 = sum(float(np.sum(np.square(a))) for a in (rot.z, rot.ri, rot.ir, rot.ii))
        err = float(np.abs(back.z - m.z).max() / np.abs(m.z).max())
        _check(err < 1e-9 and abs(norm1 / norm0 - 1.0) < 1e-9,
               f"3620 phasing: round trip error {err:.2e}, quadrature norm ratio {norm1 / norm0:.12f}")
        i1, i2 = np.unravel_index(int(np.argmax(m.z)), m.z.shape)
        say(f"twod.read_bruker_2d(3620) in {t_r:.2f} s: {m.z.shape} hypercomplex {m.nucleus} map at "
            f"{m.larmor_MHz:.3f} MHz, MAS {m.spin_rate_Hz:.0f} Hz, F2 {m.f2_ppm[0]:.0f}..{m.f2_ppm[-1]:.0f} / "
            f"F1 {m.f1_ppm[0]:.0f}..{m.f1_ppm[-1]:.0f} ppm, maximum at {m.f2_ppm[i2]:.1f} / {m.f1_ppm[i1]:.1f} ppm; "
            f"exact phase round trip ({err:.0e}), F1 rotation keeps the quadrature norm")
    finally:
        twod.clear_kernel_cache()


def _sstr_correlate(say) -> None:
    """larmor.correlate on a known decomposition: target = 2 A + 0.8 C, a
    reference = A + D on its own (descending) axis."""
    from larmor import correlate

    g = lambda x, c: np.exp(-4.0 * np.log(2.0) * ((x - c) / 6.0) ** 2)  # noqa: E731
    x1 = np.linspace(-60.0, 60.0, 1201)
    x2 = np.linspace(80.0, -40.0, 601)                     # descending, Bruker order
    target = correlate.Observable("1D", x1, 2.0 * g(x1, 10.0) + 0.8 * g(x1, -25.0))
    ref = correlate.Observable("projection", x2, g(x2, 10.0) + g(x2, 40.0), kind="projection")
    gu = correlate.common_grid([target, ref])
    gi = correlate.common_grid([target, ref], mode="intersection")
    _check(gu[0] == -60.0 and gu[-1] == 80.0 and gu.size == 1201 and gi[0] == -40.0 and gi[-1] == 60.0,
           f"common_grid: union {gu[0]}..{gu[-1]} ({gu.size}), intersection {gi[0]}..{gi[-1]}")
    at = lambda grid, amp, p: float(amp[int(np.argmin(np.abs(grid - p)))])  # noqa: E731
    grid, diff = correlate.difference(target, [ref], region=(0.0, 20.0))
    mask = (grid >= 0.0) & (grid <= 20.0)
    s_reg = correlate.scale_to(target.on(grid), ref.on(grid), mask)
    s_all = correlate.scale_to(target.on(grid), ref.on(grid))
    # inside (0, 20) ppm only the shared line lives: the least-squares scale is
    # its true ratio 2, and what is left is 0.8 C - 2 D; over everything the
    # reference's own line dilutes it to <t, r> / <r, r> = 2|A|^2 / 2|A|^2 = 1
    # (interpolation onto the union grid costs < 1e-3)
    _check(abs(s_reg - 2.0) < 0.01 and abs(at(grid, diff, 10.0)) < 0.02 and abs(at(grid, diff, -25.0) - 0.8) < 0.02
           and abs(at(grid, diff, 40.0) + 2.0) < 0.02 and abs(s_all - 1.0) < 0.01
           and correlate.scale_to(target.on(grid), -ref.on(grid)) == 0.0,
           f"difference: scale {s_reg:.4f} in the region (truth 2), {s_all:.4f} overall (truth 1); residue at "
           f"10 / -25 / 40 ppm = {at(grid, diff, 10.0):.3f} / {at(grid, diff, -25.0):.3f} / {at(grid, diff, 40.0):.3f}")
    gx, inter = correlate.intersection([target, ref])
    err = float(np.abs(inter - g(gx, 10.0)).max())
    _check(err < 0.01, f"intersection differs from the shared line by {err:.4f}")
    say(f"correlate: least-squares scale {s_reg:.4f} in the shared region (truth 2), {s_all:.4f} overall (truth 1), "
        f"un-correlated residue keeps C (0.8 -> {at(grid, diff, -25.0):.3f}) and removes A ({at(grid, diff, 10.0):+.4f}); "
        f"intersection = the shared line to {err:.1e}; descending reference axis handled")


def _sstr_convert_nuclei(say) -> None:
    """larmor.convert / larmor.nuclei beyond engine-physics: the EFG constant
    from CODATA, nu_Q for three spins, IUPAC receptivities, gamma(1H), B0, and
    the NMR table's feasibility partition."""
    from larmor import convert, nuclei

    k_au = _SSTR_E * (_SSTR_EH / (_SSTR_E * _SSTR_A0 ** 2)) * 1e-28 / _SSTR_H / 1e6   # MHz per barn x a.u.
    got = convert.cq_from_efg(1.0, 1.0)
    _check(abs(got / k_au - 1.0) < 1e-5, f"cq_from_efg(1 a.u., 1 b) = {got} MHz, CODATA e(Eh/e a0^2)(1e-28)/h = {k_au:.6f}")
    # nu_Q = 3 C_Q / [2I(2I - 1)]: C_Q/2 at I = 3/2, 3C_Q/20 at 5/2, C_Q/14 at 7/2
    for spin, frac in ((1.5, 0.5), (2.5, 0.15), (3.5, 1.0 / 14.0)):
        nq = convert.nu_q(2.6, spin)
        _check(abs(nq - 2.6 * frac) < 1e-12 and abs(convert.cq_from_nu_q(nq, spin) - 2.6) < 1e-12,
               f"nu_q(2.6 MHz, I={spin}) = {nq} (expected {2.6 * frac})")
    iso = {i.symbol: i for i in nuclei.all_isotopes()}
    # receptivity relative to 1H, IUPAC 2001 (Harris et al., Pure Appl. Chem. 73, 1795)
    rec = {}
    for sym, ref in (("27Al", 0.207), ("29Si", 3.68e-4), ("23Na", 9.27e-2)):
        rec[sym] = iso[sym].receptivity_1H
        _check(abs(rec[sym] / ref - 1.0) < 0.02, f"{sym} receptivity {rec[sym]:.4g}, IUPAC {ref:.4g}")
    b0 = nuclei.b0_from_1H(600.0)
    _check(abs(nuclei.GAMMA_1H - 42.577478518) < 1e-4 and abs(iso["1H"].larmor_MHz(b0) - 600.0) < 1e-9,
           f"gamma(1H) {nuclei.GAMMA_1H} MHz/T (CODATA 42.577478518), 600 MHz -> {b0:.4f} T")
    grid = {s for row in nuclei.PERIODIC_ROWS for s in row if s}
    seen = [el for els in nuclei.FEASIBILITY.values() for el in els]
    spots = {"Si": "favorable", "O": "challenging", "Mg": "very_difficult", "Fe": "impractical", "Ar": "impossible"}
    _check(len(seen) == len(set(seen)) and set(seen) == grid and nuclei.feasibility("Xx") is None
           and all(nuclei.feasibility(el) == cat for el, cat in spots.items()),
           f"feasibility: {len(seen)} entries for {len(grid)} table cells, spots "
           f"{ {el: nuclei.feasibility(el) for el in spots} }")
    workable = [el for cat in ("favorable", "challenging", "very_difficult") for el in nuclei.FEASIBILITY[cat]]
    missing = [el for el in workable if nuclei.primary_isotope(el) is None]
    _check(not missing, f"elements rated workable with no NMR-active isotope in the shipped table: {missing}")
    say(f"convert/nuclei: EFG constant {got:.4f} MHz/(b a.u.) (CODATA {k_au:.4f}); nu_Q(2.6 MHz) "
        f"{convert.nu_q(2.6, 1.5):.3f} / {convert.nu_q(2.6, 2.5):.3f} / {convert.nu_q(2.6, 3.5):.4f} MHz at I 3/2 5/2 7/2; "
        f"receptivity 27Al {rec['27Al']:.4f}, 29Si {rec['29Si']:.3e}, 23Na {rec['23Na']:.4f} (IUPAC); 600 MHz = "
        f"{b0:.4f} T; feasibility partitions {len(grid)} elements, {len(workable)} workable all with an isotope")


def _sstr_staticct(say) -> None:
    """larmor.staticct beyond engine-physics: the frequency surface's powder
    average against convert's d2, the eta = 0 singularities in closed form,
    and the refusal of an edge that sits on a horn."""
    from larmor import convert, staticct

    n = 400
    u = (np.arange(n) + 0.5) / n                         # cos(theta) midpoints: an area-faithful powder
    ph = (np.arange(n) + 0.5) / n * (np.pi / 2.0)       # phi midpoints over the symmetry quadrant
    th, phi = np.meshgrid(np.arccos(u), ph, indexing="ij")
    avg = []
    for cq, eta, spin, lar in ((4.0, 0.3, 2.5, 130.3), (30.0, 0.7, 1.5, 216.0)):
        mean = float(staticct.shift_ppm(th, phi, cq, eta, spin, lar).mean())
        want = convert.ct_second_order_shift_ppm(convert.pq_from_cq_eta(cq, eta), spin, lar)
        _check(abs(mean / want - 1.0) < 1e-4, f"static CT powder average {mean:.6f} ppm vs d2 {want:.6f} "
                                              f"(C_Q {cq}, eta {eta}, I {spin})")
        avg.append(f"{mean:.4f} / {want:.4f}")
    # eta = 0: nu = -(nu_Q^2 [I(I+1) - 3/4] / 16 nu0) (1 - x)(9x - 1), x = cos^2 theta,
    # spans -(16/9) to +1 in those units with both singularities at the limits
    cq, spin, lar = 10.0, 2.5, 130.3
    f = staticct.ct_static_features(cq, 0.0, spin, lar)
    nuq = 3.0 * cq * 1e6 / (2.0 * spin * (2.0 * spin - 1.0))
    a = nuq ** 2 * (spin * (spin + 1.0) - 0.75) / (lar * 1e6) / lar        # ppm
    lo, hi = -a / 9.0, a / 16.0
    span = hi - lo
    # the edges are the exact support; the horns are smoothed-histogram maxima
    # just inside a divergence at the support (measured 1.5 % of the span in)
    _check(abs(f.edges[0] - lo) < 1e-3 * span and abs(f.edges[1] - hi) < 1e-3 * span and len(f.horns) == 2
           and all(abs(hz - ref) < 0.03 * span for hz, ref in zip(f.horns, (lo, hi))),
           f"eta = 0 static CT: edges {f.edges}, horns {f.horns}; closed form {lo:.3f} / {hi:.3f} ppm")
    bad = staticct.read_cq_eta(f.horns[0], f.horns[1], f.edges[0], spin, lar)
    _check(not bad.ok and "opposite" in bad.message, f"read_cq_eta accepted an edge sitting on a horn: {bad}")
    say(f"staticct: powder average of the frequency surface = d2 ({'; '.join(avg)} ppm); eta = 0 edges "
        f"{f.edges[0]:.2f} / {f.edges[1]:.2f} ppm = closed form {lo:.2f} / {hi:.2f}, horns {f.horns[0]:.2f} / "
        f"{f.horns[1]:.2f}; a horn-coincident edge is refused")


def _sstr_herzfeld_berger(say) -> None:
    """larmor.herzfeld_berger beyond engine-physics: the moments of the
    sideband manifold against the static powder (M1 = 0, Maricq-Waugh M2,
    and M3 whose sign pins the shielding convention)."""
    from larmor import herzfeld_berger as hb

    _check(abs(3.0 * np.cos(hb.MAGIC_RAD) ** 2 - 1.0) < 1e-12, f"MAGIC_RAD {hb.MAGIC_RAD} is not the zero of P2")
    nu0, nur = 162.0, 3000.0
    step = nur / nu0
    parts = []
    for zeta, eta in ((60.0, 0.4), (-45.0, 0.0)):
        inten = hb.sideband_intensities(zeta, eta, nu0, nur)
        need = hb.n_sidebands_needed(zeta, nu0, nur)
        mom = [sum(v * (n * step) ** k for n, v in inten.items()) for k in range(4)]
        # static shift anisotropy is -zeta: <d^2> = zeta^2 (1 + eta^2/3) / 5,
        # <d^3> = (2/35)(-zeta)^3 (1 - eta^2); the moments of the gamma-averaged
        # sideband spectrum equal the static ones
        m2 = zeta ** 2 * (1.0 + eta ** 2 / 3.0) / 5.0
        m3 = -(2.0 / 35.0) * zeta ** 3 * (1.0 - eta ** 2)
        _check(max(abs(n) for n in inten) == need and mom[0] > 0.9999 and abs(mom[1]) < 1e-3 * np.sqrt(m2)
               and abs(mom[2] / m2 - 1.0) < 1e-3 and abs(mom[3] / m3 - 1.0) < 1e-2,
               f"sidebands zeta {zeta} eta {eta}: orders +-{max(abs(n) for n in inten)} (need {need}), "
               f"sum {mom[0]:.6f}, M1 {mom[1]:.2e}, M2 {mom[2]:.3f} vs {m2:.3f}, M3 {mom[3]:.1f} vs {m3:.1f} ppm^3")
        parts.append(f"zeta {zeta:+g}: M2 {mom[2]:.2f} (static {m2:.2f}), M3 {mom[3]:+.0f} ({m3:+.0f})")
    say(f"herzfeld_berger at {nur:.0f} Hz / {nu0:g} MHz: sum 1 over +-n_sidebands_needed, M1 = 0, {'; '.join(parts)} ppm^n")


def _sstr_czjzek(say) -> None:
    """larmor.czjzek_dist against mrsimulator: the d = 5 weights against the
    fit engine's kernel weights (mrsimulator's analytic density) and the
    closed-form / marginal read-outs against its random-tensor Monte-Carlo."""
    from mrsimulator.models import CzjzekDistribution

    from larmor import czjzek_dist
    from larmor.engine import CzjzekKernel

    sigma = 1.6
    cq = np.linspace(0.05, 25.0, 60)
    eta = np.linspace(0.0, 1.0, 11)
    kern = CzjzekKernel(x_ppm=np.zeros(2), K=np.zeros((cq.size * eta.size, 2), np.float32),
                        cq_grid_MHz=cq, eta_grid=eta)
    ref = kern.weights(sigma)
    w5 = czjzek_dist.czjzek_weights(sigma, 5.0, cq, eta)
    w3 = czjzek_dist.czjzek_weights(sigma, 3.0, cq, eta)
    d5, d3 = float(np.abs(w5 - ref).max() / ref.max()), float(np.abs(w3 - ref).max() / ref.max())
    _check(d5 < 1e-9 and d3 > 1e-2, f"czjzek_weights vs the kernel's mrsimulator weights: d=5 {d5:.2e}, d=3 {d3:.2e}")
    state = np.random.get_state()
    try:
        np.random.seed(20261005)            # mrsimulator samples the random tensors with np.random
        zeta, eta_mc = CzjzekDistribution(sigma, cache=False).rvs(size=40000)
    finally:
        np.random.set_state(state)
    rms_mc = float(np.sqrt(np.mean(zeta ** 2 * (1.0 + eta_mc ** 2 / 3.0))))
    ax = czjzek_dist.suggested_cq_axis(sigma, 2000)
    mean_cq = float((ax * czjzek_dist.marginal_cq(sigma, ax)).sum())
    eg = np.linspace(0.0, 1.0, 401)
    mean_eta = float((czjzek_dist.czjzek_pdf(sigma, ax, eg).sum(axis=1) * eg).sum())
    # 40000 samples: the three means are good to ~0.2 % (1 sigma); 1.5 % is > 6 sigma
    dev = (rms_mc / czjzek_dist.rms_pq(sigma) - 1.0, float(np.abs(zeta).mean()) / mean_cq - 1.0,
           float(eta_mc.mean()) / mean_eta - 1.0)
    _check(all(abs(v) < 0.015 for v in dev),
           f"Czjzek sigma {sigma}: Monte-Carlo rms P_Q {rms_mc:.4f} vs rms_pq {czjzek_dist.rms_pq(sigma):.4f}, "
           f"<|C_Q|> {np.abs(zeta).mean():.4f} vs marginal {mean_cq:.4f}, <eta> {eta_mc.mean():.4f} vs pdf {mean_eta:.4f}")
    say(f"czjzek_dist: d=5 weights = the fit kernel's mrsimulator weights to {d5:.0e} (d=3 differs by {d3:.2f}); "
        f"40000 random GIM tensors (sigma {sigma}): rms P_Q {rms_mc:.3f} vs 2 sqrt(5) sigma {czjzek_dist.rms_pq(sigma):.3f}, "
        f"<|C_Q|> {np.abs(zeta).mean():.3f} vs marginal {mean_cq:.3f}, <eta> {eta_mc.mean():.3f} vs pdf {mean_eta:.3f}")


def _sstr_estimate(say) -> None:
    """larmor.estimate beyond engine-physics: a Gaussian's start width, a MAS
    quad_ct start value against its truth, and the cq_for_width inversion."""
    from larmor import engine, estimate
    from larmor.recipe import Param, Recipe, SiteModel

    x = np.linspace(-100.0, 100.0, 2001)
    gauss = 7.0 * np.exp(-4.0 * np.log(2.0) * ((x - 12.0) / 12.0) ** 2)
    sv = estimate.start_values("gauss_lor", x, gauss, "27Al", 130.3)
    _check(abs(sv.get("shift_fwhm_ppm", 0.0) - 12.0) < 0.2, f"start_values(gauss_lor) on a 12 ppm FWHM Gaussian: {sv}")

    def quad(cq, fwhm):
        return Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, spin_rate_Hz=20000.0, sites=[
            SiteModel(model="quad_ct", label="Q", params={
                "isotropic_chemical_shift_ppm": Param(20.0), "Cq_MHz": Param(cq), "eta": Param(0.6),
                "shift_fwhm_ppm": Param(fwhm), "amplitude": Param(1.0)})])

    gx, ys, _ = engine.simulate(quad(6.0, 1.0), exp_ppm=np.linspace(-150.0, 150.0, 1201))
    mas = estimate.start_values("quad_ct", gx, ys, "27Al", 130.3, spin_rate_Hz=20000.0)
    # a seed only has to land in the basin: 15 % (measured -2 %); the old
    # static-probe bug seeded this MAS pattern ~1.9x low
    _check(abs(mas.get("Cq_MHz", 0.0) / 6.0 - 1.0) < 0.15, f"start_values(quad_ct, 20 kHz MAS, C_Q 6 MHz): {mas}")
    cq = estimate.cq_for_width(30.0, "27Al", 130.3, eta=0.6, spin_rate_Hz=20000.0)
    bx, by, _ = engine.simulate(quad(cq, 1.0), exp_ppm=np.linspace(-250.0, 250.0, 2501))
    width = estimate.band_width_ppm(bx, by, frac=estimate.CALIB_FRAC)[1]
    # the bisection stops within 5 % of the target width; 6 % plus the axis step
    _check(cq > 0 and abs(width / 30.0 - 1.0) < 0.06, f"cq_for_width(30 ppm) = {cq:.3f} MHz renders {width:.2f} ppm wide")
    say(f"estimate: gauss_lor start width {sv['shift_fwhm_ppm']:.2f} ppm (truth 12); 20 kHz MAS quad_ct seeded at "
        f"C_Q {mas['Cq_MHz']:.2f} MHz (truth 6); cq_for_width(30 ppm) = {cq:.3f} MHz renders {width:.2f} ppm at "
        f"{estimate.CALIB_FRAC:.0%} height")


def _sstr_measure(say) -> None:
    """larmor.measure beyond engine-physics: an absolute integral against the
    analytic Gaussian area on both axis orders, and S/N against the known
    height and noise."""
    from larmor import measure

    x = np.linspace(-50.0, 50.0, 2001)
    amp, fw = 3.0, 8.0
    y = amp * np.exp(-4.0 * np.log(2.0) * ((x - 5.0) / fw) ** 2)
    area = amp * fw * np.sqrt(np.pi / (4.0 * np.log(2.0)))
    a_up = measure.integrate(x, y, (30.0, -20.0))
    a_dn = measure.integrate(x[::-1], y[::-1], (30.0, -20.0))
    _check(abs(a_up / area - 1.0) < 1e-4 and abs(a_dn - a_up) < 1e-9 * area,
           f"integrate: {a_up:.6f} ascending, {a_dn:.6f} descending, analytic {area:.6f}")
    xs = np.linspace(200.0, -200.0, 4000)
    ys = 50.0 * np.exp(-4.0 * np.log(2.0) * (xs / 4.0) ** 2) + np.random.default_rng(8).normal(0.0, 0.5, xs.size)
    s = measure.snr(xs, ys, (10.0, -10.0))
    # height 50 over noise 0.5 = 100; 800 noise points estimate the RMS to ~2.5 %
    _check(abs(s / 100.0 - 1.0) < 0.1, f"snr {s:.1f} for a height-50 line over 0.5 RMS noise (truth 100)")
    say(f"measure: integral {a_up:.5f} (analytic {area:.5f}) on ascending and descending axes; S/N {s:.1f} (truth 100)")


def _sstr_refranges(say) -> None:
    """larmor.refranges beyond engine-physics: every band / position resolves
    a citation, the narrowest-band rule, the reported-position fallback, no
    answer far from everything, families, and the auto-label grammar."""
    from larmor import refranges

    bad = [(nuc, e["label"]) for nuc, rows in refranges.REF_RANGES.items() for e in rows
           if e.get("ref") not in refranges.REFS or not e["lo_ppm"] < e["hi_ppm"]]
    bad += [(nuc, e["label"]) for nuc, rows in refranges.REF_POSITIONS.items() for e in rows
            if e.get("ref") not in refranges.REFS]
    n_bands = sum(len(v) for v in refranges.REF_RANGES.values())
    n_pos = sum(len(v) for v in refranges.REF_POSITIONS.values())
    _check(not bad, f"literature entries without a resolvable citation or with lo >= hi: {bad}")
    o = refranges.assign("17O", 50.0)          # in the 40-60 BO band AND the 30-75 NBO band
    f = refranges.assign("19F", -20.0)         # no 19F band there; beta-PbF2 is reported at -20
    _check(o is not None and (o["lo_ppm"], o["hi_ppm"]) == (40.0, 60.0) and o["kind"] == "range"
           and f is not None and f["label"] == "PbF2" and f["kind"] == "position"
           and refranges.assign("19F", 150.0) is None and refranges.assign("27Al", 150.0) is None,
           f"assign: 17O 50 ppm -> {o and o['label']}, 19F -20 ppm -> {f and f['label']}")
    _check(refranges.family_for("27Al", 37.0) == "Al(V)" and refranges.family_for("29Si", -110.0) == "",
           f"family_for: 27Al 37 ppm {refranges.family_for('27Al', 37.0)!r}, 29Si {refranges.family_for('29Si', -110.0)!r}")
    autos = ("Czjzek-3", "pk-0", "line-copy", "A+1sb", "")
    _check(all(refranges.is_auto_label(s) for s in autos) and not refranges.is_auto_label("Al[4]"),
           f"is_auto_label: {[refranges.is_auto_label(s) for s in autos]}, Al[4] {refranges.is_auto_label('Al[4]')}")
    say(f"refranges: {n_bands} bands + {n_pos} positions all cited; 17O 50 ppm -> the narrower {o['label']} band, "
        f"19F -20 ppm -> {f['label']} (position), nothing at +150 ppm; 27Al 37 ppm -> family Al(V); auto labels recognised")


def _sstr_sanity(say) -> None:
    """larmor.sanity beyond engine-error-tools: one injected defect at a time
    must give exactly one warning naming that parameter."""
    from larmor import sanity
    from larmor.recipe import Param, Recipe, SiteModel

    def base():
        return Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, spin_rate_Hz=20000.0, sites=[
            SiteModel(model="quad_ct", label="Q", params={
                "isotropic_chemical_shift_ppm": Param(60.0), "Cq_MHz": Param(5.0), "eta": Param(0.4),
                "shift_fwhm_ppm": Param(2.0), "amplitude": Param(10.0)}),
            SiteModel(model="czjzek", label="C", params={
                "isotropic_chemical_shift_ppm": Param(10.0), "sigma_Cq_MHz": Param(1.5), "shift_fwhm_ppm": Param(5.0),
                "line_fwhm_ppm": Param(0.5, vary=False), "amplitude": Param(5.0)}),
            SiteModel(model="gauss_lor", label="G", params={
                "isotropic_chemical_shift_ppm": Param(-20.0), "shift_fwhm_ppm": Param(4.0), "amplitude": Param(3.0),
                "gl": Param(1.0, vary=False)})])

    window = (120.0, -60.0)
    clean = sanity.check_recipe(base(), window=window)
    _check(clean == [], f"sanity on a valid recipe: {clean}")
    cases = [(0, "eta", 1.3, "outside"), (0, "Cq_MHz", -2.0, "negative"), (0, "Cq_MHz", 150.0, "above"),
             (1, "sigma_Cq_MHz", 0.0, "non-physical width"), (2, "amplitude", -5.0, "negative"),
             (2, "shift_fwhm_ppm", float("nan"), "not finite"),
             (2, "isotropic_chemical_shift_ppm", 200.0, "outside the fit window")]
    for site, name, value, words in cases:
        rec = base()
        rec.sites[site].params[name] = Param(value)
        w = sanity.check_recipe(rec, window=window)
        _check(len(w) == 1 and w[0]["site"] == site and w[0]["param"] == name and words in w[0]["message"]
               and rec.sites[site].label in sanity.summarize(w),
               f"sanity: {name} = {value} on site {site} gave {w}")
    say(f"sanity: a valid 3-line recipe is clean; each of {len(cases)} injected defects (eta 1.3, C_Q -2 / 150 MHz, "
        f"sigma 0, amplitude -5, NaN width, a line at 200 ppm outside the window) gives exactly its one warning")


def _sstr_diagnostics(say) -> None:
    """larmor.diagnostics beyond engine-error-tools: the runs test on
    sequences with a known run count, and lag-1 autocorrelation of an AR(1)
    series with a known coefficient."""
    import math

    from larmor import diagnostics

    block = np.r_[np.ones(20), -np.ones(20)]             # 2 runs
    alt = np.tile([1.0, -1.0], 20)                       # 40 runs
    n1 = n2 = 20
    n = n1 + n2
    mu = 1.0 + 2.0 * n1 * n2 / n                         # Wald-Wolfowitz mean and variance
    sd = math.sqrt(2.0 * n1 * n2 * (2.0 * n1 * n2 - n) / (n * n * (n - 1.0)))
    rb, ra = diagnostics.runs_test(block), diagnostics.runs_test(alt)
    _check(rb["runs"] == 2 and ra["runs"] == 40 and abs(rb["z"] - (2 - mu) / sd) < 1e-9
           and abs(ra["z"] - (40 - mu) / sd) < 1e-9 and rb["structured"] and ra["structured"]
           and abs(rb["p"] - math.erfc(abs(rb["z"]) / math.sqrt(2.0))) < 1e-12,
           f"runs_test: block {rb}, alternating {ra}; expected z {(2 - mu) / sd:.3f} / {(40 - mu) / sd:.3f}")
    eps = np.random.default_rng(21).normal(0.0, 1.0, 5000)
    ar = np.empty_like(eps)
    ar[0] = eps[0]
    for i in range(1, eps.size):
        ar[i] = 0.7 * ar[i - 1] + eps[i]
    lag = diagnostics.lag1_autocorr(ar)
    verdict = diagnostics.residual_structure(ar)
    # the lag-1 estimate of an AR(1) series has SD sqrt((1 - phi^2)/N) = 0.01
    _check(abs(lag - 0.7) < 0.05 and verdict["structured"], f"AR(1) phi 0.7: lag-1 {lag:.3f}, verdict {verdict}")
    say(f"diagnostics: runs test z {rb['z']:.3f} for 2 runs / {ra['z']:+.3f} for 40 runs in 40 signs (Wald-Wolfowitz "
        f"{(2 - mu) / sd:.3f} / {(40 - mu) / sd:+.3f}), p {rb['p']:.1e}; AR(1) phi 0.7 -> lag-1 {lag:.3f}, structured")


def _sstr_identifiability(say) -> None:
    """larmor.identifiability on a truly degenerate pair from a real fit: two
    co-centred Gaussians with only their amplitudes free. The model is linear,
    so corr(a1, a2) = -<g1, g2> / (|g1| |g2|) = -sqrt(2 w1 w2 / (w1^2 + w2^2))."""
    from larmor import fit as fitmod, identifiability
    from larmor.recipe import Param, Recipe, SiteModel

    w1, w2 = 10.0, 13.0
    x = np.linspace(100.0, -100.0, 1001)
    g = lambda w: np.exp(-4.0 * np.log(2.0) * (x / w) ** 2)  # noqa: E731
    y = 60.0 * g(w1) + 40.0 * g(w2) + np.random.default_rng(4).normal(0.0, 0.5, x.size)

    def line(label, w):
        return SiteModel(model="gauss_lor", label=label, params={
            "isotropic_chemical_shift_ppm": Param(0.0, vary=False), "shift_fwhm_ppm": Param(w, vary=False),
            "amplitude": Param(50.0, min=0.0), "gl": Param(1.0, vary=False)})

    rec = Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, spin_rate_Hz=20000.0, sites=[line("A", w1), line("B", w2)])
    res = fitmod.fit(rec, x, y, window_ppm=(100.0, -100.0))
    names, corr = identifiability.corr_matrix(res.lmfit_result)
    pairs = identifiability.unidentifiable_pairs(res.lmfit_result)
    r_true = -np.sqrt(2.0 * w1 * w2 / (w1 ** 2 + w2 ** 2))
    r = float(corr[0, 1]) if corr is not None and corr.shape == (2, 2) else float("nan")
    _check(names == ["s0_amp", "s1_amp"] and abs(r - r_true) < 5e-3 and len(pairs) == 1
           and set(pairs[0][:2]) == {"s0_amp", "s1_amp"} and abs(pairs[0][2] - r) < 1e-12,
           f"identifiability: {names}, r {r:.5f} (closed form {r_true:.5f}), pairs {pairs}")
    amps = [s.params["amplitude"].value for s in rec.sites]
    say(f"identifiability: co-centred Gaussians {w1:g} / {w2:g} ppm, amplitudes free -> r {r:.4f} "
        f"(closed form {r_true:.4f}), flagged {pairs[0][0]} ~ {pairs[0][1]}; fitted {amps[0]:.1f} + {amps[1]:.1f} "
        f"(truth 60 + 40)")
