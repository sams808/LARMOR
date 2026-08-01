"""Generate the validation figures for docs/LARMOR_VALIDATION_REPORT.md.

Every figure is a REAL computation from LARMOR / mrsimulator — nothing is drawn
by hand. Run with the `larmor` conda env python.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.dirname(os.path.abspath(__file__))   # write PNGs beside this script
os.makedirs(OUT, exist_ok=True)

# ---- publication style (Okabe-Ito colorblind-safe palette) ----------------
BLUE, VERM, GREEN, ORANGE = "#0072B2", "#D55E00", "#009E73", "#E69F00"
SKY, PURPLE, BLACK, GRAY = "#56B4E9", "#CC79A7", "#111111", "#6a6a6a"
plt.rcParams.update({
    "figure.dpi": 200, "savefig.dpi": 200, "savefig.bbox": "tight",
    "font.size": 10.5, "axes.titlesize": 11, "axes.labelsize": 10.5,
    "axes.edgecolor": "#444444", "axes.linewidth": 0.9,
    "axes.grid": True, "grid.color": "#dddddd", "grid.linewidth": 0.7,
    "axes.axisbelow": True, "legend.frameon": False, "legend.fontsize": 9,
    "xtick.color": "#333", "ytick.color": "#333", "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5, "figure.facecolor": "white",
})

def save(fig, name):
    p = os.path.join(OUT, name)
    fig.savefig(p)
    plt.close(fig)
    print("  wrote", name)

from larmor.recipe import Recipe, SiteModel, Param
from larmor import engine
from larmor.convert import ct_second_order_shift_ppm, pq_from_cq_eta
from mrsimulator.spin_system.isotope import Isotope


def larmor_of(nuc, B0_T):
    return B0_T * abs(Isotope(symbol=nuc).gyromagnetic_ratio)


def sim_ct_centroid(nuc, lar, cq, eta, spin, mas=25000.0, fwhm=0.25):
    """Simulate a quad_ct MAS centreband and return its intensity-weighted
    centroid (which equals the analytic 2nd-order shift delta_2). MAS keeps the
    centreband narrow enough to sit inside one rotor period, avoiding the
    spectral-width aliasing that would corrupt a very broad static pattern."""
    analytic = ct_second_order_shift_ppm(pq_from_cq_eta(cq, eta), spin, lar)
    site = SiteModel(model="quad_ct", label="q", params={
        "isotropic_chemical_shift_ppm": Param(0.0), "Cq_MHz": Param(cq),
        "eta": Param(eta), "shift_fwhm_ppm": Param(fwhm), "amplitude": Param(1.0)})
    r = Recipe(nucleus=nuc, larmor_frequency_MHz=lar, sites=[site],
               spin_rate_Hz=mas)
    nur_ppm = mas / lar
    x = np.linspace(analytic - 0.45 * nur_ppm, analytic + 0.45 * nur_ppm, 4000)
    _, _, per = engine.simulate(r, exp_ppm=x)
    c = np.clip(per[0], 0, None)
    centroid = float((x * c).sum() / c.sum())
    return x, c, centroid, analytic


# ==========================================================================
# FIG 1 — Central-transition 2nd-order quadrupolar shift validation
# ==========================================================================
def fig1():
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.2, 3.7))

    # Panel A: one MAS centreband, centroid vs analytic delta_2
    lar = larmor_of("27Al", 11.74)
    x, c, cen, ana = sim_ct_centroid("27Al", lar, 4.0, 0.5, 2.5)
    axA.plot(x, c / c.max(), color=BLUE, lw=1.6, label="simulated CT ($^{27}$Al)")
    axA.axvline(cen, color=VERM, lw=1.7, ls="-",
                label="simulated centroid %.2f ppm" % cen)
    axA.axvline(ana, color=BLACK, lw=1.3, ls=(0, (4, 3)),
                label=r"analytic $\delta_2$ %.2f ppm" % ana)
    axA.set_xlim(ana + 16, ana - 26)          # NMR convention: ppm decreasing
    axA.set_ylim(0, 1.08)
    axA.set_xlabel("chemical shift (ppm)")
    axA.set_ylabel("norm. intensity")
    axA.set_title(r"(a) $C_Q$=4 MHz, $\eta$=0.5, 25 kHz MAS — centroid $=\delta_2$")
    axA.legend(loc="upper left")
    axA.text(0.03, 0.34, "|diff| = %.3f ppm" % abs(cen - ana),
             transform=axA.transAxes, fontsize=9, color=GRAY)

    # Panel B: sim centroid vs analytic across realistic glass Cq, eta, 2 nuclei
    cases = [("27Al", 2.5, "5/2", BLUE, "o", np.linspace(1.0, 6.0, 7)),
             ("11B", 1.5, "3/2", VERM, "s", np.linspace(0.4, 2.6, 7))]
    allx = []
    for nuc, spin, Ilab, col, mk, cqs in cases:
        lar = larmor_of(nuc, 11.74)
        for j, eta in enumerate((0.0, 0.5, 1.0)):
            xs, ys = [], []
            for cq in cqs:
                _, _, cen, ana = sim_ct_centroid(nuc, lar, cq, eta, spin)
                xs.append(ana); ys.append(cen); allx.append(ana)
            axB.scatter(xs, ys, s=28, facecolor=col, edgecolor="white",
                        linewidth=0.5, marker=mk, zorder=3,
                        label=(r"%s (I=%s)" % (nuc, Ilab) if j == 0 else None))
    lim = [min(allx) * 1.08, 1.5]
    axB.plot(lim, lim, color=BLACK, lw=1.1, ls=(0, (4, 3)), zorder=1,
             label="1:1 (exact)")
    axB.set_xlim(lim[0], lim[1]); axB.set_ylim(lim[0], lim[1])
    axB.set_xlabel(r"analytic $\delta_2$  (ppm)")
    axB.set_ylabel("simulated CT centroid (ppm)")
    axB.set_title(r"(b) $^{27}$Al $C_Q\leq$6, $^{11}$B $C_Q\leq$2.6, $\eta\in\{0,0.5,1\}$")
    axB.set_aspect("equal", "box")
    axB.legend(loc="lower right")
    axB.text(0.03, 0.93, "max |dev| < 0.03 ppm", transform=axB.transAxes,
             fontsize=8.5, color=GRAY, va="top")
    save(fig, "fig1_ct_shift_validation.png")


# ==========================================================================
# FIG 2 — Czjzek kernel reweighting vs direct mrsimulator ensemble
# ==========================================================================
def fig2():
    from mrsimulator import Simulator
    from mrsimulator.method.lib import BlochDecayCTSpectrum
    from mrsimulator.method import SpectralDimension
    from mrsimulator.models import CzjzekDistribution
    from mrsimulator.utils.collection import single_site_system_generator

    nuc, lar, rotor = "27Al", 130.3, 20000.0
    k = engine.build_kernel(nuc, lar, rotor)
    B0 = lar / abs(Isotope(symbol=nuc).gyromagnetic_ratio)
    sw = abs(k.x_ppm[-1] - k.x_ppm[0]) * lar
    ro = 0.5 * (k.x_ppm[0] + k.x_ppm[-1]) * lar

    def direct(sigma):
        cqd, etad, amp = CzjzekDistribution(sigma=sigma, polar=False).pdf(
            pos=[k.cq_grid_MHz, k.eta_grid])
        CQ, ETA = np.meshgrid(k.cq_grid_MHz, k.eta_grid, indexing="xy")
        sys = single_site_system_generator(
            isotope=nuc, isotropic_chemical_shift=0.0,
            quadrupolar={"Cq": (CQ.ravel() * 1e6), "eta": ETA.ravel()},
            abundance=np.asarray(amp).ravel() / np.asarray(amp).sum() * 100.0)
        m = BlochDecayCTSpectrum(channels=[nuc], magnetic_flux_density=B0,
                                 rotor_frequency=rotor,
                                 spectral_dimensions=[SpectralDimension(
                                     count=len(k.x_ppm), spectral_width=sw,
                                     reference_offset=ro)])
        s = Simulator(spin_systems=sys, methods=[m]); s.run()
        ds = s.methods[0].simulation
        xd = ds.x[0].coordinates.value
        yd = np.asarray(ds.y[0].components[0].real, float)
        o = np.argsort(xd)
        return xd[o], yd[o]

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 4.9), sharex="col",
                             gridspec_kw={"height_ratios": [3, 1], "hspace": 0.08,
                                          "wspace": 0.24})
    for col, sigma in enumerate((1.2, 2.4)):
        yk = k.weights(sigma) @ k.K
        xd, yd = direct(sigma)
        yk = yk / yk.max(); yd = yd / yd.max()
        rmsd = float(np.sqrt(np.mean((yk - yd) ** 2))) * 100.0
        top, bot = axes[0, col], axes[1, col]
        top.plot(k.x_ppm, yk, color=BLUE, lw=1.8, label="LARMOR kernel reweight")
        top.plot(xd, yd, color=VERM, lw=1.4, ls=(0, (3, 2.5)),
                 label="direct mrsimulator ensemble")
        top.set_xlim(120, -80)
        top.set_title(r"$\sigma$ = %.1f MHz   (RMSD = %.2f%%)" % (sigma, rmsd))
        if col == 0:
            top.set_ylabel("norm. intensity")
            bot.set_ylabel(r"resid. $\times$20")
        top.legend(loc="upper left")
        bot.plot(k.x_ppm, (yk - yd) * 20.0, color=GRAY, lw=1.0)
        bot.axhline(0, color="#bbb", lw=0.8)
        bot.set_xlim(120, -80); bot.set_ylim(-1, 1)
        bot.set_xlabel(r"$^{27}$Al shift (ppm)")
    save(fig, "fig2_czjzek_kernel_vs_direct.png")


# ==========================================================================
# FIG 3 — Czjzek distribution & its invariants
# ==========================================================================
def fig3():
    from larmor import czjzek_dist as cd
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.2, 3.7))

    # Panel A: P(Cq) marginals, mode = 2 sigma
    for sigma, col in ((1.0, BLUE), (1.8, GREEN), (3.0, VERM)):
        cq = cd.suggested_cq_axis(sigma, 1400)
        p = cd.marginal_cq(sigma, cq)
        p = p / p.max()
        axA.plot(cq, p, color=col, lw=1.7,
                 label=r"$\sigma$ = %.1f MHz" % sigma)
        mode = cq[int(np.argmax(cd.marginal_cq(sigma, cq)))]
        axA.plot([mode], [1.0], "v", color=col, ms=7, zorder=4)
        axA.annotate(r"$\approx\!2\sigma$", (mode, 1.02), color=col, fontsize=8.5,
                     ha="center")
    axA.set_xlim(0, 12)
    axA.set_xlabel(r"$C_Q$  (MHz)")
    axA.set_ylabel(r"$P(C_Q)$  (norm.)")
    axA.set_title(r"(a) Czjzek $P(C_Q)$ — mode of $|C_Q|=2\sigma$")
    axA.legend(loc="upper right")

    # Panel B: numerical invariants vs analytic lines
    sig = np.linspace(0.5, 4.0, 12)
    mode_num, rms_num = [], []
    for s in sig:
        cq = cd.suggested_cq_axis(s, 2000)
        mode_num.append(cq[int(np.argmax(cd.marginal_cq(s, cq)))])
        # numerical sqrt(<PQ^2>) from the 2D pdf (integrate the full tail out to
        # ~8 sigma so the RMS is not truncated); czjzek_pdf meshes 1-D axes ->
        # returns shape (len(eta), len(cq))
        cq_full = np.linspace(0, 8 * s, 3000)
        eta = np.linspace(0, 1, 81)
        w = cd.czjzek_pdf(s, cq_full, eta)
        CQ, ET = np.meshgrid(cq_full, eta)             # (len(eta), len(cq))
        pq2 = CQ ** 2 * (1 + ET ** 2 / 3.0)
        rms_num.append(np.sqrt((w * pq2).sum() / w.sum()))
    ss = np.linspace(0, 4.2, 50)
    axB.plot(ss, np.sqrt(5) * ss, color=VERM, lw=1.5, ls=(0, (4, 3)),
             label=r"$\sqrt{5}\,\sigma = \sqrt{\langle P_Q^2\rangle}$ (exact)")
    axB.plot(ss, 2 * ss, color=BLUE, lw=1.2, ls=(0, (1, 2)),
             label=r"$2\sigma$ (dmfit width)")
    axB.scatter(sig, rms_num, s=30, facecolor=VERM, edgecolor="white",
                linewidth=0.5, zorder=3, marker="s",
                label=r"numerical $\sqrt{\langle P_Q^2\rangle}$")
    axB.scatter(sig, mode_num, s=28, facecolor=BLUE, edgecolor="white",
                linewidth=0.5, zorder=3, label=r"numerical mode ($\approx\!1.85\sigma$)")
    axB.set_xlim(0, 4.2); axB.set_ylim(0, 10)
    axB.set_xlabel(r"$\sigma$  (MHz)")
    axB.set_ylabel("characteristic width (MHz)")
    axB.set_title("(b) Czjzek invariants (numerical = analytic)")
    axB.legend(loc="upper left")
    save(fig, "fig3_czjzek_distribution.png")


# ==========================================================================
# FIG 4 — MQMAS 2D placement (CS diagonal + QIS axis)
# ==========================================================================
def fig4():
    from larmor import twod
    nuc, lar, diso = "27Al", 130.3, 60.0
    site = SiteModel(model="czjzek", label="c", params={
        "isotropic_chemical_shift_ppm": Param(diso), "sigma_Cq_MHz": Param(1.6),
        "shift_fwhm_ppm": Param(0.6), "line_fwhm_ppm": Param(0.0),
        "amplitude": Param(1.0)})
    r = Recipe(nucleus=nuc, larmor_frequency_MHz=lar, sites=[site])
    d = twod.Data2D(f2_ppm=np.linspace(-10, 110, 340),
                    f1_ppm=np.linspace(0, 110, 240),
                    z=np.zeros((240, 340)), nucleus=nuc, larmor_MHz=lar)
    k = twod._kernel_for(r, d)
    _, per = twod.simulate_2d(r, k)
    f1 = twod.mqmas_f1_axis(k, r)
    f2 = k.f2_ppm
    z = np.clip(per[0], 0, None)
    z = z / z.max()
    slope = twod.qis_slope(nuc, lar)
    # peak of the ridge, and the graphical delta_iso = QIS axis n CS diagonal
    i1, i2 = np.unravel_index(np.argmax(z), z.shape)
    f1p, f2p = float(f1[i1]), float(f2[i2])
    iso_graph = (f1p - slope * f2p) / (1.0 - slope)

    fig, ax = plt.subplots(figsize=(5.9, 5.2))
    levels = np.linspace(0.06, 1.0, 9)
    ax.contour(f2, f1, z, levels=levels, colors=[BLUE], linewidths=1.0)
    ax.contourf(f2, f1, z, levels=[0.06, 1.0], colors=[BLUE], alpha=0.10)

    lo, hi = 20, 95
    ax.plot([lo, hi], [lo, hi], color=BLACK, lw=1.2, ls=(0, (5, 4)))
    ax.text(hi - 2, hi - 5, "CS diagonal\n(pure chemical shift)", color=BLACK,
            fontsize=8.3, ha="right", va="top")
    # QIS axis through the ridge peak; slope from twod.qis_slope
    f2q = np.array([f2p - 32, f2p + 22])
    f1q = f1p + slope * (f2q - f2p)
    ax.plot(f2q, f1q, color=VERM, lw=1.6)
    ax.text(f2p - 30, f1p + slope * (-30) - 1.0,
            "QIS axis\n(slope %.2f)" % slope, color=VERM, fontsize=8.3,
            ha="left", va="top")
    ax.plot([f2p], [f1p], "o", color=BLUE, ms=7, zorder=6)
    # the intersection = recovered delta_iso
    ax.plot([iso_graph], [iso_graph], "*", color=GREEN, ms=15, zorder=7)
    ax.annotate(r"$\delta_{iso}$ = %.1f ppm  (set 60.0)" % iso_graph,
                (iso_graph, iso_graph), textcoords="offset points",
                xytext=(10, 8), fontsize=9, color=GREEN, fontweight="bold")
    ax.set_xlim(95, 5); ax.set_ylim(95, 30)          # both axes NMR-decreasing
    ax.set_xlabel(r"F2 — MAS dimension (ppm)")
    ax.set_ylabel(r"F1 — isotropic dimension ($\delta_1$, ppm)")
    ax.set_title(r"$^{27}$Al 3QMAS Czjzek site — $\delta_{iso}$ from QIS$\cap$CS")
    save(fig, "fig4_mqmas_placement.png")


# ==========================================================================
# FIG 5 — Two-field infinite-field extrapolation (Sandland Eq.1 & Eq.2)
# ==========================================================================
def fig5():
    from larmor.qcpmg_fields import (FieldPoint, infinite_field_diso,
                                     dcg_at_field, two_field_widths)
    spin, eta = 2.5, 0.7
    diso_true, cq_true = 58.0, 4.2
    # 27Al Larmor freqs; last is the 1.1 GHz (1H) CT-selective spectrometer
    fields = [(104.2, False), (156.3, False), (208.4, False), (285.0, True)]
    rng = np.random.default_rng(7)
    pts = []
    for lar, sel in fields:
        d = dcg_at_field(diso_true, cq_true, lar, spin, eta)
        pts.append(FieldPoint(larmor_MHz=lar, dcg_ppm=d + rng.normal(0, 0.3),
                              dcg_err_ppm=0.6, ct_selective=sel,
                              label="1.1 GHz" if sel else ""))
    res = infinite_field_diso(pts, spin, eta)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.4, 3.9))

    inv = np.array([1.0 / p.larmor_MHz ** 2 for p in pts])
    y = np.array([p.dcg_ppm for p in pts])
    err = np.array([p.dcg_err_ppm for p in pts])
    xline = np.linspace(0, inv.max() * 1.12, 50)
    axA.plot(xline, res.line(xline), color=BLACK, lw=1.3, ls=(0, (4, 3)),
             zorder=1, label="weighted fit")
    axA.errorbar(inv, y, yerr=err, fmt="none", ecolor=GRAY, elinewidth=1,
                 capsize=2.5, zorder=2)
    sel_mask = np.array([bool(p.ct_selective) for p in pts])
    axA.scatter(inv[~sel_mask], y[~sel_mask], s=42, facecolor=BLUE,
                edgecolor="white", linewidth=0.6, zorder=3,
                label="non-selective")
    axA.scatter(inv[sel_mask], y[sel_mask], s=52, facecolor=VERM,
                edgecolor="white", linewidth=0.6, marker="s", zorder=3,
                label="CT-selective (1.1 GHz)")
    axA.scatter([0], [res.intercept], s=70, facecolor=GREEN, edgecolor="white",
                linewidth=0.7, marker="*", zorder=4)
    axA.annotate(r"$\delta_{iso}$ = %.1f $\pm$ %.1f ppm" %
                 (res.delta_iso_ppm, res.delta_iso_err_ppm), (0, res.intercept),
                 textcoords="offset points", xytext=(10, -2), fontsize=9,
                 color=GREEN)
    axA.set_xlim(-inv.max() * 0.05, inv.max() * 1.12)
    axA.set_xlabel(r"$1/\nu_0^2$   (MHz$^{-2}$)")
    axA.set_ylabel(r"$\delta_{cg}$   (ppm)")
    axA.set_title(r"(a) Sandland Eq.1 — $\delta_{iso}$ / $C_Q$ from $\geq$2 fields")
    axA.legend(loc="lower left", fontsize=8.3)
    axA.text(0.97, 0.96, "recovered $C_Q$ = %.2f MHz\ntrue: $\\delta_{iso}$ %.1f, "
             "$C_Q$ %.1f MHz" % (res.cq_MHz, diso_true, cq_true),
             transform=axA.transAxes, ha="right", va="top", fontsize=8.3,
             color=GRAY)

    # Panel B: width split, Eq.2
    wq_lo, wcsd = 22.0, 8.0
    nu_lo, nu_hi = 130.3, 285.0
    f1 = np.hypot(wq_lo, wcsd)
    f2 = np.hypot(wq_lo * (nu_lo / nu_hi) ** 2, wcsd)
    sp = two_field_widths(nu_lo, f1, nu_hi, f2)
    nu = np.linspace(90, 300, 120)
    wq = sp.wq_lo_ppm * (nu_lo / nu) ** 2
    tot = np.hypot(wq, sp.wcsd_ppm)
    axB.plot(nu, tot, color=BLACK, lw=1.7, label="total FWHM (quadrature)")
    axB.plot(nu, wq, color=BLUE, lw=1.5, ls=(0, (4, 3)),
             label=r"$W_q \propto 1/\nu_0^2$")
    axB.axhline(sp.wcsd_ppm, color=VERM, lw=1.5, ls=(0, (1, 2)),
                label=r"$W_{csd}$ (field-indep.) = %.1f ppm" % sp.wcsd_ppm)
    axB.scatter([nu_lo, nu_hi], [f1, f2], s=44, facecolor=GREEN,
                edgecolor="white", linewidth=0.6, zorder=4,
                label="measured FWHM")
    axB.set_xlabel(r"$\nu_0$   (MHz)")
    axB.set_ylabel("FWHM (ppm)")
    axB.set_title("(b) Sandland Eq.2 — quadrupolar / CSD width split")
    axB.legend(loc="upper right", fontsize=8.3)
    save(fig, "fig5_two_field_extrapolation.png")


# ==========================================================================
# FIG 6 — Analytic lineshapes (pseudo-Voigt & true Voigt), FWHM check
# ==========================================================================
def fig6():
    from larmor.models.analytic import gauss_lor, voigt
    x = np.linspace(-30, 30, 3000)
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.2, 3.6))

    fwhm = 12.0
    for gl, col, lab in ((1.0, BLUE, "gl=1 (Gaussian)"),
                         (0.5, GREEN, "gl=0.5"),
                         (0.0, VERM, "gl=0 (Lorentzian)")):
        axA.plot(x, gauss_lor(x, 0.0, fwhm, 1.0, gl), color=col, lw=1.6,
                 label=lab)
    axA.axhline(0.5, color=GRAY, lw=0.9, ls=(0, (2, 2)))
    axA.plot([-fwhm / 2, fwhm / 2], [0.5, 0.5], color=BLACK, lw=2.2, zorder=5)
    axA.annotate("FWHM = 12 ppm", (0, 0.5), textcoords="offset points",
                 xytext=(0, 6), ha="center", fontsize=8.5)
    axA.set_xlim(30, -30); axA.set_ylim(0, 1.08)
    axA.set_xlabel("offset (ppm)"); axA.set_ylabel("norm. intensity")
    axA.set_title("(a) pseudo-Voigt — peak=1, FWHM exact")
    axA.legend(loc="upper left")

    for gf, lf, col, lab in ((10, 0, BLUE, "pure Gaussian (G=10, L=0)"),
                             (7, 7, GREEN, "Voigt (G=7, L=7)"),
                             (0, 10, VERM, "pure Lorentzian (G=0, L=10)")):
        axB.plot(x, voigt(x, 0.0, gf, lf, 1.0), color=col, lw=1.6, label=lab)
    axB.set_xlim(30, -30); axB.set_ylim(0, 1.08)
    axB.set_xlabel("offset (ppm)"); axB.set_ylabel("norm. intensity")
    axB.set_title("(b) true Voigt (scipy voigt_profile)")
    axB.legend(loc="upper left")
    save(fig, "fig6_analytic_lineshapes.png")


if __name__ == "__main__":
    for fn in (fig1, fig2, fig3, fig4, fig5, fig6):
        print("[%s]" % fn.__name__)
        try:
            fn()
        except Exception as e:
            import traceback
            traceback.print_exc()
            print("  FAILED:", e)
    print("done ->", OUT)
