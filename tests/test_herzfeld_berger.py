"""Herzfeld-Berger sideband analysis: the exact intensities, their measurement
off a spectrum, the inversion to (zeta, eta), the tensor-convention converters,
and the Tools dialog that wires them together."""
import os

import numpy as np
import pytest

from larmor import convert as C
from larmor import herzfeld_berger as hb

NU0, NUR = 162.0, 5000.0          # 31P at 9.4 T, 5 kHz MAS


def _m2(intens):
    return sum(v * (n * NUR) ** 2 for n, v in intens.items())


@pytest.mark.parametrize("zeta,eta", [(80.0, 0.5), (-45.0, 0.0), (120.0, 1.0),
                                      (30.0, 0.8)])
def test_intensities_sum_to_one_and_carry_the_static_second_moment(zeta, eta):
    """Parseval fixes the sum; Maricq-Waugh fixes the second moment of the
    sideband pattern to the static one, (zeta nu0)^2 (1 + eta^2/3) / 5."""
    I = hb.sideband_intensities(zeta, eta, NU0, NUR)
    assert sum(I.values()) == pytest.approx(1.0, abs=1e-5)
    theory = (zeta * NU0) ** 2 * (1.0 + eta ** 2 / 3.0) / 5.0
    assert _m2(I) / theory == pytest.approx(1.0, abs=5e-4)
    assert all(v >= 0 for v in I.values())


def test_no_sidebands_when_static_or_isotropic():
    assert hb.sideband_intensities(80.0, 0.5, NU0, 0.0) == {0: 1.0}
    assert hb.sideband_intensities(0.0, 0.5, NU0, NUR) == {0: 1.0}
    assert hb.n_sidebands_needed(80.0, NU0, 0.0) == 0


def test_sign_of_zeta_mirrors_the_manifold():
    Ip = hb.sideband_intensities(80.0, 0.3, NU0, NUR, n_max=4)
    Im = hb.sideband_intensities(-80.0, 0.3, NU0, NUR, n_max=4)
    for n in Ip:
        assert Ip[n] == pytest.approx(Im[-n], abs=1e-9)
    # a positive shielding zeta (negative shift anisotropy) piles intensity on
    # the +N (higher-shift) side, matching the csa_mas model's convention
    assert Ip[1] > Ip[-1]


def test_matches_the_csa_mas_model():
    """The same physics mrsimulator computes for csa_mas: integrate its
    sidebands and compare, order by order."""
    from larmor.models._singlesite import simulate_single_site

    xs, ys = simulate_single_site("31P", int(NU0 * 1000), int(NUR), 0, 0,
                                  8000, 500, False, 32, -300.0, 300.0, 2048)
    xs, ys = np.array(xs), np.array(ys)
    step = NUR / NU0
    ref = {n: ys[np.abs(xs - n * step) < 0.3 * step].sum() for n in range(-6, 7)}
    tot = sum(ref.values())
    I = hb.sideband_intensities(80.0, 0.5, NU0, NUR, n_max=6)
    s = sum(I.values())
    for n in range(-4, 5):
        assert I[n] / s == pytest.approx(ref[n] / tot, abs=2e-3), n


def test_inversion_recovers_zeta_eta_and_sign():
    rng = np.random.default_rng(1)
    for zeta, eta in [(80.0, 0.5), (-45.0, 0.2)]:
        I = hb.sideband_intensities(zeta, eta, NU0, NUR)
        meas = {n: max(I[n] * (1 + 0.02 * rng.standard_normal()), 0.0)
                for n in range(-4, 5)}
        out = hb.fit_sidebands(meas, NU0, NUR)
        assert out.ok
        assert out.zeta_ppm == pytest.approx(zeta, abs=0.03 * abs(zeta))
        assert out.eta == pytest.approx(eta, abs=0.05)
        assert out.rms < 5e-3
    bad = hb.fit_sidebands({0: 1.0, 1: 0.0}, NU0, NUR)
    assert not bad.ok and "two sidebands" in bad.message


def _synthetic_manifold(zeta, eta, centre=-12.0, fwhm=1.5):
    """Gaussians of known integral at the sideband positions plus a floor."""
    x = np.linspace(-200.0, 200.0, 6001)
    I = hb.sideband_intensities(zeta, eta, NU0, NUR, n_max=6)
    pos = hb.sideband_positions_ppm(centre, NU0, NUR, I)
    sig = fwhm / (2 * np.sqrt(2 * np.log(2)))
    y = np.full_like(x, 0.02)
    for n, inten in I.items():
        y += inten / (sig * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - pos[n]) / sig) ** 2)
    return x, y, I


def test_measure_manifold_then_fit_reads_the_tensor_off_a_spectrum():
    x, y, I = _synthetic_manifold(60.0, 0.4)
    raw = hb.measure_manifold(x, y, -12.0, NU0, NUR, 4)
    tot = sum(raw.values())
    sub = sum(I[k] for k in range(-4, 5))
    for n in range(-4, 5):
        assert raw[n] / tot == pytest.approx(I[n] / sub, abs=2e-3)
    out = hb.fit_sidebands(raw, NU0, NUR)
    assert out.zeta_ppm == pytest.approx(60.0, abs=1.5)
    assert out.eta == pytest.approx(0.4, abs=0.05)


def test_csa_conventions_round_trip():
    d11, d22, d33 = C.csa_principal_from_haeberlen(10.0, 50.0, 0.3)
    assert d11 >= d22 >= d33
    assert (d11 + d22 + d33) / 3 == pytest.approx(10.0)
    assert min(d11, d22, d33) == pytest.approx(10.0 - 50.0)     # delta_zz = iso - zeta
    diso, zeta, eta = C.csa_haeberlen_from_principal(d11, d22, d33)
    assert (diso, zeta, eta) == pytest.approx((10.0, 50.0, 0.3))
    span, skew = C.csa_span_skew(d11, d22, d33)
    assert span == pytest.approx(d11 - d33) and -1.0 <= skew <= 1.0
    back = C.csa_principal_from_span_skew(diso, span, skew)
    assert back == pytest.approx((d11, d22, d33))
    # negative zeta: delta_zz is the most positive component
    p = C.csa_principal_from_haeberlen(0.0, -40.0, 0.0)
    assert p[0] == pytest.approx(40.0) and p[1] == pytest.approx(p[2])
    assert C.csa_haeberlen_from_principal(*p)[1] == pytest.approx(-40.0)
    # axially symmetric: skew is +-1
    assert abs(C.csa_span_skew(*p)[1]) == pytest.approx(1.0)


def test_csa_mas_sideband_count_grows_with_the_static_span():
    from larmor.models.base import SimContext
    from larmor.models.csa import _n_ssb_for

    ctx = SimContext("31P", NU0, NUR, np.linspace(-300, 300, 100))
    assert _n_ssb_for(20.0, ctx) == 8                          # floor
    assert _n_ssb_for(200.0, ctx) == 32
    assert _n_ssb_for(600.0, ctx) >= 64
    assert _n_ssb_for(2000.0, ctx) == 256                       # ceiling
    static = SimContext("31P", NU0, 0.0, np.linspace(-300, 300, 100))
    assert _n_ssb_for(200.0, static) == 8                       # floor


# ---------------------------------------------------------------- Qt
def test_dialog_measures_fits_and_seeds():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    from larmor.desktop.herzfeld_berger_dialog import HerzfeldBergerDialog

    x, y, _ = _synthetic_manifold(60.0, 0.4, centre=-12.0)
    d = HerzfeldBergerDialog(None, x, y, "31P", NU0, NUR, 0.5)
    try:
        d.centre.setValue(-12.0)
        d.n_side.setValue(4)
        assert d.table.rowCount() == 9 and len(d._regions) == 9
        d._fit()
        assert d.fit is not None and d.btnSeed.isEnabled()
        assert d.fit.zeta_ppm == pytest.approx(60.0, abs=1.5)
        assert "Ω" in d.res.text() and "δ11" in d.res.text()
        got = []
        d.seed_site.connect(lambda *a: got.append(a))
        d._seed()
        assert got and got[0][0] == pytest.approx(60.0, abs=1.5)
        assert got[0][2] == pytest.approx(-12.0)
    finally:
        d.close()


def test_convert_dialog_csa_group_links_the_three_conventions():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    from larmor.desktop.utilities import ConvertDialog

    d = ConvertDialog(None, 162.0)
    try:
        d.cs_iso.setValue(10.0); d.cs_zeta.setValue(50.0); d.cs_eta.setValue(0.3)
        d11, d22, d33 = C.csa_principal_from_haeberlen(10.0, 50.0, 0.3)
        assert d.cs_d11.value() == pytest.approx(d11, abs=0.01)
        assert d.cs_d33.value() == pytest.approx(d33, abs=0.01)
        span, skew = C.csa_span_skew(d11, d22, d33)
        assert d.cs_span.value() == pytest.approx(span, abs=0.01)
        assert d.cs_skew.value() == pytest.approx(skew, abs=0.002)
        # editing a principal component drives the Haeberlen row back
        d.cs_d11.setValue(d11 + 10.0)
        assert d.cs_iso.value() == pytest.approx(10.0 + 10.0 / 3, abs=0.01)
    finally:
        d.close()


def test_hb_menu_and_seed_in_the_app():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    os.environ["LARMOR_NO_SESSION"] = "1"
    from PySide6.QtWidgets import QApplication, QMenu
    QApplication.instance() or QApplication([])
    from larmor.desktop.app import MainWindow

    win = MainWindow()
    try:
        labels = [a.text() for m in win.menuBar().findChildren(QMenu)
                  for a in m.actions()]
        assert any("Herzfeld" in t for t in labels)
        x, y, _ = _synthetic_manifold(60.0, 0.4)
        win._display_1d(x, y, "31P", NU0, NUR, "t", "x")
        win._hb_seed(60.0, 0.4, -12.0, 1.0)
        s = win.recipe["sites"][-1]
        assert s["model"] == "csa_mas"
        assert s["params"]["zeta_ppm"]["value"] == 60.0
        assert s["params"]["eta"]["value"] == 0.4
        assert s["params"]["isotropic_chemical_shift_ppm"]["value"] == -12.0
    finally:
        win.close()
