"""larmor.sidebands -- the +-nu_rot repeat found in the data itself.

Synthetic manifolds come from the registered ``sidebands`` model through
engine.simulate (spacing = spin_rate_Hz / larmor_MHz, geometric ratio), so the
truth is known to the digit; the counter-examples (static pattern, singlet,
two lines, J-multiplet, noise) pin what must NOT trigger. The two in-repo
example EXPNOs run as consistency tests that never skip; one real 27Al 1r
smoke skips cleanly when the local data tree is absent.
"""
from pathlib import Path

import numpy as np
import pytest

from conftest import BRUKER_1R, require
from larmor import engine
from larmor import sidebands as sb
from larmor.recipe import Param, Recipe, SiteModel

REPO = Path(__file__).resolve().parents[1]
X = np.linspace(-50.0, 50.0, 4000)          # 1H at 100 MHz: 1 kHz = 10 ppm


def _manifold(x=X, larmor=100.0, nu=1000.0, ratio=0.4, n=3, fwhm=1.0,
              centre=0.0):
    """y of a ``sidebands``-model manifold on x (the test_models_v2 idiom)."""
    r = Recipe(nucleus="1H", larmor_frequency_MHz=larmor, spin_rate_Hz=nu,
               sites=[SiteModel(model="sidebands", label="S", params={
                   "isotropic_chemical_shift_ppm": Param(centre),
                   "shift_fwhm_ppm": Param(fwhm), "amplitude": Param(1.0),
                   "ssb_ratio": Param(ratio), "n_ssb": Param(float(n)),
                   "gl": Param(1.0)})])
    xc, y, _ = engine.simulate(r, exp_ppm=np.asarray(x, float))
    assert np.allclose(xc, x)
    return np.asarray(y, float)


def _gauss(x, centre, fwhm, height):
    return height * np.exp(-4.0 * np.log(2.0) * ((x - centre) / fwhm) ** 2)


def _matched_k(det):
    return sorted(o.k for o in det.orders if o.matched)


# ---------------------------------------------------------------- profile
def test_correlation_profile_is_normalised_and_matches_direct_shift():
    y = _manifold()
    lags, c = sb.correlation_profile(X, y)
    assert lags[0] == 0.0 and np.all(np.diff(lags) > 0)
    assert c[0] == pytest.approx(1.0, abs=1e-9)
    assert np.all(np.abs(c) <= 1.0 + 1e-9)
    assert lags.max() <= 50.01                          # cut at half the span
    fft_at_10 = float(np.interp(10.0, lags, c))
    assert fft_at_10 == pytest.approx(sb.shift_correlation(X, y, 10.0), abs=0.01)
    assert fft_at_10 > 0.2
    i = int(np.argmin(np.abs(lags - 10.0)))
    assert c[i] >= c[i - 1] and c[i] >= c[i + 1]        # a local maximum
    # measured on this manifold (ratio 0.4, n 3): 0.688 / 0.386 / 0.188
    assert fft_at_10 == pytest.approx(0.688, abs=0.02)
    assert float(np.interp(20.0, lags, c)) == pytest.approx(0.386, abs=0.02)
    assert float(np.interp(30.0, lags, c)) == pytest.approx(0.188, abs=0.02)


# ----------------------------------------------------------------- detect
def test_detects_manifold_at_the_recorded_rate():
    y = _manifold()
    det = sb.detect(X, y, 100.0, nu_rot_Hz=1000.0)
    assert det.ok and det.from_meta and not det.scanned
    assert det.nu_rot_Hz == pytest.approx(1000.0, rel=2e-3)
    assert det.spacing_ppm == pytest.approx(10.0, rel=2e-3)
    assert det.centre_ppm == pytest.approx(0.0, abs=0.1)
    assert det.centre_fwhm_ppm == pytest.approx(1.0, abs=0.1)
    assert det.centre_height == pytest.approx(1.0, abs=0.02)
    assert _matched_k(det) == [-3, -2, -1, 1, 2, 3]
    assert det.n_orders == 3
    for o in det.matched():
        assert o.height / det.centre_height == pytest.approx(0.4 ** abs(o.k),
                                                             abs=0.05)
        assert o.ppm == pytest.approx(o.k * 10.0, abs=0.05)
    assert 0.0 < det.score <= 1.0
    assert det.confidence >= 0.7
    assert det.sides() == (1, -1)
    assert [o.k for o in det.matched(1)] == [1, 2, 3]
    assert [o.k for o in det.matched(-1)] == [-1, -2, -3]
    # orders are reported |k| ascending, both signs
    assert [abs(o.k) for o in det.orders] == sorted(abs(o.k) for o in det.orders)
    assert "1 000" in det.message


def test_recovers_a_rate_off_by_1p5_percent():
    """The +-2 % window re-centres on the DATA rate, not the candidate."""
    y = _manifold()
    for guess in (1015.0, 985.0):
        det = sb.detect(X, y, 100.0, nu_rot_Hz=guess)
        assert det.ok and det.from_meta, det.message
        assert det.nu_rot_Hz == pytest.approx(1000.0, rel=3e-3)


def test_far_off_rate_needs_the_scan():
    y = _manifold()
    det = sb.detect(X, y, 100.0, nu_rot_Hz=1100.0, scan=False)
    assert not det.ok and "2 %" in det.message and "1 100" in det.message
    det = sb.detect(X, y, 100.0, nu_rot_Hz=1100.0, scan=True)
    assert det.ok and det.scanned and not det.from_meta
    assert det.nu_rot_Hz == pytest.approx(1000.0, rel=3e-3)
    assert det.confidence == pytest.approx(sb.SCAN_PENALTY, abs=1e-6)


def test_scan_finds_an_unknown_rate_and_the_fundamental():
    det = sb.detect(X, _manifold(), 100.0, nu_rot_Hz=0.0, scan=True)
    assert det.ok and det.scanned
    assert det.nu_rot_Hz == pytest.approx(1000.0, rel=3e-3)
    # strong higher orders: c(2d) is comparable to c(d) -- the smallest
    # candidate within 90 % of the best prominence is the fundamental
    det = sb.detect(X, _manifold(ratio=0.9), 100.0, nu_rot_Hz=0.0, scan=True)
    assert det.ok
    assert det.spacing_ppm == pytest.approx(10.0, rel=3e-3)
    assert det.prominence > 0.8
    det = sb.detect(X, _manifold(n=4), 100.0, nu_rot_Hz=0.0, scan=True)
    assert det.ok and det.n_orders == 4
    assert _matched_k(det) == [-4, -3, -2, -1, 1, 2, 3, 4]


def test_static_rate_zero_returns_not_ok_without_scan():
    """Auto mode never scans a certain static recipe, even on data that DOES
    repeat."""
    det = sb.detect(X, _manifold(), 100.0, nu_rot_Hz=0.0, scan=False)
    assert not det.ok and "static" in det.message
    assert det.orders == () and det.matched() == [] and det.sides() == ()


def test_static_pattern_yields_nothing():
    """An 81Br static second-order CT pattern (the test_staticct recipe) with
    1 % noise: a two-horn pattern correlates with itself along a broad lobe,
    never at an interior maximum, and its horns have no symmetric partner."""
    site = SiteModel(model="quad_ct", label="s", params={
        "isotropic_chemical_shift_ppm": Param(0.0), "Cq_MHz": Param(30.0),
        "eta": Param(0.5), "shift_fwhm_ppm": Param(8.0),
        "amplitude": Param(1.0)})
    r = Recipe(nucleus="81Br", larmor_frequency_MHz=216.0, spin_rate_Hz=0.0,
               sites=[site])
    x = np.linspace(-5000.0, 3000.0, 6001)
    _, y, _ = engine.simulate(r, exp_ppm=x)
    y = y + np.random.default_rng(3).normal(0.0, 0.01 * y.max(), x.size)
    for nu in (0.0, 14000.0):
        det = sb.detect(x, y, 216.0, nu_rot_Hz=nu, scan=True)
        assert not det.ok and det.message
        assert det.confidence < sb.MIN_CONFIDENCE or det.sides() != (1, -1)
    det = sb.detect(x, y, 216.0, nu_rot_Hz=14000.0, scan=False)
    assert not det.ok and "2 %" in det.message


def test_singlet_two_lines_noise_and_jmultiplet_do_not_trigger():
    # (a) a lone line
    ys = _gauss(X, 0.0, 1.0, 1.0)
    assert not sb.detect(X, ys, 100.0, nu_rot_Hz=1000.0).ok
    assert not sb.detect(X, ys, 100.0, nu_rot_Hz=0.0, scan=True).ok
    # (b) two unrelated lines 12 ppm apart: the even autocorrelation peaks at
    # 1200 Hz (c = a b / (a^2 + b^2) = 0.40) -- only the symmetric +-1
    # criterion tells this from a manifold
    y2 = _gauss(X, 0.0, 1.0, 1.0) + _gauss(X, 12.0, 1.0, 0.5)
    det = sb.detect(X, y2, 100.0, nu_rot_Hz=1200.0)
    assert not det.ok and det.score == pytest.approx(0.40, abs=0.02)
    assert _matched_k(det) == [1] and det.sides() == (1,)
    assert "±1" in det.message and "−1" in det.message
    assert not sb.detect(X, y2, 100.0, nu_rot_Hz=0.0, scan=True).ok
    # (c) pure noise, several seeds
    for seed in range(6):
        yn = np.random.default_rng(seed).normal(0.0, 1.0, X.size)
        assert not sb.detect(X, yn, 100.0, nu_rot_Hz=1000.0).ok
        assert not sb.detect(X, yn, 100.0, nu_rot_Hz=0.0, scan=True).ok
    # (d) a 1:3:3:1 quartet (J 1 kHz) recorded at 20 kHz: the +-2 % window
    # around the recorded rate never sees the J comb
    from larmor import models as model_registry
    pj = {p.name: Param(p.default)
          for p in model_registry.get("jmultiplet").params}
    pj["j_hz"] = Param(1000.0)
    pj["n_j"] = Param(3.0)
    pj["shift_fwhm_ppm"] = Param(1.0)
    rj = Recipe(nucleus="31P", larmor_frequency_MHz=160.46, spin_rate_Hz=20000.0,
                sites=[SiteModel(model="jmultiplet", label="J", params=pj)])
    xj = np.linspace(-200.0, 200.0, 8000)
    _, yj, _ = engine.simulate(rj, exp_ppm=xj)
    det = sb.detect(xj, yj, 160.46, nu_rot_Hz=20000.0, scan=False)
    assert not det.ok and "2 %" in det.message


def test_noise_and_baseline_offset_are_tolerated():
    rng = np.random.default_rng(1)
    y = _manifold() + 0.3 + rng.normal(0.0, 0.02, X.size)
    det = sb.detect(X, y, 100.0, nu_rot_Hz=1000.0)
    assert det.ok, det.message
    assert det.nu_rot_Hz == pytest.approx(1000.0, rel=5e-3)
    assert {-2, -1, 1, 2} <= set(_matched_k(det))
    assert det.centre_height == pytest.approx(1.0, abs=0.1)   # floor removed
    assert det.confidence >= 0.9                              # measured 1.00
    # heavier noise still confirms the recorded rate
    y = _manifold() + rng.normal(0.0, 0.1, X.size)
    det = sb.detect(X, y, 100.0, nu_rot_Hz=1000.0)
    assert det.ok and det.nu_rot_Hz == pytest.approx(1000.0, rel=0.02)


def test_descending_and_nonuniform_axes():
    y = _manifold()
    ref = sb.detect(X, y, 100.0, nu_rot_Hz=1000.0)
    desc = sb.detect(X[::-1], y[::-1], 100.0, nu_rot_Hz=1000.0)
    assert desc.ok
    assert desc.nu_rot_Hz == pytest.approx(ref.nu_rot_Hz, rel=1e-6)
    assert _matched_k(desc) == _matched_k(ref)
    u = np.linspace(-1.0, 1.0, 4000)
    xn = 50.0 * np.sign(u) * np.abs(u) ** 1.03          # monotone, non-uniform
    yn = _manifold(x=xn)
    det = sb.detect(xn, yn, 100.0, nu_rot_Hz=1000.0)
    assert det.ok and det.nu_rot_Hz == pytest.approx(1000.0, rel=0.01)
    # oversized traces are block-averaged, not point-sampled
    xb = np.linspace(-50.0, 50.0, 40000)
    det = sb.detect(xb, _manifold(x=xb), 100.0, nu_rot_Hz=1000.0)
    assert det.ok and det.nu_rot_Hz == pytest.approx(1000.0, rel=3e-3)
    assert _matched_k(det) == [-3, -2, -1, 1, 2, 3]


def test_recentre_moves_orders_but_not_the_rate():
    y = _manifold()
    ref = sb.detect(X, y, 100.0, nu_rot_Hz=1000.0)
    det = sb.detect(X, y, 100.0, nu_rot_Hz=1000.0, centre_ppm=10.0)
    assert det.ok
    assert det.nu_rot_Hz == pytest.approx(ref.nu_rot_Hz, rel=1e-6)
    assert det.centre_ppm == pytest.approx(10.0, abs=0.1)
    minus_one = next(o for o in det.orders if o.k == -1)
    assert minus_one.matched and minus_one.ppm == pytest.approx(0.0, abs=0.1)
    assert _matched_k(det) == [-4, -3, -2, -1, 1, 2]


def test_degenerate_input_never_raises():
    y = _manifold()
    assert not sb.detect(X[:30], y[:30], 100.0, 1000.0).ok
    assert "Larmor" in sb.detect(X, y, 0.0, 1000.0).message
    assert "flat" in sb.detect(X, np.zeros_like(y), 100.0, 1000.0).message
    assert "flat" in sb.detect(X, np.full_like(y, 3.0), 100.0, 1000.0).message
    det = sb.detect(X, y, 100.0, nu_rot_Hz=8000.0)          # 80 ppm > half span
    assert not det.ok and "half the axis" in det.message
    det = sb.detect(X, y, 100.0, nu_rot_Hz=8000.0, scan=True)
    assert det.ok and det.scanned                            # falls back to a scan
    yn = y.copy()
    yn[::7] = np.nan
    assert sb.detect(X, yn, 100.0, 1000.0).ok                # non-finite dropped
    assert not sb.detect(X, y, 100.0, 1000.0, min_confidence=1.01).ok
    assert not sb.detect(["a"], ["b"], 100.0, 1000.0).ok


# ------------------------------------------------------- seeds and wording
def test_seed_ratio_and_describe():
    det = sb.detect(X, _manifold(), 100.0, nu_rot_Hz=1000.0)
    assert sb.seed_ratio(det) == pytest.approx(0.4, abs=0.05)
    one = sb.SidebandDetection(
        True, "", nu_rot_Hz=1000.0, spacing_ppm=10.0, centre_height=1.0,
        orders=(sb.SidebandOrder(1, 10.0, 5.0, True),
                sb.SidebandOrder(-1, -10.0, 0.001, True)))
    assert 0.02 <= sb.seed_ratio(one) <= 0.98
    assert sb.seed_ratio(sb.SidebandDetection(False, "x")) == 0.3   # nothing
    text = sb.describe(det)
    assert text.startswith("Spinning sidebands: νrot 1 000 Hz")
    assert "10.0 ppm" in text and "100.00 MHz" in text
    assert "3 orders each side" in text and "confidence 1.00" in text
    assert "⚠" not in text and "⚠" not in sb.describe(det, 1000.0)
    warned = sb.describe(det, recipe_rate_Hz=1100.0)
    assert "⚠" in warned and ("1 100" in warned or "1100" in warned)
    assert "⚠" not in sb.describe(det, 1004.0)          # inside RATE_MISMATCH
    assert sb.format_hz(20000) == "20 000" and sb.format_hz(950.4) == "950"


# ------------------------------------------------------------- real data
def test_example_27Al_26kHz_is_consistent_with_acqus():
    """examples/pCABS2-4/3616: 27Al MAS at 26 kHz (acqus and title agree).
    The visible sidebands are the SATELLITE-transition manifold, displaced
    from the central-transition peak by the satellites' own second-order
    shift (about +10.6 ppm here) and only ~6 % tall, so the CT-anchored
    comb has no +-1 pair exactly at +-nu_rot: the detector must NOT claim
    a repeat at the recorded rate, and must say why."""
    from larmor import loader

    ppm, amp, rec, _meta, _warn = loader.load_any(REPO / "examples/pCABS2-4/3616")
    assert rec["spin_rate_Hz"] == 26000.0 and not rec.get("mas_uncertain")
    det = sb.detect(ppm, amp, 130.323, nu_rot_Hz=26000.0)
    assert isinstance(det, sb.SidebandDetection) and isinstance(det.message, str)
    if det.ok:                                     # pin only measured truths
        assert det.nu_rot_Hz == pytest.approx(26000.0, rel=0.02)
        for o in det.matched():
            assert ppm.min() <= o.ppm <= ppm.max()
    else:
        assert "2 %" in det.message or "±1" in det.message
    # the scan never invents a rate far from the acquisition's
    scan = sb.detect(ppm, amp, 130.323, nu_rot_Hz=0.0, scan=True)
    assert not scan.ok or scan.nu_rot_Hz == pytest.approx(26000.0, rel=0.06)


def test_example_11B_never_invents_a_third_rate():
    """examples/pCABS2-4/1118: acqus MASR 30 000 vs title 'MAS=20 kHz'; the
    loader keeps 30 000 flagged uncertain ('highest wins'). At 30 kHz the
    spacing (187 ppm) exceeds half the 371 ppm axis; the scan's repeat sits
    on the -1 tooth at -124.6 ppm, i.e. the operator-typed 20 kHz, but the
    +1 side is a negative-going lobe so nothing is confirmed. Whatever the
    verdict, the rate it reports must be one of the two sources' rates."""
    from larmor import loader

    ppm, amp, rec, _meta, _warn = loader.load_any(REPO / "examples/pCABS2-4/1118")
    assert rec["spin_rate_Hz"] == 30000.0 and rec.get("mas_uncertain") is True
    det = sb.detect(ppm, amp, 160.4616, nu_rot_Hz=0.0, scan=True)
    assert isinstance(det.message, str) and det.message
    if det.nu_rot_Hz > 0:
        assert (det.nu_rot_Hz == pytest.approx(20000.0, rel=0.02)
                or det.nu_rot_Hz == pytest.approx(30000.0, rel=0.02))
    at30 = sb.detect(ppm, amp, 160.4616, nu_rot_Hz=30000.0)
    assert not at30.ok and "half the axis" in at30.message


def test_real_27Al_1r_smoke():
    from larmor import loader

    ppm, amp, rec, _meta, _warn = loader.load_any(str(require(BRUKER_1R)))
    det = sb.detect(ppm, amp, rec["larmor_frequency_MHz"],
                    nu_rot_Hz=rec.get("spin_rate_Hz", 0.0), scan=False)
    assert isinstance(det, sb.SidebandDetection)
    if det.ok:
        assert det.sides() == (1, -1)
        assert det.confidence >= sb.MIN_CONFIDENCE
    text = sb.describe(det) if det.spacing_ppm else det.message
    assert isinstance(text, str) and text


# ---------------------------------------------------------------- follow νrot
# WA (2026-09): "make sure that if there are bands identified as ssb, they
# move if I change the vrot used for fitting". A linked copy carries
# site["sideband"] = {"parent": i, "k": order}; refresh_linked recomputes its
# position constraint from the recipe's spin_rate_Hz.

def _linked_sites(spacing=153.465):
    def site(pos, expr=None, mark=None):
        d = {"model": "gauss_lor", "label": "x", "params": {
            "isotropic_chemical_shift_ppm": {"value": pos, "stderr": None,
                                             "vary": True, "min": None,
                                             "max": None, "expr": expr},
            "amplitude": {"value": 1.0, "stderr": None, "vary": True,
                          "min": 0.0, "max": None, "expr": None}}}
        if mark:
            d["sideband"] = dict(mark)
        return d
    return [
        site(60.0),
        site(60.0 + spacing, sb.linked_position_expr(0, 1, spacing), {"parent": 0, "k": 1}),
        site(60.0 - 2 * spacing, sb.linked_position_expr(0, -2, spacing), {"parent": 0, "k": -2}),
        # made before the marker existed: a constant offset that stays put
        site(60.0 - spacing, "s0.isotropic_chemical_shift_ppm - 153.465"),
        # unlinked by hand since: the marker is dropped, the value kept
        site(5.0, None, {"parent": 0, "k": 3}),
    ]


def test_linked_position_expr_is_the_table_form():
    assert sb.linked_position_expr(0, 1, 124.642) == \
        "s0.isotropic_chemical_shift_ppm + 124.642"
    assert sb.linked_position_expr(2, -2, 124.642) == \
        "s2.isotropic_chemical_shift_ppm - 249.284"
    from larmor import cellparse
    assert cellparse.format_link(sb.linked_position_expr(0, 1, 124.642),
                                 "isotropic_chemical_shift_ppm") == "A+124.642"


def test_refresh_linked_moves_marked_copies_to_the_new_rate():
    lar = 130.3
    sites = _linked_sites(20000.0 / lar)
    moved = sb.refresh_linked(sites, lar, 22000.0)
    assert moved == 2
    d = 22000.0 / lar
    p1 = sites[1]["params"]["isotropic_chemical_shift_ppm"]
    p2 = sites[2]["params"]["isotropic_chemical_shift_ppm"]
    assert p1["expr"] == sb.linked_position_expr(0, 1, d)
    assert p2["expr"] == sb.linked_position_expr(0, -2, d)
    assert p1["value"] == pytest.approx(60.0 + d)
    assert p2["value"] == pytest.approx(60.0 - 2 * d)
    # the pre-marker constant offset is untouched
    assert sites[3]["params"]["isotropic_chemical_shift_ppm"]["expr"] == \
        "s0.isotropic_chemical_shift_ppm - 153.465"
    # the hand-unlinked copy lost its marker and kept its value
    assert "sideband" not in sites[4]
    assert sites[4]["params"]["isotropic_chemical_shift_ppm"]["value"] == 5.0
    # the same rate again moves nothing; a static / unknown rate moves nothing
    assert sb.refresh_linked(sites, lar, 22000.0) == 0
    assert sb.refresh_linked(sites, lar, 0.0) == 0
    assert sb.refresh_linked(sites, 0.0, 22000.0) == 0
    # the constraint evaluates in the fit engine
    from larmor import fit as fitmod
    from larmor.recipe import Recipe
    rec = Recipe.from_dict({"nucleus": "27Al", "larmor_frequency_MHz": lar,
                            "spin_rate_Hz": 22000.0, "sites": sites})
    params = fitmod._make_params(rec)
    assert params["s1_pos"].value == pytest.approx(60.0 + d)
    assert params["s2_pos"].value == pytest.approx(60.0 - 2 * d)


def test_sideband_marker_survives_the_recipe_round_trip_and_site_edits():
    from larmor.constraints_util import (remap_exprs_after_delete,
                                         remap_exprs_after_move)
    from larmor.recipe import Recipe

    sites = _linked_sites()
    rec = Recipe.from_dict({"nucleus": "27Al", "larmor_frequency_MHz": 130.3,
                            "spin_rate_Hz": 20000.0, "sites": sites})
    assert rec.sites[1].sideband == {"parent": 0, "k": 1}
    assert rec.sites[0].sideband is None
    d = rec.to_dict()
    assert d["sites"][1]["sideband"] == {"parent": 0, "k": 1}
    assert "sideband" not in d["sites"][0]         # omitted when None
    assert not any("unknown site fields" in n for n in rec.notes)

    # delete the line before the parent: the marker's parent index shifts
    sites = [{"model": "gauss_lor", "label": "first", "params": {
        "isotropic_chemical_shift_ppm": {"value": 0.0, "expr": None}}}] + _linked_sites()
    for s in sites[1:]:
        p = s["params"]["isotropic_chemical_shift_ppm"]
        if p.get("expr"):
            p["expr"] = p["expr"].replace("s0.", "s1.")
        if s.get("sideband"):
            s["sideband"]["parent"] = 1
    sites.pop(0)                                # the app removes, then remaps
    remap_exprs_after_delete(sites, 0)
    assert sites[1]["sideband"] == {"parent": 0, "k": 1}
    assert sites[2]["sideband"] == {"parent": 0, "k": -2}
    assert sites[1]["params"]["isotropic_chemical_shift_ppm"]["expr"].startswith("s0.")
    # delete the parent itself: markers (and the constraints) go
    sites = _linked_sites()
    sites.pop(0)
    remap_exprs_after_delete(sites, 0)
    assert all("sideband" not in s for s in sites)
    assert all(not s["params"]["isotropic_chemical_shift_ppm"].get("expr")
               for s in sites)
    # a reorder keeps the marker pointing at the parent
    sites = _linked_sites()
    order = [1, 0, 2, 3, 4]                     # parent moves to index 1
    sites = [sites[i] for i in order]
    remap_exprs_after_move(sites, {old: new for new, old in enumerate(order)})
    assert sites[0]["sideband"] == {"parent": 1, "k": 1}
    assert sites[2]["sideband"] == {"parent": 1, "k": -2}
