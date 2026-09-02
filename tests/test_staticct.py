"""The static-CT singularity reader: the forward map must agree with both
the analytic centroid and an actual mrsimulator static simulation, and the
inversion must recover known parameters from feature positions alone."""
import os

import numpy as np
import pytest

os.environ.setdefault("LARMOR_NO_SESSION", "1")

from larmor.staticct import (StaticReading, ct_static_features, read_cq_eta,
                             shift_ppm)


def test_powder_average_matches_canonical_delta2():
    """⟨ν(θ,φ)⟩ over the sphere must equal convert.ct_second_order_shift_ppm
    — the same δ₂ every other module uses. This pins the coefficient set."""
    from larmor.convert import ct_second_order_shift_ppm, pq_from_cq_eta

    rng = np.random.default_rng(7)
    n = 400_000
    # uniform over the sphere: cosθ uniform in [-1, 1]
    theta = np.arccos(rng.uniform(-1.0, 1.0, n))
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    for cq, eta, spin, lar in ((4.0, 0.3, 2.5, 130.32),
                               (30.0, 0.7, 1.5, 216.0)):
        mean = float(np.mean(shift_ppm(theta, phi, cq, eta, spin, lar)))
        want = ct_second_order_shift_ppm(pq_from_cq_eta(cq, eta), spin, lar)
        assert mean == pytest.approx(want, rel=0.01), (cq, eta)


@pytest.mark.slow
def test_features_match_mrsimulator_static_pattern():
    """The predicted horns must sit on the simulated pattern's maxima and
    the predicted edges must bracket its support."""
    from larmor import engine
    from larmor.recipe import Param, Recipe, SiteModel

    cq, eta, spin, lar = 30.0, 0.5, 1.5, 216.0            # an 81Br-like case
    site = SiteModel(model="quad_ct", label="s", params={
        "isotropic_chemical_shift_ppm": Param(0.0), "Cq_MHz": Param(cq),
        "eta": Param(eta), "shift_fwhm_ppm": Param(8.0),
        "amplitude": Param(1.0)})
    r = Recipe(nucleus="81Br", larmor_frequency_MHz=lar, spin_rate_Hz=0.0,
               sites=[site])
    x = np.linspace(-5000, 3000, 16001)
    gx, y, _ = engine.simulate(r, exp_ppm=x)
    y = np.clip(y, 0, None)

    f = ct_static_features(cq, eta, spin, lar)
    assert len(f.horns) == 2
    # each predicted horn lands within a linewidth of a REAL tall feature
    # (these are the two tallest maxima of the ideal lineshape)
    for h in f.horns:
        i = int(np.argmin(np.abs(gx - h)))
        win = y[max(0, i - 40): i + 40]
        assert win.max() > 0.5 * y.max(), h
    # support: (almost) all intensity lives inside the predicted edges
    inside = (gx >= f.edges[0] - 25) & (gx <= f.edges[1] + 25)
    assert y[~inside].max() < 0.02 * y.max()
    total = np.trapezoid(y, gx)
    assert np.trapezoid(y[inside], gx[inside]) > 0.98 * total


def test_read_cq_eta_round_trip():
    """Features generated from known (Cq, η, δiso) must invert back."""
    spin, lar = 1.5, 216.0
    for cq, eta, diso in ((30.0, 0.5, -300.0), (12.0, 0.2, 150.0),
                          (45.0, 0.85, 0.0)):
        f = ct_static_features(cq, eta, spin, lar)
        assert len(f.horns) == 2, (cq, eta)
        span = f.edges[1] - f.edges[0]
        # the informative edge: the one NOT sitting on a horn
        dists = [min(abs(e - h) for h in f.horns) for e in f.edges]
        edge = f.edges[int(np.argmax(dists))]
        out = read_cq_eta(f.horns[0] + diso, f.horns[1] + diso,
                          edge + diso, spin, lar)
        assert out.ok, out.message
        assert out.Cq_MHz == pytest.approx(cq, rel=0.03), (cq, eta)
        assert out.eta == pytest.approx(eta, abs=0.03), (cq, eta)
        assert out.delta_iso_ppm == pytest.approx(diso, abs=0.02 * span)
        assert out.Pq_MHz == pytest.approx(
            cq * np.sqrt(1 + eta * eta / 3), rel=0.03)
        assert out.resid_ppm < 0.03 * span
        # the degenerate (horn-coincident) edge is refused with guidance
        bad = f.edges[int(np.argmin(dists))]
        if min(abs(bad - h) for h in f.horns) < 0.05 * (f.horns[1] - f.horns[0]):
            out2 = read_cq_eta(f.horns[0] + diso, f.horns[1] + diso,
                               bad + diso, spin, lar)
            assert not out2.ok and "opposite" in out2.message


def test_read_cq_eta_rejects_nonsense():
    out = read_cq_eta(10.0, 10.0, 50.0, 2.5, 130.32)
    assert isinstance(out, StaticReading) and not out.ok
    out = read_cq_eta(0.0, 100.0, 200.0, 0.5, 130.32)   # spin-1/2
    assert not out.ok
    out = read_cq_eta(0.0, 100.0, 200.0, 2.5, 0.0)      # no field
    assert not out.ok
