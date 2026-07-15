import numpy as np
import pytest

from larmor import engine
from larmor.io import fxmla
from larmor import fit as fitmod

from conftest import CAALGLASS, require


def test_gauss_lor_shapes():
    x = np.linspace(-50, 50, 2001)
    g = engine.gauss_lor(x, 0.0, 10.0, 2.0, gl=1.0)
    l = engine.gauss_lor(x, 0.0, 10.0, 2.0, gl=0.0)
    assert g.max() == pytest.approx(2.0, rel=1e-6)
    assert l.max() == pytest.approx(2.0, rel=1e-6)
    # both peak-normalized with the same FWHM...
    half = np.where(g >= 1.0)[0]
    assert x[half[-1]] - x[half[0]] == pytest.approx(10.0, abs=0.2)
    # ...but the Lorentzian has heavier tails
    tail = np.abs(x) > 25
    assert l[tail].sum() > 5 * g[tail].sum()


@pytest.mark.slow
def test_kernel_and_fit_caalglass():
    """End-to-end: import the dmfit fit, refine with lmfit, beat Phase 0's RMSD."""
    dm = fxmla.read(require(CAALGLASS))
    recipe, _ = fxmla.to_recipe(dm)

    kernel = engine.build_kernel(
        recipe.nucleus, recipe.larmor_frequency_MHz, recipe.spin_rate_Hz)
    assert kernel.K.shape == (80 * 11, 2048)
    assert np.all(np.diff(kernel.x_ppm) > 0)

    # sigma reweighting sanity: larger sigma -> broader lineshape
    y_narrow = kernel.spectrum(1.0, 0.0, 1.0, 1.0)
    y_broad = kernel.spectrum(4.0, 0.0, 1.0, 1.0)
    width = lambda y: np.sum(y > y.max() / 2)
    assert width(y_broad) > width(y_narrow)

    exp_ppm, exp_amp = dm.spectrum.ppm, dm.spectrum.amplitude
    result = fitmod.fit(recipe, exp_ppm, exp_amp,
                        window_ppm=(150.0, -80.0), kernel=kernel)

    # refined fit should beat the fixed-parameter Phase 0 replay (RMSD 0.027)
    assert result.rmsd < 0.01

    # dmfit's ad-hoc Gauss/Lor sideband lines (236/208 ppm) sit outside the
    # window and must be frozen, not left to wander
    assert len(result.frozen_sites) == 2

    site1 = recipe.sites[0]
    assert site1.params["isotropic_chemical_shift_ppm"].value == pytest.approx(
        66.2, abs=5.0)
    assert site1.params["sigma_Cq_MHz"].value == pytest.approx(2.0, abs=0.7)
    # the whole point: uncertainties exist now
    assert result.lmfit_result.errorbars
    assert site1.params["sigma_Cq_MHz"].stderr is not None
    assert 0 < site1.params["sigma_Cq_MHz"].stderr < 0.5
