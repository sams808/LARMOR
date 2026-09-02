"""VOCS (frequency-stepped) stitching: sub-spectra at stepped offsets must
reassemble into the true wide pattern, without double-counting overlaps."""
import numpy as np
import pytest

from larmor.vocs import stitch


def _true_pattern(x):
    """A broad asymmetric two-horn pattern, roughly a static CT."""
    return (np.exp(-0.5 * ((x + 350) / 120) ** 2)
            + 1.4 * np.exp(-0.5 * ((x - 150) / 60) ** 2))


def _acquire(centre, half_width, n=801, roll_frac=0.25):
    """Simulate one offset acquisition: the true pattern inside a window
    around the carrier, with a cosine excitation roll-off at the edges."""
    x = np.linspace(centre - half_width, centre + half_width, n)
    y = _true_pattern(x)
    edge = roll_frac * half_width
    for sign in (-1, 1):
        border = centre + sign * half_width
        d = np.abs(x - border)
        m = d < edge
        y[m] *= np.sin(0.5 * np.pi * d[m] / edge) ** 2
    return x, y


#: offsets 350 ppm apart with 300 ppm half-width: 250 ppm of raw overlap.
#: roll-off eats the outer 10% of each piece, and the stitch tests trim 15%
#: -- trimming past the roll-off is exactly what a real average needs, and
#: 15% still leaves ~70 ppm of clean overlap per junction.
PIECES = [_acquire(c, 300.0, roll_frac=0.10) for c in (-500.0, -150.0, 200.0)]


def test_skyline_recovers_the_true_pattern():
    res = stitch(PIECES, method="skyline", trim_frac=0.15)
    truth = _true_pattern(res.ppm)
    covered = res.coverage > 0
    err = np.max(np.abs(res.amp[covered] - truth[covered]))
    assert err < 0.02 * truth.max(), err
    # ascending, spanning the union of the trimmed pieces
    assert np.all(np.diff(res.ppm) > 0)
    assert res.ppm[0] >= -710 and res.ppm[-1] <= 410


def test_average_does_not_double_count_overlaps():
    res = stitch(PIECES, method="average", trim_frac=0.15)
    truth = _true_pattern(res.ppm)
    covered = res.coverage > 0
    # overlap regions are covered by 2 pieces; the mean must stay on the
    # truth (a plain sum would sit at 2x there)
    assert np.any(res.coverage >= 2)
    err = np.max(np.abs(res.amp[covered] - truth[covered]))
    assert err < 0.02 * truth.max(), err


def test_gap_between_pieces_is_warned_not_invented():
    far = [_acquire(-500.0, 200.0), _acquire(500.0, 200.0)]
    res = stitch(far, method="skyline", trim_frac=0.1)
    assert any("gap" in n.lower() for n in res.notes)
    mid = (res.ppm > -200) & (res.ppm < 200)
    assert np.all(res.amp[mid] == 0.0) and np.all(res.coverage[mid] == 0)


def test_stitch_validates_inputs():
    with pytest.raises(ValueError, match="no spectra"):
        stitch([])
    with pytest.raises(ValueError, match="unknown stitch method"):
        stitch(PIECES, method="sum")
    with pytest.raises(ValueError, match="matching 1D"):
        stitch([(np.arange(5.0), np.arange(4.0))])


def test_finest_piece_sets_the_grid():
    coarse = _acquire(-200.0, 300.0, n=101)
    fine = _acquire(200.0, 300.0, n=3001)
    res = stitch([coarse, fine], trim_frac=0.0)
    dx = np.median(np.diff(res.ppm))
    assert dx <= np.median(np.diff(fine[0])) * 1.05
