"""larmor.display -- the Qt-free arithmetic behind View > Y axis, the
Datasets dock's per-overlay scale / shift / offset and the Full view's
extents and snap-back rule."""
import numpy as np
import pytest

from larmor import display


def _spectrum(n=801):
    x = np.linspace(200.0, -100.0, n)                  # descending, NMR order
    y = 5000.0 * np.exp(-((x - 60.0) / 15.0) ** 2) + 800.0 * np.exp(
        -((x + 20.0) / 8.0) ** 2)
    return x, y


# ------------------------------------------------------------- y_factor / area
def test_area_is_positive_on_a_descending_axis_and_restricts_to_a_region():
    x, y = _spectrum()
    total = display.area(x, y)
    assert total > 0
    # the analytic Gaussian areas: A sqrt(pi) w
    assert total == pytest.approx(5000.0 * np.sqrt(np.pi) * 15.0
                                  + 800.0 * np.sqrt(np.pi) * 8.0, rel=1e-3)
    part = display.area(x, y, region=(120.0, 0.0))     # the big line only
    assert part == pytest.approx(5000.0 * np.sqrt(np.pi) * 15.0, rel=1e-3)
    assert display.area(x, y, region=(0.0, -100.0)) == pytest.approx(
        800.0 * np.sqrt(np.pi) * 8.0, rel=1e-3)
    # a region cutting a line at its centre keeps half of it
    assert display.area(x, y, region=(-20.0, -100.0)) == pytest.approx(
        0.5 * 800.0 * np.sqrt(np.pi) * 8.0, rel=2e-2)
    # region given lo-first is the same region
    assert display.area(x, y, region=(0.0, 120.0)) == part
    # degenerate inputs
    assert display.area([], []) == 0.0
    assert display.area([1.0], [2.0]) == 0.0
    assert display.area(x, y, region=(1000.0, 900.0)) == 0.0
    assert display.area([0, 1, np.nan], [1, 1, 1]) == pytest.approx(1.0)


def test_y_factor_per_mode():
    x, y = _spectrum()
    assert display.y_factor("raw", x, y) == 1.0
    assert display.y_factor("max", x, y) == pytest.approx(1.0 / y.max())
    assert (y * display.y_factor("max", x, y)).max() == pytest.approx(1.0)
    fa = display.y_factor("area", x, y)
    assert display.area(x, y * fa) == pytest.approx(1.0)
    fr = display.y_factor("region", x, y, region=(120.0, 0.0))
    assert display.area(x, y * fr, region=(120.0, 0.0)) == pytest.approx(1.0)
    # the whole-axis area is then larger than 1
    assert display.area(x, y * fr) > 1.0


def test_y_factor_falls_back_to_raw_on_degenerate_data():
    x, y = _spectrum()
    assert display.y_factor("max", x, np.zeros_like(y)) == 1.0
    assert display.y_factor("max", x, -y) == 1.0                # inverted trace
    assert display.y_factor("max", None, None) == 1.0
    assert display.y_factor("max", x, []) == 1.0
    assert display.y_factor("area", x, np.zeros_like(y)) == 1.0
    assert display.y_factor("region", x, y, region=None) == 1.0
    assert display.y_factor("region", x, y, region=(1000.0, 900.0)) == 1.0
    assert display.y_factor("nonsense", x, y) == 1.0
    assert display.y_factor("max", x, np.full_like(y, np.nan)) == 1.0


def test_labels_and_region_text():
    assert display.y_axis_label("raw") == "intensity"
    assert display.y_axis_label("max") == "intensity (normalised to max)"
    assert display.y_axis_label("area") == "intensity (normalised to area)"
    assert display.y_axis_label("region", (100.0, -20.0)) == \
        "intensity (normalised to area 100…−20 ppm)"
    assert display.y_axis_label("region", None) == "intensity"
    assert display.format_region((-20.0, 100.0)) == "100…−20 ppm"
    assert display.parse_region("100,-20") == (100.0, -20.0)
    assert display.parse_region("-20,100") == (100.0, -20.0)
    assert display.parse_region("") is None
    assert display.parse_region(None) is None
    assert display.parse_region("5,5") is None
    assert display.parse_region("a,b") is None
    assert display.parse_region("1") is None
    assert tuple(display.Y_MODE_LABELS) == display.Y_MODES


# ------------------------------------------------------------- overlays
def test_overlay_display_leaves_the_arrays_alone_and_applies_every_knob():
    ppm = np.linspace(80.0, -20.0, 501)
    amp = 40.0 * np.exp(-((ppm - 10.0) / 4.0) ** 2)
    ppm0, amp0 = ppm.copy(), amp.copy()
    x, y = display.overlay_display(ppm, amp)
    assert np.array_equal(x, ppm) and np.array_equal(y, amp)
    assert x is not ppm and y is not amp                 # never the stored arrays
    x, y = display.overlay_display(ppm, amp, scale=2.5, shift=1.2, yoff=0.3,
                                   stack=0.4, span=100.0)
    assert np.allclose(x, ppm + 1.2)
    assert np.allclose(y, amp * 2.5 + 70.0)
    assert np.array_equal(ppm, ppm0) and np.array_equal(amp, amp0)
    # own normalisation under max: reaches `scale`, offsets on top
    x, y = display.overlay_display(ppm, amp, mode="max")
    assert y.max() == pytest.approx(1.0)
    x, y = display.overlay_display(ppm, amp, mode="max", scale=0.5, yoff=0.25,
                                   span=2.0)
    assert y.max() == pytest.approx(0.5 + 0.5)
    # own area normalisation: unit area, whatever the intensity
    x, y = display.overlay_display(ppm, 1e6 * amp, mode="area")
    assert display.area(x, y) == pytest.approx(1.0)
    # a region follows the SHIFTED trace (the user lined a peak up first)
    x, y = display.overlay_display(ppm, amp, shift=30.0, mode="region",
                                   region=(50.0, 30.0))
    assert display.area(x, y, region=(50.0, 30.0)) == pytest.approx(1.0)


def test_match_scale_and_badge():
    a = np.array([0.0, 100.0, 20.0])
    b = np.array([0.0, 25.0, 5.0])
    assert display.match_scale(a, b) == pytest.approx(4.0)
    assert (b * display.match_scale(a, b)).max() == pytest.approx(a.max())
    assert display.match_scale(a, np.zeros(3)) == 1.0
    assert display.match_scale(-a, b) == 1.0
    assert display.match_scale([], b) == 1.0
    assert display.overlay_badge() == ""
    assert display.overlay_badge(1.0, 0.0, 0.0) == ""
    assert display.overlay_badge(2.5, 1.2, 0.0) == "×2.5 · +1.2 ppm"
    assert display.overlay_badge(1.0, -3.0, 0.3) == "−3 ppm · ↑0.3" \
        or display.overlay_badge(1.0, -3.0, 0.3) == "-3 ppm · ↑0.3"
    assert display.overlay_badge(1.0, 0.0, -0.5) == "↓0.5"


# ------------------------------------------------------------- the Full view
def test_full_extents_come_from_the_data_with_margin_and_residual_room():
    x, y = _spectrum()
    (xlo, xhi), (ylo, yhi) = display.full_extents(x, y)
    m = display.FULL_MARGIN_FRAC * 300.0
    assert xlo == pytest.approx(-100.0 - m) and xhi == pytest.approx(200.0 + m)
    top = float(y.max())
    lo = min(float(y.min()), -display.RESID_ROOM_FRAC * top)
    pad = display.Y_PAD_FRAC * (top - lo)
    assert yhi == pytest.approx(top + pad) and ylo == pytest.approx(lo - pad)
    # a time axis (ms, ascending) works the same way
    t = np.linspace(0.0, 50.0, 256)
    (tlo, thi), _ = display.full_extents(t, np.exp(-t / 10.0))
    assert tlo == pytest.approx(-1.0) and thi == pytest.approx(51.0)
    assert display.full_extents([], []) is None
    assert display.full_extents([np.nan], [1.0]) is None
    # NaN points are ignored, not propagated
    (xlo2, xhi2), _ = display.full_extents(np.r_[x, np.nan], np.r_[y, 1.0])
    assert (xlo2, xhi2) == (xlo, xhi)


def test_snap_back_rule():
    data = (-100.0, 200.0)                       # span 300
    assert display.snap_back((250.0, 400.0), data)         # entirely right of the data
    assert display.snap_back((400.0, 250.0), data)         # order does not matter
    assert display.snap_back((-500.0, -101.0), data)       # entirely left
    assert display.snap_back((-600.0, 700.0), data)        # 1300 > 3 x 300
    assert not display.snap_back((-400.0, 500.0), data)    # exactly 3 spans: kept
    assert not display.snap_back((10.0, 60.0), data)       # a zoom inside
    assert not display.snap_back((-150.0, 20.0), data)     # half outside: kept
    assert not display.snap_back((10.0, 60.0), (5.0, 5.0)) # degenerate data
    assert not display.snap_back((np.nan, 1.0), data)
    assert display.snap_back((-100.0, 200.0), data, width_factor=0.5)


def test_view_limits_envelope():
    lim = display.view_limits((-100.0, 200.0), (0.0, 1000.0))
    assert lim == {"xMin": -400.0, "xMax": 500.0, "yMin": -2000.0, "yMax": 3000.0}
    lim = display.view_limits((200.0, -100.0), (1000.0, 0.0), x_spans=0.5, y_spans=1.0)
    assert lim == {"xMin": -250.0, "xMax": 350.0, "yMin": -1000.0, "yMax": 2000.0}
    lim = display.view_limits((5.0, 5.0), (3.0, 3.0))   # no span: a nominal one
    assert lim["xMin"] < 5.0 < lim["xMax"] and lim["yMin"] < 3.0 < lim["yMax"]
