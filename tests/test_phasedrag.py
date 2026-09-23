"""Pure-Python unit tests for larmor/phasedrag.py (no Qt, no data).

The drag-to-phase arithmetic: the panel's p0 wrap rule, the pixel -> degree
gesture with its dominant-axis lock and Shift fine control, and the pivot
re-expression, proved exact against processing.op_phase.
"""
import numpy as np
import pytest

from larmor.phasedrag import (FINE, LOCK_PX, P0_DEG_PER_PX, P1_DEG_PER_PX,
                              DragGesture, compensate_p0, wrap_p0)


def test_constants_are_the_documented_sensitivities():
    assert P0_DEG_PER_PX == 0.25 and P1_DEG_PER_PX == 1.0
    assert FINE == 0.1 and LOCK_PX == 4.0


def test_wrap_p0_matches_the_panel_rule():
    assert wrap_p0(190) == -170
    assert wrap_p0(-190) == 170
    assert wrap_p0(180) == 180
    assert wrap_p0(-180) == -180
    assert wrap_p0(540) == 180
    assert wrap_p0(0.25) == 0.25


def test_horizontal_drag_locks_to_p0_at_quarter_degree_per_pixel():
    g = DragGesture()
    assert g.move(40, 2) == (10.0, 0.0)      # 2 px wobble never leaks into p1
    assert g.lock == "p0"
    for _ in range(5):
        res = g.move(40, 0)
    assert res == (60.0, 0.0)
    assert g.move(0, 80) == (60.0, 0.0)      # vertical travel ignored once locked
    assert g.changed is True


def test_vertical_drag_up_is_positive_p1_at_one_degree_per_pixel():
    g = DragGesture()
    assert g.move(1, 50) == (0.0, 50.0)
    assert g.lock == "p1"
    assert g.move(0, -10) == (0.0, 40.0)
    assert g.move(30, 0) == (0.0, 40.0)      # horizontal travel ignored once locked


def test_lock_withholds_then_releases_the_first_pixels():
    g = DragGesture()
    assert g.move(3, 0) is None and g.lock is None
    assert g.move(5, 0) == (2.0, 0.0)        # all 8 px applied on lock, nothing lost

    g = DragGesture()
    assert g.move(2, 3) is None
    assert g.move(0, 2) == (0.0, 5.0)
    assert g.lock == "p1"

    g = DragGesture()
    assert g.move(1, 1) is None
    assert g.changed is False


def test_shift_is_tenfold_finer_and_can_toggle_mid_gesture():
    g = DragGesture()
    assert g.move(100, 0) == (25.0, 0.0)
    dp0, _ = g.move(100, 0, fine=True)
    assert dp0 == pytest.approx(27.5)
    dp0, _ = g.move(-100, 0)                 # incremental: no jump on release
    assert dp0 == pytest.approx(2.5)

    g = DragGesture()                        # fine from the very first event
    dp0, dp1 = g.move(40, 0, fine=True)
    assert dp0 == pytest.approx(1.0) and dp1 == 0.0


def test_pivot_move_compensation_is_exact_against_op_phase():
    from larmor import processing as proc

    rng = np.random.default_rng(7)
    n = 512
    y = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    x = np.linspace(-100.0, 100.0, n)
    p0, p1 = 37.0, 210.0
    s1 = proc.op_phase(proc.from_processed(x, y.copy(), 100.0), p0, p1,
                       pivot_frac=0.30)
    s2 = proc.op_phase(proc.from_processed(x, y.copy(), 100.0),
                       compensate_p0(p0, p1, 0.30, 0.62), p1, pivot_frac=0.62)
    assert np.allclose(s1.y, s2.y, atol=1e-9)
    assert compensate_p0(170, 100, 0.0, 0.5) == -140     # wrapped
