"""Universal Bruker reader (files/folders, 1D/2D, raw/processed) and the
ssNake-style 1D/2D Fourier processing."""
import numpy as np
import pytest

from larmor import fourier
from larmor.io import bruker

from conftest import (BRUKER_1R, BRUKER_2RR_PSEUDO, BRUKER_2RR_MQMAS,
                      BRUKER_FID, BRUKER_SER, require)


# ---------------------------------------------------------------- resolve
def test_resolve_1r_file():
    ref = bruker.resolve(require(BRUKER_1R))
    assert ref.target == "1r" and ref.ndim == 1 and ref.procno == 1
    assert (ref.expno / "acqus").exists()


def test_resolve_2rr_and_fid_and_folder():
    assert bruker.resolve(require(BRUKER_2RR_MQMAS)).target == "2rr"
    assert bruker.resolve(require(BRUKER_FID)).target == "fid"
    # EXPNO folder resolves to its processed 2rr when present
    expno = require(BRUKER_2RR_MQMAS).parent.parent.parent
    assert bruker.resolve(expno).ndim == 2


def test_resolve_rejects_junk(tmp_path):
    with pytest.raises((ValueError, FileNotFoundError)):
        bruker.resolve(tmp_path)


# ---------------------------------------------------------------- read 1D
def test_read_1r_gives_ppm_spectrum():
    d = bruker.read(require(BRUKER_1R))
    assert d.ndim == 1 and d.domain == "freq"
    assert d.nucleus == "27Al"
    ppm = d.axes[0].values
    assert np.all(np.diff(ppm) > 0)
    assert d.data.shape == ppm.shape
    assert d.data.dtype.kind == "f"          # real spectrum


def test_read_fid_is_time_domain_complex():
    d = bruker.read(require(BRUKER_FID))
    assert d.ndim == 1 and d.domain == "time"
    assert np.iscomplexobj(d.data)
    assert any("digital filter" in w for w in d.warnings)


# ---------------------------------------------------------------- read 2D
def test_read_real_mqmas_2d_has_two_ppm_axes():
    d = bruker.read(require(BRUKER_2RR_MQMAS))
    assert d.ndim == 2 and d.domain == "freq" and not d.is_pseudo2d
    f1, f2 = d.axes
    assert f1.unit == "ppm" and f2.unit == "ppm"
    assert d.data.shape == (len(f1.values), len(f2.values))


def test_read_pseudo2d_keeps_arrayed_axis():
    d = bruker.read(require(BRUKER_2RR_PSEUDO))
    assert d.ndim == 2 and d.is_pseudo2d
    # indirect axis is NOT a chemical shift
    assert d.axes[0].unit in ("point", "s")
    assert any("not a chemical shift" in w for w in d.warnings)


def test_read_ser_reshapes_and_uses_vdlist():
    d = bruker.read(require(BRUKER_SER))
    assert d.ndim == 2 and d.domain == "time" and d.is_pseudo2d
    assert d.data.ndim == 2 and np.iscomplexobj(d.data)
    # a saturation-recovery ser carries a vdlist -> delay axis in seconds
    assert d.axes[0].unit == "s"


def test_read_is_strictly_readonly():
    ref = bruker.resolve(require(BRUKER_1R))
    before = bruker.snapshot(ref.expno)
    bruker.read(require(BRUKER_1R), verify=True)
    bruker.verify_untouched(ref.expno, before)      # raises if anything changed


# ---------------------------------------------------------------- 1D FT
def test_ft1d_synthetic_peak():
    sw, n, freq = 10000.0, 1024, 1500.0
    t = np.arange(n) / sw
    fid = np.exp(2j * np.pi * freq * t) * np.exp(-t / 0.02)
    ppm, spec = fourier.ft1d(fid, sw, 100.0, ops=[{"op": "em", "lb_hz": 20}])
    peak = ppm[np.argmax(np.abs(spec))]
    assert abs(abs(peak) - 15.0) < 0.3       # 1500 Hz / 100 MHz = 15 ppm


def test_ft1d_appends_ft_if_missing():
    fid = np.ones(256, complex)
    ppm, spec = fourier.ft1d(fid, 10000.0, 100.0, ops=[])
    assert spec.size == 256 and ppm is not None


# ---------------------------------------------------------------- 2D FT
def test_states_recombine_shapes():
    ser = np.arange(8 * 4).reshape(8, 4).astype(complex)
    hyper = fourier.states_recombine(ser, "States")
    assert hyper.shape == (4, 4)             # cos+i*sin halves the rows
    assert np.iscomplexobj(hyper)
    # QF passes through unchanged
    assert fourier.states_recombine(ser, "QF").shape == (8, 4)


def test_ft2d_places_a_peak_at_known_frequencies():
    """A single 2D oscillation must land at the right (F1, F2) ppm."""
    n1, n2, sw1, sw2 = 32, 256, 4000.0, 20000.0
    sfo1 = sfo2 = 100.0
    f2_hz, f1_hz = 2000.0, 500.0
    t2 = np.arange(n2) / sw2
    t1 = np.arange(n1) / sw1
    # cosine- and sine-modulated States pair
    cos = np.exp(-t1[:, None] / 0.05) * np.cos(2 * np.pi * f1_hz * t1)[:, None]
    sin = np.exp(-t1[:, None] / 0.05) * np.sin(2 * np.pi * f1_hz * t1)[:, None]
    direct = (np.exp(2j * np.pi * f2_hz * t2) * np.exp(-t2 / 0.02))[None, :]
    ser = np.empty((2 * n1, n2), complex)
    ser[0::2] = cos * direct
    ser[1::2] = sin * direct
    p = fourier.FT2DParams(mode="States", f2_ops=[{"op": "em", "lb_hz": 20}],
                           f1_ops=[{"op": "em", "lb_hz": 20}])
    f2, f1, z = fourier.ft2d(ser, sw2, sw1, sfo2, sfo1, params=p)
    i1, i2 = np.unravel_index(int(np.argmax(np.abs(z))), z.shape)
    assert abs(f2[i2] - 20.0) < 0.5          # 2000/100 = 20 ppm
    assert abs(f1[i1] - 5.0) < 0.5           # 500/100 = 5 ppm


@pytest.mark.slow
def test_raw_fid_ft_matches_topspin_1r_peak():
    """The clean handedness+correctness proof: FT a real 1D fid myself and it
    must peak at the SAME ppm as TopSpin's own processed 1r for that EXPNO,
    up to the SR calibration (a few ppm). A mirror-flip would land the peak at
    the opposite sign."""
    fid = bruker.read(require(BRUKER_FID))               # 27Al zg
    ppm, spec = fourier.ft1d(fid.data, fid.axes[0].sw_Hz,
                             fid.meta["larmor_MHz"],
                             ops=[{"op": "fcor", "factor": 0.5},
                                  {"op": "em", "lb_hz": 100}])
    my_peak = ppm[np.argmax(np.abs(spec))]

    r1r = bruker.read(require(BRUKER_1R))
    topspin_peak = r1r.axes[0].values[np.argmax(np.abs(r1r.data))]
    # same resonance to within the SR offset TopSpin applied at processing
    assert abs(my_peak - topspin_peak) < 4.0
    assert np.sign(my_peak) == np.sign(topspin_peak)     # not mirror-flipped


@pytest.mark.slow
def test_ft2d_from_real_ser_is_sane():
    """FT the real pseudo-2D ser end to end: the result must be a finite 2D
    spectrum with the right shape and a real 27Al resonance somewhere in the
    plausible range (structural check; exact ppm depends on SR + baseline)."""
    raw = bruker.read(require(BRUKER_SER))
    d = fourier.ft2d_from_nmrdata(raw, fourier.FT2DParams(
        f2_ops=[{"op": "fcor", "factor": 0.5}, {"op": "em", "lb_hz": 200}]))
    assert d.z.shape[0] == raw.data.shape[0]
    assert np.isfinite(d.z).all() and d.z.max() > 0
    # every retained slice carries real signal (no leaked blank row)
    per_row = np.abs(d.z).max(axis=1)
    assert (per_row > 0).all()
