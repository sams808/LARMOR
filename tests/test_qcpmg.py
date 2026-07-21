"""QCPMG processing: period detection, spikelets, and echo coaddition."""
import numpy as np
import pytest

from larmor import qcpmg


def _synthetic_train(period=80, n_echoes=40, decay=0.02):
    """A CPMG-like echo train: repeated symmetric echoes with T2 decay."""
    half = period // 2
    t = np.arange(period) - half
    echo = np.exp(-(t / 6.0) ** 2)                # one symmetric echo
    fid = np.concatenate([echo * np.exp(-decay * k) for k in range(n_echoes)])
    return fid.astype(complex), period


def test_detect_period_matches_construction():
    fid, period = _synthetic_train(period=80, n_echoes=40)
    assert abs(qcpmg.detect_period(fid) - period) <= 1


def test_coadd_removes_spikelets():
    """The spikelet spectrum is a comb; the coadded envelope is continuous."""
    from scipy.signal import find_peaks

    fid, period = _synthetic_train(period=64, n_echoes=48)
    sw, sfo = 100000.0, 100.0
    ppm_s, spec_s = qcpmg.spikelet_spectrum(fid, sw, sfo, lb_Hz=20)
    ppm_c, env = qcpmg.coadd_spectrum(fid, period, sw, sfo, lb_Hz=50)
    n_spikes, _ = find_peaks(np.abs(spec_s), height=np.abs(spec_s).max() * 0.2,
                             distance=period // 3)
    n_env, _ = find_peaks(env, height=env.max() * 0.3, distance=period // 3)
    assert len(n_spikes) > 3                       # a genuine spikelet manifold
    assert len(n_env) < len(n_spikes)              # coaddition removed spikelets


def test_spikelet_spacing():
    assert qcpmg.spikelet_spacing_ppm(100, 200000.0, 107.8) == pytest.approx(
        (200000.0 / 100) / 107.8, rel=1e-6)


def test_coadd_echoes_shape():
    fid, period = _synthetic_train(period=50, n_echoes=20)
    echo = qcpmg.coadd_echoes(fid, period)
    assert echo.shape == (period,)
    # coadding boosts the aligned echo top above a single echo's
    assert np.abs(echo).max() > 5.0


def test_split_trims_trailing_zeros():
    fid, period = _synthetic_train(period=40, n_echoes=15)
    fid = np.concatenate([fid, np.zeros(40 * 5, complex)])   # 5 blank slots
    ech = qcpmg.split_echoes(fid, period)
    assert ech.shape == (15, period)                        # blanks dropped


def test_fit_t2_recovers_decay():
    """The echo-top decay yields the transverse relaxation time."""
    period, n = 64, 60
    tau = 5e-4                       # 0.5 ms echo spacing
    T2_true = 12e-3                  # 12 ms
    half = period // 2
    t = np.arange(period) - half
    base = np.exp(-(t / 5.0) ** 2)
    echoes = np.array([base * np.exp(-k * tau / T2_true) for k in range(n)],
                      dtype=complex)
    top = qcpmg.echo_top_point(echoes)
    T2, curve = qcpmg.fit_t2(tau, qcpmg.echo_decay(echoes, top))
    assert T2 == pytest.approx(T2_true, rel=0.05)


def test_sum_echo_spectrum_is_phasable_absorption():
    fid, period = _synthetic_train(period=64, n_echoes=40)
    ppm, spec = qcpmg.sum_echo_spectrum(fid, period, 100000.0, 100.0, gb_Hz=50)
    assert np.iscomplexobj(spec)
    p0 = qcpmg.autophase0(spec)
    real = np.real(spec * np.exp(-1j * np.deg2rad(p0)))
    # the autophased real spectrum has a clear positive main peak
    assert real.max() > 3 * abs(real.min())
