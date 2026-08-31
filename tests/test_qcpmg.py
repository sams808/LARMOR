"""QCPMG processing: period detection, spikelets, and echo coaddition."""
import pathlib

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


def test_split_echoes_honours_an_explicit_echo_count():
    """The number of echoes is stated (ssNake's 'trim to an exact multiple'),
    not inferred from an all-zero test -- real trains end in NOISE, never in
    exact zeros, so the old trim silently kept the noise slots."""
    fid, period = _synthetic_train(period=40, n_echoes=15)
    fid = np.concatenate([fid, np.zeros(40 * 5, complex)])   # 5 blank slots
    assert qcpmg.split_echoes(fid, period).shape == (20, period)
    assert qcpmg.split_echoes(fid, period, n_echoes=15).shape == (15, period)


def test_split_echoes_rejects_a_period_longer_than_the_train():
    with pytest.raises(ValueError, match="exceeds"):
        qcpmg.split_echoes(np.ones(50, complex), 80)


def test_n_usable_echoes_stops_at_the_noise_floor():
    period, n_real = 64, 20
    fid, _ = _synthetic_train(period=period, n_echoes=n_real, decay=0.05)
    rng = np.random.default_rng(0)
    noise = rng.normal(0, 2e-3, period * 15) + 1j * rng.normal(0, 2e-3, period * 15)
    ech = qcpmg.split_echoes(np.concatenate([fid, noise]), period)
    assert ech.shape[0] == 35                       # all slots present...
    assert abs(qcpmg.n_usable_echoes(ech) - n_real) <= 3   # ...but only ~20 real


def test_fit_t2_recovers_decay_and_reports_uncertainty():
    """The echo-top decay yields T2, with an error bar and both time axes."""
    period, n = 64, 60
    tau = 5e-4                       # 0.5 ms echo spacing
    T2_true = 12e-3                  # 12 ms
    half = period // 2
    t = np.arange(period) - half
    base = np.exp(-(t / 5.0) ** 2)
    echoes = np.array([base * np.exp(-k * tau / T2_true) for k in range(n)],
                      dtype=complex)
    top = qcpmg.echo_top_point(echoes)
    f = qcpmg.fit_t2(tau, qcpmg.echo_decay(echoes, top), period=period)
    assert f.ok and not f.pinned
    assert f.T2_s == pytest.approx(T2_true, rel=0.05)
    assert np.isfinite(f.T2_err_s) and f.T2_err_s > 0     # was discarded before
    assert f.r2 > 0.99
    # the two axes differ by exactly the echo length -- the ssNake trap
    assert f.T2_ssnake_s == pytest.approx(f.T2_s / period, rel=1e-9)
    assert f.lb_Hz == pytest.approx(1.0 / (np.pi * f.T2_s), rel=1e-9)
    assert f.lb_ssnake_Hz == pytest.approx(f.lb_Hz * period, rel=1e-6)


def test_fit_t2_flags_failure_instead_of_fabricating_a_number():
    """A fit that cannot be made must say so -- it used to return a plausible
    bold number (the bound itself) for a flat or NaN decay."""
    bad = qcpmg.fit_t2(1e-4, np.array([1.0, np.nan, 0.5, 0.2, 0.1]))
    assert not bad.ok and bad.T2_s == 0.0
    assert not qcpmg.fit_t2(1e-4, np.zeros(20)).ok
    assert not qcpmg.fit_t2(1e-4, np.array([1.0, 1.0, 1.0])).ok   # too few


def test_matched_lb_relation():
    assert qcpmg.matched_lb_Hz(4.241e-3) == pytest.approx(75.07, rel=1e-3)
    assert qcpmg.matched_lb_Hz(0.0) == 0.0


def test_echo_top_is_robust_to_noise_in_any_single_echo():
    """The top is taken from the MEAN of all echoes; one noisy row must not
    move it (a single mid-train row picked 14 instead of 147 on real data)."""
    period, n, true_top = 200, 40, 100
    t = np.arange(period) - true_top
    base = np.exp(-(t / 4.0) ** 2)
    rng = np.random.default_rng(1)
    ech = np.array([base * np.exp(-k * 0.12)
                    + rng.normal(0, 0.25, period) for k in range(n)], complex)
    assert abs(qcpmg.echo_top_point(ech, "mean") - true_top) <= 2


def test_split_alignment_falls_when_the_period_is_wrong():
    fid, period = _synthetic_train(period=80, n_echoes=30)
    good = qcpmg.split_alignment(qcpmg.split_echoes(fid, period))
    bad = qcpmg.split_alignment(qcpmg.split_echoes(fid, period - 7))
    assert good > 0.9 and bad < good


def test_first_last_echo_returns_both_ends():
    fid, period = _synthetic_train(period=64, n_echoes=25, decay=0.05)
    ech = qcpmg.split_echoes(fid, period)
    first, last = qcpmg.first_last_echo(ech)
    assert first.shape == last.shape == (period,)
    assert first.max() > last.max()                  # it decayed


def test_echo_period_from_meta_reads_the_pulse_program():
    cnst = [0.0] * 10
    cnst[7] = 533.3333                                # spikelet spacing, Hz
    per, src = qcpmg.echo_period_from_meta({"sw_Hz": 156250.0, "cnst": cnst})
    assert per == pytest.approx(292.9688, rel=1e-4) and src == "CNST7"
    cnst2 = [0.0] * 10
    cnst2[8] = 293.0
    assert qcpmg.echo_period_from_meta({"sw_Hz": 0.0, "cnst": cnst2}) == (293.0, "CNST8")
    # nothing recorded -> say so, never invent a number
    assert qcpmg.echo_period_from_meta({"sw_Hz": 156250.0}) == (0.0, "none")


def test_detect_period_reports_failure_rather_than_a_wrong_integer():
    rng = np.random.default_rng(2)
    assert qcpmg.detect_period(rng.normal(0, 1, 4000) + 0j) == 0
    assert qcpmg.detect_period(np.zeros(500, complex)) == 0


def test_detect_period_survives_a_wide_echo():
    """A fixed min_period=8 left the autocorrelation's central lobe intact for
    a wide echo, so the lobe won the argmax (real 35Cl train: 8 vs true 293)."""
    period = 300
    t = np.arange(period) - period // 2
    echo = np.exp(-(t / 40.0) ** 2)                  # much wider than 8 points
    fid = np.concatenate([echo * np.exp(-0.05 * k) for k in range(25)]).astype(complex)
    assert abs(qcpmg.detect_period(fid) - period) <= 2


def test_carrier_ppm_uses_the_processing_reference():
    """(SFO1-SF)*1e6/SF, not O1/BF1 -- the two differ by 50.8 ppm on a real
    referenced dataset, which shifts every reported chemical shift."""
    meta = {"larmor_MHz": 78.35411487, "sf_MHz": 78.3621718681335,
            "o1_Hz": -4074.13, "bf1_MHz": 78.358189}
    ppm, referenced = qcpmg.carrier_ppm(meta)
    assert referenced and ppm == pytest.approx(-102.82, abs=0.01)
    no_procs = dict(meta); no_procs.pop("sf_MHz")
    ppm2, referenced2 = qcpmg.carrier_ppm(no_procs)
    assert not referenced2 and ppm2 == pytest.approx(-51.99, abs=0.01)


def test_sum_echo_spectrum_is_phasable_absorption():
    fid, period = _synthetic_train(period=64, n_echoes=40)
    ppm, spec = qcpmg.sum_echo_spectrum(fid, period, 100000.0, 100.0, gb_Hz=50)
    assert np.iscomplexobj(spec)
    p0 = qcpmg.autophase0(spec)
    real = qcpmg.phase_spectrum(spec, p0).real
    # the autophased real spectrum has a clear positive main peak
    assert real.max() > 3 * abs(real.min())


def _whole_echo(period=293, sw=156250.0, top=146, f0=300.0):
    """One symmetric (whole) echo whose exact transform is pure absorption.

    Deliberately NOT strongly damped: with a fast decay the samples that a
    mis-placed zero-fill boundary moves are numerically zero, so the test
    passes for almost any split index and cannot see the bug it guards.
    Truncated to the SYMMETRIC extent about the top (|t| <= min(top, m-1-top))
    so it is a genuine whole echo for any top -- otherwise the surplus
    one-sided samples are themselves a truncation artefact."""
    k = np.arange(period) - top
    t = k / sw
    y = np.exp(-np.abs(t) * f0 * np.pi) * np.exp(2j * np.pi * f0 * t)
    reach = min(top, period - 1 - top)
    return np.where(np.abs(k) <= reach, y, 0.0)


@pytest.mark.parametrize("zf", [1, 2, 4, 16])
@pytest.mark.parametrize("top", [146, 147, 100])
def test_whole_echo_ft_is_pure_absorption_at_every_zerofill(zf, top):
    """THE regression: zero-fill must go MID-ARRAY, split at (m - top).

    Appending zeros after the negative-time samples re-reads them as large
    positive times and turns a pure absorption spectrum into a ~60%-dispersive
    one. Splitting at m//2 instead of m-top is subtler but still wrong for any
    top away from the centre."""
    period, sw = 293, 156250.0
    ppm, spec = qcpmg.whole_echo_ft(_whole_echo(period, sw, top), top,
                                    sw, 100.0, zf=zf)
    assert np.abs(spec.imag).max() / np.abs(spec.real).max() < 0.02


def test_whole_echo_ft_puts_the_line_at_the_right_frequency():
    """A mis-placed zero-fill split also SHIFTS the peak (measured 38 Hz = 2
    bins at top=100), which would move every reported chemical shift."""
    period, sw, sfo, f0 = 293, 156250.0, 78.0, 300.0
    for top in (146, 147, 100, 80):
        ppm, spec = qcpmg.whole_echo_ft(_whole_echo(period, sw, top, f0), top,
                                        sw, sfo, zf=16)
        peak_hz = ppm[int(np.argmax(spec.real))] * sfo
        assert peak_hz == pytest.approx(f0, abs=15.0), top


def test_fit_t2_honours_an_explicit_time_axis_after_exclusions():
    """THE click-to-exclude regression: dropping a point must not re-time the
    survivors onto 0, tau, 2*tau... Excluding the second echo of a real train
    shifted T2 by -36% when the axis was rebuilt contiguously."""
    tau, T2_true, n = 1.875e-3, 4.24e-3, 38
    t = np.arange(n) * tau
    decay = 0.02 + 1.0 * np.exp(-t / T2_true)
    keep = np.ones(n, bool)
    keep[1] = False                      # exclude an EARLY point: worst case
    good = qcpmg.fit_t2(tau, decay[keep], t_s=t[keep], period=293)
    assert good.T2_s == pytest.approx(T2_true, rel=0.02)
    # without the explicit axis the same points give a materially different T2
    naive = qcpmg.fit_t2(tau, decay[keep], period=293)
    assert abs(naive.T2_s - T2_true) > 5 * abs(good.T2_s - T2_true)
    with pytest.raises(ValueError, match="t_s has"):
        qcpmg.fit_t2(tau, decay[keep], t_s=t)          # length mismatch


def test_fit_t2_refuses_pure_noise():
    """A noise decay converges happily (R2 ~ 0.02) and would otherwise hand
    the user a matched filter computed from nothing."""
    rng = np.random.default_rng(0)
    f = qcpmg.fit_t2(1.875e-3, rng.normal(0, 1, 38), period=293)
    assert not f.ok
    assert f.r2 < 0.5


def test_echo_top_uses_the_coherent_average():
    """Averaging MAGNITUDES adds the noise floor and can land a point off the
    true top; the echoes add in phase, so average them coherently. One point
    matters: it changed a real measured T2 by 110%."""
    period, n, true_top = 293, 38, 147
    t = np.arange(period) - true_top
    # a SHARP echo: adjacent points must differ by more than the noise, or
    # neither estimator could resolve 146 from 147 even in principle
    base = np.exp(-(t / 2.0) ** 2).astype(complex)
    rng = np.random.default_rng(3)
    ech = np.array([base * np.exp(-k * 0.45)
                    + rng.normal(0, 0.02, period)
                    + 1j * rng.normal(0, 0.02, period) for k in range(n)])
    assert qcpmg.echo_top_point(ech) == true_top
    cands = qcpmg.echo_top_candidates(ech)
    assert set(cands) == {"coherent", "magnitude", "first"}
    assert cands["coherent"] == true_top


def test_echo_period_from_meta_rejects_an_unset_constant():
    """CNST defaults to 1.0 on a Bruker dataset, so an untouched CNST7 would
    'read' a period of the whole sweep width and blow up the split."""
    cnst = [1.0] * 10                                  # pulse-program default
    assert qcpmg.echo_period_from_meta(
        {"sw_Hz": 156250.0, "cnst": cnst}, n_points=11194) == (0.0, "none")
    # and a rotor-synchronised guess that is far too short is rejected too
    assert qcpmg.echo_period_from_meta(
        {"sw_Hz": 156250.0, "masr_Hz": 16000.0}, n_points=11194)[1] == "MASR"


def test_split_echoes_rejects_dropping_every_echo():
    with pytest.raises(ValueError, match="leaves no echoes"):
        qcpmg.split_echoes(np.arange(1000) + 0j, 250, drop_first=4)


def test_measurement_helpers_refuse_degenerate_input():
    """The most degenerate cases must not come back looking the most
    confident (a zero-width window used to report sigma = 0.000)."""
    x = np.linspace(-200, 200, 4001)
    # a peak pinned to the very edge -> no window can be walked
    hi, lo = qcpmg.cg_window(x, np.exp(-((x + 200) / 10.0) ** 2))
    assert (not np.isfinite(hi)) or hi > lo
    cg, sig = qcpmg.centre_of_gravity(x, np.zeros_like(x), (float("nan"),) * 2)
    assert not np.isfinite(cg) and not np.isfinite(sig)
    # equal +/- lobes: the weights cancel, the ratio is meaningless
    y = np.exp(-((x - 50) / 10.0) ** 2) - np.exp(-((x + 50) / 10.0) ** 2)
    cg2, _ = qcpmg.centre_of_gravity(x, y, (200.0, -200.0))
    assert not np.isfinite(cg2)
    assert qcpmg.fwhm_hz(x, np.zeros_like(x), 100.0) == 0.0


def test_coadd_spectrum_agrees_with_the_shared_whole_echo_transform():
    """One whole-echo convention in the module: the legacy coadd path used to
    keep the old one-sided window and end-appended zero-fill, giving a 53%
    wider line than sum_echo_spectrum for the same data."""
    fid, period = _synthetic_train(period=64, n_echoes=40, decay=0.03)
    sw, sfo = 100000.0, 100.0
    ppm_c, env = qcpmg.coadd_spectrum(fid, period, sw, sfo, lb_Hz=0.0, zf=8)
    ech = qcpmg.split_echoes(fid, period)
    top = qcpmg.echo_top_point(ech)
    ppm_s, spec = qcpmg.whole_echo_ft(qcpmg.sum_echoes(ech, period / sw), top,
                                      sw, sfo, zf=8)
    a = env / (env.max() or 1.0)
    b = np.abs(spec) / (np.abs(spec).max() or 1.0)
    assert ppm_c.shape == ppm_s.shape
    assert float(np.abs(a - b).max()) < 0.05


def test_whole_echo_apodization_window_is_symmetric_about_the_top():
    """After the swap the echo's LEFT half sits at the END of the array; a
    one-sided window multiplied it by ~0, discarding 44% of the echo."""
    n, sw = 512, 50000.0
    w = qcpmg._gaussian_apod(n, sw, 2000.0, whole_echo=True)
    assert w[0] == pytest.approx(1.0)
    assert w[1] == pytest.approx(w[-1], rel=1e-9)     # symmetric
    assert w[n // 2] == pytest.approx(w.min())        # minimum in the middle
    one_sided = qcpmg._gaussian_apod(n, sw, 2000.0, whole_echo=False)
    assert one_sided[-1] < 1e-6 < w[-1]               # the old behaviour


def test_whole_echo_apodization_keeps_the_spectrum_absorptive():
    period, sw, top = 512, 50000.0, 147
    echo = _whole_echo(period, sw, top)
    for gb in (500.0, 2000.0):
        ppm, spec = qcpmg.whole_echo_ft(echo, top, sw, 100.0, gb_Hz=gb, zf=4)
        assert np.abs(spec.imag).max() / np.abs(spec.real).max() < 0.05


def test_whole_echo_ft_rejects_an_out_of_range_top():
    with pytest.raises(ValueError, match="outside the echo"):
        qcpmg.whole_echo_ft(np.ones(50, complex), 180, 50000.0, 100.0)


def test_autophase0_sign_matches_how_callers_apply_it():
    """autophase0 scored one sign and every consumer applied the other, so
    applying the returned angle INVERTED the spectrum (max(real) went negative)."""
    n = 2048
    x = np.arange(n) - n / 2
    line = np.exp(-(x / 30.0) ** 2).astype(complex)      # positive absorption
    for true_p0 in (60.0, -100.0, 155.0):
        dephased = line * np.exp(1j * np.deg2rad(true_p0))
        p0 = qcpmg.autophase0(dephased)
        real = qcpmg.phase_spectrum(dephased, p0).real
        assert real.max() > 0.99 * np.abs(line).max()
        assert real.max() > 20 * abs(real.min())


def test_autophase_recovers_a_positive_lineshape():
    n = 2048
    x = np.arange(n) - n / 2
    line = np.exp(-(x / 40.0) ** 2).astype(complex)
    for true_p0 in (60.0, 120.0, -100.0):
        p0, p1 = qcpmg.autophase(line * np.exp(1j * np.deg2rad(true_p0)))
        real = qcpmg.phase_spectrum(line * np.exp(1j * np.deg2rad(true_p0)), p0, p1).real
        assert real.max() > 0.95 * np.abs(line).max()


def test_autophase_handles_a_line_far_from_the_pivot():
    """A genuine linear phase across a line at ~25% of the axis needs
    |p1| ~ 180-360 deg; the old 1e-6 quadratic penalty charged that as much
    as a visible negative lobe, so the optimiser returned a smaller p1 and a
    partly dispersive spectrum."""
    n = 2048
    x = np.arange(n)
    line = np.exp(-((x - 0.25 * n) / 30.0) ** 2).astype(complex)
    for true_p0, true_p1 in ((30.0, 240.0), (-45.0, -300.0)):
        ramp = x / (n - 1) - 0.5
        dephased = line * np.exp(1j * np.deg2rad(true_p0 + true_p1 * ramp))
        p0, p1 = qcpmg.autophase(dephased)
        real = qcpmg.phase_spectrum(dephased, p0, p1).real
        assert real.max() > 0.95 * np.abs(line).max(), (true_p0, true_p1)
        assert np.abs(real[real < 0]).sum() < 0.05 * real.sum()


def test_fwhm_is_the_outermost_crossing_at_any_sampling():
    """FWHM = the outermost half-max crossings, by convention: a two-horned
    pattern whose saddle dips below half spans BOTH horns, and the value must
    not depend on how finely the axis is sampled (a sample-count 'noise
    pruning' rule made the same lineshape report 1367 vs 1465 Hz at zf 2 vs 4
    -- resolution-dependent, so it was removed). A spike in the window is
    handled by the WINDOW, which is draggable."""
    sfo = 100.0
    for npts in (401, 4001):                          # coarse and fine axes
        x = np.linspace(-200.0, 200.0, npts)
        horn = lambda c: np.exp(-((x - c) / 8.0) ** 2)
        w2 = qcpmg.fwhm_hz(x, horn(-30.0) + horn(30.0), sfo)
        assert 60.0 * sfo < w2 < 90.0 * sfo, npts      # spans both horns
    x = np.linspace(-200.0, 200.0, 4001)
    horn = np.exp(-(x / 8.0) ** 2)
    spiked = horn.copy()
    spiked[3800] = 0.9                                 # one noisy sample far out
    clean = qcpmg.fwhm_hz(x, horn, sfo)
    # the designed remedy: a window excluding the spike recovers the line
    assert qcpmg.fwhm_hz(x, spiked, sfo, (100.0, -100.0)) == pytest.approx(
        clean, rel=1e-6)


def test_sum_echoes_normalised_is_invariant_to_t2_weighting():
    """Toggling the matched filter must not rescale the spectrum (it changed
    the peak by 0.56x before), or 'before/after' plots are not comparable."""
    fid, period = _synthetic_train(period=64, n_echoes=40, decay=0.05)
    ech = qcpmg.split_echoes(fid, period)
    tau = 1e-4
    plain = qcpmg.sum_echoes(ech, tau, None, normalise=True)
    weighted = qcpmg.sum_echoes(ech, tau, 2e-3, normalise=True)
    # weighting still favours the early (larger) echoes, so the peak may rise
    # somewhat -- what must not happen is the old wholesale rescaling
    ratio = float(np.abs(weighted).max() / np.abs(plain).max())
    assert 0.7 < ratio < 1.6, ratio
    raw = qcpmg.sum_echoes(ech, tau, None, normalise=False)
    assert np.abs(raw).max() > 5 * np.abs(plain).max()     # opt-out still works


def test_apodize_echoes_returns_the_weighted_matrix():
    fid, period = _synthetic_train(period=32, n_echoes=12)
    ech = qcpmg.split_echoes(fid, period)
    tau, t2 = 1e-4, 5e-4
    out = qcpmg.apodize_echoes(ech, tau, t2)
    assert out.shape == ech.shape
    assert out[3] == pytest.approx(ech[3] * np.exp(-3 * tau / t2))
    w = qcpmg.echo_weights(12, tau, t2)
    assert w[0] == 1.0 and w[-1] < w[0]


def test_centre_of_gravity_and_fwhm_on_a_known_lineshape():
    ppm = np.linspace(-200, 200, 4001)
    y = np.exp(-((ppm + 50.0) / 20.0) ** 2)          # centred at -50 ppm
    cg, sigma = qcpmg.centre_of_gravity(ppm, y, (100.0, -200.0))
    assert cg == pytest.approx(-50.0, abs=0.5)
    assert sigma >= 0.0
    # FWHM of a Gaussian exp(-(x/a)^2) is 2*a*sqrt(ln2)
    expected_ppm = 2 * 20.0 * np.sqrt(np.log(2))
    assert qcpmg.fwhm_hz(ppm, y, 100.0, (100.0, -200.0)) == pytest.approx(
        expected_ppm * 100.0, rel=0.02)


def test_cg_sigma_grows_when_the_window_is_sloppy():
    """sigma is the QUALITY FLAG that tells the user to place the edges by
    hand: it must react to a window whose edges sit on a sloping tail."""
    ppm = np.linspace(-400, 200, 6001)
    y = np.exp(-((ppm + 50.0) / 25.0) ** 2) + 0.25 * np.exp(-((ppm + 250.0) / 90.0) ** 2)
    _, tight = qcpmg.centre_of_gravity(ppm, y, (10.0, -110.0))
    _, sloppy = qcpmg.centre_of_gravity(ppm, y, (100.0, -400.0))
    assert sloppy > tight


def test_overlay_pair_puts_both_traces_on_one_scale():
    a = np.array([0.0, 2.0, 1.0]); b = np.array([0.0, 0.5, 0.25])
    sa, sb = qcpmg.overlay_pair(a, b)
    assert sa.max() == pytest.approx(1.0) and sb.max() == pytest.approx(1.0)
    za, zb = qcpmg.overlay_pair(np.zeros(3), b)      # no divide-by-zero
    assert np.all(np.isfinite(za)) and np.all(np.isfinite(zb))


# --------------------------------------------------------------------------
# Acceptance against the real, published MagLab 35Cl dataset. Skipped when the
# data is not on this machine, but the expected values stay here so the
# protocol is documented either way.

_DATA = pathlib.Path(
    r"C:\Users\samso\Desktop\WSU_work\NMR\MagLab\DATA\35Cl_2025-12")

#: EXPNO -> (sample, published ssNake T in s)  [xlsx 'T2' sheet]
_PUBLISHED = {1: ("LAW3CL0CA", 1.446e-05), 3: ("LAW3CL4CA", 6.381e-05),
              4: ("LAW3CL1CA", 1.783e-05), 5: ("LAW3CL3CA", 4.078e-05),
              6: ("NS3", 1.216e-05), 7: ("Ab-3Cl", 8.232e-06),
              8: ("SV6-Cl", 2.548e-05), 9: ("CN1-3Cl", 3.290e-05),
              11: ("LAW3CL2CA", 2.910e-05), 14: ("SV-1Cl", 2.685e-05),
              15: ("CS-3Cl", 4.694e-05), 16: ("An-3Cl", 7.477e-05)}

_needs_data = pytest.mark.skipif(not _DATA.exists(),
                                 reason="MagLab 35Cl dataset not on this machine")


def _load_real(expno: int):
    from larmor.io import bruker
    d = bruker.read(str(_DATA / str(expno) / "fid"))
    return np.asarray(d.data, complex), dict(d.meta)


@_needs_data
def test_real_data_metadata_is_read_not_guessed():
    fid, meta = _load_real(1)
    per, src = qcpmg.echo_period_from_meta(meta, n_points=fid.size)
    assert src == "CNST7" and per == pytest.approx(292.9688, rel=1e-5)
    ppm, referenced = qcpmg.carrier_ppm(meta)
    assert referenced
    # TopSpin's own OFFSET - SW_p/(2*SF) for this dataset
    assert ppm == pytest.approx(-102.817, abs=0.01)


@_needs_data
@pytest.mark.parametrize("expno", sorted(_PUBLISHED))
def test_real_data_reproduces_the_published_t2(expno):
    """The acceptance test: LARMOR must reproduce the ssNake T value the
    published protocol recorded, for every sample in the set."""
    name, t_pub = _PUBLISHED[expno]
    fid, meta = _load_real(expno)
    per = int(round(qcpmg.echo_period_from_meta(meta, n_points=fid.size)[0]))
    assert per == 293
    ech = qcpmg.split_echoes(fid, per)
    top = qcpmg.echo_top_point(ech)
    assert top == 147, f"{name}: echo top moved (ssNake used 147)"
    f = qcpmg.fit_t2(per / meta["sw_Hz"], qcpmg.echo_decay(ech, top, "real"),
                     offset=True, period=per)
    assert f.T2_ssnake_s == pytest.approx(t_pub, rel=0.07), name
    assert f.lb_ssnake_Hz == pytest.approx(1.0 / (np.pi * t_pub), rel=0.07)


@_needs_data
def test_real_data_sum_echo_is_absorptive_and_measurable():
    fid, meta = _load_real(1)
    per = int(round(qcpmg.echo_period_from_meta(meta, n_points=fid.size)[0]))
    ech = qcpmg.split_echoes(fid, per)
    top = qcpmg.echo_top_point(ech)
    f = qcpmg.fit_t2(per / meta["sw_Hz"], qcpmg.echo_decay(ech, top, "real"),
                     offset=True, period=per)
    carrier, _ = qcpmg.carrier_ppm(meta)
    ppm, spec = qcpmg.sum_echo_spectrum(
        fid, per, meta["sw_Hz"], meta["larmor_MHz"], carrier, top=top,
        t2_weight_s=f.T2_s, zf=16)
    real = qcpmg.phase_spectrum(spec, qcpmg.autophase0(spec)).real
    # a correct whole-echo transform needs only a zero-order phase
    assert real.min() / real.max() > -0.10
    fwhm = qcpmg.fwhm_hz(ppm, real, meta["larmor_MHz"], qcpmg.cg_window(ppm, real))
    assert fwhm == pytest.approx(4638.0, rel=0.05)      # published 4638 Hz


def test_find_period_by_correlation_recovers_the_period():
    """The data-measured period: candidates from the autocorrelation,
    verified by echo-anchored correlation, refined by a three-periods-apart
    block (a one-point error drifts three points there)."""
    per, true = 250, 250
    t = np.arange(per) - 125
    echo = np.exp(-np.abs(t) / 20.0) * np.exp(2j * np.pi * 0.03 * t)
    rng = np.random.default_rng(1)
    # a long healthy train
    tr = np.concatenate([echo * np.exp(-k / 8.0) for k in range(40)])
    tr = tr + rng.normal(0, 0.01, tr.size) + 1j * rng.normal(0, 0.01, tr.size)
    p, sc = qcpmg.find_period_by_correlation(tr)
    assert p == true and sc > 0.9
    # the 81Br situation: T2' of ~1.5 echoes, the tail is pure noise
    tr2 = np.concatenate([echo * np.exp(-k / 1.5) for k in range(24)])
    tr2 = tr2 + rng.normal(0, 0.02, tr2.size) + 1j * rng.normal(0, 0.02, tr2.size)
    p2, sc2 = qcpmg.find_period_by_correlation(tr2)
    assert p2 == true and sc2 > 0.8
    # pure noise: refuse, never invent
    assert qcpmg.find_period_by_correlation(
        rng.normal(0, 1, 6000) + 0j) == (0, 0.0)


def test_period_correlation_discriminates_and_anchors():
    """The score must collapse at a wrong period, and must NOT be fooled by
    the smooth pre-echo ramp (blocks are anchored on the strongest echo --
    unanchored first-blocks correlate at ANY small lag)."""
    per = 300
    t = np.arange(per) - 220                    # echo far into the block
    echo = np.exp(-np.abs(t) / 15.0) * np.exp(2j * np.pi * 0.05 * t)
    tr = np.concatenate([echo * np.exp(-k / 6.0) for k in range(20)])
    good = qcpmg.period_correlation(tr, per)
    assert good > 0.95
    assert qcpmg.period_correlation(tr, int(per * 0.63)) < 0.5 * good
    # NOTE: tiny periods still correlate (any smooth echo is coherent over a
    # few samples) -- that degeneracy is handled by find_period_by_correlation
    # only scoring autocorrelation-proposed candidates, never raw small lags


@_needs_data
def test_real_data_period_finder():
    """81Br (CNST7 unset, MASR guess 179 pts, true 475) and 35Cl (true 293):
    the finder must land both from the data alone."""
    from larmor.io import bruker
    d = bruker.read(str(_DATA / "1" / "fid"))
    p, sc = qcpmg.find_period_by_correlation(np.asarray(d.data, complex))
    assert p == 293 and sc > 0.9
    br = _DATA.parent / "81Br_2026-08" / "30" / "fid"
    if br.exists():
        d2 = bruker.read(str(br))
        p2, sc2 = qcpmg.find_period_by_correlation(np.asarray(d2.data, complex))
        assert p2 == 475 and sc2 > 0.9


def test_magnitude_costs_no_width_on_a_whole_echo():
    """THE correction: sqrt(3) is the widening of a CAUSAL (one-sided) FID,
    where absorption and dispersion are Hilbert partners. A whole echo is
    symmetric about t=0, so its spectrum is REAL and |spec| recovers the
    absorption lineshape unchanged. The dialog told users the opposite."""
    sw, sfo, per = 100000.0, 100.0, 512
    t = np.arange(per) - per // 2
    echo = np.exp(-np.abs(t) / 40.0).astype(complex)      # symmetric = whole
    f, sp = qcpmg.whole_echo_ft(echo, per // 2, sw, sfo, zf=16)
    hz = f * sfo

    def fwhm(y):
        i = np.where(y >= 0.5 * y.max())[0]
        return abs(float(hz[i[-1]] - hz[i[0]]))

    assert np.abs(sp.imag).max() / np.abs(sp.real).max() < 1e-3   # real
    assert fwhm(np.abs(sp)) / fwhm(sp.real) == pytest.approx(1.0, abs=0.01)

    # ... while the sqrt(3) IS right for the one-sided case the constant names
    causal = np.fft.fftshift(np.fft.fft(np.exp(-np.arange(per) / 40.0) + 0j,
                                        n=per * 16))
    hz = np.fft.fftshift(np.fft.fftfreq(per * 16, 1.0 / sw))
    assert fwhm(np.abs(causal)) / fwhm(causal.real) == pytest.approx(
        qcpmg.MAGNITUDE_LORENTZ_WIDENING, rel=0.05)


def test_magnitude_spectrum_removes_the_rectified_floor():
    """Rectification turns zero-mean noise into a POSITIVE pedestal, which
    drags an integrated centre of gravity toward the window centre. The floor
    is estimated from the EDGES: the median of the whole trace stops being a
    floor once a wide pattern fills the window (it over-estimated a real
    81Br floor 6.2x and moved delta_CG by 44 ppm)."""
    rng = np.random.default_rng(0)
    noise = rng.normal(0, 1, 20000) + 1j * rng.normal(0, 1, 20000)
    assert np.abs(noise).mean() > 1.0                     # the pedestal
    assert abs(qcpmg.magnitude_spectrum(noise).mean()) < 0.15
    assert qcpmg.magnitude_spectrum(noise, subtract_floor=False).mean() > 1.0
    # a real line survives the subtraction
    x = np.linspace(-100.0, 100.0, 4001)
    line = np.exp(-((x - 10.0) / 3.0) ** 2).astype(complex)
    m = qcpmg.magnitude_spectrum(line)
    assert x[int(np.argmax(m))] == pytest.approx(10.0, abs=0.2)
    assert m.max() == pytest.approx(1.0, abs=0.01)
    # and a pattern filling most of the window is NOT eaten: the edge floor
    # stays at the baseline where the global median would sit inside the line
    wide = np.exp(-(x / 30.0) ** 2) + 0.02
    assert qcpmg.noise_floor(wide) == pytest.approx(0.02, abs=0.005)
    assert float(np.median(wide)) > 3 * qcpmg.noise_floor(wide)
    kept = qcpmg.magnitude_spectrum(wide.astype(complex))
    assert kept.max() == pytest.approx(1.0, abs=0.03)


def test_magnitude_is_phase_independent():
    """The whole point of mc: any p0/p1 leaves it unchanged."""
    rng = np.random.default_rng(3)
    spec = rng.normal(0, 1, 2000) + 1j * rng.normal(0, 1, 2000)
    base = qcpmg.magnitude_spectrum(spec)
    for p0, p1 in ((37.0, 0.0), (-120.0, 250.0), (0.0, -1000.0)):
        rot = qcpmg.phase_spectrum(spec, p0, p1)
        assert np.allclose(qcpmg.magnitude_spectrum(rot), base, atol=1e-9)


@_needs_data
def test_magnitude_rescues_an_unphaseable_real_pattern():
    """The 81Br LAW4Ca line is wider than the refocusing bandwidth, so its
    phase error is non-linear and no p0/p1 flattens it (the absorption
    spectrum keeps ~-47 % dips, which wreck the delta_CG integral).
    Magnitude needs no phase and gives a usable number."""
    from larmor.io import bruker
    src = _DATA.parent / "81Br_2026-08" / "34" / "fid"
    if not src.exists():
        pytest.skip("81Br dataset not on this machine")
    d = bruker.read(str(src))
    fid = np.asarray(d.data, complex)
    sw, sfo = float(d.meta["sw_Hz"]), float(d.meta["larmor_MHz"])
    per, _ = qcpmg.find_period_by_correlation(fid)
    ech = qcpmg.split_echoes(fid, per)
    top = qcpmg.echo_top_point(ech)
    tau = per / sw
    f = qcpmg.fit_t2(tau, qcpmg.echo_decay(ech, top, "real"), offset=True,
                     period=per)
    _, sp = qcpmg.whole_echo_ft(qcpmg.sum_echoes(ech, tau, f.T2_s), top,
                                sw, sfo, zf=16)
    absn = qcpmg.phase_spectrum(sp, *qcpmg.autophase(sp)).real
    mag = qcpmg.magnitude_spectrum(sp)
    assert absn.min() / absn.max() < -0.3        # unphaseable, by a mile
    assert mag.min() / mag.max() > -0.2          # magnitude is far cleaner


def test_second_order_phase_rescues_a_swept_pulse_chirp():
    """A frequency-swept (WURST/chirp) refocusing pulse imprints a QUADRATIC
    phase that p0/p1 cannot remove. autophase_best must notice and add p2 --
    and must NOT add it to an ordinary echo."""
    n = 4096
    x = np.linspace(-1.0, 1.0, n)
    pattern = (np.exp(-((x + 0.15) / 0.22) ** 2)
               + 0.8 * np.exp(-((x - 0.2) / 0.18) ** 2))
    ramp = np.arange(n) / (n - 1) - 0.5

    # (a) a genuine chirp phase: p0/p1 cannot flatten it, p0/p1/p2 can
    chirped = pattern * np.exp(1j * np.deg2rad(
        25.0 + 60.0 * ramp + 900.0 * (2.0 * ramp) ** 2))
    p0, p1 = qcpmg.autophase(chirped, order=1)
    lin = qcpmg.phase_spectrum(chirped, p0, p1).real
    q0, q1, q2 = qcpmg.autophase_best(chirped)
    quad = qcpmg.phase_spectrum(chirped, q0, q1, p2_deg=q2).real
    assert abs(q2) > 100.0                       # the term was needed...
    assert quad.min() / quad.max() > 3 * (lin.min() / lin.max())   # ...and helped
    assert quad.min() / quad.max() > -0.10

    # (b) an ordinary p0/p1 error: p2 must stay at exactly zero
    plain = pattern * np.exp(1j * np.deg2rad(40.0 + 30.0 * ramp))
    r0, r1, r2 = qcpmg.autophase_best(plain)
    assert r2 == 0.0
    got = qcpmg.phase_spectrum(plain, r0, r1, p2_deg=r2).real
    assert got.min() / got.max() > -0.05


def test_phase_spectrum_p2_is_a_pure_quadratic():
    """p2 is defined in degrees at the ENDS of the axis, zero at the pivot,
    so it never moves the centre of a centred line."""
    n = 1001
    spec = np.ones(n, complex)
    out = qcpmg.phase_spectrum(spec, 0.0, 0.0, p2_deg=90.0)
    ang = np.degrees(np.angle(out))
    assert ang[n // 2] == pytest.approx(0.0, abs=1e-9)     # zero at the pivot
    assert ang[0] == pytest.approx(-90.0, abs=1e-6)        # -p2 at both ends
    assert ang[-1] == pytest.approx(-90.0, abs=0.5)
    assert np.allclose(qcpmg.phase_spectrum(spec, 0.0, 0.0, p2_deg=0.0), spec)


@_needs_data
def test_p2_and_magnitude_agree_on_the_real_wcpmg_dataset():
    """The cross-validation that makes both trustworthy: on a real WURST-CPMG
    81Br sample, p0/p1 alone gives a meaningless delta_CG, while the
    p2-phased absorption and the magnitude spectrum agree to a few ppm."""
    from larmor.io import bruker
    src = _DATA.parent / "81Br_2026-08" / "34" / "fid"
    if not src.exists():
        pytest.skip("81Br dataset not on this machine")
    d = bruker.read(str(src))
    fid = np.asarray(d.data, complex)
    sw, sfo = float(d.meta["sw_Hz"]), float(d.meta["larmor_MHz"])
    per, _ = qcpmg.find_period_by_correlation(fid)
    ech = qcpmg.split_echoes(fid, per)
    top = qcpmg.echo_top_point(ech)
    tau = per / sw
    f = qcpmg.fit_t2(tau, qcpmg.echo_decay(ech, top, "real"), offset=True,
                     period=per)
    ppm, sp = qcpmg.whole_echo_ft(qcpmg.sum_echoes(ech, tau, f.T2_s), top,
                                  sw, sfo, zf=16)
    p0, p1, p2 = qcpmg.autophase_best(sp)
    assert abs(p2) > 100.0                      # a swept pulse: p2 is needed
    phased = qcpmg.phase_spectrum(sp, p0, p1, p2_deg=p2).real
    assert phased.min() / phased.max() > -0.10  # was -0.47 with p0/p1 only

    def cg(y):
        hi, lo = qcpmg.cg_window(ppm, y)
        return qcpmg.centre_of_gravity(ppm, y, (hi, lo))[0]

    assert cg(phased) == pytest.approx(cg(qcpmg.magnitude_spectrum(sp)), abs=15.0)


def test_period_is_read_from_the_nmrfam_constants_too():
    """Not every QCPMG sequence uses Bruker's CNST7. The NMRFAM/Perras
    qcpmg.av4.nmrfam writes CNST11 (spikelet Hz), CNST14 (points) and CNST15
    (period us) and leaves CNST7 at its default 1.0 -- which LARMOR then
    rejected, fell back to a MASR guess, and split the train wrongly."""
    def meta(**kw):
        c = [1.0] * 16
        for i, v in kw.items():
            c[int(i[1:])] = v
        return {"sw_Hz": 200000.0, "cnst": c}

    assert qcpmg.echo_period_from_meta(              # spikelet spacing, Hz
        meta(c11=500.0), n_points=20978) == (400.0, "CNST11")
    assert qcpmg.echo_period_from_meta(              # echo period, us
        meta(c15=2000.0), n_points=20978) == (400.0, "CNST15")
    assert qcpmg.echo_period_from_meta(              # points per echo
        meta(c14=400.0), n_points=20978) == (400.0, "CNST14")
    # Bruker's own still win when both are present
    assert qcpmg.echo_period_from_meta(
        meta(c7=533.3333, c11=500.0), n_points=20978)[1] == "CNST7"
    # and an untouched CNST11 (its default 1.0) is still rejected
    assert qcpmg.echo_period_from_meta(meta(), n_points=20978) == (0.0, "none")


def test_centre_offset_finds_a_train_that_starts_on_an_echo_top():
    """Some sequences begin acquiring AT a top, so the natural blocks hold
    the right half of one echo and the left half of the NEXT -- two different
    echoes glued into a fake one. Measured consequences on real 35Cl data:
    T2 4.0 ms instead of 10.3, p1 = 493 deg and p2 = 331 deg needed to phase,
    FWHM inflated 24 %."""
    period, n = 400, 30
    t = np.arange(period) - period // 2
    echo = np.exp(-np.abs(t) / 25.0) * np.exp(2j * np.pi * 0.02 * t)
    train = np.concatenate([echo * np.exp(-k / 8.0) for k in range(n)])
    # conventional acquisition: already centred, so nothing is changed
    assert qcpmg.centre_offset(train, period) == 0
    # start recording at the first top instead
    shifted = train[period // 2:]
    off = qcpmg.centre_offset(shifted, period)
    assert abs(off - period // 2) <= 3
    tops = [int(np.argmax(np.abs(e)))
            for e in qcpmg.split_echoes(shifted, period, first=off)[:5]]
    assert all(abs(tp - period // 2) <= 3 for tp in tops)
    # a too-short train cannot be judged: say 0 rather than guess
    assert qcpmg.centre_offset(train[:100], period) == 0
