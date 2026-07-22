"""Variable-temperature analysis: Arrhenius and Vogel–Fulcher–Tammann fits of a
rate (or 1/T1, correlation rate…) versus temperature."""
from __future__ import annotations

import numpy as np

R = 8.314462618          # J·mol⁻¹·K⁻¹


def fit_arrhenius(T_K, rate) -> dict:
    """k(T) = A·exp(−Ea/RT). Linear fit of ln k vs 1/T. Returns Ea (kJ/mol),
    A, and a callable curve(T)."""
    T = np.asarray(T_K, float)
    k = np.asarray(rate, float)
    x = 1.0 / T
    ylog = np.log(k)
    # ylog = lnA − (Ea/R)·x
    Amat = np.vstack([np.ones_like(x), -x]).T
    (lnA, Ea_over_R), *_ = np.linalg.lstsq(Amat, ylog, rcond=None)
    Ea = Ea_over_R * R
    resid = ylog - (lnA - Ea_over_R * x)
    return {
        "Ea_kJmol": Ea / 1000.0, "A": float(np.exp(lnA)), "lnA": float(lnA),
        "rmsd_lnk": float(np.sqrt(np.mean(resid ** 2))),
        "curve": lambda TT, lnA=lnA, Ea_over_R=Ea_over_R:
            np.exp(lnA - Ea_over_R / np.asarray(TT, float)),
    }


def fit_vft(T_K, rate) -> dict:
    """k(T) = A·exp(−B/(T−T0)) (Vogel–Fulcher–Tammann). Returns A, B (K), T0 (K)
    and a callable."""
    from scipy.optimize import curve_fit

    T = np.asarray(T_K, float)
    k = np.asarray(rate, float)

    def model(TT, lnA, B, T0):
        return lnA - B / (TT - T0)

    p0 = [float(np.log(k.max())), 1000.0, float(T.min()) - 50.0]
    popt, _ = curve_fit(model, T, np.log(k), p0=p0,
                        bounds=([-50, 1.0, -1e4], [50, 1e6, float(T.min()) - 1.0]),
                        maxfev=20000)
    lnA, B, T0 = popt
    return {
        "A": float(np.exp(lnA)), "B_K": float(B), "T0_K": float(T0),
        "curve": lambda TT, lnA=lnA, B=B, T0=T0:
            np.exp(lnA - B / (np.asarray(TT, float) - T0)),
    }
