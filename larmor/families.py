"""Site families: tagged lines summed into structural species, with named
ratios and honestly-stated uncertainties.

A glass paper rarely quotes the population of every fitted line -- it quotes
the BO4 fraction N4 = BO4/(BO3+BO4), the mean aluminium coordination ⟨CN⟩,
the mean Qn connectivity ⟨n⟩ or a grouped fraction ("bonded" phosphorus).
Each is a sum of site populations, and the error of a sum is NOT the
quadrature of the site errors: two overlapping lines have strongly
anticorrelated amplitudes, so their sum is known far better than either
member. This module carries that propagation on three stated bases:

``covariance``
    lmfit's ``MinimizerResult.uvars`` (``uncertainties`` ufloats built from
    the fit covariance, expr-linked parameters included): the family integral
    is ``F = Σ k_i·a_i`` with ``k_i = |integral_i| / amplitude_i`` fixed at
    the best fit, so fractions and ratios evaluated in ufloat arithmetic carry
    the full first-order propagation over the amplitude covariance (shape
    covariance neglected, as quantify() states for its own rows).
``montecarlo``
    per-trial site integrals from ``autofit.monte_carlo_errors`` (every
    synthetic refit is re-integrated, which also carries the shape
    covariance): the family/ratio value is recomputed per trial and the
    spread (population std, like ``MCParam``) is its error. It has to be
    integrals, not amplitudes -- amplitude is the peak height for
    gauss_lor / Czjzek / quadrupolar models and an area only for gl_norm.
``independent``
    neither available: quadrature of the per-site integral errors for a
    family and the delta method for a ratio, FLAGGED as independent so no
    table passes an independence assumption off as a propagated error.

Reported values are always the best-fit numbers, never a trial mean. Family
percentages are of the FULL total (untagged lines included, listed in the
note); ratio denominators use the tagged families only. Qt-free.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Sequence

import numpy as np

from larmor import cellparse

__all__ = ["PRESETS", "RatioDef", "RATIOS", "DERIVED_PARAMS", "GROUP_PARAMS",
           "presets_for", "ratios_for", "guess_from_label", "summarize",
           "trial_series"]

#: per-nucleus family vocabulary offered by the Fit-parameters table's Family
#: submenu (free text is always allowed). 31P carries both the Qn scheme and
#: the PBi grouping (isolated / bonded / chains / Bi-contact).
PRESETS: dict[str, tuple[str, ...]] = {
    "11B": ("BO3", "BO4"),
    "27Al": ("Al(IV)", "Al(V)", "Al(VI)"),
    "29Si": ("Q0", "Q1", "Q2", "Q3", "Q4"),
    "31P": ("Q0", "Q1", "Q2", "Q3", "isolated", "bonded", "chains", "Bi-contact"),
    "17O": ("BO", "NBO"),
}


@dataclass(frozen=True)
class RatioDef:
    """A named ratio ``Σ w_n·F_n / Σ w_d·F_d`` over family integrals ``F``.

    ``numerator`` / ``denominator`` map family -> weight; ``fmt`` formats the
    value (ratios are fractions, 0-1, or mean numbers); ``description`` is
    the formula a table prints next to the name."""
    name: str
    numerator: dict
    denominator: dict
    fmt: str = "{:.3f}"
    description: str = ""


def _mean_n(name: str, families: Sequence[str], description: str) -> RatioDef:
    """⟨n⟩ = Σ n·Qn / Σ Qn over the Qn families listed."""
    return RatioDef(name, {f: float(f[1:]) for f in families if f != "Q0"},
                    {f: 1.0 for f in families}, "{:.2f}", description)


#: the named-ratio catalogue per nucleus. A ratio is a genuinely derived
#: number (never a repeat of a family row): N4, ⟨CN⟩ Al, ⟨n⟩.
RATIOS: dict[str, tuple[RatioDef, ...]] = {
    "11B": (RatioDef("N4", {"BO4": 1.0}, {"BO3": 1.0, "BO4": 1.0},
                     "{:.3f}", "BO4/(BO3+BO4)"),),
    "27Al": (RatioDef("⟨CN⟩ Al", {"Al(IV)": 4.0, "Al(V)": 5.0, "Al(VI)": 6.0},
                      {"Al(IV)": 1.0, "Al(V)": 1.0, "Al(VI)": 1.0},
                      "{:.2f}", "mean Al coordination number"),),
    "29Si": (_mean_n("⟨n⟩ Si", ("Q0", "Q1", "Q2", "Q3", "Q4"),
                     "mean Qn connectivity"),),
    "31P": (_mean_n("⟨n⟩ P", ("Q0", "Q1", "Q2", "Q3"),
                    "mean Qn connectivity"),),
}

#: batch long-table params that are DERIVED from the fit (never a site's
#: fittable parameter): the per-site population and the two family kinds
DERIVED_PARAMS = frozenset({"population_pct", "family_pct", "ratio"})
#: the subset that describes a group of sites (row ids f<j> / r<j>), sorted
#: after every site column in the wide batch table
GROUP_PARAMS = frozenset({"family_pct", "ratio"})


def _nuc(nucleus) -> str:
    return (nucleus or "").strip()


def presets_for(nucleus) -> list[str]:
    """The preset family tags of a nucleus ([] when none are compiled)."""
    return list(PRESETS.get(_nuc(nucleus), ()))


def ratios_for(nucleus) -> tuple[RatioDef, ...]:
    """The named ratios defined for a nucleus (() when none)."""
    return RATIOS.get(_nuc(nucleus), ())


#: conservative label -> family aliases; anything else stays untagged. The
#: alias must stand alone (no letter/digit glued to it): lookarounds rather
#: than \b, because "]" / ")" at the end of "B[3]" is not a word boundary.
_GUESS = {
    "11B": ((re.compile(r"(?<!\w)(?:bo3|b3|b\[3\]|b\(3\))(?!\w)", re.I), "BO3"),
            (re.compile(r"(?<!\w)(?:bo4|b4|b\[4\]|b\(4\))(?!\w)", re.I), "BO4")),
    "27Al": ((re.compile(r"(?<!\w)(?:al4|aliv|al\(iv\)|al\[4\])(?!\w)", re.I), "Al(IV)"),
             (re.compile(r"(?<!\w)(?:al5|alv|al\(v\)|al\[5\])(?!\w)", re.I), "Al(V)"),
             (re.compile(r"(?<!\w)(?:al6|alvi|al\(vi\)|al\[6\])(?!\w)", re.I), "Al(VI)")),
    "29Si": ((re.compile(r"(?<!\w)q([0-4])(?!\w)", re.I), "Q{}"),),
    "31P": ((re.compile(r"(?<!\w)q([0-3])(?!\w)", re.I), "Q{}"),),
}


def guess_from_label(nucleus, label) -> str:
    """The preset family a site LABEL names ('' unless an alias matches):
    'B[4] (BO4)' -> 'BO4', 'AlIV' -> 'Al(IV)', 'Q3 site' -> 'Q3'; 'P2 chain'
    and 'Czjzek-0' -> ''. Used to order the presets submenu, never to tag a
    line silently."""
    text = str(label or "")
    for rx, fam in _GUESS.get(_nuc(nucleus), ()):
        m = rx.search(text)
        if m:
            return fam.format(*m.groups()) if m.groups() else fam
    return ""


# ------------------------------------------------------------------ summary
def _letters(idx: Sequence[int]) -> str:
    return "+".join(cellparse.index_to_letter(i) for i in idx)


def _family_order(tags: Sequence[str], nucleus) -> tuple[list[str], dict[str, list[int]]]:
    """Family names in display order (nucleus presets first, then other tags
    as first seen) and their member site indices."""
    order = [f for f in presets_for(nucleus) if f in tags]
    for f in tags:
        if f and f not in order:
            order.append(f)
    members = {f: [i for i, t in enumerate(tags) if t == f] for f in order}
    return order, members


def trial_series(samples, families: Sequence[str], nucleus) -> dict[str, np.ndarray]:
    """Per-trial derived quantities from Monte-Carlo per-site integrals
    (``MonteCarloResult.site_integrals``, n_trials x n_sites): the population
    % of every site (``s<i>.population_pct``), of every family
    (``family.<name>``) and every defined ratio (``ratio.<name>``), each as
    an array over the trials whose total (resp. denominator) is positive.
    The Monte-Carlo dialog histograms these; summarize() takes their spread.
    Empty for a malformed block."""
    S = np.abs(np.asarray(samples, float))
    tags = [str(f or "").strip() for f in families]
    if S.ndim != 2 or S.shape[1] != len(tags) or S.shape[0] == 0:
        return {}
    T = S.sum(axis=1)
    ok = T > 0
    if not np.any(ok):
        return {}
    S, T = S[ok], T[ok]
    out: dict[str, np.ndarray] = {}
    for i in range(S.shape[1]):
        out[f"s{i}.population_pct"] = 100.0 * S[:, i] / T
    order, members = _family_order(tags, nucleus)
    Ft = {f: S[:, idx].sum(axis=1) for f, idx in members.items()}
    for f in order:
        out[f"family.{f}"] = 100.0 * Ft[f] / T
    zero = np.zeros(S.shape[0])
    for rd in ratios_for(nucleus):
        present = [f for f in rd.denominator if f in members]
        if len(present) < 2:
            continue
        nt = sum((w * Ft[f] for f, w in rd.numerator.items() if f in members), zero)
        dt = sum((w * Ft[f] for f, w in rd.denominator.items() if f in members), zero)
        good = dt > 0
        if np.any(good):
            out[f"ratio.{rd.name}"] = nt[good] / dt[good]
    return out


def _finite(v) -> bool:
    return v is not None and np.isfinite(v)


def _std(v) -> float | None:
    return float(v.std_dev) if v is not None else None


def summarize(integrals: Sequence[float], amplitudes: Sequence[float],
              families: Sequence[str], nucleus: str, *,
              sigma: Sequence[float | None] | None = None,
              uvars: dict | None = None,
              samples: np.ndarray | None = None) -> dict:
    """Family populations and named ratios from per-site window integrals.

    ``integrals`` / ``amplitudes`` / ``families`` are per site at the best
    fit; ``sigma`` the per-site absolute integral errors (quantify's
    ``integral_err``, None where unknown). One error input picks the basis:
    ``samples`` (n_trials x n_sites per-trial integrals) -> ``montecarlo``;
    else ``uvars`` ({site index: ufloat amplitude}) -> ``covariance``; else
    ``independent``. Returns ``{families, ratios, untagged, basis, note}``;
    an untagged recipe gives ``families == [] and ratios == []``.
    """
    I = np.abs(np.asarray(integrals, float))
    A = np.asarray(amplitudes, float)
    tags = [str(f or "").strip() for f in families]
    n = len(tags)
    untagged = [i for i in range(n) if not tags[i]]
    total = float(I.sum())

    # family order: nucleus presets first, then other tags as first seen
    order, members = _family_order(tags, nucleus)
    if not order:
        return {"families": [], "ratios": [], "untagged": untagged,
                "basis": "independent", "note": ""}

    # ---- which basis
    S = None
    if samples is not None:
        S = np.abs(np.asarray(samples, float))
        if S.ndim != 2 or S.shape[1] != n or S.shape[0] < 2:
            S = None
    U = None
    if S is None and uvars is not None and total > 0:
        # one ufloat (uncertainties, lmfit's dependency) per site, no more
        # checks than the interface the arithmetic below needs
        try:
            U = {i: uvars[i] for i in range(n)}
            if not all(hasattr(u, "std_dev") and hasattr(u, "nominal_value")
                       for u in U.values()):
                U = None
        except (KeyError, TypeError):
            U = None

    # ---- best-fit family integrals and fractions
    F = {f: float(I[idx].sum()) for f, idx in members.items()}
    frac = {f: (100.0 * F[f] / total if total > 0 else 0.0) for f in order}
    sig = list(sigma) if sigma is not None else [None] * n

    fam_err: dict[str, float | None] = {}
    ratio_rows = []
    basis = "independent"
    note_bits: list[str] = []

    series: dict = {}
    if S is not None:
        basis = "montecarlo"
        series = trial_series(S, tags, nucleus)
        for f in order:
            fr = series.get(f"family.{f}")
            fam_err[f] = (float(np.std(fr)) if fr is not None and fr.size >= 2
                          else None)
        note_bits.append(f"Monte-Carlo, {S.shape[0]} trials (per-trial "
                         "re-integrated sums; carries amplitude and shape "
                         "covariance)")
    elif U is not None:
        basis = "covariance"
        # k_i = |I_i| / a_i at the best fit (0 for a zeroed amplitude), so
        # k_i·u_i reproduces |I_i| at the nominal value
        k = np.where(A != 0, I / np.where(A != 0, A, 1.0), 0.0)
        Tu = sum(float(k[i]) * U[i] for i in range(n))
        Fu = {f: sum(float(k[i]) * U[i] for i in idx) for f, idx in members.items()}
        for f in order:
            try:
                fam_err[f] = _std(100.0 * Fu[f] / Tu)
            except (ZeroDivisionError, AttributeError):
                fam_err[f] = None
        note_bits.append("covariance (lmfit uvars: full first-order "
                         "propagation over the amplitude covariance, "
                         "integral/amplitude fixed at the best fit; "
                         "lineshape covariance neglected)")
    else:
        for f, idx in members.items():
            s = [sig[i] for i in idx]
            if any(not _finite(v) for v in s):
                fam_err[f] = None
            else:
                fam_err[f] = (100.0 * float(np.sqrt(sum(v * v for v in s))) / total
                              if total > 0 else None)
        note_bits.append("independent — lines treated as independent "
                         "(quadrature of the site errors); run a fit or "
                         "Monte-Carlo errors for a propagated value")

    fam_rows = [{"family": f, "sites": members[f], "letters": _letters(members[f]),
                 "integral": F[f], "fraction_pct": frac[f],
                 "fraction_err_pct": fam_err[f]} for f in order]

    # ---- named ratios
    for rd in ratios_for(nucleus):
        present = [f for f in rd.denominator if f in members]
        row = {"name": rd.name, "description": rd.description, "fmt": rd.fmt,
               "value": None, "err": None, "defined": False, "note": ""}
        if len(present) < 2:
            missing = [f for f in rd.denominator if f not in members]
            row["note"] = (f"{rd.name} needs " + ", ".join(missing)
                           + " tagged as well")
            ratio_rows.append(row)
            continue
        row["defined"] = True
        num = sum(w * F[f] for f, w in rd.numerator.items() if f in members)
        den = sum(w * F[f] for f, w in rd.denominator.items() if f in members)
        if den <= 0:
            row["note"] = f"{rd.name}: every tagged family integrates to zero"
            ratio_rows.append(row)
            continue
        row["value"] = num / den
        if S is not None:
            rt = series.get(f"ratio.{rd.name}")
            row["err"] = (float(np.std(rt)) if rt is not None and rt.size >= 2
                          else None)
        elif U is not None:
            nu = sum(w * Fu[f] for f, w in rd.numerator.items() if f in members)
            du = sum(w * Fu[f] for f, w in rd.denominator.items() if f in members)
            try:
                row["err"] = _std(nu / du)
            except (ZeroDivisionError, AttributeError):
                row["err"] = None
        else:
            # delta method, families independent: dr/dF_f = (w_n D - N w_d)/D²
            var = 0.0
            for f in present:
                s_f = fam_err[f]
                if s_f is None:
                    var = None
                    break
                s_abs = s_f * total / 100.0           # back to integral units
                d = (rd.numerator.get(f, 0.0) * den - num * rd.denominator.get(f, 0.0)) / den ** 2
                var += (d * s_abs) ** 2
            row["err"] = float(np.sqrt(var)) if var is not None else None
        ratio_rows.append(row)

    if untagged:
        note_bits.append("untagged: " + ", ".join(
            cellparse.index_to_letter(i) for i in untagged)
            + " (counted in the total, in no family)")
    return {"families": fam_rows, "ratios": ratio_rows, "untagged": untagged,
            "basis": basis, "note": "; ".join(note_bits)}
