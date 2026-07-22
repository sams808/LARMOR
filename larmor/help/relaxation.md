# Relaxation series (T1 / T2)

Tools ▸ Relaxation opens the guided workflow on a pseudo-2D (`ser` + vdlist).

## Steps
1. **Process slices** — every slice is EM/ZF/FT'd and phased on the
   most-relaxed slice; the relaxed spectrum is shown.
2. **Slices [first]…[last]** — restrict the range: trim early-stopped or
   already-relaxed acquisitions.
3. **Integration zones** — drag one or more regions over your peak(s); each
   becomes a build-up curve.
4. **Fit T1/T2** — mono-exponential by default (TopSpin's 3-parameter form);
   tick **stretched β** for a distribution. Click a point to **exclude an
   outlier** and it refits. Toggle **log delay** and use **Fit view** to
   rescale to the points.
5. **Copy CSV** for the delays, integrals and fitted τ.

Kinds: satrec / invrec / cpmg / t1rho, auto-detected from the pulse program.
