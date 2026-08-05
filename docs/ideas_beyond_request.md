# Further ideas for LARMOR — status

Focused on the paramount goal — **physically relevant lineshapes and fits** — and
on making the batch/series workflow something you can lean on for real work.

**Status legend:** ✅ implemented · ⏳ roadmap (not yet built) · ✋ intentionally
left to the human (per your call: too sample/system/model-dependent to automate).

## Physically-honest fits (the core)
1. ✅ **Global constraint library.** Save/name reusable constraint sets and apply
   them across models — *Decomposition ▸ Advanced ▸ Save / Apply constraints*
   (`larmor/constraint_library.py`).
2. ✋ **Prior/penalty terms.** Soft priors from a known crystal — the choice of
   prior and its strength is a physical judgement; kept manual.
3. ✋ **Model-selection scores (AIC/BIC, F-test).** Whether an extra site is
   "justified" depends on the system and the question; kept a human decision.
4. ✋ **Residual autocorrelation flag.** What counts as structured residual vs a
   fine you accept is model-dependent; kept manual. *(A basic "residual N× noise"
   note already appears after each fit.)*
5. ✅ **Parameter identifiability map.** Pairs with |r| ≥ 0.95 flagged as
   unidentifiable after every fit and in the Correlations dialog
   (`larmor/identifiability.py`).
6. ✅ **Lineshape sanity checks.** η∈[0,1], positive widths, non-negative
   amplitudes, sites inside the window — flagged (never auto-corrected) after each
   fit (`larmor/sanity.py`).

## Batch / series power tools
7. ⏳ **Series-aware bounds** (monotonic/smooth penalty). *Partly delivered:* the
   sequential fit's between-pass **trajectory smoothing** regularises toward
   smooth series; an explicit monotonic penalty is still roadmap.
8. ✅ **Per-spectrum quality gates.** Batch fit flags RMSD outliers and low-S/N
   spectra (red cell + reason tooltip).
9. ⏳ **Two-way batch table ↔ live plots** (click a row → highlight the spectrum;
   drag a point → re-fit).
10. ✅ **Batch templates.** Save/load a named batch setup (release set, ±, baseline,
    threshold, toggles).
11. ⏳ **Population from integrals, in-batch** (the F6 quantification as the batch
    population column, with error propagation).
12. ⏳ **Bootstrapped series errors** (a Monte-Carlo pass across the whole batch for
    honest error bands on the evolution plot).

**NEW — ✅ Sequential (forward–backward) fit.** Your idea: fit an end-member, carry
its parameters to the next spectrum, sweep to the far end and back to smooth, with
automatic *N*-pass sweeps (1/2/4/8/16), live RMSD evolution and trajectory
smoothing (`larmor/seqfit.py`, *Tools ▸ Sequential fit*, `larmor seqfit` CLI).

## Plotting / publication
13. ✅ **Figure style presets.** `nature` / `acs` / `rsc` (plus article/thesis/
    presentation) in the Plotting studio.
14. ⏳ **Annotation layer** (draggable text/arrows/site labels persisted with the
    spec). *(Spec-level `annotations` already render; interactive editing is
    roadmap.)*
15. ⏳ **Reference overlays** (known-compound spectrum / literature C_Q·η markers).
    *Partly delivered:* the studio already overlays any spectrum/fit as a trace and
    draws iso/quad reference lines in 2D.
16. ✅ **Difference/stack normalisation modes.** `max` / `area` / `noise`
    normalisation and *difference vs first trace* in the studio.

## Data & workflow
17. ⏳ **Project bundle** (`.larproj` capturing data refs, fits, baselines, batches,
    figures). *(A `.larproj` loader exists for fits; a full bundle is roadmap.)*
18. ⏳ **Watch-folder / auto-fit queue.**
19. ⏳ **Undo history for fits** (a timeline of fit states, branchable).
20. ✅ **Headless CLI / notebook API.** `larmor batchfit …` and `larmor seqfit …`
    reuse the exact GUI engines for scripted/large studies and CI.

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
