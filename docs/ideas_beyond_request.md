# 20 further ideas for LARMOR (beyond what was asked)

Focused on the paramount goal — **physically relevant lineshapes and fits** — and
on making the batch/series workflow something you can lean on for real work.

## Physically-honest fits (the core)
1. **Global constraint library.** Save/name reusable constraint sets ("tie all
   δ_iso to site 0", "C_Q(field2) = C_Q(field1)") and apply them across batches so
   the physics you trust is one click, not re-typed per project.
2. **Prior/penalty terms.** Optional soft priors on parameters (e.g. C_Q ~ N(μ,σ)
   from a known crystal) added to the residual, so a glassy fit is nudged toward
   physically sensible values instead of overfitting noise.
3. **Model-selection scores.** Report AIC/BIC and reduced-χ² per fit and per batch,
   with an "is the extra site justified?" F-test, so you don't add lines the data
   can't support.
4. **Residual autocorrelation flag.** Detect structured (non-white) residuals — the
   tell-tale of a missing/oversimplified lineshape — and warn on the plot.
5. **Parameter identifiability map.** From the fit covariance, flag pairs that are
   >0.95 correlated (e.g. width vs C_Q at one field) so you know what the data
   truly constrains vs what is degenerate.
6. **Lineshape sanity checks.** Guard against unphysical results (η>1, negative
   width, populations summing to ≠100 %, a site fully outside the window) with a
   pre-fit and post-fit checklist.

## Batch / series power tools
7. **Series-aware bounds.** Let parameters vary *monotonically* or *smoothly* along
   the series (a light penalty on jumps), which is often the real physics of a
   composition/temperature series.
8. **Per-spectrum quality gates.** Auto-flag spectra with poor SNR, truncated
   windows, or failed convergence; exclude/re-fit them without redoing the batch.
9. **Two-way batch table** ↔ live plots: click a row → the spectrum's plot
   highlights; drag a parameter in the series plot → re-fit that point.
10. **Batch templates.** Save a whole batch setup (model, window, baseline type,
    release set, threshold) as a named template for the next series.
11. **Population from integrals, in-batch.** Offer the same integral-over-window
    quantification as Report (F6) as the population column, not just amplitude
    fraction, with error propagation.
12. **Bootstrapped series errors.** A Monte-Carlo pass across the whole batch so the
    series-evolution plot carries honest error bands, not just covariance bars.

## Plotting / publication
13. **Figure style presets** (Nature/ACS/thesis) — fonts, sizes, tick density,
    palette — one click to a journal's spec, saved with the project.
14. **Annotation layer.** Draggable text, arrows, site labels, and ppm brackets on
    any plot, persisted with the figure spec.
15. **Reference overlays.** Drop a known-compound spectrum or literature C_Q/η
    markers onto a plot for visual assignment.
16. **Difference/stack normalisation modes.** Area-, max-, or noise-normalised
    stacks and difference plots for series comparison.

## Data & workflow
17. **Project bundle.** One `.larproj` capturing data references, every fit,
    baselines, batches and figures, so a paper is fully reproducible and shareable.
18. **Watch-folder / queue.** Point LARMOR at an acquisition folder; new EXPNOs
    appear in the Explorer and can be auto-fit against a template as they land.
19. **Undo history for fits.** A timeline of fit states (like processing undo) so
    you can step back to a previous solution and branch.
20. **Headless CLI / notebook API.** `larmor batch-fit …` and a thin Python API so
    a whole series can be fit and tabulated in a script or Jupyter, reusing the
    exact engines behind the GUI (great for large studies and CI).

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
