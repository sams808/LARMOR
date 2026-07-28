# Further improvement ideas (from field use, 2026-07)

Concrete directions to make LARMOR better, gathered after real fitting sessions.
Ranked loosely by value-to-effort for solid-state NMR of glasses/materials.

## Fitting & physics

1. **Derived-quantity panel for quadrupolar sites.** Extend the new Czjzek
   Cq/νQ read-out into a full derived table per site: **P_Q**, √⟨P_Q²⟩, the
   second-order isotropic shift δ₂, and Cq/η — shown live and written into the
   report/CSV. For a disordered site, report only the invariants (P_Q), clearly.

2. **Plot the fitted Czjzek distribution.** Alongside the lineshape, draw the
   actual P(ν_Q, η) the fitted σ implies (the Eq.-6 PDF from the manual). The
   user *sees* the physical distribution, not just a width — pedagogically strong
   and it validates the model choice.

3. **Correlation & error visualisation.** Show the lmfit covariance as a
   correlation heat-map, and 2D χ² maps for parameter pairs (δiso–Cq, σ–dCS).
   Makes "which parameters the data actually determines" visual, not just a ±.

4. **Residual diagnostics.** Overlay the residual with a ±noise band, show the
   reduced χ², and flag *residual structure* (autocorrelation / runs test) so a
   systematically-wrong model is caught even when the RMSD looks small.

5. **Spinning-sideband auto-detect.** Find peaks at ±k·νrot and offer to add
   them as a linked `sidebands` manifold in one click (dmfit "ss band" comfort).

6. **MQMAS auto-referencing presets** per spin/method (3Q/5Q for ²⁷Al/¹¹B/²³Na):
   pick nucleus + coherence and the F1 shear/scale is set correctly, removing the
   manual F1-reference step for standard cases.

## Workflow & throughput

7. **Batch fitting over a series.** Apply one saved recipe to a folder of
   spectra (a Ca/Na composition series) and get a table of populations and
   parameters vs sample, with error bars — directly serves the glass work.

8. **Project/session files.** Save an entire session (all workspaces, overlays,
   fits, processing) as one file that reopens exactly, not just per-spectrum
   recipes.

9. **Compare against a reference fit.** Overlay a previously-saved fit or a
   literature dmfit `.fxmla` on the current data, with a side-by-side parameter
   diff — quick sanity-check against published values.

10. **Copy-ready outputs.** "Copy fit summary" → a formatted (site, δiso, P_Q,
    %, ±err) table on the clipboard, plus a **LaTeX table** and a short methods
    paragraph for papers/notebooks (ties into the figure studio).

11. **Crash-safe autosave.** Snapshot the session every few minutes and offer
    recovery on restart, so a crash never loses an hour of fitting.

## Interaction polish (TopSpin/dmfit parity)

12. **Mouse phasing gesture.** On top of the new pivot line, add TopSpin's
    drag-to-phase: horizontal drag = p0, vertical = p1, about the pivot — the
    fastest way to phase by hand.

13. **Keyboard-driven fitting.** Shortcuts to nudge the selected parameter,
    cycle sites, toggle pin, and fit — power-user speed like dmfit/TopSpin.

14. **Peak-pick + assignment table.** A pickable, labelable peak list (AlIV /
    AlV / AlVI…) that exports and doubles as fit starting points — bridges
    "look at the spectrum" and "place lines".

15. **Full undo for processing & referencing.** Undo currently covers the model
    only; extend it to the experiment axis (SR/calibrate) and the processing
    pipeline so every action is reversible.

## Reach & robustness

16. **More input formats.** JCAMP-DX, Magritek/Spinsolve, Varian/Agilent FIDs —
    so collaborators aren't Bruker-only.

17. **Per-nucleus smart defaults.** Opening a ²⁷Al glass pre-populates sensible
    Czjzek starting points and kernel resolution; remember the user's typical
    Cq ranges and preferred models per nucleus.

18. **Live fit animation.** Draw intermediate model curves during the fit (not
    just the progress bar) so convergence — or divergence — is visible as it
    happens; great for teaching and for stopping a runaway fit early.

19. **dmfit amplitude calibration.** A one-time "calibrate export to dmfit":
    paste one dmfit amp for a known exported line and LARMOR locks the exact
    Czjzek seed-amplitude factor, removing the reconstruction guesswork.

20. **Automated figure + report.** One-click "publication figure + parameter
    table + methods sentence" for a fitted spectrum.

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
