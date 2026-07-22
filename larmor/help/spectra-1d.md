# 1D spectra — processing & fitting

Everything you do to a normal 1D spectrum: open it, process it, fit it.

## Open
File ▸ Open (or the Explorer on the left) accepts dmfit `.fxmla`, a LARMOR
recipe `.json`, a Bruker `1r`/`fid`, or a two-column `.csv`/`.txt`. A raw `fid`
opens with Process ▸ Open FID for apodization/phasing before the transform.

## Process (right-hand panel — live)
- **WDW / LB / GB / SSB** — apodization (EM, GM, sine bells, TRAF).
- **Phase** — p0/p1 sliders with **−90 / +90 / 180°** quick steps, and Autophase.
- **Baseline** — automatic (polynomial order) or manual PCHIP anchors.
- **SR / Calibrate** — reference the axis: type an SR, or Process ▸ Calibrate,
  click a peak, and set its ppm. Double-click the experiment strip to edit
  nucleus / field / νrot / SR (or copy SR from another spectrum).
The pipeline is stored in the recipe and replayed on reload, so a fit is
reproducible from the raw data.

## Fit
1. Pick a model from the **Models** menu (or the toolbar), then **click the
   spectrum** to drop lines — the mode stays on, so place as many as you want;
   click the model again (or Esc) to stop.
2. Drag the **paddles** (position/amplitude square, width side-handles) or edit
   the **Fit-Parameters** table: a value, a bound `[0..100]`, a link to another
   line by its letter (`A`, `A+20`, `A+20kHz`, `0.5B`), or pin ☑ to fix.
3. **Fit** (F5). Read **RMSD** and **S/N** next to the buttons; **Report** (F6)
   gives populations (% ± error).
4. **Zones** restrict the fit to chosen regions; **Auto Fit** does multi-start;
   **Errors Analysis** profiles χ² around a parameter.

## Models (which lineshape)
- **Gauss/Lorentz** — pseudo-Voigt (gl = 1 Gaussian … 0 Lorentzian).
- **Gauss/Lorentz (area)** — same, but amplitude is the *integral* (populations
  read straight off the amplitudes).
- **Voigt (true)** — Gaussian ⊗ Lorentzian with independent widths.
- **J-multiplet** — binomial multiplet split by J (Hz).
- **Czjzek / ext. Czjzek** — amorphous quadrupolar (glasses).
- **Quad CT / Quad 1st / Quad CT+CSA** — crystalline quadrupolar.
- **CSA powder** — shielding anisotropy with spinning sidebands.
- **Spectrum (background)** — fit another measured spectrum's amplitude/shift.

## Measure & export
- Tools ▸ **Integrals & measurements** — drag regions → integral, %, centre,
  FWHM; Copy CSV.
- File ▸ **Copy plot** (with all lines) / **Save plot image** for slides.
- **Save fit as** (txt/csv/json/dmfit) writes data + model + residual + each
  component; **Save spectrum as** writes a reopenable CSV.
