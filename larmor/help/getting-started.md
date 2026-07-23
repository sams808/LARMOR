# Getting started with LARMOR

> **LARMOR** — *Lineshape Analysis & Refinement for Magnetic-resonance Of solids*.
> An open desktop successor to dmfit for solid-state NMR: dmfit's fitting
> fluency, ssNake's processing depth, and the **mrsimulator** physics stack, with
> uncertainties on every number and reproducible, diffable recipes.

This guide is the front door. It explains the ideas the whole program is built
on — the **workspace model**, the **Explorer**, and the *open-anything* import —
and then points you to the manual for whatever experiment is in front of you.

---

## The three ideas

| Idea | What it means for you |
|---|---|
| **Open anything** | 1D, 2D, a raw FID or `ser`, a dmfit `.fxmla`, a plain CSV — LARMOR reads it, shows it, and *then* offers the right processing. Nothing is rejected at the door. |
| **Workspaces** | Every dataset you open lives in its own switchable workspace (like TopSpin windows). Extract a trace from a 2D and it becomes a new workspace; the map stays where it was. |
| **Reproducible recipes** | Processing and fitting are stored as a **recipe** replayed from the raw data. A saved recipe *is* the analysis — it reopens exactly, and it diffs cleanly in version control. |

Every fitted number is reported **with its uncertainty**, or with the reason
there is none. Instrument folders are opened **read-only** — LARMOR never writes
into your Bruker data.

---

## 1 · Opening data

**File ▸ Open**, the **Explorer** dock on the left, or **drag-and-drop** onto the
plot. LARMOR's universal reader accepts:

| Source | Notes |
|---|---|
| Bruker `1r` / `2rr` | processed 1D / 2D spectrum (or point at the `pdata` folder, or the EXPNO) |
| Bruker `fid` / `ser` | raw time-domain — opens with a processing preview so you apodize/phase, then transform |
| dmfit `.fxmla` | a dmfit model + data (1D and MQMAS) imported directly |
| LARMOR `.json` | a saved recipe (data reference + full processing + model) |
| `.csv` / `.txt` / `.dat` | a two-column *(shift, intensity)* spectrum, with an optional `# key=value` metadata header |

The reader detects **1D vs 2D**, **time vs frequency** domain, real vs
**hypercomplex** (it loads the `2rr`/`2ri`/`2ir`/`2ii` quadrants when present),
and it distinguishes a genuine spectroscopic 2D from a **pseudo-2D** relaxation
series (a `vdlist`/`vclist` with a placeholder F1 axis). What you opened decides
where it lands:

- a **1D frequency** spectrum → the workbench (process & fit);
- a **2D** dataset → the contour view;
- a raw **FID** → a magnitude-FT preview so you can see it before you phase;
- a raw **`ser`** → a 2D preview (or the guided relaxation tool if it's a series).

## 2 · The Explorer

The Explorer browses a sample folder and **auto-identifies** each experiment —
nucleus, 1D/2D, and kind (single pulse, MQMAS, satrec, QCPMG…) — by reading
`acqus` and the pulse program, so you can see at a glance what an EXPNO holds
before opening it. **File ▸ Open sample** points it at a folder.

## 3 · Workspaces

The **Workspaces** dock lists everything you have open. Each entry carries an
icon for its kind (bare spectrum, spectrum-with-fit, 2D map). Click to switch,
**Close** to free it, **Save** to write its recipe. Opening a 2D, or pulling a
trace off a map, spawns a **new** workspace so your existing fit is never
disturbed; **Back to 2D map** (Ctrl + 2) returns to the parent map. Snapshots are
lightweight (arrays + a recipe reference), and closing a workspace frees it.

## 4 · Processing, fitting, tools

Once your data is open, the workflow is the same shape everywhere — process,
then fit, then measure/export — but the details depend on the experiment. The
**Process** panel is live (edits re-apply as you type); the **Fit-Parameters**
table at the bottom is a dmfit-style spreadsheet with paddles on the plot; the
**Tools** and **Decomposition** menus hold the experiment-specific machinery.

Each experiment has its own manual with worked steps and the science behind it:

| If you have… | Read |
|---|---|
| a normal 1D spectrum | **1D spectra — processing & fitting** |
| a raw FID / echo you must transform | **1D spectra** + **Processing reference** |
| a QCPMG echo train (broad lines) | **QCPMG** |
| a T₁/T₂ relaxation series | **Relaxation (T1/T2)** |
| a 2D MQMAS map | **MQMAS (2D)** |
| an HMQC / DQ-SQ correlation | **HMQC & correlation** |
| any 2D to phase / measure | **2D processing** |
| several datasets to compare or co-fit | **Multi-dataset & co-fitting** |
| a question about a lineshape model | **Lineshapes — models & physics** |
| a question about a processing step | **Processing reference** |

All are under **? ▸ User manuals**; the reference documents are direct **?**-menu
items. Every tool with a **Help** button opens the matching section.

---

## Conventions used everywhere

- **Chemical shift** increases to the **left** (decreasing frequency, IUPAC δ
  scale). In 2D, F2 (direct) high-shift is left, F1 (indirect) high-shift is top.
- **Quadrupolar coupling** is reported as $C_Q = e^2qQ/h$ (MHz) and asymmetry
  $\eta_Q\in[0,1]$; for disordered sites the **quadrupolar product**
  $P_Q = C_Q\sqrt{1+\eta_Q^2/3}$ is the invariant that keeps its meaning.
- **Shielding anisotropy** uses the **Haeberlen** convention (δ_iso, ζ, η_CS).
- **Widths** may be entered in **ppm or Hz** (`300Hz`, `1.5kHz`, `2ppm`);
  LARMOR converts with the Larmor frequency.

---

## References

- D. Massiot, F. Fayon, M. Capron, I. King, S. Le Calvé, B. Alonso, J.-O. Durand,
  B. Bujoli, Z. Gan, G. Hoatson, "Modelling one- and two-dimensional solid-state
  NMR spectra", *Magn. Reson. Chem.* **40**, 70 (2002). *(dmfit)*
- S. G. J. van Meerten, W. M. J. Franssen, A. P. M. Kentgens, "ssNake: A
  cross-platform open-source NMR data processing and fitting application",
  *J. Magn. Reson.* **301**, 56 (2019). *(ssNake)*
- D. J. Srivastava, P. J. Grandinetti *et al.*, **mrsimulator**,
  github.com/deepanshs/mrsimulator. *(the physics engine)*
- M. H. Levitt, *Spin Dynamics: Basics of Nuclear Magnetic Resonance*, 2nd ed.,
  Wiley (2008). *(general reference)*

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
