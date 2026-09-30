# Reading a pulse program

> Every Bruker EXPNO carries the sequence TopSpin ran for it, as a text file
> called `pulseprogram`. **Tools ▸ Pulse program…** (or right-click an EXPNO in
> the Explorer ▸ **Pulse program…**) shows it three ways: a TopSpin-like
> **timing diagram**, the **text** itself, and a **table of the parameters** it
> uses with the values acqus recorded. This page explains the syntax you will
> meet, how to read the diagram, and what the usual solid-state sequences do,
> pulse by pulse. Nothing here writes into the EXPNO.

---

## 1 · What the file is and where TopSpin keeps it

- The **source** of a pulse program lives in TopSpin's library,
  `exp/stan/nmr/lists/pp/` (Bruker's own sequences: `zg`, `cp`, `mp3qdfsz`, …)
  or `exp/stan/nmr/lists/pp/user/` (your facility's, e.g. `hahnecho.nmrfam`).
  The `PULPROG` parameter in `acqus` names the one that was used.
- At acquisition TopSpin copies the **compiled** version next to `acqus` as
  `pulseprogram` (TopSpin 4 may write `pulseprogram.precomp` as well). It is the
  source after the C preprocessor: `#include` files are pasted in, `#ifdef`
  branches are resolved from `ZGOPTNS` (so only the branch that actually ran
  remains), the `mc` macro is spelled out, and every pasted piece is preceded by
  a line such as `# 66 "/opt/topspin3.6.2/exp/stan/nmr/lists/pp/zg"` — a
  **line marker**, not a statement. LARMOR strips those markers from the Text
  tab and keeps everything else.
- The file has four regions: a **header** of `;` comments (what the sequence is,
  who wrote it, and a **legend** of the parameters: `;p1 : 90 degree pulse`),
  the **definitions** (`"p2=p1*2"`, `define delay del7`), the **body** from the
  first label (`1 ze`) to `exit`, and the **phase programs** after `exit`
  (`ph1=0 2 2 0`).

---

## 2 · Syntax primer

### Delays

| Written | Meaning | Value from |
|---|---|---|
| `d1` … `d31` | a delay parameter | `D[n]` in acqus, in **seconds** |
| `10u`, `30m`, `1s` | a literal: 10 µs, 30 ms, 1 s | the number itself |
| `del7`, `zfilter`, `MCWRK` | a delay declared with `define delay` and given by a quoted definition | evaluated from the definition |
| `vd` | the next entry of the variable-delay list (`vdlist`, `VDLIST` in acqus) | the list, one entry per slice |
| `d1 pl1:f1` | a delay during which the power level of channel f1 is set to pl1 | — |
| `MCWRK  * 2` | a delay multiplied | — |

A delay is a gap on the timeline. Tokens written after it on the same line
(`wr #0`, `zd`, `ivd`, `ip1`, `id0`) are housekeeping that happens during that
delay — see the MC macros below.

### Pulses and phases

- `p1 ph1` — a rectangular pulse of length **p1** (`P[1]` in acqus, in **µs**)
  with phase program **ph1**, on channel f1 by default.
- `(p1 pl1 ph1):f1` — the same, spelled fully: power level pl1, phase ph1,
  channel f1. The parentheses group everything that belongs to that pulse.
- `(p1 ph1 d0 p1 ph2):f1` — a sequence of elements on one channel; several
  groups on one line, `(p15 pl1 ph2):f1 (p15:sp0 ph10):f2`, run
  **simultaneously** (a CP contact). `(center … )` centres the groups on each
  other (an HMQC's 180° in the middle of the t1 pair).
- `p1*0.5 ph1` — a pulse of half the length; `vp` — the next entry of a
  variable-pulse list.
- A phase written `ph1:r` also resets the receiver phase reference.

### Power levels

- `pl0` … `pl31` are the power levels (`PLW[n]` in acqus, in **watts**; TopSpin
  also shows them in dB: PLdB = −10·log₁₀(PLW / 1 W), so **a larger dB number is
  less power** — 250 W is −24 dB, 7.5 W is −8.8 dB, 0.24 W is +6 dB). `sp0` …
  `sp31` are the powers of shaped pulses (`SPW[n]`).
- A level is *set* on a channel with `pl1:f1` (usually on the recycle delay:
  `d1 pl1:f1 pl2:f2`) or carried by the pulse itself (`(p1 pl1 ph1):f1`). Once set
  it stays until changed, which is why a decoupling program started with
  `cpds2:f2` runs at whatever level was last set on f2 (`pl12:f2` a line above).
- The legend at the top of the file says what each level is for (`;pl21 : power
  level for selective pulse, ca. plw11 / 1024`).

### Channels

`f1` is the observe channel (`NUC1`, frequency `SFO1`), `f2` and `f3` the
decoupler / second and third nuclei (`NUC2`, `NUC3`). A `:fN` suffix selects the
channel: `(p3 pl2 ph1):f2`, `cpds2:f2`, `pl12:f2`. The Diagram tab labels each
row with the nucleus acqus recorded for it.

### Shaped pulses

`(p2:sp1 ph2):f1` is a shaped pulse: length p2, **shape file `spnam1`**
(`gauss`, `ramp.100`, `dfs`, `wurst`, …), power `spw1`, frequency offset
`spoffs1`. The diagram draws a shaped pulse with an envelope (a ramp for
`ramp.*`, a swept oscillation for a DFS / WURST sweep, a bell otherwise) and
names the shape under it. Sweeps take their frequencies from constants the
legend documents (`cnst1`, `cnst2` in kHz for a DFS).

### Decoupling

- `cpds2:f2` (or `cpd2:f2`) starts the composite-pulse decoupling program
  **`cpdprg2`** (`spinal64`, `tppm15`, `cw`, …) on f2, with the pulse length
  **`pcpd2`** and the power last set on f2. The `s` form restarts the program
  from its first pulse each time. `cw:f2` is continuous-wave decoupling.
- `do:f2` (**d**ecoupler **o**ff) stops it. The diagram draws decoupling as a
  hatched bar from the element it starts on to the `do`; a `cpd2:f2` written
  on the `go=` line runs during the acquisition only.

### Loops, labels and counters

- A **label** is a number or a name at the start of a line (`1 ze`, `2 d1`,
  `LBLF1, MCWRK`, `Passaq, 1m`); it marks a place to jump back to.
- `lo to 3 times l20` repeats everything since label 3, **l20** times (`L[20]`
  in acqus): a saturation comb, an echo train, a REDOR loop.
- `go=2` loops to label 2 by itself, **NS** times (after DS dummy scans);
  `rcyc=2` does the same for a scan started with `adc`.
- `lo to LBLF1 times td1` and `lo to LBLAV times tdav` are the 2D increment
  and averaging loops the MC macro inserts (below).
- Variable lists: `vd` / `vp` / `vc` take the next entry of the vdlist /
  vplist / vclist; `ivd` / `ivp` / `ivc` advance the pointer (once per slice,
  in the housekeeping delay), `rvd` resets it. `lo to 3 times vc` uses the
  counter list as a loop count (a CPMG T₂ series).
- `goto Passaq`, `if "aq < 50.1m" goto Passaq`, `print "…"`: the `aq_prot`
  guard many solids programs include — it refuses to run with an acquisition
  time above 50 ms (duty cycle on the probe).

### Acquisition and the receiver phase

- `go=2 ph31` — acquire one FID with **receiver phase ph31**, add it to
  memory, step every phase program, loop to label 2; NS scans in all.
  `gosc` / `goscnp` acquire a single scan without the loop (`np`: without
  stepping the phases). The acquisition lasts **aq** = TD / (2·SW_h) after the
  pre-scan delay `de`.
- `adc ph31` starts the digitiser explicitly; the delays that follow are
  acquisition windows (a QCPMG echo train records through its 180° pulses
  this way), `eosc` / `eoscnp` close the scan and `rcyc=2` loops for the next.
- `acqt0` is the time between the effective start of the signal and the first
  point, used by *baseopt* digitisation to make the first-order phase zero:
  `"acqt0=-p1*2/3.1416"` in `zg` places time zero 2/π of the way into the 90°
  pulse (its centre of effect); an echo sequence sets `acqt0=0` and aims the
  echo top at the first point instead.

### Phase cycling — what `ph31` means

- A phase program lists the phase of a pulse scan by scan, in **units of 90°**:
  `0 1 2 3` is x, y, −x, −y. `ph1=(12) 0 2 4 6 8 10` changes the unit to
  360°/12 = 30°, so the list is 0°, 60°, 120°, … ; `(360)` writes degrees.
  `{0}*24 {1}*24` repeats an entry, `{0 2}^1^2^3` appends the list shifted by
  1, 2 and 3 units.
- `ph31` (sometimes `ph30`) is the **receiver** phase. A pulse whose phase steps
  through N values of 360°/N, with the receiver following −Δp·φ, keeps only the
  signal whose coherence order changed by Δp (modulo N) across that pulse —
  the coherence-pathway rule (Bodenhausen, Kogler & Ernst 1984). The six 60°
  steps of `ph1` in `mp3qdfsz` select Δp = ±3: triple-quantum excitation.
  The 4-step `0 2 2 0 1 3 3 1` of `zg` is CYCLOPS plus a 180° alternation,
  cancelling receiver imbalance and DC offsets.
- The **cycle length** is the longest list; shorter lists repeat. NS should be
  a multiple of it (the legend usually says so: `;ns : 96 * n`).
- `ip1` advances phase program ph1 by one step for the next increment — how
  States and TPPI acquisitions build the second dimension — and `rp1` resets
  it.

### `define`, quoted definitions and the MC macros

- `define delay del7`, `define loopcounter ST1CNT`, `define pulse pX` declare
  names; `"del7=((1s*l1)/cnst31)-(p2/2)"` gives them (or a standard parameter)
  a value at compile time. Inside the quotes, numbers carry their units
  (`1s`, `30m`, `10u`), `cnst31` is `CNST[31]`, `l1` is `L[1]`, and the usual
  arithmetic applies; `trunc(td1 / 2)` rounds down. LARMOR evaluates these
  definitions from the acqus values and, when a definition overrides a
  parameter that acqus also stores (`"p2=p1*cnst2"`), shows the computed value
  and notes the acqus one if the two differ.
- `"in0=inf1"` ties the t1 increment to the F1 spectral width (`INF[1]`, in µs).
- `#include <Avance.incl>` pulls in the standard macro definitions; `prosol
  relations=<solids_default>` names the prosol table `getprosol` reads to
  fill pulses and powers from the probe file; `dccorr` switches on the
  digitiser's DC correction.
- The **`mc` macro** — `30m mc #0 to 2 F1PH(ip1, id0)` in the source — is
  TopSpin's shorthand for the housekeeping of a 2D or arrayed acquisition. In
  the compiled file it is spelled out as delays **`MCWRK`** (a fraction of the
  30 ms) and **`MCREST`** (the remainder), labels `LBLF1` (the F1 increment
  loop, `td1` times), `LBLSTS1` (the States pair: two scans with `ip1` in
  between, hence `ST1CNT = trunc(td1 / 2)`), `LBLAV` (`tdav` averaging), `LBLF0`
  (`td0` repetitions of a 1D), and the tokens `wr #0` (write the FID),
  `if #0` (move to the next FID of the ser), `zd` (zero the buffer), `id0` (add
  `in0` to d0 — the t1 step), `ivd` (next vd entry), `rf #0` (rewind). None of
  it is part of the spin physics: the diagram draws these delays small and the
  summary line leaves them out.

---

## 3 · How to read the diagram

- **One row per rf channel** in use, labelled with the nucleus from acqus, and
  a bottom row **acq** for the receiver. Time runs left to right, **not to
  scale**: a hard pulse is one nominal unit wide, a shaped pulse or a named
  delay two, the acquisition four — the way TopSpin's own pulse-program display
  draws it. A 5 s recycle delay and a 20 µs z-filter are the same width; read
  the numbers, not the lengths.
- A filled **box** is a pulse: its name above (`p1`), its phase program below
  (`ph1`), then the resolved length (`2.1 µs`) and the power level it carries
  (`pl11 250 W`). A box with a **curved top** is a shaped pulse, its shape file
  named under it.
- A **gap** with a name and a value is a delay; a small tick with `pl1` is a
  power setting; a note under a delay (`wr #0 zd`, `+in0 per increment`) is
  the housekeeping or the increment that rides on it.
- A **hatched bar** is decoupling, running until its `do`; the bar names the
  cpd program and its pulse length.
- The **decaying oscillation** on the acq row is the acquisition, labelled with
  the command and the receiver phase (`go=2`, `ph31`) and the acquisition time.
- **Brackets** under the rows are loops, each with its count: `× l20 = 200`,
  `× ns = 512`, `× td1 = 22`. The shortest loop sits nearest the rows.
- A grey note under the brackets says when a value could not be resolved
  (a definition that uses a constant acqus does not carry) and which symbols it
  needs; the Parameters tab lists every symbol with its source (acqus,
  definition, literal, list).
- **Export…** (or right-click the diagram) saves the figure as PNG or SVG,
  or copies it to the clipboard.

---

## 4 · Worked examples

### `zg` — one pulse, acquire

```
  d1
  p1 ph1
  go=2 ph31
```

Recycle delay **d1** (1–5 T₁ for a quantitative spectrum), a pulse **p1** at
power **pl1**, acquisition. Look at: `p1` and `pl1` against the calibrated
90° (the title line often records the flip angle), `d1` against T₁, `ns`
against the 8-step cycle. `acqt0=-p1*2/3.1416` shows this dataset was
digitised with *baseopt*.

### Hahn echo — `hahnecho.nmrfam`

```
  d1
  (p1 pl1 ph1):f1
  d6
  (p2 ph2):f1
  d7
  go=2 ph31
```

90° – τ – 180° – τ' – acquire, rotor-synchronised. The definitions carry the
physics: `"p2=p1*cnst2"` makes the second pulse a 180° (`cnst2 = 2`) or a
second 90° for a solid echo (`cnst2 = 1`); `"d6=((1s*l1)/cnst31)-(p1/2)-(p2/2)"`
is **l1** rotor periods at the spinning frequency **cnst31**, minus half of each
pulse, so the echo forms on a rotor echo; `del7` shows the ideal second delay
and `d7` is the value actually used (tuned so the first point sits on the echo
top: "half-echo acquisition"). `ph2` runs 0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3
against `ph1 = 0 1 2 3` — a 16-step cycle that keeps only the refocused
pathway. Check: `l1 / cnst31` against the MAS rate, and `d7` against `del7`.

### Saturation recovery — `satrect1`

```
3 d20
  (p1 pl1 ph4):f1
  lo to 3 times l20
  vd
  (p1 pl1 ph1):f1
  go=2 ph31
```

A comb of **l20** pulses spaced **d20** destroys the magnetisation, the
variable delay **vd** lets it recover, a 90° reads it out. The `ivd` in the
housekeeping steps the vdlist once per slice and `lo to LBLF1 times td1` runs
the **td1** slices — a pseudo-2D whose `ser` the Relaxation tool
(Tools ▸ Relaxation) fits. Check: `l20 × d20` comfortably longer than T₂, the
vdlist span against the expected T₁, `d1` (which need only be short here).

### QCPMG and WCPMG

```
  (p1 ph1):f1
  d3
  (p2 ph2):f1
  d3
  1u adc ph31 syrec
  d6
3 (p2 ph2):f1
  d6
  lo to 3 times l22
  rcyc=2
```

An echo is formed, the digitiser is started once (`adc`), and a train of
**l22** 180° pulses refocuses the signal again and again while the receiver
records through the whole train (`d6` windows either side of each pulse).
The spectrum is a comb of spikelets spaced by 1 / (echo period); the
Tools ▸ QCPMG page co-adds or Fourier-transforms the echoes. WCPMG replaces the
rectangular pulses with WURST sweeps (`(p1:sp1 ph1):f1`, `(p2:sp2 ph2):f1`)
so that a pattern hundreds of kHz wide is excited uniformly. Check: `l22` and
`d6` against `aq` (the train must fit in TD points), `d3` against `d6` (the
definition usually ties them), the sweep width of `sp1` against the pattern.

### CP/MAS — `cp`

```
  (p3 pl2 ph1):f2
  (p15 pl1 ph2):f1 (p15:sp0 ph10):f2
  1u cpds2:f2
  go=2 ph31
  1m do:f2
```

A 1H 90° (**p3** at **pl2** on f2), then the **contact**: for **p15** the X
channel is on at **pl1** while 1H is on with the ramp shape **sp0**
(`ramp.100`) — the two groups on one line are simultaneous, and the diagram
aligns their boxes. `cpds2:f2` starts 1H decoupling (`cpdprg2`, e.g. SPINAL-64,
at **pl12** with **pcpd2**) just before the acquisition and `do:f2` stops it
after. `ph2 = 0 0 2 2` alternates the spin temperature. Check: the
Hartmann–Hahn match (pl1 against sp0), `p15` against the build-up time, `d1`
against the **1H** T₁ (not the X one), `pcpd2` against the 1H 180°.

### 3QMAS with z-filter — `mp3qdfsz`

```
  d1
  (p1 pl11 ph1):f1
  d0
  (p2:sp1 ph2):f1
  d10
  (p3 pl21 ph3):f1
  d4
  (p3 ph4):f1
  go=1 ph31
```

Four pulses on one channel. **p1** at the high power **pl11** excites
triple-quantum coherence; **d0** is t1, incremented by **in0** (one rotor period
here: `"in0=inf1"`, rotor-synchronised); the shaped **p2:sp1** is a
double-frequency sweep (`spnam1 = dfs`, from `cnst1` to `cnst2` kHz) that
converts 3Q into observable coherence with more signal than a hard pulse —
its length is defined as `1s/(cnst31*cnst0)`, a fraction 1/cnst0 of the rotor
period; **d10** is the split-t1 delay (0 for spin 5/2); **p3** at **pl21**,
about 30 dB lower (`plw11 / 1024`), is a central-transition-selective 90°;
**d4** is the z-filter (20 µs); the second **p3** reads out. `ph1=(12) 0 2 4 6
8 10` selects ±3Q, `ph4` runs 24 steps, and the 48-step `ph31` follows; the
`ip1` in the housekeeping and the `lo to LBLSTS1 times 2` pair implement
States, so `td1` rows give `td1 / 2` complex t1 points. Check: `p1` against
the ν₁ needed for 3Q excitation, `pl21` ≈ `pl11` / 1024, `in0 × cnst31 = 1`,
`ns` a multiple of 96. Shearing and referencing of the result are on the
MQMAS page.

---

## 5 · Where the numbers come from

| Symbol | acqus array | Unit in the file | Shown as |
|---|---|---|---|
| `p0`–`p31` | `P` | µs | µs |
| `d0`–`d31` | `D` | s | s / ms / µs |
| `pl0`–`pl31` | `PLW` (older data: `PL`, dB only) | W | W, with dB in the comment |
| `sp0`–`sp31`, `spnam` | `SPW`, `SPNAM` | W, file name | W · shape |
| `l0`–`l31` | `L` | count | count |
| `cnst0`–`cnst31` | `CNST` | as documented in the legend | number |
| `in0`–`in31`, `inf1` | `IN` (s), `INF` (µs) | | s / µs |
| `pcpd2`, `cpdprg2` | `PCPD`, `CPDPRG` | µs, program name | |
| `ns`, `ds`, `td`, `td0`, `td1` | `NS`, `DS`, `TD`, `TD0`, `acqu2s TD` | count | |
| `aq`, `dw` | derived: TD / (2·SW_h), 1 / (2·SW_h) | | s |
| `vd`, `vp`, `vc` | `VDLIST`, `VPLIST`, `VCLIST` | list file | the list's name |

A quoted definition wins over the stored array (TopSpin recomputes it when the
experiment starts); the comment column then shows the definition and, if it
disagrees with the array, the array's value too. A definition that needs a
constant this acqus does not carry is left blank and named in the note under
the diagram — never guessed.

---

## References

- Bruker BioSpin, *TopSpin Pulse Programming Manual* (User Guide), the
  reference for the statements, the `mc` macro and the phase-program syntax.
- J. Keeler, *Understanding NMR Spectroscopy*, 2nd ed., Wiley (2010) — pulse
  sequences, coherence orders and phase cycling.
- G. Bodenhausen, H. Kogler, R. R. Ernst, "Selection of coherence-transfer
  pathways in NMR pulse experiments", *J. Magn. Reson.* **58**, 370 (1984).
- A. Paterson, the documented NMRFAM solid-state sequences (`hahnecho.nmrfam`
  and its siblings in the facility's `pp/user` library, National Magnetic
  Resonance Facility at Madison) — an example of a pulse program whose header
  states its instructions, parameters and calculations.
- E. L. Hahn, "Spin echoes", *Phys. Rev.* **80**, 580 (1950).
- H. Y. Carr, E. M. Purcell, *Phys. Rev.* **94**, 630 (1954); S. Meiboom,
  D. Gill, *Rev. Sci. Instrum.* **29**, 688 (1958). *(CPMG)*
- F. H. Larsen, H. J. Jakobsen, P. D. Ellis, N. C. Nielsen,
  "Sensitivity-enhanced quadrupolar-echo NMR of half-integer quadrupolar
  nuclei", *J. Phys. Chem. A* **101**, 8597 (1997). *(QCPMG)*
- L. A. O'Dell, R. W. Schurko, "QCPMG using adiabatic pulses for faster
  acquisition of ultra-wideline NMR spectra", *Chem. Phys. Lett.* **464**, 97
  (2008). *(WCPMG)*
- A. Pines, M. G. Gibby, J. S. Waugh, "Proton-enhanced NMR of dilute spins in
  solids", *J. Chem. Phys.* **59**, 569 (1973). *(cross polarisation)*
- B. M. Fung, A. K. Khitrin, K. Ermolaev, "An improved broadband decoupling
  sequence for liquid crystals and solids", *J. Magn. Reson.* **142**, 97
  (2000). *(SPINAL-64)*
- L. Frydman, J. S. Harwood, "Isotropic spectra of half-integer quadrupolar
  spins from bidimensional magic-angle spinning NMR", *J. Am. Chem. Soc.*
  **117**, 5367 (1995). *(MQMAS)*
- J.-P. Amoureux, C. Fernandez, S. Steuernagel, "Z filtering in MQMAS NMR",
  *J. Magn. Reson. A* **123**, 116 (1996).
- D. Iuga, H. Schäfer, R. Verhagen, A. P. M. Kentgens, "Population and
  coherence transfer induced by double frequency sweeps in half-integer
  quadrupolar spin systems", *J. Magn. Reson.* **147**, 192 (2000). *(DFS)*
- D. J. States, R. A. Haberkorn, D. J. Ruben, "A two-dimensional nuclear
  Overhauser experiment with pure absorption phase in four quadrants",
  *J. Magn. Reson.* **48**, 286 (1982). *(States acquisition)*

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
