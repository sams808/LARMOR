"""Engine stage ``engine-records`` of the distribution check: the records a
fit leaves and the tables a user publishes, each checked against something
outside the code under test -- the raw acqus / procs / title / pulse-program
files read with a plain regex, the closed-form areas of Gaussian lines, the
first-order error formulas, numpy, IUPAC / CODATA constants, round trips.

Two halves, one stage: ``_records_acquisition`` (larmor.acquisition,
comparability, inventory, referencing, provenance, pulseprog,
quantitativity on the bundled EXPNOs; synthetic EXPNOs with the same values
when the examples are absent) and ``_records_tables`` (quantify + families,
methods, recipe + recipe_diff, paramstatus on the shipped fits, project,
fittable, series_table, series_grid, aliases). Qt-free, deterministic, every
file under ``ctx["out"]``; see ``larmor.distcheck`` for the contract.
"""
from __future__ import annotations

import contextlib
import os
import time
from pathlib import Path

import numpy as np

from larmor.distcheck.engines import (ISO, _check, _examples, _expno, _import_mrsimulator, _out,
                                      _write_csv_spectrum)

__all__ = ["records_stage"]


def records_stage(say, ctx):
    """The acquisition records (Experimental sentences, Table S1, the
    comparability verdict, the session inventory and referencing audit, the
    source hash and reload checks, the pulse program, the quantitativity
    judgement) and the published tables (populations, families and ratios,
    the Methods text, recipes and their diff, the project bundle, the
    parameter and series tables, the batch grid resolver, the alias store)."""
    _records_acquisition(say, ctx)
    _records_tables(say, ctx)


# ------------------------------------------------------------ records: acquisition
def _racq_jcamp(path) -> dict:
    """A Bruker JCAMP parameter file read with a plain regex -- independently
    of nmrglue, the reader every module under test goes through: scalars as
    written ('<zg>' keeps its brackets), a '(0..N)' array as the list of the
    tokens on the lines that follow it."""
    import re

    text = re.sub(r"\$\$[^\n]*", "", Path(path).read_text(encoding="latin-1"))      # '$$' comments
    out = {}
    for m in re.finditer(r"^##\$(\w+)=[ \t]*(.*?)\s*(?=^##|\Z)", text, re.M | re.S):
        arr = re.match(r"\(0\.\.\d+\)\s*(.*)", m.group(2), re.S)
        out[m.group(1)] = arr.group(1).split() if arr else m.group(2).strip()
    return out


def _racq_truth(e: Path) -> dict:
    """What one EXPNO's own files state, read without LARMOR: acqus (scalars,
    the D / P / PLW / SPW arrays, TE), pdata/1/procs, the title's 'MAS=<n> kHz',
    acqu2s, the TopSpin version on the acqus TITLE line, and from the pulse
    program text the power token of p1 and the number of pulse statements."""
    import re

    acq, procs = _racq_jcamp(e / "acqus"), _racq_jcamp(e / "pdata" / "1" / "procs")
    title = (e / "pdata" / "1" / "title").read_text(encoding="utf-8", errors="replace")
    head = (e / "acqus").read_text(encoding="latin-1").splitlines()[0]
    prog = (e / "pulseprogram").read_text(encoding="latin-1")

    def arr(k):
        return [float(x) for x in acq.get(k, [])]

    mas = re.search(r"\bMAS\s*=\s*([\d.]+)\s*kHz", title)
    body = prog[re.search(r"^\s*(?:\d+\s+)?ze\s*$", prog, re.M).end():]
    body = body[:re.search(r"^exit\s*$", body, re.M).start()]
    p1_line = re.search(r"^\s*\(?\s*p1\b[^\n]*", body, re.M).group(0)
    p1_pl = re.search(r"\bpl(\d+)\b", p1_line)
    t = {"nucleus": acq["NUC1"].strip("<>"), "nuc2": acq["NUC2"].strip("<>"),
         "bf1": float(acq["BF1"]), "bf2": float(acq["BF2"]), "sfo1": float(acq["SFO1"]),
         "sw": float(acq["SW_h"]), "td": int(acq["TD"]), "ns": int(acq["NS"]), "ds": int(acq["DS"]),
         "rg": float(acq["RG"]), "date": int(acq["DATE"]), "masr": float(acq["MASR"]),
         "pulprog": acq["PULPROG"].strip("<>"), "probe": acq["PROBHD"].strip("<>"),
         "te": float(acq["TE"]) if "TE" in acq else None,
         "d": arr("D"), "p": arr("P"), "plw": arr("PLW"), "spw": arr("SPW"),
         "topspin": re.search(r"TopSpin\s+([\d.]+)", head).group(1),
         "sf": float(procs["SF"]), "lb": float(procs["LB"]), "wdw": int(procs["WDW"]),
         "si": int(procs["SI"]), "tdeff": int(procs["TDeff"]), "ph_mod": int(procs["PH_mod"]),
         "absg": int(procs["ABSG"]), "phc0": float(procs["PHC0"]), "phc1": float(procs["PHC1"]),
         "title": title, "title_mas": float(mas.group(1)) * 1e3 if mas else None,
         "flip_claim": bool(re.search(r"P1\(90\)|\d\s*(?:deg|°)", title, re.I)),
         "twod": (e / "acqu2s").is_file(), "has_1r": (e / "pdata" / "1" / "1r").is_file(),
         "auditp": (e / "pdata" / "1" / "auditp.txt").is_file(), "prog": prog,
         "p1_power": int(p1_pl.group(1)) if p1_pl else 1,          # no token: the default pl1
         "n_pulses": len(re.findall(r"^\s*\(?\s*p\d+\b", body, re.M))}
    t["sr"] = (t["sf"] - t["bf1"]) * 1e6          # TopSpin: SR = (SF - BF1) * 1e6 Hz
    t["aq"] = t["td"] / (2.0 * t["sw"])           # TD real points at a dwell of 1 / (2 SW_h)
    t["td1"] = int(_racq_jcamp(e / "acqu2s")["TD"]) if t["twod"] else None
    return t


def _racq_write_jcamp(path: Path, params: dict) -> Path:
    """A TopSpin-style JCAMP parameter file; a dict value is a sparse
    '(0..63)' array (numbers, or '<name>' strings)."""
    lines = ["##TITLE= Parameter file, TopSpin 3.6.3", "##JCAMPDX= 5.0",
             "##DATATYPE= Parameter Values", "##ORIGIN= Bruker BioSpin GmbH", "##OWNER= nmr"]
    for k, v in params.items():
        if isinstance(v, dict):
            blank = "<>" if any(isinstance(x, str) for x in v.values()) else 0
            lines += [f"##${k}= (0..63)", " ".join(str(v.get(i, blank)) for i in range(64))]
        else:
            lines.append(f"##${k}= {v}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines + ["##END="]) + "\n", encoding="ascii")
    return path


def _racq_copy(src: Path, dst: Path, sf: float | None = None) -> Path:
    """A copy of an EXPNO's acqus and pdata/1 procs / title / 1r under the
    stage folder (the instrument files are never touched); ``sf`` rewrites
    the copy's procs SF."""
    import re
    import shutil

    for rel in ("acqus", "pdata/1/procs", "pdata/1/title", "pdata/1/1r"):
        if (src / rel).is_file():
            (dst / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src / rel, dst / rel)
    if sf is not None:
        p = dst / "pdata" / "1" / "procs"
        p.write_text(re.sub(r"(?m)^##\$SF= .*$", f"##$SF= {sf!r}", p.read_text(encoding="latin-1")),
                     encoding="latin-1")
    return dst


def _racq_synthetic(root: Path) -> Path:
    """``root/pCABS2-4``: three EXPNOs with the bundled data's acquisition,
    processing, title and pulse-program values -- an 11B zg whose acqus MASR
    (30 kHz) contradicts its title (20 kHz), a 27Al zg, a 27Al mp3qdfsz 3QMAS
    (acqu2s, ser, 2rr) -- for an installation without the examples."""
    zg = (";zg\n;avance-version (12/01/11)\n;1D sequence\n;$CLASS=HighRes\n;$DIM=1D\n\n"
          '"acqt0=-p1*2/3.1416"\n; dimension 1D; AQ_mode\ndefine delay MCWRK\ndefine delay MCREST\n'
          '"MCREST = 30m - 30m"\n"MCWRK = 0.333333*30m"\n\n    dccorr\n1 ze\n2 MCWRK  * 2\n'
          "LBLF0, MCWRK\n  MCREST\n  d1\n  p1 ph1\n  go=2 ph31\n  MCWRK wr #0\n  MCWRK zd\n"
          "  lo to LBLF0 times td0\n\nexit\n\nph1=0 2 2 0 1 3 3 1\nph31=0 2 2 0 1 3 3 1\n"
          ";pl1 : f1 channel - power level for pulse (default)\n;p1 : f1 channel -  high power pulse\n"
          ";d1 : relaxation delay; 1-5 * T1\n;ns: 1 * n, total number of scans: NS * TD0\n")
    mq = (";mp3qdfsz\n;3Q MAS pulse program for odd half integer spin nuclei\n;pl1 : =0W, not used\n"
          ";pl11 : power level for excitation pulse\n;pl21 : power level for selective pulse\n"
          ";spw1 : power level for frequency sweep\n;p1 : excitation pulse at plw11\n"
          ";p2 : =1s/(cnst31*l0) duration of sweep at spw1\n;p3 : 90 degree selective pulse at plw21\n"
          ";d1 : recycle delay\n;$CLASS=Solids\n;$DIM=2D\n;$SUBTYPE=MQMAS\n\nprosol relations=<solids_mqmas>\n\n"
          '"p2=1s/(cnst31*cnst0)"\n"acqt0=1u*cnst11"\n"in0=inf1"\n\n; dimension 2D; AQ_mode  (F1) States\n'
          'define delay MCWRK\ndefine delay MCREST\ndefine loopcounter ST1CNT\n"ST1CNT = trunc(td1 / 2)"\n'
          '"MCREST = 10m - 10m"\n"MCWRK = 0.250000*10m"\n\n    dccorr\n  ze\nLBLAV, MCWRK\n1 MCWRK \n'
          "LBLSTS1, MCWRK  * 2\nLBLF1, MCWRK\n  MCREST\n  d1\n  (p1 pl11 ph1):f1\n  d0\n  (p2:sp1 ph2):f1\n"
          "  d10\n  (p3 pl21 ph3):f1\n  d4\n  (p3 ph4):f1\n  go=1 ph31\n  MCWRK  wr #0 if #0 zd ip1\n"
          "  lo to LBLSTS1 times 2\n  MCWRK  rp1 id0  MCWRK  id10\n  lo to LBLF1 times ST1CNT\n"
          "  MCWRK rf #0\n  lo to LBLAV times tdav\nexit\n\nph1=(12) 0 2 4 6 8 10\n"
          "ph2= {0}*24 {1}*24 {2}*24 {3}*24\nph3= 0\nph4= 0 0 0 0 0 0 1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3\n"
          "ph31= 0 2 0 2 0 2 1 3 1 3 1 3 2 0 2 0 2 0 3 1 3 1 3 1\n"
          "      2 0 2 0 2 0 3 1 3 1 3 1 0 2 0 2 0 2 1 3 1 3 1 3\n")
    common = {"NUC2": "<1H>", "BF2": 500.13, "INSTRUM": "<spect>", "RG": 287, "TE": 298.8605,
              "PROBHD": "<H8984_0013 (PH MASDVT 500W2 BL2.5 N-P/F-H)>", "BYTORDA": 0, "DTYPA": 0,
              "TD0": 1, "TDav": 1, "PARMODE": 0}
    al = {"NUC1": "<27Al>", "BF1": 130.318169, "SFO1": 130.32296462, "O1": 4795.62, "MASR": 26000}
    spec = {
        1118: ({"NUC1": "<11B>", "BF1": 160.461579, "SFO1": 160.461579, "O1": 0, "SW_h": 59523.8095238095,
                "TD": 2048, "NS": 512, "DS": 0, "DATE": 1643982460, "MASR": 30000, "PULPROG": "<zg>",
                "D": {1: 5, 2: 0.00345}, "P": {1: 2.1}, "PLW": {0: 30, 1: 5.3}},
               {"SF": 160.460649890015, "SI": 4096, "OFFSET": 191.2682, "SW_p": 59523.8095238095,
                "LB": 50, "TDeff": 1536, "ABSG": 2, "PHC0": 130.1749, "PHC1": -109.95},
               "02/12/21 : conditions OK\nd1=5s MAS=20 kHz\npi/12 a 20 kHz OK", zg),
        3616: ({**al, "SW_h": 100000, "TD": 2048, "NS": 2048, "DS": 0, "DATE": 1643970823,
                "PULPROG": "<zg>", "D": {1: 1, 2: 0.00345}, "P": {1: 2}, "PLW": {0: 30, 1: 7.5}},
               {"SF": 130.317561599976, "SI": 2048, "OFFSET": 425.1385, "SW_p": 100000, "LB": 100,
                "TDeff": 800, "ABSG": 0, "PHC0": 22.76813, "PHC1": -197.05},
               "02/12/21 : conditions OK\nMAS=26 kHz\npure Ca - 2GPa - sans iode\nexc pi/12 a 20 kHz", zg),
        3620: ({**al, "SW_h": 119047.619047619, "TD": 1024, "NS": 3600, "DS": 96, "DATE": 1644264246,
                "PULPROG": "<mp3qdfsz>", "PARMODE": 1, "D": {0: 1e-06, 1: 0.5, 4: 2e-05, 10: 1e-07},
                "P": {1: 3.1, 2: 9.615385, 3: 4}, "PLW": {0: 30, 11: 250, 21: 7.5}, "SPW": {0: 60, 1: 220},
                "SPNAM": {0: "<ramp.100>", 1: "<dfs>"}, "CNST": {0: 4, 11: 1, 31: 26000}},
               {"SF": 130.317561599976, "SI": 1024, "OFFSET": 498.2201, "SW_p": 119047.619047619,
                "LB": 200, "TDeff": 2048, "ABSG": 5, "PHC0": -89.8125, "PHC1": 0},
               "02/12/21 : conditions OK\n2D MAS=26 kHz", mq)}
    folder = root / "pCABS2-4"
    for n, (acq, procs, title, prog) in spec.items():
        e = folder / str(n)
        _racq_write_jcamp(e / "acqus", {**common, **acq})
        _racq_write_jcamp(e / "pdata" / "1" / "procs", {"BYTORDP": 0, "DTYPP": 0, "NC_proc": 0, "XDIM": 0,
                                                         "WDW": 1, "GB": 0, "SSB": 0, "PH_mod": 1, "FT_mod": 6,
                                                         **procs})
        (e / "pdata" / "1" / "title").write_text(title, encoding="utf-8")
        (e / "pulseprogram").write_text(prog, encoding="ascii")
        if acq.get("PARMODE"):
            _racq_write_jcamp(e / "acqu2s", {"NUC1": acq["NUC1"], "TD": 22, "SW_h": 26041.6666666667,
                                             "FnMODE": 4})
            (e / "ser").write_bytes(bytes(64))
            (e / "pdata" / "1" / "2rr").write_bytes(bytes(64))
        else:
            x = np.arange(procs["SI"])
            y = 1e6 * np.exp(-0.5 * ((x - procs["SI"] / 3.0) / 12.0) ** 2)
            (e / "pdata" / "1" / "1r").write_bytes(y.astype("<i4").tobytes())
            (e / "fid").write_bytes(np.zeros(acq["TD"], "<i4").tobytes())
    return folder


def _racq_acquisition(say, E, T) -> dict:
    """larmor.acquisition: the record, the sentences, Table S1, the paragraph."""
    import csv
    import datetime as dt
    import json
    import re

    from larmor import acquisition as A

    gamma_1h = 42.577478518                 # MHz/T: the proton's gamma / 2 pi (CODATA 2018)
    wdw = {0: "none", 1: "EM", 2: "GM", 3: "SINE", 4: "QSINE"}
    ph_mod = {0: "no", 1: "pk", 2: "mc", 3: "ps"}
    blocks, rates = {}, {}
    for n, e in E.items():
        t = T[n]
        b = blocks[n] = A.read_block(e)
        json.dumps(b)                       # every value a JSON scalar
        _check(set(b) == set(A.BLOCK_KEYS), f"read_block({n}) returned " + (
            "{} (no readable acqus)" if not b else f"keys off BLOCK_KEYS: {sorted(set(b) ^ set(A.BLOCK_KEYS))}"))
        want = {"nucleus": t["nucleus"], "bf1_MHz": t["bf1"], "sfo1_MHz": t["sfo1"], "sf_MHz": t["sf"],
                "ns": t["ns"], "ds": t["ds"], "d1_s": t["d"][1], "p1_us": t["p"][1], "plw1_W": t["plw"][1],
                "sw_Hz": t["sw"], "td": t["td"], "pulprog": t["pulprog"], "probe": t["probe"],
                "topspin": t["topspin"], "mas_acqus_Hz": t["masr"], "mas_title_Hz": t["title_mas"],
                "lb_hz": t["lb"], "si": t["si"], "tdeff": t["tdeff"], "absg": t["absg"],
                "wdw": wdw.get(t["wdw"]), "ph_mod": ph_mod.get(t["ph_mod"]),
                "date_utc": dt.datetime.fromtimestamp(t["date"], dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
        bad = {k: (b.get(k), v) for k, v in want.items() if b.get(k) != v}
        _check(not bad, f"read_block({n}) disagrees with its own acqus / procs / title (record, file): {bad}")
        # the magnet: the module goes through BF1 / Xi(nucleus); acqus' BF2 is the 1H channel itself
        b0 = t["bf2"] / gamma_1h
        _check(abs(b["sr_hz"] - t["sr"]) < 1e-6 and abs(b["aq_s"] - t["aq"]) < 1e-12 and t["nuc2"] == "1H"
               and abs(b["magnet_1h_MHz"] - t["bf2"]) < 1e-6 * t["bf2"] and abs(b["b0_T"] - b0) < 1e-4,
               f"read_block({n}): SR {b['sr_hz']} vs (SF - BF1) 1e6 = {t['sr']} Hz, AQ {b['aq_s']} vs "
               f"{t['aq']} s, 1H {b['magnet_1h_MHz']} vs BF2 {t['bf2']} MHz, B0 {b['b0_T']} vs {b0} T")
        # the MAS rule: two sources within 2 % settle it, two disagreeing keep the highest, flagged
        srcs = [v for v in (t["masr"], t["title_mas"]) if v]
        agree = len(srcs) < 2 or abs(srcs[0] - srcs[1]) <= 0.02 * max(srcs)
        rates[n] = (srcs[0] if agree else max(srcs), not agree)
        _check((b["spin_rate_Hz"], b["mas_uncertain"]) == rates[n],
               f"read_block({n}): MAS {b['spin_rate_Hz']} Hz uncertain {b['mas_uncertain']}, but acqus says "
               f"{t['masr']:.0f} Hz and the title {t['title_mas']} Hz")
        if not t["flip_claim"]:
            _check(b["flip_deg"] is None and b["flip_source"] == "",
                   f"read_block({n}): flip {b['flip_deg']} ({b['flip_source']}) without any title evidence")
        say(f"read_block {n}: {t['nucleus']} {t['pulprog']} SFO1 {t['sfo1']} MHz, D1 {t['d'][1]:g} s, "
            f"P1 {t['p'][1]:g} µs at PLW1 {t['plw'][1]:g} W, NS {t['ns']}, TD {t['td']}, SR {b['sr_hz']:.2f} Hz "
            f"= (SF - BF1)·1e6, B0 {b['b0_T']:.4f} T (BF2 {t['bf2']} MHz -> {b0:.4f}); MAS acqus "
            f"{t['masr']:.0f} / title {t['title_mas']:.0f} Hz -> {b['spin_rate_Hz']:.0f} Hz"
            + (" flagged" if b["mas_uncertain"] else "") + "; all = the files read independently")
    # ---- the Experimental sentences: what the files support, nothing else
    mas_said = {}
    for n, b in blocks.items():
        t, (rate, unc) = T[n], rates[n]
        sents = A.sentences(b)
        txt = " ".join(sents)
        mas_said[n] = sents[1] if len(sents) > 1 else ""
        musts = [f"The {t['nucleus']} MAS NMR spectrum was acquired", f"B₀ = {t['bf2'] / gamma_1h:.1f} T",
                 f"{t['sf']:.2f} MHz for {t['nucleus']}", f"a {t['p'][1]:g} µs pulse",
                 f"a recycle delay of {t['d'][1]:g} s", f"{t['ns']} transients",
                 f"The spectral width was {t['sw'] / 1e3:g} kHz (TD = {t['td']} points)",
                 "SR = " + f"{t['sr']:.2f}".replace("-", "−") + " Hz [state the reference standard]",
                 f"TopSpin {t['topspin']}", f"zero-filling to {t['si']} points (SI)",
                 f"[confirm MAS rate: acqus {t['masr'] / 1e3:.3g} kHz, title {t['title_mas'] / 1e3:.3g} kHz]"
                 if unc else f"The magic-angle spinning rate was {rate / 1e3:.1f} kHz."]
        if t["wdw"] == 1:
            musts.append(f"exponential apodization (LB = {t['lb']:g} Hz)")
        if 0 < t["tdeff"] < t["td"]:
            musts.append(f"the first {t['tdeff']} points of the fid (TDeff)")
        if t["p1_power"] == 1:              # only where p1 plays at the default pl1 is PLW1 its power
            musts.append(f"a {t['p'][1]:g} µs pulse at {t['plw'][1]:g} W")
        missing = [m for m in musts if m not in txt]
        _check(not missing, f"sentences({n}) lack {missing}: {txt}")
        # no title words for these, no flip evidence, no 1H reference; acqus TE is never a temperature
        never = ["decoupl", "ambient temperature", "adamantane"]
        never += [] if t["flip_claim"] else ["flip angle", "°"]
        never += ["spinning rate was"] if unc else []
        never += [str(t["te"]), f"{t['te']:g}", f"{t['te']:.1f}"] if t["te"] else []
        found = [w for w in never if w in txt]
        _check(not found and not re.search(r"\b\d+(?:\.\d+)?\s*K\b", txt),
               f"sentences({n}) state what the files do not support: {found} / a temperature: {txt}")
        lines = A.summary_lines(b)
        mas_line = next((ln for ln in lines if ln.startswith("MAS:")), "")
        _check(f"NS {t['ns']}" in lines[0] and f"D1 {t['d'][1]:g} s" in lines[0]
               and f"acqus {t['masr']:.0f} Hz" in mas_line and ("(confirm)" in mas_line) == unc,
               f"summary_lines({n}): {lines}")
    say(f"sentences {sorted(E)}: nucleus, B0 / SF, P1, D1, NS, SW / TD, SR + '[state the reference standard]', "
        f"TopSpin processing as in the files; MAS {'; '.join(f'{n} {mas_said[n]!r}' for n in sorted(E))}; "
        "no flip angle / decoupling / temperature (acqus TE) / adamantane claimed")
    # ---- Table S1: rows in (nucleus, date) order, the varying columns marked
    tab = A.table(list(blocks.values()))
    order = sorted(E, key=lambda n: (T[n]["nucleus"], T[n]["date"]))
    _check([int(r["expno"]) for r in tab.rows] == order, f"table rows {[r['expno'] for r in tab.rows]} != {order}")
    keymap = {"nucleus": "nucleus", "sfo1_MHz": "sfo1", "probe": "probe", "pulprog": "pulprog",
              "p1_us": ("p", 1), "plw1_W": ("plw", 1), "d1_s": ("d", 1), "ns": "ns", "sw_Hz": "sw",
              "td": "td", "tdeff": "tdeff", "si": "si", "lb_hz": "lb", "absg": "absg", "wdw": "wdw",
              "ph_mod": "ph_mod"}

    def fv(n, spec):
        return T[n][spec[0]][spec[1]] if isinstance(spec, tuple) else T[n][spec]

    vary = {k for k, s in keymap.items() if len({fv(n, s) for n in E}) > 1}
    probe = T[order[0]]["probe"] if "probe" not in vary else None
    _check(tab.varying & set(keymap) == vary and tab.constants.get("probe") == probe,
           f"Table S1 varying {sorted(tab.varying & set(keymap))} vs the files {sorted(vary)}; "
           f"constants {tab.constants}")
    rows = list(csv.reader(A.table_csv(tab).strip().splitlines()))
    head = rows[0]
    c_d1, c_ns, c_pr = head.index("D1 (s)"), head.index("NS"), head.index("Probe")
    _check(len(rows) == len(E) + 2 and rows[-1][0] == "# varies"
           and [r[c_d1] for r in rows[1:-1]] == [f"{T[n]['d'][1]:g}" for n in order]
           and [r[c_ns] for r in rows[1:-1]] == [str(T[n]["ns"]) for n in order]
           and rows[-1][c_d1] == ("yes" if "d1_s" in vary else "") and rows[-1][c_pr] == "",
           f"table_csv: {rows}")
    tex = A.table_latex(tab, caption="Acquisition parameters.")
    ncol = len({"sample", "expno", "nucleus"} | tab.varying)
    _check(max(map(ord, tex)) < 128 and "\\begin{tabular}{" + "l" * ncol + "}" in tex
           and tex.count("\\\\\n") >= len(E) + 1 and (probe is None or "Probe = " + probe.replace("_", "\\_") in tex),
           "table_latex: not pure ASCII, or the column count / rows / constant caption are wrong")
    md = A.table_markdown(tab).splitlines()
    _check(len(md) == len(E) + 2 and ("**D1 (s)**" in md[0]) == ("d1_s" in vary), f"table_markdown: {md[:2]}")
    say(f"Table S1: rows {order} by (nucleus, date); varying {sorted(vary)} = the files' own differences; "
        f"constant probe {probe!r}; CSV D1 column {[r[c_d1] for r in rows[1:-1]]}, LaTeX {ncol} columns "
        f"pure ASCII, Markdown {len(md)} lines")
    # ---- the paragraph: one per nucleus, a parameter that varies as a range
    paras = A.paragraph(list(blocks.values())).split("\n\n")
    nucs = sorted({T[n]["nucleus"] for n in E})
    _check(len(paras) == len(nucs), f"paragraph: {len(paras)} paragraphs for nuclei {nucs}")
    said = []
    for p, nuc in zip(paras, nucs):
        grp = sorted(n for n in E if T[n]["nucleus"] == nuc)
        d1s = sorted({T[n]["d"][1] for n in grp})
        rng = (f"recycle delays of {d1s[0]:g}–{d1s[-1]:g} s (Table S1)" if len(d1s) > 1
               else f"a recycle delay of {d1s[0]:g} s")
        were = "spectra were acquired" if len(grp) > 1 else "spectrum was acquired"
        _check(p.startswith(f"The {nuc} ") and rng in p and were in p, f"paragraph {nuc}: no {rng!r} / {were!r}: {p}")
        said.append(f"{nuc} {grp}: '{rng}'")
    say(f"paragraph: {len(paras)} paragraphs, one per nucleus; " + "; ".join(said))
    return blocks


def _racq_comparability(say, E, T, out):
    """larmor.comparability on two EXPNOs and on one EXPNO against itself."""
    from larmor import comparability as C

    keys = {"NUC1": "nucleus", "SFO1": "sfo1", "PULPROG": "pulprog", "SW_h": "sw", "TD": "td",
            "D1": ("d", 1), "P1": ("p", 1), "PLW1": ("plw", 1), "NS": "ns", "RG": "rg", "PROBHD": "probe",
            "WDW": "wdw", "LB": "lb", "TDeff": "tdeff", "SI": "si", "PH_mod": "ph_mod", "ABSG": "absg",
            "SF": "sf"}

    def fv(n, spec):
        return T[n][spec[0]][spec[1]] if isinstance(spec, tuple) else T[n][spec]

    pa, pb = C.read_params(E[1118]), C.read_params(E[3616])
    for p, n in ((pa, 1118), (pb, 3616)):
        bad = {k: (p.value(k), fv(n, s)) for k, s in keys.items() if p.value(k) != fv(n, s)}
        _check(not bad and p.commands == [] and not T[n]["auditp"],
               f"read_params({n}) disagrees with its files (read, file): {bad}; commands {p.commands}")
    cmp = C.compare([pa, pb], ["11B", "27Al"])
    # every differing value here differs far beyond its rule's tolerance (SFO1 / D1 5 %, others 1e-6)
    want = {k for k, s in keys.items() if fv(1118, s) != fv(3616, s)}
    got = {r.key for r in cmp.reports if r.differs} & set(keys)
    nuc = cmp.report("NUC1")
    _check(got == want and cmp.level == "bad" and nuc is not None and nuc.level == "bad"
           and "mixed nuclei" in nuc.message,
           f"compare(1118, 3616): differing {sorted(got)} vs the files {sorted(want)}, level {cmp.level}")
    cmp.to_csv(out / "compare_1118_3616.csv")
    n_csv = len((out / "compare_1118_3616.csv").read_text(encoding="utf-8").splitlines())
    _check(n_csv == 1 + len(cmp.reports), f"Comparison.to_csv wrote {n_csv} lines for {len(cmp.reports)} keys")
    t = T[3616]
    same = C.compare([pb, C.read_params(E[3616])])
    summ = same.summary()
    _check(same.level == "ok" and not same.flagged() and not any(r.differs for r in same.reports)
           and summ.startswith("✓ 2 spectra") and f"D1 {t['d'][1]:g} s" in summ and f"NS {t['ns']}" in summ,
           f"compare(3616, 3616): level {same.level}, flagged {[r.key for r in same.flagged()]}: {summ}")
    # TopSpin's procs as LARMOR ops: TDeff counts real points (complex fid -> TDeff // 2), the stored
    # phase replays as phase(-PHC0, PHC1) pivoted at 1.0, the carrier sits at (SFO1 - SF) / SF
    ops = {o["op"]: o for o in C.ops_for(pb)}
    want_ops = {"em": {"op": "em", "lb_hz": t["lb"]}, "zf": {"op": "zf", "si": t["si"]},
                "phase": {"op": "phase", "p0": -t["phc0"], "p1": t["phc1"], "pivot_frac": 1.0}}
    if 0 < t["tdeff"] < t["td"]:
        want_ops["tdeff"] = {"op": "tdeff", "points": t["tdeff"] // 2}
    carrier = (t["sfo1"] - t["sf"]) / t["sf"] * 1e6
    _check(all(ops.get(k) == v for k, v in want_ops.items()) and abs(ops["ft"]["offset_ppm"] - carrier) < 1e-4,
           f"ops_for(3616): {list(ops.values())} vs {list(want_ops.values())}, carrier {carrier:.4f} ppm")
    say(f"comparability: 1118 vs 3616 level {cmp.level}, differing {sorted(got)} = the files' differences "
        f"(flagged {len(cmp.flagged())}); 3616 vs itself '{summ[:60]}...'; ops_for -> "
        f"{[o['op'] for o in C.ops_for(pb)]}, carrier {ops['ft']['offset_ppm']:.4f} ppm = (SFO1 - SF)/SF")


def _racq_inventory(say, E, T, root, out):
    """larmor.inventory on the session folder that holds the three EXPNOs."""
    import re

    from larmor import inventory

    inv = inventory.build(root)
    rows = {int(r.info.expno): r for r in inv.rows}
    _check(set(rows) == set(E), f"inventory.build({root.name}): EXPNOs {sorted(rows)} != {sorted(E)}")
    bad = {n: (r.nucleus, r.info.ndim, r.ns, r.d1_s) for n, r in rows.items()
           if (r.nucleus, r.info.ndim, r.ns, r.d1_s) != (T[n]["nucleus"], 2 if T[n]["twod"] else 1,
                                                         T[n]["ns"], T[n]["d"][1])}
    _check(not bad, f"inventory rows (nucleus, ndim, NS, D1) disagree with the files: {bad}")
    picks = sorted(r.expno for r in inv.picks())
    # one sample, one 1D EXPNO with a processed 1r per nucleus: each is its block's production spectrum
    want = sorted(n for n in E if not T[n]["twod"] and T[n]["has_1r"])
    one_d = sorted({T[n]["nucleus"] for n in E if not T[n]["twod"]},       # mass order, 1H last
                   key=lambda s: (s == "1H", int(re.match(r"\d+", s).group()), s))
    _check(picks == want and all(rows[n].role == "2D" for n in E if T[n]["twod"]) and inv.nuclei() == one_d,
           f"inventory: picks {picks} (files: {want}), roles {[(n, r.role) for n, r in rows.items()]}, "
           f"nuclei {inv.nuclei()}")
    no_1h = not any(T[n]["nucleus"] == "1H" for n in E)
    _check(not no_1h or all(r.sr_verdict == "no reference" for r in inv.rows),
           f"inventory SR verdicts {[(n, r.sr_verdict) for n, r in rows.items()]} in a session without 1H")
    csv_lines = inventory.to_csv(inv, out / "inventory.csv").read_text(encoding="utf-8").splitlines()
    _check(len(csv_lines) == 1 + len(rows), f"inventory.to_csv: {len(csv_lines)} lines")
    say(f"inventory.build({root.name}): {len(rows)} EXPNOs "
        f"{ {n: (r.nucleus, f'{r.info.ndim}D', r.role) for n, r in sorted(rows.items())} }, picks {picks}, "
        f"nuclei {inv.nuclei()}, SR verdict '{rows[3616].sr_verdict}' (no 1H reference in the session)")


def _racq_referencing(say, E, T, root, out, blocks):
    """larmor.referencing: the session scan, the Xi ratios, the audit -- on the
    session as it is (no 1H) and on a synthetic adamantane-referenced copy."""
    from larmor import acquisition as A, inventory, referencing as R

    xi = {"1H": 1.0, "11B": 0.32083974, "27Al": 0.26056859}   # IUPAC 2001 (Harris et al.) Xi / 100
    acqs = R.scan_session(root)
    by = {a.expno: a for a in acqs}
    _check([a.expno for a in acqs] == sorted(E, key=lambda n: T[n]["date"]),
           f"scan_session: {[a.expno for a in acqs]}, the files' date order {sorted(E, key=lambda n: T[n]['date'])}")
    bad = {n: (a.nucleus, a.sr_hz, a.ns, a.d1_s) for n, a in by.items()
           if a.nucleus != T[n]["nucleus"] or abs(a.sr_hz - T[n]["sr"]) > 1e-6 or a.ns != T[n]["ns"]
           or a.d1_s != T[n]["d"][1] or abs(a.magnet_1h_MHz - T[n]["bf2"]) > 1e-6 * T[n]["bf2"]}
    _check(not bad, f"scan_session acquisitions disagree with the files: {bad}")
    got_xi = {nuc: R.xi_ratio(nuc) for nuc in sorted({T[n]["nucleus"] for n in E} | {"1H"})}
    # the IUPAC table gives Xi to 8 decimals: agreement to half a unit of the last one
    _check(all(abs(v - xi[k]) <= 5e-9 for k, v in got_xi.items()), f"xi_ratio {got_xi} vs IUPAC {xi}")
    verdicts = {a.acq.expno: a.verdict for a in R.audit(acqs)}
    _check(set(verdicts.values()) == {"no reference"}, f"audit without a 1H reference: {verdicts}")
    say(f"scan_session: {len(acqs)} EXPNOs in date order, SR {', '.join(f'{n} {by[n].sr_hz:.2f}' for n in sorted(by))} "
        f"Hz = (SF - BF1)·1e6; xi_ratio {got_xi} = IUPAC; audit: all 'no reference' (no 1H in the session)")
    # a session with a 1H adamantane reference whose IUPAC-Xi prediction is 3616's stored SF
    t = T[3616]
    rs = out / "refsession"
    sf_h = t["sf"] / xi[t["nucleus"]]
    ref = rs / "adamantane" / "1"
    _racq_write_jcamp(ref / "acqus", {"NUC1": "<1H>", "BF1": t["bf2"], "SFO1": t["bf2"], "O1": 0,
                                      "DATE": t["date"] - 3600, "PROBHD": f"<{t['probe']}>",
                                      "PULPROG": "<zg>", "NS": 8, "TD": 1024, "SW_h": 10000})
    _racq_write_jcamp(ref / "pdata" / "1" / "procs", {"SF": sf_h, "SI": 1024})
    (ref / "pdata" / "1" / "title").write_text("1H adamantane reference", encoding="utf-8")
    ok = _racq_copy(E[3616], rs / "pCABS2-4" / "3616")
    off = _racq_copy(E[3616], rs / "pCABS2-4" / "3617", sf=t["bf1"])      # xiref forgotten: SF = BF1
    rows = {r.acq.expno: r for r in R.audit(R.scan_session(rs))}
    want_sr = (sf_h * xi[t["nucleus"]] - t["bf1"]) * 1e6                  # = the file's own SR
    want_d = (t["bf1"] - t["sf"]) / t["sf"] * 1e6
    r_ok, r_off = rows[3616], rows[3617]
    # 0.02 ppm (2.6 Hz at 130 MHz): the 8-decimal precision of the tabulated Xi, below the 0.05 ppm verdict
    tol_hz = 0.02 * t["sf"]
    _check(rows[1].verdict == "1H reference" and r_ok.verdict == "ok" and abs(r_ok.delta_ppm) < 0.02
           and r_ok.ref is not None and r_ok.ref.expno == 1 and r_off.verdict == "unreferenced"
           and abs(r_off.expected_sr_hz - want_sr) < tol_hz and abs(r_off.delta_ppm - want_d) < 0.02,
           f"audit: {[(n, r.verdict, r.delta_ppm, r.expected_sr_hz) for n, r in sorted(rows.items())]}; "
           f"expected SR {want_sr:.2f} Hz, offset {want_d:.3f} ppm")
    ev = A.referencing_evidence([A.read_block(ok)])
    none_off, none_here = A.referencing_evidence([A.read_block(off)]), A.referencing_evidence([blocks[3616]])
    ok_txt = " ".join(A.sentences(A.read_block(ok), referencing_text=ev))
    _check(bool(ev) and "adamantane" in ev and f"{R.ADAMANTANE_1H_PPM:g} ppm" in ev and none_off is None
           and none_here is None and f"Chemical shifts were {ev} (SR = " in ok_txt
           and "[state the reference standard]" not in ok_txt,
           f"referencing_evidence: ok {ev!r}, unreferenced {none_off!r}, session without 1H {none_here!r}")
    inv = inventory.build(rs)
    cand = inventory.sr_candidates(inv)
    fix = next((f for f in inventory.fixes_for(cand[0], inv) if f.kind == "sr_override"), None) if cand else None
    _check([r.expno for r in cand] == [3617] and fix is not None
           and abs(fix.payload["new_sr_hz"] - want_sr) < tol_hz and abs(fix.payload["old_sr_hz"]) < 0.5
           and Path(fix.payload["reference"]) == ref,
           f"inventory on the referenced session: candidates {[r.expno for r in cand]}, fix "
           f"{None if fix is None else fix.payload}")
    say(f"audit with a 1H adamantane reference (SF {sf_h:.6f} MHz): 3616 copy 'ok' ({r_ok.delta_ppm:+.5f} ppm), "
        f"SF = BF1 copy 'unreferenced', expected SR {r_off.expected_sr_hz:.2f} Hz (file {t['sr']:.2f}), offset "
        f"{r_off.delta_ppm:+.3f} ppm (BF1 vs SF: {want_d:+.3f}); evidence clause only for the 'ok' one; "
        f"inventory remedy: sr {fix.payload['new_sr_hz']:.2f} on {cand[0].expno}")


def _racq_provenance(say, ctx, E, T, out, ex):
    """larmor.provenance: the source hash, the shipped recipes' hash, the
    reload checks on a modified copy, the software stamp, carry_source."""
    import datetime as dt
    import hashlib
    import json
    import re

    import larmor
    from larmor import acquisition as A, provenance as P
    from larmor.loader import load_any

    def sha(p):
        return hashlib.sha256(Path(p).read_bytes()).hexdigest()

    f1r = E[3616] / "pdata" / "1" / "1r"
    rec = {"source_path": str(E[3616]), "source_kind": "bruker"}
    h = P.source_sha256(f1r)
    _check(h == sha(f1r) == P.source_sha256(f1r) == P.source_sha256(rec) and P.data_file_for(rec) == f1r
           and P.data_file_for(rec, from_raw=True) == E[3616] / "fid",
           f"source_sha256: {h} vs hashlib {sha(f1r)}; data_file_for {P.data_file_for(rec)}")
    line = f"source_sha256(3616 1r) {h[:16]}… = hashlib.sha256 of the bytes, twice"
    if ex is not None:
        for name, n in (("pCABS2-4_27Al.recipe.json", 3616), ("pCABS2-4_11B.recipe.json", 1118)):
            d = json.loads((ex / name).read_text(encoding="utf-8"))
            d["source_path"] = str(E[n])                # the shipped path is relative to the repository
            one = E[n] / "pdata" / "1" / "1r"
            _check(d["source_sha256"] == sha(one) and P.data_file_for(d) == one,
                   f"{name}: source_sha256 {d['source_sha256'][:16]} is not sha256({n}/pdata/1/1r) {sha(one)[:16]}")
            # reopened against its own data: neither changed bytes nor a moved
            # axis (the stored SR is the data's own, as the fit saw it)
            lines = P.verify_source(d, load_any(str(E[n]))[2])
            _check(lines == [] and abs(float(d.get("sr_hz") or 0.0) - T[n]["sr"]) < P.SR_TOL_HZ,
                   f"examples/{name} (sr_hz {d.get('sr_hz')}, the data's SR {T[n]['sr']:.2f} Hz) reopens with {lines}")
        line += "; the shipped 27Al / 11B recipes hash 3616 / 1118 pdata/1/1r and reopen without a warning"
    else:
        line += "; the shipped-recipe hashes need the bundled examples"
    # the reload checks on a copy: unchanged -> silent, one byte -> hash line, SF moved -> SR line
    vc = _racq_copy(E[3616], out / "verify" / "pCABS2-4" / "3616")
    one = vc / "pdata" / "1" / "1r"
    stored = load_any(str(vc))[2]
    _check(stored["source_kind"] == "bruker" and stored["source_sha256"] == sha(one)
           and stored["acquisition"] == A.read_block(vc) and abs(stored["sr_hz"] - T[3616]["sr"]) < 1e-6,
           f"load_any(copy): sha {stored['source_sha256'][:16]}, SR {stored['sr_hz']}, block equal "
           f"{stored['acquisition'] == A.read_block(vc)}")
    quiet = P.verify_source(stored, load_any(str(vc))[2])
    raw = one.read_bytes()
    one.write_bytes(bytes([raw[0] ^ 1]) + raw[1:])
    changed = P.verify_source(stored, load_any(str(vc))[2])
    one.write_bytes(raw)
    procs = vc / "pdata" / "1" / "procs"
    moved = T[3616]["sf"] + 100e-6                      # a TopSpin 'sr' 100 Hz later
    procs.write_text(re.sub(r"(?m)^##\$SF= .*$", f"##$SF= {moved!r}", procs.read_text(encoding="latin-1")),
                     encoding="latin-1")
    rerefd = P.verify_source(stored, load_any(str(vc))[2])
    _check(quiet == [] and len(changed) == 1 and changed[0].startswith(P.SOURCE_CHANGED_PREFIX)
           and len(rerefd) == 1 and rerefd[0].startswith(P.REREFERENCED_PREFIX)
           and f"{T[3616]['sr']:.2f}" in rerefd[0] and f"{T[3616]['sr'] + 100:.2f}" in rerefd[0],
           f"verify_source: unchanged {quiet}, one byte {changed}, SF +100 Hz {rerefd}")
    say(line + f"; reload checks on a copy: unchanged {quiet}, one byte -> {changed[0][:28]}…, "
        f"SF +100 Hz -> '{rerefd[0][:48]}…'")
    # the software stamp: the running versions, read back independently
    import mrsimulator

    s, s2 = P.software_stamp(), P.software_stamp()
    when = dt.datetime.fromisoformat(s["fitted"])
    want = {"larmor": larmor.__version__, "numpy": np.__version__, "mrsimulator": mrsimulator.__version__}
    frozen = bool(ctx.get("frozen"))
    off = {k: (s.get(k), v) for k, v in want.items() if s.get(k) != v and not (frozen and k != "larmor" and s.get(k) == "")}
    long = P.git_commit()
    short = s.get("git_commit") or ""
    _check(not off and abs((dt.datetime.now() - when).total_seconds()) < 300
           and {k: v for k, v in s.items() if k != "fitted"} == {k: v for k, v in s2.items() if k != "fitted"}
           and (not short or not long or long.startswith(short)),
           f"software_stamp: {off}, fitted {s['fitted']}, git {short!r} vs file-read {long!r}")
    # carry_source: the per-spectrum copy keeps kind and block; from the raw fid it hashes the fid
    fid = E[3616] / "fid"
    cs = P.carry_source(stored, raw_expno=E[3616])
    _check(cs["source_kind"] == "bruker" and cs["source_sha256"] == sha(fid)
           and cs["acquisition"] == {**stored["acquisition"], "data_file": "fid"}
           and P.carry_source({}) == {"source_kind": "", "source_sha256": "", "acquisition": {}},
           f"carry_source: {cs['source_kind']} {cs['source_sha256'][:16]} vs fid {sha(fid)[:16]}")
    say(f"software_stamp: LARMOR {s['larmor']}, numpy {s['numpy']}, mrsimulator {s['mrsimulator']} (= the "
        f"imported modules), git {short or '-'} (file read {long[:12] or '-'}); carry_source(raw fid) hashes "
        f"the fid {cs['source_sha256'][:12]}… and marks data_file 'fid'")


def _racq_pulseprog(say, E, T):
    """larmor.pulseprog: parse, resolve against acqus / acqu2s, timeline, describe."""
    import re

    from larmor import pulseprog as pp

    for n, e in E.items():
        t = T[n]
        rec = pp.load(e)
        prog, vals = rec.program, rec.values
        _check(rec.path is not None and rec.path.name == "pulseprogram" and prog is not None
               and not rec.warnings and rec.pulprog == t["pulprog"] and prog.name == t["pulprog"]
               and rec.nuclei.get("f1") == t["nucleus"] and rec.nuclei.get("f2") == t["nuc2"],
               f"pulseprog.load({n}): {rec.path}, {rec.pulprog}, {rec.nuclei}, {rec.warnings}")
        stmts = [s for s in prog.walk() if s.kind in ("pulse", "shaped_pulse")]
        _check(len(stmts) == t["n_pulses"] and (stmts[0].name, stmts[0].power or "pl1") == ("p1", f"pl{t['p1_power']}"),
               f"pulseprog {n}: pulses {[(s.name, s.power) for s in stmts]}; the program text has "
               f"{t['n_pulses']} with p1 at pl{t['p1_power']}")
        bad = {}
        for s in prog.walk():
            if s.kind in ("pulse", "shaped_pulse", "delay") and re.fullmatch(r"[pd]\d+", s.name):
                v = vals.get(s.name)
                want = t["p"][int(s.name[1:])] * 1e-6 if s.name[0] == "p" else t["d"][int(s.name[1:])]
                # a duration the program defines ("p2=1s/(cnst31*cnst0)") is recomputed; acqus keeps
                # TopSpin's own value to 7 digits
                tol = (1e-6 if s.name in prog.definitions else 1e-12) * max(want, 1e-12)
                if v is None or v.si is None or abs(v.si - want) > tol:
                    bad[s.name] = (None if v is None else v.si, want)
            if s.kind in ("pulse", "shaped_pulse") and s.power:
                m = re.fullmatch(r"(pl|sp)(\d+)", s.power)
                want = (t["plw"] if m.group(1) == "pl" else t["spw"])[int(m.group(2))]
                if vals[s.power].si != want:
                    bad[s.power] = (vals[s.power].si, want)
        if vals["ns"].si != t["ns"] or abs(vals["aq"].si - t["aq"]) > 1e-15:
            bad["ns/aq"] = (vals["ns"].si, vals["aq"].si, t["ns"], t["aq"])
        if t["twod"] and (vals["td1"].si != t["td1"] or ("ST1CNT" in vals and vals["ST1CNT"].si != t["td1"] // 2)):
            bad["td1"] = (vals["td1"].si, t["td1"])
        _check(not bad, f"pulseprog.resolve({n}) disagrees with acqus / acqu2s (resolved, file): {bad}")
        tail = t["prog"][re.search(r"^exit\s*$", t["prog"], re.M).end():]
        ph31 = re.search(r"^ph31\s*=([^\n]*(?:\n[ \t]+[\d \t]+)*)", tail, re.M)
        n31 = len(re.findall(r"\d+", ph31.group(1)))
        tl = pp.timeline(prog, vals)
        ev = [x for x in tl.events if x.kind in ("pulse", "shaped_pulse")]
        go = next(x for x in tl.events if x.kind == "acquire")
        d1 = next(x for x in tl.events if x.kind == "delay" and x.label == "d1")
        ns_loop = next(lp for lp in tl.loops if lp.times == "ns")
        _check(len(prog.phase_programs.get("ph31", [])) == n31 and len(ev) == len(stmts)
               and [x.value for x in ev] == [vals[s.name].text for s in stmts]
               and all(x.channel == (s.channel or "f1") for x, s in zip(ev, stmts))
               and all(a.end <= b.start for a, b in zip(ev, ev[1:])) and d1.end <= ev[0].start
               and ev[-1].end <= go.start and go.end <= tl.total and tl.channels[-1] == "acq"
               and "f1" in tl.channels and not tl.notes and d1.value == vals["d1"].text
               and go.value == f"aq {vals['aq'].text}" and ns_loop.text == f"× ns = {t['ns']}",
               f"timeline({n}): {[(x.kind, x.label, x.value, x.start) for x in tl.events]}, loops "
               f"{[(lp.times, lp.text) for lp in tl.loops]}, notes {tl.notes}, ph31 {n31} steps")
        desc = pp.describe(prog, vals)
        _check(desc.startswith(f"{t['pulprog']}: d1 – p1") and desc.endswith(f"phase cycle of {prog.phase_cycle_length()} steps"),
               f"describe({n}): {desc}")
        pw = stmts[0].power
        say(f"pulseprog {n} {t['pulprog']}: {len(stmts)} pulse(s) {[s.name + (':' + s.power if s.power else '') for s in stmts]} "
            f"on {tl.channels}; p1 {vals['p1'].text} (= acqus P[1]) at "
            + (f"{pw} {vals[pw].text} (= PLW[{pw[2:]}])" if pw else "the default pl1")
            + f"; d1 {vals['d1'].text}, aq {vals['aq'].text} (= TD/(2 SW_h)), ns {t['ns']}"
            + (f", td1 {t['td1']} (acqu2s)" if t["twod"] else "")
            + f"; ph31 {n31} steps as written; timeline {len(tl.events)} events / {len(tl.loops)} loops")


def _racq_quantitativity(say, E, T):
    """larmor.quantitativity: the D1 / flip-angle judgement recomputed from acqus."""
    import math
    import re

    from larmor import quantitativity as Q

    spin = {"1H": 0.5, "11B": 1.5, "27Al": 2.5}           # nuclear spins
    # the titles say 'pi/12 a 20 kHz': a 20 kHz rf field has a 90-degree pulse of 1 / (4 x 20 kHz)
    p90 = 1e6 / (4 * 20e3)
    t1_a, t1_b = 0.5, 2.0                                 # typed T1: site A per site, the rest overall
    sites = [{"model": "gauss_lor", "label": lab, "params": {"isotropic_chemical_shift_ppm": {"value": v}}}
             for lab, v in (("A", 60.0), ("B", 5.0))]
    t1_kinds = re.compile(r"satrec|t1ir|invrec", re.I)
    parts = []
    for n, e in E.items():
        t = T[n]
        i_ = spin[t["nucleus"]]
        f = Q.facts_for(e)
        a = f.acquisition
        _check(a is not None and (a.nucleus, a.d1_s, a.p1_us, a.plw1_w, a.ns) == (t["nucleus"], t["d"][1], t["p"][1], t["plw"][1], t["ns"])
               and abs(a.aq_s - t["aq"]) < 1e-12 and f.spin == i_,
               f"facts_for({n}): {a}, spin {f.spin}")
        sibling = any(t1_kinds.search(T[m]["pulprog"]) and T[m]["nucleus"] == t["nucleus"] for m in E if m != n)
        _check(sibling or f.t1_status == "none", f"facts_for({n}): T1 status {f.t1_status} with no T1 EXPNO")
        single = t["n_pulses"] == 1
        limit = 30.0 / (i_ + 0.5)                         # Eden Eq. 38 for a half-integer quadrupole
        chk = Q.check(f, sites)
        judged = (chk.flip_assumed_90 and chk.flip_deg is None and chk.excitation_judged
                  and abs(chk.flip_limit_deg - limit) < 1e-12) if single else \
            (chk.flip_source == "echo(90)" and chk.flip_deg == 90.0 and not chk.excitation_judged)
        _check(judged and chk.recoveries == [] and chk.recovery_unknown_text() == f"recycle {t['d'][1]:g} s — T1 unknown",
               f"check({n}) with nothing typed: flip {chk.flip_deg} ({chk.flip_source}), judged "
               f"{chk.excitation_judged}, limit {chk.flip_limit_deg}, recoveries {chk.recoveries}")
        chk2 = Q.check(f, sites, {"p90_us": p90, "t1_s": t1_b, "t1_by_site": {"A": t1_a}})
        flip = 90.0 * t["p"][1] / p90 if single else 90.0
        theta = min((i_ + 0.5) * flip, 180.0) if single else 90.0     # CT nutation (I + 1/2) theta
        recycle = t["d"][1] + t["aq"]

        def ss(t1):
            x = math.exp(-recycle / t1)
            return (1.0 - x) / (1.0 - x * math.cos(math.radians(theta)))

        got = {r.label: r for r in chk2.recoveries}
        _check(abs(chk2.flip_deg - flip) < 1e-9 and set(got) == {"A", "B"}
               and (got["A"].source, got["B"].source) == ("per_site", "user")
               and abs(got["A"].recovery - ss(t1_a)) < 1e-9 and abs(got["B"].recovery - ss(t1_b)) < 1e-9
               and abs(got["A"].ratio - recycle / t1_a) < 1e-9
               and chk2.excitation_over() == (single and flip > limit)
               and chk2.recovery_ok() == (min(ss(t1_a), ss(t1_b)) >= Q.RECOVERY_MIN),
               f"check({n}) typed P1(90) {p90} µs, T1 {t1_a} / {t1_b} s: flip {chk2.flip_deg} vs {flip}, "
               f"recoveries {[(r.label, r.recovery, r.source) for r in chk2.recoveries]} vs "
               f"{ss(t1_a):.6f} / {ss(t1_b):.6f}")
        parts.append(f"{n} {'flip %.1f° vs limit %g°' % (flip, limit) if single else 'MQMAS 90° (not judged)'}, "
                     f"recycle {recycle:.4g} s -> {100 * ss(t1_a):.1f} / {100 * ss(t1_b):.1f} %")
    say("quantitativity: acqus D1 + TD/(2 SW_h), P1, PLW1 and the spin read right; nothing typed -> flip "
        "unknown (90° assumed), T1 unknown; typed P1(90) 12.5 µs (the titles' 'pi/12 a 20 kHz') and T1 "
        "0.5 s (A) / 2 s: " + "; ".join(parts) + " = (1-E)/(1-E cos (I+½)θ) recomputed")


def _records_acquisition(say, ctx):
    """The acquisition records a user publishes, each checked against the raw
    files read independently (a regex over acqus / procs / title / acqu2s /
    the pulse program): larmor.acquisition (the flat record of every EXPNO,
    the Experimental sentences, Table S1 as CSV / LaTeX / Markdown, the
    per-nucleus paragraph), larmor.comparability (an 11B / 27Al pair flagged
    on exactly the keys that differ, an EXPNO against itself alike, procs ->
    LARMOR ops), larmor.inventory (the session's EXPNOs, nuclei, 1D / 2D and
    picks), larmor.referencing (the session scan, the IUPAC Xi ratios, the
    audit verdicts on a synthetic adamantane-referenced session, the
    evidence clause, the inventory's SR remedy), larmor.provenance (the
    source hash, the shipped recipes' hashes, the reload checks on a modified
    copy, the software stamp, carry_source), larmor.pulseprog (zg and
    mp3qdfsz parsed, every pulse / delay / power resolved against acqus, the
    timeline) and larmor.quantitativity (the D1 / flip-angle judgement
    recomputed from acqus). The bundled EXPNOs 1118 / 3616 / 3620 when
    present, else synthetic EXPNOs with the same values under the stage
    folder."""
    out = _out(ctx, "records_acq")
    _import_mrsimulator(say)           # BF1 / Xi and the audit read mrsimulator's isotope table
    t0 = time.perf_counter()
    E = {n: _expno(ctx, n) for n in (1118, 3616, 3620)}
    ex = _examples(ctx)
    if ex is None or any(v is None for v in E.values()):
        folder = _racq_synthetic(out / "synthetic")
        E, ex = {n: folder / str(n) for n in E}, None
        say("no bundled example data: EXPNOs 1118 / 3616 / 3620 written as synthetic acqus / procs / title / "
            "pulseprogram / 1r files with the bundled values; the shipped-recipe hash check needs the examples")
    root = E[3616].parent.parent       # the session folder the inventory / audit scan
    T = {n: _racq_truth(e) for n, e in E.items()}
    blocks = _racq_acquisition(say, E, T)
    _racq_comparability(say, E, T, out)
    _racq_inventory(say, E, T, root, out)
    _racq_referencing(say, E, T, root, out, blocks)
    _racq_provenance(say, ctx, E, T, out, ex)
    _racq_pulseprog(say, E, T)
    _racq_quantitativity(say, E, T)
    say(f"records_acq: {time.perf_counter() - t0:.1f} s after the mrsimulator import")


def _records_tables(say, ctx):
    """The records a fit leaves and the tables a user publishes: per-site
    populations, tagged families and the named ratios against the closed-form
    areas of Gaussian lines and the first-order error formulas (independent,
    covariance and Monte-Carlo bases); the software block, the Methods text
    and the LaTeX table; every shipped recipe loaded, saved and reloaded, the
    schema-migration hook, forward compatibility and a one-value diff; the
    derived parameter statuses of the shipped fits; the project bundle
    reopened after the project moved with its data; the multi-fit parameter
    table (CSV / TSV re-read against the recipes); the series table (CSV
    join, replicate statistics, OLS, the series and DUST CSVs, species bars);
    the batch-grid resolver; the alias store and the rename log (both
    redirected under the stage folder)."""
    import shutil

    out = _out(ctx, "records_tables")
    for old in list(out.iterdir()):              # a kept --out folder from an earlier run
        if old.is_dir():
            shutil.rmtree(old)
        else:
            old.unlink()
    ex = _examples(ctx)
    fits = _rtab_fits(say, out, ex)
    tagged, quant = _rtab_quantify(say)
    _rtab_methods(say, fits, tagged, quant)
    _rtab_recipes(say, out, fits, ex)
    _rtab_statuses(say, fits)
    _rtab_fittable(say, out, fits)
    _rtab_project(say, out, fits, ex)
    _rtab_series(say, out)
    _rtab_grid(say, out)
    _rtab_aliases(say, out)


# ------------------------------------------------------------ records: inputs
def _rtab_fits(say, out: Path, ex) -> list:
    """[(path, Recipe)] for every shipped examples/*.recipe.json; without the
    bundled data, two synthetic fits saved under the stage folder (a 27Al
    three-line Czjzek fit with errors, an 11B Amorphous + Gauss/Lorentz fit
    with held gl) so every check that follows runs on saved files either way."""
    from larmor.recipe import Recipe

    if ex is not None:
        paths = sorted(ex.glob("*.recipe.json"))
        _check(len(paths) >= 2, f"examples/ holds {len(paths)} recipe file(s); the 27Al and 11B fits are expected")
        return [(p, Recipe.load(p)) for p in paths]
    say("no bundled example data: the recipe / table checks run on two synthetic saved fits")
    fits = []
    for rec in _rtab_synthetic_fits():
        p = out / f"{rec.sample}.recipe.json"
        rec.save(p)
        fits.append((p, Recipe.load(p)))
    return fits


def _rtab_synthetic_fits() -> list:
    """Stand-ins for the shipped fits: the same models, bounds, held gl and
    stderr pattern."""
    from larmor import models
    from larmor.recipe import Param, Recipe, SiteModel

    def cz(label, pos, sigma, dcs, amp, errs):
        return SiteModel(model="czjzek", label=label, params={
            ISO: Param(pos, stderr=errs[0]), "sigma_Cq_MHz": Param(sigma, stderr=errs[1], min=0.05),
            "shift_fwhm_ppm": Param(dcs, stderr=errs[2], min=0.1), "amplitude": Param(amp, stderr=errs[3], min=0.0)})

    # full-precision values, as a fit writes them, so the tables' digits are tested
    al = Recipe(sample="synthetic_27Al", source_kind="csv", source_path="synthetic_27Al.csv", nucleus="27Al",
                larmor_frequency_MHz=130.3175616, spin_rate_Hz=26000.0, fit_window_ppm=(150.0, -80.0),
                fit_rmsd=0.045798748, sites=[
                    cz("AlO4", 62.815324, 1.4698706, 14.014254, 3751815.8, (1.0703761, 0.33867882, 0.94338737, 761667.17)),
                    cz("AlO5", 28.591567, 0.64155377, 33.253842, 1264940.6, (5.2134172, 0.25507268, 7.4162388, 697700.76)),
                    cz("AlO6", -0.18255056, 0.80884753, 17.688993, 573644.39, (6.2088373, 0.68935698, 8.2007048, 559775.2))])
    amo = models.get("amorphous").defaults()
    # (eta_fwhm off its registry default 0 = its lower bound: a saved fit with a
    # value AT a bound carries the fit's at-bound note, which this one does not)
    for key, value in ((ISO, 16.369827), ("Cq_MHz", 2.0991939), ("eta", 0.073527829), ("eta_fwhm", 0.1605395),
                       ("shift_fwhm_ppm", 5.2865868), ("amplitude", 2692771.5)):
        amo[key].value = value
    b = Recipe(sample="synthetic_11B", source_kind="csv", source_path="synthetic_11B.csv", nucleus="11B",
               larmor_frequency_MHz=160.4606499, spin_rate_Hz=20000.0, fit_window_ppm=(30.0, -30.0),
               fit_rmsd=0.0037235576, sites=[SiteModel(model="amorphous", label="BO3", params=amo),
                                             SiteModel(model="gauss_lor", label="BO4", params={
                                                 ISO: Param(-0.30681179), "shift_fwhm_ppm": Param(2.8744391, min=0.1),
                                                 "amplitude": Param(8583843.7, min=0.0), "gl": Param(1.0, vary=False)})])
    return [al, b]


def _rtab_line(pos, fwhm, amp, label, family="", stderr=None):
    """A pure Gaussian (gl = 1) Gauss/Lorentz line: its area is
    amp * fwhm * sqrt(pi / (4 ln 2)) in closed form."""
    from larmor.recipe import Param, SiteModel

    return SiteModel(model="gauss_lor", label=label, family=family, params={
        ISO: Param(pos), "shift_fwhm_ppm": Param(fwhm), "amplitude": Param(amp, stderr=stderr),
        "gl": Param(1.0, vary=False)})


@contextlib.contextmanager
def _rtab_env(**values):
    """The environment variables set to ``values`` for the block, the
    previous state (set or unset) put back afterwards."""
    old = {k: os.environ.get(k) for k in values}
    os.environ.update({k: str(v) for k, v in values.items()})
    try:
        yield
    finally:
        for k, v in old.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


# ------------------------------------------------------------ quantification
def _rtab_quantify(say):
    """quantify + families on four Gaussian 11B lines (two tagged BO3, one
    BO4, one untagged) and three tagged 27Al lines. Returns the tagged 11B
    recipe and its quantify() table for the Methods checks."""
    import math

    from uncertainties import correlated_values

    from larmor import families, quantify
    from larmor.recipe import Recipe

    k_area = math.sqrt(math.pi / (4.0 * math.log(2.0)))      # area of a unit-height Gaussian per ppm of FWHM
    f2s = 1.0 / (2.0 * math.sqrt(2.0 * math.log(2.0)))       # FWHM -> sigma
    # (position ppm, FWHM ppm, amplitude, label, family, amplitude stderr)
    lines = [(17.0, 5.0, 80.0, "B3s", "BO3", 2.0), (12.0, 6.0, 50.0, "B3a", "BO3", 2.5),
             (0.5, 3.0, 120.0, "B4", "BO4", 1.5), (-12.0, 4.0, 10.0, "imp", "", 1.0)]
    rec = Recipe(nucleus="11B", larmor_frequency_MHz=160.46, spin_rate_Hz=20000.0, sample="tagged",
                 source_kind="csv", sites=[_rtab_line(*ln) for ln in lines])
    pos, fw, amp, sa = (np.array([ln[i] for ln in lines], float) for i in (0, 1, 2, 5))
    area = amp * fw * k_area                    # the truth
    sig = area * sa / amp                       # integral errors at a fixed shape
    total = float(area.sum())
    win = (60.0, -40.0)
    q = quantify.quantify(rec, win)
    rows = q["rows"]
    got = np.array([r["integral"] for r in rows])
    # every line is sampled at <= 0.23 sigma on the 2048-point default axis and
    # the window keeps > 17 sigma of each tail: the trapezoid is exact to ~1e-12
    _check(np.allclose(got, area, rtol=1e-6, atol=0.0), f"quantify integrals {got} vs the Gaussian areas {area}")
    frac = np.array([r["fraction_pct"] for r in rows])
    _check(abs(frac.sum() - 100.0) < 1e-9 and np.allclose(frac, 100.0 * area / total, rtol=0.0, atol=1e-6),
           f"populations {frac} (sum {frac.sum()}) vs {100.0 * area / total}")
    # f_i = I_i / T over independent integrals: dI_i -> (T - I_i)/T^2, dI_j -> -I_i/T^2
    want = [100.0 * math.sqrt(((total - area[i]) * sig[i]) ** 2
                              + sum((area[i] * sig[j]) ** 2 for j in range(len(lines)) if j != i)) / total ** 2
            for i in range(len(lines))]
    errs = [r["fraction_err_pct"] for r in rows]
    _check(all(e is not None for e in errs) and np.allclose(errs, want, rtol=1e-9, atol=0.0),
           f"population errors {errs} vs first-order {want}")
    one = quantify.quantify(Recipe(nucleus="11B", larmor_frequency_MHz=160.46, sites=[_rtab_line(*lines[0])]),
                            win)["rows"][0]
    _check(one["fraction_pct"] == 100.0 and one["fraction_err_pct"] == 0.0,
           f"a single line: {one['fraction_pct']} ± {one['fraction_err_pct']} % (100 ± 0 by definition)")
    say(f"quantify: 4 Gaussian 11B lines, integrals = closed-form areas (max rel. dev. "
        f"{np.max(np.abs(got / area - 1.0)):.1e}), populations {' / '.join(f'{v:.2f}' for v in frac)} % "
        f"(sum {frac.sum():.9f}), first-order errors {' / '.join(f'{e:.3f}' for e in errs)}; one line 100 ± 0 %")

    # tagged families and N4 = BO4 / (BO3 + BO4): the untagged line counts in
    # the total of the family percentages but never in a ratio
    fams = {f["family"]: f for f in q["families"]}
    f3, f4 = float(area[0] + area[1]), float(area[2])
    _check(list(fams) == ["BO3", "BO4"] and fams["BO3"]["letters"] == "A+B" and fams["BO4"]["sites"] == [2]
           and q["untagged"] == [3] and abs(fams["BO3"]["fraction_pct"] - 100.0 * f3 / total) < 1e-6
           and abs(fams["BO4"]["fraction_pct"] - 100.0 * f4 / total) < 1e-6
           and abs(fams["BO3"]["fraction_pct"] + fams["BO4"]["fraction_pct"] + frac[3] - 100.0) < 1e-9,
           f"families {[(f['family'], f['letters'], f['fraction_pct']) for f in q['families']]}, untagged "
           f"{q['untagged']} vs BO3 {100 * f3 / total:.4f} %, BO4 {100 * f4 / total:.4f} %")
    n4 = next(r for r in q["ratios"] if r["name"] == "N4")
    n4_true = f4 / (f3 + f4)
    # delta method over independent family integrals (the 'independent' basis)
    n4_err = math.sqrt((f4 * math.hypot(sig[0], sig[1])) ** 2 + (f3 * sig[2]) ** 2) / (f3 + f4) ** 2
    _check(n4["defined"] and abs(n4["value"] - n4_true) < 1e-9 and n4["err"] is not None
           and abs(n4["err"] / n4_err - 1.0) < 1e-6 and q["family_basis"] == "independent",
           f"N4 {n4['value']} ± {n4['err']} ({q['family_basis']}) vs {n4_true} ± {n4_err}")
    # the family shares: first-order propagation over independent integrals
    # with the normalisation in the derivative, sqrt(R² σ_F² + F² σ_R²) / T²
    # -- so a one-line family carries exactly its site row's error
    fam_want = {}
    for name, idx in (("BO3", [0, 1]), ("BO4", [2])):
        rest = [i for i in range(len(lines)) if i not in idx]
        f_in, f_out = float(area[idx].sum()), float(area[rest].sum())
        fam_want[name] = 100.0 * math.sqrt(f_out ** 2 * float((sig[idx] ** 2).sum())
                                           + f_in ** 2 * float((sig[rest] ** 2).sum())) / total ** 2
    _check(all(fams[k]["fraction_err_pct"] is not None and abs(fams[k]["fraction_err_pct"] / v - 1.0) < 1e-9
               for k, v in fam_want.items()) and abs(fams["BO4"]["fraction_err_pct"] / errs[2] - 1.0) < 1e-9,
           f"family errors {[(k, fams[k]['fraction_err_pct']) for k in fam_want]} vs first-order {fam_want}; "
           f"the one-line BO4 vs its own row {errs[2]}")
    # covariance basis: lmfit-style ufloats with the overlapping BO3 pair
    # anticorrelated, against the analytic Jacobian over the amplitudes
    corr = np.eye(len(lines))
    corr[0, 1] = corr[1, 0] = -0.6
    corr[1, 2] = corr[2, 1] = 0.3
    cov = corr * np.outer(sa, sa)
    qc = quantify.quantify(rec, win, uvars=dict(enumerate(correlated_values(amp.tolist(), cov))))
    kk = fw * k_area                            # integral per unit amplitude
    g_n4 = np.array([-kk[0] * f4, -kk[1] * f4, kk[2] * f3, 0.0]) / (f3 + f4) ** 2
    n4c = next(r for r in qc["ratios"] if r["name"] == "N4")
    bad = [] if (qc["family_basis"] == "covariance" and qc["rows"] == rows) else ["basis / rows"]
    if abs(n4c["err"] / math.sqrt(g_n4 @ cov @ g_n4) - 1.0) > 1e-6:
        bad.append(f"N4 ± {n4c['err']} vs {math.sqrt(g_n4 @ cov @ g_n4)}")
    fam_c = {f["family"]: f for f in qc["families"]}
    cov_errs = {}
    for name, member in (("BO3", [1.0, 1.0, 0.0, 0.0]), ("BO4", [0.0, 0.0, 1.0, 0.0])):
        m = np.array(member)
        g = 100.0 * kk * (m * total - float(m @ area)) / total ** 2
        cov_errs[name] = math.sqrt(g @ cov @ g)
        if abs(fam_c[name]["fraction_err_pct"] / cov_errs[name] - 1.0) > 1e-6:
            bad.append(f"{name} ± {fam_c[name]['fraction_err_pct']} vs {cov_errs[name]}")
    _check(not bad, f"covariance basis: {bad}")
    # Monte-Carlo basis: the error is the spread of the per-trial ratio, the
    # value stays the best fit
    trials = area + np.random.default_rng(11).normal(0.0, 1.0, (400, len(lines))) * sig
    mc = families.summarize(area, amp, [ln[4] for ln in lines], "11B", sigma=sig, samples=trials)
    n4m = next(r for r in mc["ratios"] if r["name"] == "N4")
    per = np.abs(trials[:, 2]) / np.abs(trials[:, :3]).sum(axis=1)
    _check(mc["basis"] == "montecarlo" and abs(n4m["value"] - n4_true) < 1e-12
           and abs(n4m["err"] - float(np.std(per))) < 1e-12,
           f"Monte-Carlo basis: N4 {n4m['value']} ± {n4m['err']} vs best {n4_true}, trial std {np.std(per)}")
    # a weighted ratio: <CN> Al = (4 Al(IV) + 5 Al(V) + 6 Al(VI)) / sum
    al_lines = [(60.0, 12.0, 100.0, "AlIV", "Al(IV)"), (32.0, 14.0, 30.0, "AlV", "Al(V)"),
                (5.0, 10.0, 40.0, "AlVI", "Al(VI)")]
    qa = quantify.quantify(Recipe(nucleus="27Al", larmor_frequency_MHz=130.3,
                                  sites=[_rtab_line(*ln) for ln in al_lines]), (150.0, -80.0))
    ia = np.array([a * w * k_area for _p, w, a, _l, _f in al_lines])
    cn_true = float(np.dot([4.0, 5.0, 6.0], ia) / ia.sum())
    cn_name = families.ratios_for("27Al")[0].name
    cn = next((r for r in qa["ratios"] if r["name"] == cn_name), {})
    _check(cn.get("defined") and abs(cn["value"] - cn_true) < 1e-9, f"<CN> Al {cn.get('value')} vs {cn_true}")
    say(f"families: BO3 (A+B) {fams['BO3']['fraction_pct']:.2f} ± {fams['BO3']['fraction_err_pct']:.3f} %, "
        f"BO4 (C) {fams['BO4']['fraction_pct']:.2f} ± {fams['BO4']['fraction_err_pct']:.3f} % (= first-order, BO4 = its "
        f"row), D untagged (in the total only); N4 {n4['value']:.4f} ± {n4['err']:.4f} = truth / delta method; covariance "
        f"basis (rho_AB -0.6) N4 ± {n4c['err']:.4f}, BO3 ± {cov_errs['BO3']:.3f} %, BO4 ± {cov_errs['BO4']:.3f} % = "
        f"analytic Jacobian; Monte-Carlo N4 ± {n4m['err']:.4f} = std of 400 trials; <CN> Al {cn['value']:.3f}")

    # the fit-health tail chip and the tail-covering window against erfc /
    # the normal quantile (the in-window INTEGRAL of a line the window cuts is
    # not checked against erfc: quantify._integrate keeps only the grid points
    # inside the window, so up to one grid step is lost at each cut edge)
    qt = quantify.quantify(rec, (20.0, -40.0))
    tail = qt["rows"][0]["tail_outside_pct"]
    tail_true = 50.0 * math.erfc((20.0 - pos[0]) / (fw[0] * f2s) / math.sqrt(2.0))
    # the share outside is read off the cumulative area interpolated linearly
    # between the 0.29 ppm grid points: 0.25 percentage point covers that
    _check(tail is not None and abs(tail - tail_true) < 0.25 and all(r["tail_outside_pct"] < 1e-6 for r in rows),
           f"tail_outside_pct {tail} vs erfc {tail_true}; full window {[r['tail_outside_pct'] for r in rows]}")
    hi, lo, clipped = quantify.window_containing_tails(rec)
    z = 2.807033768343804                      # one-sided normal quantile of 0.25 % (99.5 % coverage)
    hi_t, lo_t = float(np.max(pos + z * fw * f2s)), float(np.min(pos - z * fw * f2s))
    _check(abs(hi - hi_t) < 0.1 and abs(lo - lo_t) < 0.1 and clipped == {},
           f"window_containing_tails {hi:.3f}..{lo:.3f} vs Gaussian quantiles {hi_t:.3f}..{lo_t:.3f} ({clipped})")
    ints, w2 = quantify.site_integrals(rec, (-40.0, 60.0))
    _check(w2 == (60.0, -40.0) and np.array_equal(ints, got), f"site_integrals {w2} {ints} vs quantify {got}")
    say(f"tails: line A cut at 20 ppm leaves {tail:.2f} % outside (erfc {tail_true:.2f} %); 99.5 % window "
        f"{hi:.2f}..{lo:.2f} ppm (normal quantiles {hi_t:.2f}..{lo_t:.2f}); site_integrals = the table's")
    return rec, q


# ------------------------------------------------------------ Methods text
def _rtab_methods(say, fits, tagged, quant):
    """software_versions, the Methods sentence / paragraph and the LaTeX
    table, against the recipes and the populations they were given."""
    import platform

    import larmor
    from larmor import methods

    sw = methods.software_versions()
    _check(sw.get("larmor") == larmor.__version__ and sw.get("numpy") == np.__version__
           and sw.get("python") == platform.python_version() and sw.get("lmfit") and sw.get("mrsimulator"),
           f"software_versions {sw} (larmor {larmor.__version__}, numpy {np.__version__})")
    words = {"czjzek": "Czjzek", "amorphous": "Amorphous", "gauss_lor": "Gauss/Lorentz"}
    named = []
    for path, rec in fits:
        d = rec.to_dict()
        text = methods.methods_sentence(d, software=sw)
        need = [rec.nucleus, f"{rec.larmor_frequency_MHz:.1f} MHz", f"{len(rec.sites)} sites", f"LARMOR {sw['larmor']}"]
        need += sorted({words[s.model] for s in rec.sites if s.model in words})
        if any(s.model == "czjzek" for s in rec.sites):
            need += ["P_Q", "4σ"]               # the Czjzek width convention is always stated
        missing = [w for w in need if w not in text]
        _check(not missing, f"Methods text of {Path(path).name} lacks {missing}: {text[:240]}")
        named.append("+".join(sorted({words.get(s.model, s.model) for s in rec.sites})))
        # the LaTeX table prints the recipe's own positions, row by row
        lines = methods.latex_table(d).splitlines()
        body = lines[lines.index(r"\midrule") + 1:lines.index(r"\bottomrule")]
        _check(len(body) == len(rec.sites), f"latex_table of {Path(path).name}: {len(body)} rows for {len(rec.sites)} sites")
        for row, site in zip(body, rec.sites):
            p = site.params[ISO]
            want = f"{p.value:.2f}" + (f" ± {p.stderr:.2f}" if p.stderr is not None else "")
            cells = row.removesuffix(r" \\").split(" & ")
            _check(cells[0] == site.label and cells[1] == want,
                   f"latex_table of {Path(path).name}: row {cells[:2]} vs ({site.label!r}, {want!r})")
    # the tagged synthetic fit: families sentence, and the table's population
    # cells are quantify's numbers (verified above against the closed form)
    d = tagged.to_dict()
    text = methods.methods_sentence(d, quant=quant)
    need = ["11B", "4 sites", "Gauss/Lorentz", "structural families (BO3, BO4)", "N4 (BO4/(BO3+BO4)) is reported",
            "treated as independent"]
    missing = [w for w in need if w not in text]
    _check(not missing, f"Methods sentence with families lacks {missing}: {text}")
    lines = methods.latex_table(d, quant).splitlines()
    i0 = lines.index(r"\midrule") + 1
    for row, site, r in zip(lines[i0:i0 + len(tagged.sites)], tagged.sites, quant["rows"]):
        cells = row.removesuffix(r" \\").split(" & ")
        want = f"{r['fraction_pct']:.1f} ± {r['fraction_err_pct']:.1f}"
        _check(cells[0] == site.label and cells[-1] == want, f"latex population cell {cells} vs {want}")
    for f in quant["families"]:
        row = next((ln for ln in lines if ln.startswith(f"Σ {f['family']} & ")), "")
        _check(row.removesuffix(r" \\").split(" & ")[-1].split(" ± ")[0] == f"{f['fraction_pct']:.1f}",
               f"latex family row {row!r} vs {f['fraction_pct']:.1f}")
    n4 = next(r for r in quant["ratios"] if r["name"] == "N4")
    n4_line = f"= {n4['value']:.3f} ± {n4['err']:.3f}"
    _check(any(ln.startswith(r"\multicolumn") and r"N$_{4}$" in ln and n4_line in ln for ln in lines),
           f"latex N4 line ({n4_line}) missing")
    d["software"] = sw
    para = methods.methods_paragraph(d)
    _check(para == methods.methods_sentence(d, software=sw), f"methods_paragraph of a CSV source: {para[:200]}")
    say(f"methods: software larmor {sw['larmor']}, numpy {sw['numpy']}, lmfit {sw['lmfit']}, mrsimulator "
        f"{sw['mrsimulator']}, git {sw.get('git_commit') or '-'}; Methods text names {', '.join(named)} for "
        f"{len(fits)} fits; LaTeX positions = the recipes', populations / Σ families / N4 {n4['value']:.3f} = "
        f"quantify's; paragraph of a CSV source = the sentence + versions")


# ------------------------------------------------------------ recipes
def _rtab_recipes(say, out: Path, fits, ex):
    """Every shipped recipe: lossless load, save / reload, the source data
    where the recipe says; the schema-migration hook (as tests/
    test_recipe_diff drives it), forward compatibility; recipe_diff on one
    changed value."""
    import copy
    import json

    from larmor import models
    from larmor import recipe as R
    from larmor.recipe_diff import recipe_diff

    blank = json.loads(json.dumps(R.Recipe().to_dict()))
    added_all = []
    for path, rec in fits:
        name = Path(path).name
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        d = json.loads(json.dumps(rec.to_dict()))
        changed = sorted(k for k in raw if k not in d or d[k] != raw[k])
        added = sorted(k for k in d if k not in raw)
        # an older file lacks the fields added since: they come back at their defaults
        _check(not changed and all(d[k] == blank[k] for k in added),
               f"{name}: load -> to_dict changed {changed}, added non-default {[k for k in added if d[k] != blank[k]]}")
        added_all += added
        _check(rec.sites and all(s.model in models.REGISTRY for s in rec.sites) and rec.nucleus
               and rec.larmor_frequency_MHz > 0 and raw.get("larmor_recipe_version") == R.RECIPE_VERSION,
               f"{name}: {len(rec.sites)} sites {[s.model for s in rec.sites]}, nucleus {rec.nucleus!r}, "
               f"version {raw.get('larmor_recipe_version')}")
        if ex is not None:
            # shipped recipes name their data relative to the folder holding examples/
            _check((ex.parent / rec.source_path / "acqus").is_file(), f"{name}: source {rec.source_path} not shipped")
        copy_path = out / f"resaved_{name}"
        rec.save(copy_path)
        back = R.Recipe.load(copy_path)
        _check(back.to_dict() == rec.to_dict(), f"{name}: save / reload changed the recipe")
        back.save(out / f"resaved2_{name}")
        _check((out / f"resaved2_{name}").read_bytes() == copy_path.read_bytes(), f"{name}: a second save differs")
    rec0 = fits[0][1]
    # the migration hook: an old-version dict passes through _MIGRATIONS hop by
    # hop, each noted -- driven exactly as tests/test_recipe_diff does, with a
    # temporary v0 hop (the shipped table is empty while the schema is v1)
    old = rec0.to_dict()
    old["larmor_recipe_version"] = 0
    old["sample"] = "v0 file"
    bump = lambda dd: {**dd, "sample": dd.get("sample", "") + " [migrated]"}  # noqa: E731
    had, prev = 0 in R._MIGRATIONS, R._MIGRATIONS.get(0)
    R._MIGRATIONS[0] = bump
    try:
        mig = R.Recipe.from_dict(old)
    finally:
        if had:
            R._MIGRATIONS[0] = prev
        else:
            R._MIGRATIONS.pop(0, None)
    _check(mig.sample == "v0 file [migrated]" and any("migrated from schema v0 to v1" in n for n in mig.notes)
           and [s.params for s in mig.sites] == [s.params for s in rec0.sites]
           and R.Recipe.from_dict(rec0.to_dict()).notes == rec0.notes,
           f"recipe migration: sample {mig.sample!r}, notes {mig.notes}")
    _check(R.run_migrations({"sample": "x"}, 0, {0: bump}, 1, "recipe")
           == ({"sample": "x [migrated]"}, ["recipe migrated from schema v0 to v1"])
           and R.run_migrations({"a": 1}, 0, {}, 2) == ({"a": 1}, []), "run_migrations")
    # forward compatibility: fields a newer LARMOR writes are dropped with a note
    newer = rec0.to_dict()
    newer["future_field"] = {"x": 1}
    newer["sites"][0]["colour"] = "#ff0000"
    fwd = R.Recipe.from_dict(newer)
    new_notes = fwd.notes[len(rec0.notes):]
    strip = lambda dd: {k: v for k, v in dd.items() if k != "notes"}  # noqa: E731
    _check(len(new_notes) == 2 and "future_field" in new_notes[0] and "colour" in new_notes[1]
           and strip(fwd.to_dict()) == strip(rec0.to_dict()), f"forward compatibility: notes {new_notes}")
    # recipe_diff: one value changed -> exactly one difference, naming it
    ref = rec0.to_dict()
    cur = copy.deepcopy(ref)
    i_site = len(cur["sites"]) - 1
    pname = "shift_fwhm_ppm"
    v_ref = ref["sites"][i_site]["params"][pname]["value"]
    cur["sites"][i_site]["params"][pname]["value"] = v_ref + 0.5
    rows = recipe_diff(cur, ref)
    n_par = sum(1 for s in ref["sites"] for p in s["params"] if p != "gl")
    moved = [r for r in rows if r["delta"] != 0.0]
    _check(len(rows) == n_par and len(moved) == 1 and moved[0]["site"] == i_site and moved[0]["param"] == pname
           and abs(moved[0]["delta"] - 0.5) < 1e-9 and abs(moved[0]["delta_pct"] - 50.0 / v_ref) < 1e-9
           and not any(r["param"] == "gl" for r in rows),
           f"recipe_diff: {len(rows)} rows (expected {n_par}), differences {moved}")
    extra = copy.deepcopy(cur)
    extra["sites"].append(copy.deepcopy(cur["sites"][0]))
    lone = [r for r in recipe_diff(extra, ref) if r["site"] == len(ref["sites"])]
    _check(lone and all(r["reference"] is None and r["delta"] is None for r in lone), f"recipe_diff extra site: {lone}")
    say(f"recipes: {len(fits)} ({', '.join(Path(p).name for p, _ in fits)}) load losslessly"
        + (f" ({len(added_all)} field(s) an older file lacks restored at their defaults)" if added_all else "")
        + f", re-save byte-stable; v0 -> v1 migration hop applied and noted; unknown fields dropped with 2 notes; "
        f"diff of one changed value -> 1 of {len(rows)} rows (site {i_site} {pname} +0.5, "
        f"{moved[0]['delta_pct']:+.2f} %)")


def _rtab_statuses(say, fits):
    """paramstatus beyond engine-constraints: the derived statuses of the
    shipped fits against what their files store (vary / expr) and the fit's
    own at-bound note, site_constraints, and the silent 'default' kind."""
    import json

    from larmor import models, paramstatus
    from larmor.recipe import Param

    lines = []
    for path, rec in fits:
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        held_file, linked_file = set(), set()
        for i, s in enumerate(raw["sites"]):
            for p, v in s["params"].items():
                if v.get("expr"):
                    linked_file.add((i, p))
                elif v.get("vary") is False:
                    held_file.add((i, p))
        st = paramstatus.recipe_statuses(rec)
        held = {k for k, s in st.items() if s.held}
        linked = {k for k, s in st.items() if s.kind == "linked"}
        at_bound = {k for k, s in st.items() if s.kind == "at_bound"}
        fit_note = any(n.startswith(paramstatus.AT_BOUND_NOTE_PREFIX) for n in rec.notes)
        _check(held == held_file and linked == linked_file and bool(at_bound) == fit_note,
               f"{Path(path).name}: held {sorted(held)} vs file {sorted(held_file)}, linked {sorted(linked)} vs "
               f"{sorted(linked_file)}, at a bound {sorted(at_bound)} vs the fit's note {fit_note}")
        for i in range(len(rec.sites)):
            sc = paramstatus.site_constraints(rec, i)
            mine = lambda keys: {p for (j, p) in keys if j == i}  # noqa: E731
            _check(set(sc["fixed"]) == {p for p in mine(held) if st[(i, p)].kind == "fixed"}
                   and {p for p, _e in sc["linked"]} == mine(linked) and {p for p, _s, _b in sc["at_bound"]} == mine(at_bound),
                   f"{Path(path).name} site_constraints({i}) {sc}")
        lines.append(f"{Path(path).name}: {len(held)} held, {len(linked)} linked, {len(at_bound)} at a bound")
    lb = models.get("czjzek").defaults()["line_fwhm_ppm"]
    s0 = paramstatus.param_status("czjzek", "line_fwhm_ppm", lb)
    s1 = paramstatus.param_status("czjzek", "line_fwhm_ppm", Param(2.0, vary=False))
    s2 = paramstatus.param_status("czjzek", "sigma_Cq_MHz", Param(0.05, min=0.05))
    _check(s0.kind == "default" and s0.held and s0.marker == "" and s0.csv_flag == "fixed"
           and s0.word == "held at the model default" and s1.kind == "fixed" and s1.marker == "†"
           and s1.word == "held fixed" and s2.kind == "at_bound" and s2.csv_flag == "at_min"
           and s2.word == "at its lower bound 0.05",
           f"param_status kinds: {s0} / {s1} / {s2}")
    say(f"paramstatus on the saved fits: {'; '.join(lines)} (= the files' vary / expr and the fit's own at-bound "
        f"note); site_constraints agree; Czjzek lb at its registry default -> '{s0.kind}' (held, no glyph, CSV "
        f"'{s0.csv_flag}'), pinned elsewhere -> '{s1.kind}' {s1.marker}")


# ------------------------------------------------------------ the parameter table
def _rtab_fittable(say, out: Path, fits):
    """larmor.fittable on the saved fits: FitEntry facts, the long and wide
    layouts, CSV written and re-read, TSV, every value against the recipe."""
    import csv

    from larmor import fittable, paramstatus

    entries = [fittable.load_fit_file(p) for p, _ in fits]
    for e, (p, rec) in zip(entries, fits):
        _check(e.name == fittable.fit_name(p) and e.source == str(p) and e.nucleus == rec.nucleus
               and e.larmor_MHz == rec.larmor_frequency_MHz and e.spin_rate_Hz == rec.spin_rate_Hz
               and e.rmsd == rec.fit_rmsd and len(e.sites) == len(rec.sites),
               f"FitEntry {e.name}: {e.nucleus} {e.larmor_MHz} {e.spin_rate_Hz} rmsd {e.rmsd} vs the recipe")
    flat = [(rec, i, site, pname, prm) for _p, rec in fits for i, site in enumerate(rec.sites)
            for pname, prm in site.params.items()]
    statuses = [paramstatus.recipe_statuses(rec) for _p, rec in fits]
    held = {(id(rec), i, pname) for (_p, rec), st in zip(fits, statuses) for (i, pname), s in st.items() if s.held}

    def num(cell):
        return float(cell) if cell != "" else None

    def close(a, b, rel):
        return (a is None and b is None) or (a is not None and b is not None and abs(a - b) <= rel * max(abs(b), 1e-300))

    hl, rl = fittable.build_long(entries)
    fittable.to_csv(hl, rl, out / "fit_table_long.csv")
    with open(out / "fit_table_long.csv", encoding="utf-8", newline="") as f:
        back = list(csv.reader(f))
    _check(tuple(back[0]) == fittable.LONG_HEADERS and len(back) == len(flat) + 1,
           f"long CSV: header {back[0]}, {len(back) - 1} rows for {len(flat)} parameters")
    for row, (rec, i, site, pname, prm) in zip(back[1:], flat):
        is_held = (id(rec), i, pname) in held
        want_status = "linked" if prm.expr else ("fixed" if is_held else "")
        # cells are %.8g: 1e-7 relative covers the rounding
        _check(row[3] == pname and row[2] == site.model and close(num(row[4]), prm.value, 1e-7)
               and close(num(row[5]), None if is_held else prm.stderr, 1e-7) and row[6] == want_status
               and close(num(row[7]), prm.min, 1e-7) and close(num(row[8]), prm.max, 1e-7),
               f"long CSV row {row} vs {site.label} {pname} = {prm}")
    tsv = [ln.split("\t") for ln in fittable.to_tsv(hl, rl).split("\n")]
    _check(tsv == back, "the TSV (clipboard) cells differ from the CSV's")
    hw, rw = fittable.build_wide(entries)
    fittable.to_csv(hw, rw, out / "fit_table_wide.csv")
    with open(out / "fit_table_wide.csv", encoding="utf-8", newline="") as f:
        wide = list(csv.reader(f))
    sites = [(rec, i, site) for _p, rec in fits for i, site in enumerate(rec.sites)]
    _check(wide[0] == list(hw) and wide[1:] == [[str(c) for c in r] for r in rw] and hw[:3] == ["fit", "site", "model"]
           and hw[-5:] == ["RMSD", "nucleus", "ν0 (MHz)", "νrot (Hz)", "source"] and len(rw) == len(sites),
           f"wide table: header {hw}, {len(rw)} rows for {len(sites)} sites")
    # the parameter columns: those some site uses, in PARAM_COLUMNS order, then
    # any key the curated list lacks, headed from the model's own definition
    curated = dict(fittable.PARAM_COLUMNS)
    col_keys = [k for k, _h in fittable.PARAM_COLUMNS if any(k in s.params for _r, _i, s in sites)]
    heads = [fittable.header_text(curated[k]) for k in col_keys]
    for _r, _i, s in sites:
        for k in s.params:
            if k not in curated and k not in col_keys:
                col_keys.append(k)
                heads.append(fittable.header_text(fittable.auto_header(s.model, k)))
    _check(hw[3:-5] == heads, f"wide parameter columns {hw[3:-5]} vs {heads}")
    n_cells = n_held = 0
    for row, (rec, i, site) in zip(rw, sites):
        _check(row[2] == site.model and close(float(row[-5]), rec.fit_rmsd, 1e-4), f"wide row {row[:3]} / RMSD {row[-5]}")
        for col, key in enumerate(col_keys, start=3):
            prm = site.params.get(key)
            cell = row[col]
            if prm is None:
                _check(cell == "", f"wide cell {key} of {site.label}: {cell!r} for a parameter the model lacks")
                continue
            # value %.5g: 1e-4 relative; a held value never carries an error bar
            is_held = (id(rec), i, key) in held
            _check(close(float(cell.split(" ")[0]), prm.value, 1e-4) and not (is_held and "±" in cell),
                   f"wide cell {key} of {site.label}: {cell!r} vs {prm}")
            n_cells += 1
            n_held += is_held
    say(f"fittable: {len(entries)} fits -> long {len(rl)} rows / wide {len(rw)} x {len(hw)}; CSV + TSV re-read: every "
        f"value, error, status word and bound back to 8 digits; {n_cells} wide cells = the recipes ({n_held} held, "
        f"no error bar)")


# ------------------------------------------------------------ project bundle
def _rtab_project(say, out: Path, fits, ex):
    """project.build_bundle on every workspace kind, then the restore half
    after the project folder moved with its data: load_bundle, the recipe and
    spectrum back, relocate (2D by reference), relocate_batch_state; the v1
    -> v2 migration of an old bundle."""
    import json
    import shutil

    import larmor
    from larmor import project
    from larmor.recipe import Recipe

    proj = out / "project"
    data = proj / "data"
    two = data / "S1" / "35" / "pdata" / "1" / "2rr"
    two.parent.mkdir(parents=True, exist_ok=True)
    two.write_bytes(b"\0" * 8)                   # a 2D map is saved BY REFERENCE: its bytes are never read
    x = np.linspace(-20.0, 40.0, 121)
    y = np.exp(-0.5 * ((x - 15.0) / 3.0) ** 2) + 0.01 * np.sin(x)
    ov = _write_csv_spectrum(data / "overlay.csv", x, 0.5 * y)
    if ex is not None:
        spec = json.loads((ex / "pCABS2-4_fits.figure.json").read_text(encoding="utf-8"))
    else:
        spec = {"kind": "batch_grid", "panels": [{"recipe": str(fits[0][0]), "title": "fit"}], "figsize": (7.0, 3.0)}
    rec = fits[0][1]
    state = {"paths": [str(ov)], "per_spectrum": {str(ov): {"excluded": [1]}},
             "result": {"recipes": [{"source_path": str(ov)}]},
             "series": {"version": 1, "rows": [{"source_path": str(ov), "display_name": "ov"}], "columns": []}}
    workspaces = [
        {"kind": "1d", "title": "glass A", "series": {"id": "s1", "order": 0, "name": "glass A"},
         "snap": {"source_path": str(data / "a.csv"), "recipe": rec.to_dict(), "hidden": {2, 0},
                  "exp_ppm": x, "exp_amp": y,
                  "overlays": [{"label": "ov", "source": str(ov), "color": "#1f77b4", "visible": False, "scale": 2.0},
                               {"label": "typed in", "source": ""}]}},
        {"kind": "1d", "title": "no spectrum", "snap": {"exp_ppm": None}},
        {"kind": "2d", "title": "3QMAS", "snap": {"source_path": str(two), "recipe": None, "hidden": set(),
                                                  "fittable": True,
                                                  "view2d": {"ops": [{"op": "shear", "factor": 0.5}], "nlevels": 9}}},
        {"kind": "figure", "title": "figure", "snap": {"spec": spec}},
        {"kind": "batch", "title": "batch", "snap": {"state": state}}]
    bundle, dropped = project.build_bundle(workspaces, 3, proj)
    kinds = [w["kind"] for w in bundle["workspaces"]]
    _check(dropped == 1 and bundle["active"] == 2 and kinds == ["1d", "2d", "figure", "batch"]
           and bundle["larmor_project_version"] == project.PROJECT_BUNDLE_VERSION
           and bundle["software"].get("larmor") == larmor.__version__
           and project.summary(bundle) == "1 spectrum, 1 2D map, 1 figure, 1 batch session",
           f"build_bundle: dropped {dropped}, active {bundle['active']}, kinds {kinds}, {project.summary(bundle)!r}")
    (proj / "session.larproj.json").write_text(json.dumps(bundle), encoding="utf-8")
    # the project moves together with its data; the old folder is gone
    moved = out / "project_moved"
    shutil.copytree(proj, moved)
    shutil.rmtree(proj)
    back, notes = project.load_bundle(moved / "session.larproj.json")
    w1, w2, wf, wb = back["workspaces"]
    _check(notes == [] and back["active"] == 2
           and Recipe.from_dict(w1["recipe"]).to_dict() == rec.to_dict()
           and np.array_equal(w1["exp_ppm"], x) and np.array_equal(w1["exp_amp"], y) and w1["hidden"] == [0, 2]
           and w1["overlays"] == [{"label": "ov", "color": "#1f77b4", "visible": False, "source": str(ov),
                                   "scale": 2.0, "shift": 0.0, "yoff": 0.0}]
           and w1["series"] == workspaces[0]["series"],
           f"reopened 1D entry: notes {notes}, hidden {w1['hidden']}, overlays {w1['overlays']}")
    new_two = moved / "data" / "S1" / "35" / "pdata" / "1" / "2rr"
    got_two = project.relocate(w2["source_path"], moved, w2["source_rel"])
    _check(not Path(w2["source_path"]).exists() and got_two is not None and Path(got_two) == new_two
           and w2["view"]["ops"] == [{"op": "shear", "factor": 0.5}] and w2["view"]["nlevels"] == 9.0
           and w2["view"]["flip"] == {"f2": True, "f1": False} and w2["fittable"] is True,
           f"reopened 2D entry: relocated to {got_two}, view {w2['view']}")
    _check(wf["spec"] == json.loads(json.dumps(spec)), "the figure spec changed through the bundle")
    st, missing = project.relocate_batch_state(wb["state"], moved)
    new_ov = str(moved / "data" / "overlay.csv")
    _check(missing == [] and st["paths"] == [new_ov] and list(st["per_spectrum"]) == [new_ov]
           and st["result"]["recipes"][0]["source_path"] == new_ov
           and st["series"]["rows"][0]["source_path"] == new_ov and wb["state"]["paths"] == [str(ov)],
           f"relocate_batch_state: {st['paths']} missing {missing}")
    v1 = {"larmor_project_version": 1, "active": 0,
          "workspaces": [{"title": "old", "recipe": w1["recipe"], "exp_ppm": [1.0, 2.0], "exp_amp": [0.0, 1.0]}]}
    (moved / "v1.larproj.json").write_text(json.dumps(v1), encoding="utf-8")
    b1, notes1 = project.load_bundle(moved / "v1.larproj.json")
    _, notes9 = project.migrate_bundle({"larmor_project_version": 99, "workspaces": []})
    try:
        project.load_bundle(fits[0][0])
        refused = ""
    except ValueError as exc:
        refused = str(exc)
    _check(b1["larmor_project_version"] == project.PROJECT_BUNDLE_VERSION and [w["kind"] for w in b1["workspaces"]] == ["1d"]
           and notes1[:1] == ["project migrated from schema v1 to v2"] and b1["workspaces"][0]["recipe"] == w1["recipe"]
           and len(notes9) == 1 and "newer LARMOR" in notes9[0] and "not a LARMOR project bundle" in refused,
           f"bundle migration: {b1.get('larmor_project_version')} {notes1} / newer {notes9} / recipe file {refused!r}")
    say(f"project bundle: {project.summary(bundle)} (1 empty workspace dropped, active 3 -> 2); reopened after the "
        f"folder moved: recipe, {x.size}-point spectrum, hidden lines, overlay, figure spec identical; 2D map and "
        f"batch spectra relocated through their project-relative paths; v1 bundle -> v{b1['larmor_project_version']} "
        f"({notes1[0]})")


# ------------------------------------------------------------ series table
def _rtab_series(say, out: Path):
    """larmor.series_table: names from Bruker-shaped paths, a joined
    composition CSV, replicate statistics, OLS, the series and DUST CSVs,
    the species-bar spec and the .series.json round trip."""
    import csv

    from larmor import series_table as st

    root = out / "series_paths"                  # Bruker-shaped paths that do not exist: names come from the path
    paths = [root / "2026-03b" / "03232026_P5-Bi8-12_SS_ALP" / "3104" / "pdata" / "1" / "1r",
             root / "2026-05" / "04272026_P5-Bi8-12_SS_ALP" / "3114" / "pdata" / "1" / "1r",
             root / "2026-05" / "05082026_P5-Bi8-12_SS_ALP" / "3102",
             root / "2026-03b" / "03232026_P5-Bi0_SS_ALP" / "3104"]
    t = st.SeriesTable.from_paths([str(p) for p in paths])
    names = ["P5-Bi8-12 (03232026)", "P5-Bi8-12 (04272026)", "P5-Bi8-12 (05082026)", "P5-Bi0"]
    _check(t.labels() == names and t.groups() == {"P5-Bi8-12": [0, 1, 2], "P5-Bi0": [3]} and t.has_replicates()
           and t.rows[3].folder == "03232026_P5-Bi0_SS_ALP", f"from_paths: {t.labels()} groups {t.groups()}")
    comp = out / "series_composition.csv"
    comp.write_text("sample;P2O5_mol;P2O5_mol_sd;Bi2O3 (mol%);Vm;method\n"
                    "P5-Bi8-12;5.2;0.1;7.9;27.3;EPMA\nP5Bi0;5.0;0.1;0.0;26.7;nominal\n", encoding="utf-8")
    headers, rows, sample_h = st.read_wide_csv(comp)
    matches = st.propose_mapping(t, rows, sample_h)
    how = [(m.how, m.csv_row) for m in matches]
    _check(sample_h == "sample" and how == [("exact", 0)] * 3 + [("normalised", 1)], f"propose_mapping: {how}")
    cols = st.apply_join(t, rows, sample_h, matches, ["P2O5_mol", "Bi2O3 (mol%)", "Vm", "method"], tag="analysed",
                         source=comp.name)
    xv, xe = t.x_values("P2O5_mol")
    _check([c.key for c in cols] == ["P2O5_mol", "Bi2O3 (mol%)", "Vm", "method"] and cols[0].err_key == "P2O5_mol_sd"
           and xv.tolist() == [5.2, 5.2, 5.2, 5.0] and xe.tolist() == [0.1] * 4
           and all(r.values["method"] is None for r in t.rows), f"apply_join: {xv} ± {xe}")
    n4 = np.array([0.40, 0.42, 0.45, 0.30])
    n4e = np.array([0.010, 0.012, 0.011, 0.020])
    rs = st.replicate_stats(t.groups(), n4, n4e, x=xv, xerr=xe)
    _check(rs["labels"] == ["P5-Bi8-12", "P5-Bi0"] and np.allclose(rs["y"], [n4[:3].mean(), n4[3]], rtol=0, atol=1e-12)
           and np.allclose(rs["yerr"], [np.std(n4[:3], ddof=1), n4e[3]], rtol=0, atol=1e-12) and rs["n"].tolist() == [3, 1]
           and np.allclose(rs["x"], [5.2, 5.0]) and np.allclose(rs["xerr"], [0.1, 0.1]),
           f"replicate_stats {rs}")
    xs = np.array([0.0, 5.0, 10.0, 15.0, np.nan, 20.0])
    ys = 0.0123 * xs + 0.381
    ys[4] = 9.9                                  # paired with the NaN x: skipped
    fit = st.ols(xs, ys)
    _check(abs(fit["slope"] - 0.0123) < 1e-12 and abs(fit["intercept"] - 0.381) < 1e-12 and abs(fit["r"] - 1.0) < 1e-12
           and fit["n"] == 5 and fit["slope_err"] < 1e-9, f"ols on exact linear data: {fit}")
    p = st.write_series_csv(t, out / "series_table_series.csv", extra={"N4": n4})
    with open(p, encoding="utf-8", newline="") as f:
        back = list(csv.reader(f))
    want_head = ["position", "name", "group", "folder", "title", "source_path", "analysed P2O5_mol",
                 "analysed P2O5_mol ±", "analysed Bi2O3 (mol%)", "analysed Vm", "analysed method", "N4"]
    _check(back[0] == want_head and [r[1] for r in back[1:]] == names and [r[0] for r in back[1:]] == ["1", "2", "3", "4"]
           and [float(r[6]) for r in back[1:]] == xv.tolist() and [float(r[-1]) for r in back[1:]] == n4.tolist()
           and back[1][5] == str(paths[0]) and back[1][10] == "", f"write_series_csv: {back[:2]}")
    head, drows, skipped = st.dust_rows(t, n4, n4e, average=False)
    head_a, arows, _ = st.dust_rows(t, n4, n4e, average=True)
    _check(head == head_a == ["Sample", "P2O5", "Bi2O3", "N4_measured", "N4_measured_err"]
           and skipped == ["analysed Vm", "analysed method"] and [r[0] for r in drows] == names
           and drows[0][1:] == ["5.2", "7.9", "0.4", "0.01"] and [r[0] for r in arows] == ["P5-Bi8-12", "P5-Bi0"]
           and abs(float(arows[0][3]) - n4[:3].mean()) < 1e-6 and abs(float(arows[0][4]) - np.std(n4[:3], ddof=1)) < 1e-6
           and arows[1][1:] == ["5", "0", "0.3", "0.02"], f"dust_rows: {drows} / {arows}")
    pops = np.array([[60.0, 40.0], [55.0, np.nan], [50.0, 50.0], [70.0, 30.0]])
    spec = st.species_bar_spec(t.labels(), ["A", "B"], pops, xlabel="analysed P2O5_mol", colors=["#1f77b4", None])
    _check(spec["kind"] == "species_bar" and spec["categories"] == names
           and spec["series"] == [{"label": "A", "values": [60.0, 55.0, 50.0, 70.0], "color": "#1f77b4"},
                                  {"label": "B", "values": [40.0, 0.0, 50.0, 30.0]}]
           and spec["xlabel"] == "analysed P2O5_mol" and spec["ylabel"] == "population (%)", f"species_bar_spec {spec}")
    t.save(out / "series.series.json")
    _check(st.SeriesTable.load(out / "series.series.json") == t, "the .series.json round trip changed the table")
    say(f"series table: {', '.join(names)}; joined P2O5 {xv.tolist()} ± {xe[0]} ('normalised' P5Bi0 match); "
        f"replicates N4 {rs['y'][0]:.4f} ± {rs['yerr'][0]:.4f} (numpy, ddof 1, n 3); OLS slope {fit['slope']:.4f}, "
        f"intercept {fit['intercept']:.3f}, r {fit['r']:.3f}; series CSV {len(back[0])} columns and DUST rows "
        f"{arows} re-read; species bars; .series.json round trip")


# ------------------------------------------------------------ batch grid
def _rtab_grid(say, out: Path):
    """larmor.series_grid: saved fits found from a list, a folder and a batch
    CSV (exact scope before substring), panels, the fit rebuilt from the
    CSV's own rows, an excluded site, a manual locate; the CSV pivoted per
    scope by series_table.read_batch_csv_wide."""
    import csv

    from larmor import series_grid, series_table
    from larmor.recipe import Param, Recipe, SiteModel

    gdir = out / "grid"
    gdir.mkdir(parents=True, exist_ok=True)
    x = np.linspace(-30.0, 30.0, 121)
    saved = {}
    for k, (scope, stem, with_data) in enumerate((("AbNa", "AbNa", True), ("Ab", "Ab", True),
                                                   ("NS3", "NS3_batch", False))):
        rec = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sample=scope, sites=[
            SiteModel(model="gauss_lor", label="A", params={
                ISO: Param(10.0 + 0.25 * k, stderr=0.1), "shift_fwhm_ppm": Param(4.0 + 0.5 * k, stderr=0.2),
                "amplitude": Param(100.0 - 10.0 * k, stderr=2.0), "gl": Param(1.0, vary=False)}),
            SiteModel(model="gauss_lor", label="B", params={
                ISO: Param(-8.0), "shift_fwhm_ppm": Param(3.0), "amplitude": Param(40.0 + 5.0 * k),
                "gl": Param(1.0, vary=False)})])
        if with_data:
            yk = (100.0 - 10.0 * k) * np.exp(-0.5 * ((x - 10.0) / 2.0) ** 2)
            rec.source_path = str(_write_csv_spectrum(gdir / f"{stem}_raw.csv", x, yk, nucleus="11B", larmor=160.0))
            rec.source_kind = "csv"
        rec.save(gdir / f"{stem}.recipe.json")
        saved[scope] = (gdir / f"{stem}.recipe.json", rec)
    # a batch export in the long layout (shared ladder + per-scope rows, the
    # derived population_pct rows included), written in the CSV's row order
    table = gdir / "batch_table.csv"
    with open(table, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["scope", "site", "label", "param", "value", "stderr", "model", "source_path"])
        for j, lab in enumerate("AB"):
            w.writerow(["shared", f"s{j}", lab, "gl", "1", "", "gauss_lor", ""])
        for scope, (_p, rec) in saved.items():
            for j, site in enumerate(rec.sites):
                for pname in (ISO, "shift_fwhm_ppm", "amplitude"):
                    prm = site.params[pname]
                    w.writerow([scope, f"s{j}", site.label, pname, repr(prm.value),
                                "" if prm.stderr is None else repr(prm.stderr), "gauss_lor", rec.source_path])
                w.writerow([scope, f"s{j}", site.label, "population_pct", "50.0", "", "gauss_lor", rec.source_path])
    files = [str(saved[s][0]) for s in ("AbNa", "Ab", "NS3")]
    empty = gdir / "empty"
    empty.mkdir(exist_ok=True)
    r_list, r_dir, r_empty = (series_grid.resolve_paths(files), series_grid.resolve_paths(str(gdir)),
                              series_grid.resolve_paths(str(empty)))
    near = series_grid.find_recipes_near_csv(table)
    hints = series_grid.read_csv_hints(table)
    _check(r_list == (files, []) and sorted(r_dir[0]) == sorted(files) and not r_dir[1] and r_empty[0] == []
           and len(r_empty[1]) == 1 and near == files
           and hints == {s: {"source_path": rec.source_path, "models": ("gauss_lor",)} for s, (_p, rec) in saved.items()},
           f"resolve_paths / find_recipes_near_csv: {near} (expected {files}), hints {hints}")
    panels, warns = series_grid.load_panels(str(gdir))
    by = {p.sample: p for p in panels}
    via_csv, warns_csv = series_grid.load_panels(str(table))
    _check(not warns and sorted(by) == ["Ab", "AbNa", "NS3"] and by["Ab"].has_data and by["AbNa"].has_data
           and not by["NS3"].has_data and all(p.models == ("gauss_lor", "gauss_lor") and p.n_sites == 2 for p in panels)
           and not warns_csv and sorted(p.sample for p in via_csv) == ["Ab", "AbNa", "NS3"]
           and not any(p.reconstructed or p.needs_manual for p in via_csv),
           f"load_panels: {[(p.sample, p.has_data, p.models) for p in panels]} {warns}; via CSV "
           f"{[(p.sample, p.reconstructed, p.needs_manual) for p in via_csv]} {warns_csv}")
    by_scope = series_grid.csv_rows_by_scope(table)
    _p, ab = saved["Ab"]
    rebuilt = series_grid.recipe_from_csv_rows(by_scope["shared"], by_scope["Ab"], ab.source_path)
    same = all(rebuilt.sites[j].params[pn].value == site.params[pn].value
               and rebuilt.sites[j].params[pn].stderr == site.params[pn].stderr
               for j, site in enumerate(ab.sites) for pn in site.params)
    _check(same and [s.label for s in rebuilt.sites] == ["A", "B"] and rebuilt.nucleus == "11B"
           and rebuilt.larmor_frequency_MHz == 160.0 and all("population_pct" not in s.params for s in rebuilt.sites),
           f"recipe_from_csv_rows: {[(s.label, {k: v.value for k, v in s.params.items()}) for s in rebuilt.sites]}, "
           f"{rebuilt.nucleus} {rebuilt.larmor_frequency_MHz}")
    only_a = [r for r in by_scope["AbNa"] if r["site"] == "s0"]
    excl = series_grid.recipe_from_csv_rows(by_scope["shared"], only_a)
    amp_b = excl.sites[1].params["amplitude"] if len(excl.sites) == 2 else None
    _check(amp_b is not None and excl.sites[1].label == "B" and amp_b.value == 0.0 and amp_b.vary is False
           and amp_b.min == 0.0 and amp_b.max == 0.0 and excl.sites[0].params["amplitude"].value == 100.0,
           f"an excluded site is not zeroed at its index: {[(s.label, s.params['amplitude']) for s in excl.sites]}")
    lost = series_grid.Panel(path="", sample="g9", nucleus="", models=(), has_data=False, needs_manual=True)
    found = series_grid.resolve_manual(lost, gdir / "Ab_raw.csv")
    _check(found.has_data and not found.needs_manual and found.data_path == str(gdir / "Ab_raw.csv") and lost.needs_manual,
           f"resolve_manual: {found}")
    headers, wrows, scope_h = series_table.read_batch_csv_wide(table)
    amps = {r["scope"]: float(r["A amplitude"]) for r in wrows}
    errs = {r["scope"]: float(r["A amplitude ±"]) for r in wrows}
    _check(scope_h == "scope" and amps == {s: rec.sites[0].params["amplitude"].value for s, (_p, rec) in saved.items()}
           and set(errs.values()) == {2.0} and series_table.is_long_batch_csv(series_table.read_wide_csv(table)[0]),
           f"read_batch_csv_wide: {amps} ± {errs}")
    say(f"series_grid: list / folder / batch CSV -> {', '.join(Path(f).name for f in near)} (exact scope before "
        f"substring); panels has_data {[(p.sample, p.has_data) for p in panels]}; the Ab fit rebuilt from the CSV "
        f"rows = the saved values ({rebuilt.nucleus} at {rebuilt.larmor_frequency_MHz:g} MHz from its spectrum); "
        f"excluded site B zeroed at index 1; manual locate; CSV pivoted per scope {amps}")


# ------------------------------------------------------------ aliases
def _rtab_aliases(say, out: Path):
    """larmor.aliases with both stores pointed under the stage folder (never
    the user's, even with --real-settings): an alias reaches scan, the series
    labels and the window title; the refused renames; a rename on disk that
    moves the alias and writes one rename-log record."""
    import json

    import larmor
    from larmor import aliases, series_table
    from larmor.io import scan

    store = out / "alias_store"
    before = (os.environ.get("LARMOR_ALIASES"), os.environ.get("LARMOR_RENAME_LOG"))
    with _rtab_env(LARMOR_ALIASES=store / "aliases.json", LARMOR_RENAME_LOG=store / "rename_log.jsonl"):
        _check(aliases.aliases_path() == store / "aliases.json" and aliases.rename_log_path() == store / "rename_log.jsonl",
               f"alias stores at {aliases.aliases_path()} / {aliases.rename_log_path()}")
        folder = "01192026_SR31649_Base0Ca_SS_ALP"
        sample = out / "alias_data" / folder
        acqus = "##TITLE= Parameter file\n##JCAMPDX= 5.0\n##$NUC1= <27Al>\n##$SFO1= 130.3\n##END=\n"
        for e in ("10", "12"):
            (sample / e / "pdata" / "1").mkdir(parents=True, exist_ok=True)
            (sample / e / "acqus").write_text(acqus, encoding="ascii")
            (sample / e / "pdata" / "1" / "1r").write_bytes(b"\0" * 8)
            (sample / e / "pdata" / "1" / "title").write_text("27Al MAS\n", encoding="ascii")
        expno = sample / "10"
        r1 = str(expno / "pdata" / "1" / "1r")

        def label():
            return series_table.SeriesTable.from_paths([r1]).labels()[0]

        start = (aliases.alias_for(sample), aliases.display_name(sample), label(), aliases.window_label(r1))
        _check(start == ("", folder, "Base0Ca", "1r"), f"before any alias: {start}")
        aliases.set_alias(sample, "Glass A")
        nm = scan.sample_name(expno, "27Al MAS")
        probe = str(sample).upper() if os.name == "nt" else str(sample)     # Windows keys ignore case
        _check((store / "aliases.json").is_file() and aliases.alias_for(probe) == "Glass A"
               and (nm.key, nm.source) == ("Glass A", "alias") and label() == "Glass A"
               and aliases.window_label(r1) == "Glass A · 10",
               f"sample alias: {nm}, label {label()!r}, title {aliases.window_label(r1)!r}")
        aliases.set_alias(expno, "run 10")
        _check(label() == "run 10" and aliases.window_label(r1) == "run 10",
               f"EXPNO alias: label {label()!r}, title {aliases.window_label(r1)!r}")
        refusals = {}
        for name, kw in (("echo", {}), ("12", {}), ("   ", {}), ("13", {"open_paths": [r1]})):
            try:
                aliases.check_rename(expno, name, **kw)
                refusals[name] = ""
            except aliases.RenameError as exc:
                refusals[name] = str(exc)
        _check("must stay a number" in refusals["echo"] and "already exists" in refusals["12"]
               and "empty" in refusals["   "] and "open in LARMOR" in refusals["13"] and expno.is_dir()
               and not aliases.rename_log_path().exists(), f"check_rename refusals: {refusals}")
        new = aliases.rename_folder(expno, "11")
        log = [json.loads(ln) for ln in aliases.rename_log_path().read_text(encoding="utf-8").splitlines()]
        _check(new == sample / "11" and new.is_dir() and not expno.exists() and (new / "pdata" / "1" / "1r").is_file()
               and aliases.alias_for(new) == "run 10" and aliases.alias_for(expno) == "" and len(log) == 1
               and log[0]["kind"] == "expno" and Path(log[0]["old"]) == expno and Path(log[0]["new"]) == new
               and log[0]["alias_moved"] is True and log[0]["larmor"] == larmor.__version__,
               f"rename_folder: {new}, alias {aliases.alias_for(new)!r}, log {log}")
        aliases.set_alias(new, "")
        aliases.set_alias(sample, folder)            # the folder's own name removes the alias
        _check(aliases.load() == {}, f"aliases left after removal: {aliases.load()}")
    _check((os.environ.get("LARMOR_ALIASES"), os.environ.get("LARMOR_RENAME_LOG")) == before,
           "the alias store variables were not restored")
    say(f"aliases: 'Glass A' then 'run 10' reach scan.sample_name, the series label and the window title; "
        f"{len(refusals)} renames refused ({', '.join(repr(k.strip() or 'empty') for k in refusals)}); EXPNO 10 -> 11 "
        f"moved its alias and wrote 1 rename-log record; stores under {store.name}/ only")
