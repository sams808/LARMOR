"""larmor.pulseprog: the Bruker pulse-program parser, the acqus resolver,
the not-to-scale timeline and the one-line summary.

Fixtures: the two shipped programs (examples/pCABS2-4/1118 is TopSpin's
``zg`` as compiled at acquisition, examples/pCABS2-4/3620 is ``mp3qdfsz``,
a 4-pulse z-filtered 3QMAS with a shaped DFS pulse and nested 2D loops)
plus inline snippets in the style of the sequences met on real
instruments: a rotor-synchronised Hahn echo with an ``aq_prot`` include, a
saturation-recovery pseudo-2D with a ``vd`` list, a QCPMG echo train with
explicit ``adc`` acquisition, a CP/MAS with simultaneous contact pulses
and decoupling, an HMQC with a multi-line ``(center ...)`` block and an
uncompiled ``zg`` with ``#include`` and the ``mc`` macro."""
from pathlib import Path

import pytest

from larmor import pulseprog as pp

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples" / "pCABS2-4"
ZG = EX / "1118"
MQ = EX / "3620"

pytestmark = pytest.mark.skipif(not (ZG / "pulseprogram").exists(),
                                reason="shipped example data missing")


def _acqus(expno: Path, name: str = "acqus") -> dict:
    from larmor.referencing import _jcamp
    return _jcamp(expno / name) or {}


# ------------------------------------------------------------------ snippets
HAHN = """\
# 1 "/opt/topspin3.6.2/exp/stan/nmr/lists/pp/user/hahnecho.nmrfam"
;hahnecho.nmrfam
;hahn echo for half-echo aquisition
; Author: Alexander Paterson
;$COMMENT=Hahn echo with half-echo acquisition
;$CLASS=Solids
;PARAMETERS:
;ns : 16*n as implemented
;d1 : recycle delay
;pl1 : power level for 90 and 180
;p1 : 90 degree pulse
;p2 : p1*2, 180 degree pulse
;d6 : first delay period, calculated as d6=((1s*l1)/cnst31)-(p1/2)-(p2/2)
;d7 : second delay period; may require optimization due to hardware limits
;del7 : calculates approximate d7 value
;cnst31 : MAS spin rate (or = 1e6 for static)
;l1 : Number of rotor cycles for MAS (or d6 in microseconds for static)
;cnst2 : to select between solid echo (1) or Hahn echo (2)

"p2=p1*cnst2"
"d6=((1s*l1)/cnst31)-(p1/2)-(p2/2)"
define delay del7
"del7=((1s*l1)/cnst31)-(p2/2)"
"acqt0=0"

prosol relations=<solids_default>

# 1 "mc_line 66 file /opt/topspin3.6.2/exp/stan/nmr/lists/pp/user/hahnecho.nmrfam exp. def. part of mc cmd. before ze"
; dimension 1D; AQ_mode
define delay MCWRK
define delay MCREST
"MCREST = 10m - 10m"
"MCWRK = 0.333333*10m"

    dccorr
# 66 "/opt/topspin3.6.2/exp/stan/nmr/lists/pp/user/hahnecho.nmrfam"
1 ze
2 MCWRK  * 2
LBLF0, MCWRK
  MCREST
  del7 ; displays del7 on acquisition panel
  d1

  (p1 pl1 ph1):f1

  d6

  (p2 ph2):f1

  d7

  go=2 ph31

  MCWRK wr #0
  MCWRK zd
  lo to LBLF0 times td0
exit

ph1=0 1 2 3
ph2=0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3
ph31=0 3 2 1 2 1 0 3
"""

HAHN_ACQUS = {"P": [0, 3.25, 6.5] + [0] * 61, "D": [0, 30.0, 0, 0, 0, 0, 0, 24.75e-6] + [0] * 56,
              "PLW": [0, 28.8] + [0] * 62, "L": [1] * 32,
              "CNST": [1, 1, 2] + [1] * 28 + [35714.0] + [1] * 32,
              "NS": 16, "TD": 4096, "SW_h": 413907.0, "TD0": 1, "NUC1": "19F", "NUC2": "off"}

SATREC = """\
# 1 "/opt/topspin3.6.2/exp/stan/nmr/lists/pp/satrect1"
;satrect1
;saturation recovery T1 experiment
;pl1 : X power level
;p1 : X 90 degree pulse
;d1 : recycle delay
;d20 : delay in saturation pulse train
;l20 : number of pulses in saturation pulse train, 0 if undesired
;vdlist : list containing tau delays
;$DIM=pseudo 2D

prosol relations=<solids_default>

"acqt0=-p1/2"
; dimension 2D; AQ_mode  (F1) QF
define delay MCWRK
define delay MCREST
"MCREST = 10m - 10m"
"MCWRK = 0.500000*10m"

    dccorr
1 ze
LBLAV, MCWRK
2 MCWRK
LBLF1, MCWRK
  MCREST
  d1

# 1 "/opt/topspin3.6.2/exp/stan/nmr/lists/pp/aq_prot.incl" 1
;aq_prot.incl
;protect against too long acquisition time
1m
if "aq < 50.1m" goto Passaq
2u
print "acquisition time exceeds 50m limit!"
goto HaltAcqu
Passaq, 1m
# 46 "/opt/topspin3.6.2/exp/stan/nmr/lists/pp/satrect1" 2

3 d20
  (p1 pl1 ph4):f1
  lo to 3 times l20
  vd			;recovery delay
  (p1 pl1 ph1):f1
  go=2 ph31
  MCWRK  wr #0 if #0 zd ivd
  lo to LBLF1 times td1
  MCWRK rf #0
  lo to LBLAV times tdav
HaltAcqu, 1m
exit

ph1= 0 0 2 2 1 1 3 3
ph4= 0
ph31= 0 0 2 2 1 1 3 3
"""

SATREC_ACQUS = {"P": [0, 0.5] + [0] * 62, "D": [0, 0.1] + [0] * 18 + [1e-3] + [0] * 43,
                "PLW": [0, 150.0] + [0] * 62, "L": [1] * 20 + [200] + [1] * 11,
                "NS": 8, "TD": 8192, "SW_h": 100000.0, "TD0": 1, "TDav": 1,
                "VDLIST": "<T1_standard_array>", "NUC1": "23Na"}

QCPMG = """\
;qcpmg
;quadrupolar Carr-Purcell Meiboom-Gill echo train
;p1 : 90 degree pulse
;p2 : 180 degree pulse
;d3 : half echo delay
;d6 : acquisition window
;l22 : number of echoes
;pl12 : decoupling power
;pcpd2 : decoupling pulse length

"d3=d6/2-p2/2"
"acqt0=0"

1 ze
2 d1 pl1:f1
  1u pl12:f2
  (p1 ph1):f1
  d3
  1u cpds2:f2
  (p2 ph2):f1
  d3
  1u adc ph31 syrec
  1u ph30:r
  d6
3 (p2 ph2):f1
  d6
  lo to 3 times l22
  1u eoscnp
  1m do:f2
  rcyc=2
  10m wr #0
exit

ph1=0 2 1 3
ph2=1 1 2 2 3 3 0 0
ph30=0
ph31=0 2 1 3
"""

QCPMG_ACQUS = {"P": [0, 2.0, 4.0] + [0] * 61, "D": [0, 0.5, 0, 0, 0, 0, 200e-6] + [0] * 57,
               "PLW": [0, 100.0] + [0] * 10 + [20.0] + [0] * 51, "L": [1] * 22 + [64] + [1] * 9,
               "PCPD": [100, 100, 6.8] + [100] * 7, "CPDPRG": ["", "", "spinal64"] + [""] * 6,
               "NS": 128, "TD": 32768, "SW_h": 2e6, "TD0": 1, "NUC1": "35Cl", "NUC2": "1H"}

CPMAS = """\
;cp
;basic cp experiment
;pl1 : f1 channel - power level for contact
;pl2 : f2 channel - power level for 90 (1H)
;pl12 : f2 channel - decoupling power
;sp0 : f2 channel - ramp for contact
;p3 : f2 channel - 90 degree pulse
;p15 : contact time
;pcpd2 : f2 channel - pulse length in decoupling sequence

1 ze
2 d1 do:f2
  (p3 pl2 ph1):f2
  (p15 pl1 ph2):f1 (p15:sp0 ph10):f2
  1u cpds2:f2
  go=2 ph31
  1m do:f2
  30m mc #0 to 2 F0(zd)
exit

ph1= 1 3
ph2= 0 0 2 2 1 1 3 3
ph10= 0
ph31= 0 2 2 0 1 3 3 1
"""

HMQC = """\
;pf_hmqc
;p1 		: pi/2 canal F1 pl1
;p11 		: pi/2 canal F2 pl2
"p5=2*p1"
"d0=p5"
"l0=1"
"d6=(l6*l0)*(1s/cnst31)-(p1/2)-(p5/2)-p11"
define delay zfilter
"zfilter=5*(1s/cnst31)"
	ze
1	1m
2 MCWRK pl1:f1 pl2:f2
	(p1 pl1 ph1):f1
	d6
(center
	(p11 pl2 ph11 d0 p11 pl2 ph12):f2
	(p5 pl1 ph2):f1
)
	d6
 	go=2 ph31
  MCWRK  wr #0
exit
ph1  = 0 0
ph2  = 0 0 0 0  1 1 1 1  2 2 2 2  3 3 3 3
ph11 = 0 0 2 2
ph12 = 0 2 2 0
ph31 = 0 2 0 2  2 0 2 0
"""

ZG_SOURCE = """\
;zg
;avance-version (12/01/11)
;1D sequence

#include <Avance.incl>

"acqt0=-p1*2/3.1416"

1 ze
2 30m
  d1
  p1 ph1
  go=2 ph31
  30m mc #0 to 2 F0(zd)
exit

ph1=0 2 2 0 1 3 3 1
ph31=0 2 2 0 1 3 3 1

;pl1 : f1 channel - power level for pulse (default)
;p1 : f1 channel -  high power pulse
;d1 : relaxation delay; 1-5 * T1
"""


def _kinds(prog, *skip):
    return [s.kind for s in prog.statements if s.kind not in skip]


def _find(prog, kind, **attrs):
    out = []
    for s in prog.walk():
        if s.kind == kind and all(getattr(s, k) == v for k, v in attrs.items()):
            out.append(s)
    return out


# ------------------------------------------------------------------- files
def test_find_pulseprogram_shipped_pdata_and_precomp(tmp_path):
    assert pp.find_pulseprogram(ZG) == ZG / "pulseprogram"
    assert pp.find_pulseprogram(ZG / "pdata" / "1") == ZG / "pulseprogram"
    assert pp.find_pulseprogram(ZG / "fid") == ZG / "pulseprogram"
    e = tmp_path / "7"
    e.mkdir()
    assert pp.find_pulseprogram(e) is None
    (e / "acqus").write_text("##TITLE= x\n")
    (e / "pulseprogram.precomp").write_text(";zg\n1 ze\nexit\n")
    assert pp.find_pulseprogram(e) == e / "pulseprogram.precomp"
    (e / "pulseprogram").write_text(";zg\n1 ze\nexit\n")
    assert pp.find_pulseprogram(e) == e / "pulseprogram"


def test_clean_text_drops_cpp_markers_and_keeps_the_rest():
    text = (ZG / "pulseprogram").read_text()
    assert text.startswith('# 1 "/opt/topspin')
    clean = pp.clean_text(text)
    assert not any(line.startswith("# ") for line in clean.splitlines())
    assert "1 ze" in clean and "  go=2 ph31" in clean
    assert ";p1 : f1 channel -  high power pulse" in clean
    assert "\n\n\n" not in clean


# ------------------------------------------------------------- shipped zg
def test_parse_shipped_zg():
    prog = pp.parse((ZG / "pulseprogram").read_text())
    assert prog.name == "zg"
    assert prog.includes == ["Avance.incl"]
    assert prog.definitions == {"acqt0": "-p1*2/3.1416", "MCREST": "30m - 30m",
                                "MCWRK": "0.333333*30m"}
    assert prog.defined_delays == ["MCWRK", "MCREST"]
    assert prog.meta["CLASS"].startswith("HighRes") and prog.dimension == "1D"
    labels = [s.name for s in prog.statements if s.kind == "label"]
    assert labels == ["1", "2", "LBLF0"]
    d1, = _find(prog, "delay", name="d1")
    assert d1.duration_expr == "d1" and d1.channel == ""
    p1, = _find(prog, "pulse")
    assert (p1.name, p1.phase, p1.channel, p1.power) == ("p1", "ph1", "f1", "")
    go, = _find(prog, "acquire")
    assert (go.name, go.phase, go.target_label) == ("go=2", "ph31", "2")
    lo, = _find(prog, "loop")
    assert (lo.target_label, lo.times_expr) == ("LBLF0", "td0")
    wr = [s for s in prog.statements if s.kind == "delay" and s.extras]
    assert [s.extras for s in wr] == [["wr #0"], ["zd"]]
    m2, = [s for s in prog.statements if s.duration_expr == "MCWRK*2"]
    assert m2.line_no == 46
    assert prog.phase_programs == {"ph1": [0, 2, 2, 0, 1, 3, 3, 1], "ph31": [0, 2, 2, 0, 1, 3, 3, 1]}
    assert prog.legend["p1"] == "f1 channel -  high power pulse"
    assert prog.legend["d1"] == "relaxation delay; 1-5 * T1"
    assert prog.legend["ns"].startswith("1 * n")
    assert prog.legend_for("PL1") == prog.legend["pl1"]
    assert not prog.warnings
    assert prog.channels() == ["f1"]
    assert prog.symbols_used()[:4] == ["MCWRK", "MCREST", "d1", "p1"]
    assert "ns" in prog.symbols_used() and "ph31" in prog.symbols_used()


def test_resolve_shipped_zg_takes_p1_and_d1_from_acqus():
    prog = pp.parse((ZG / "pulseprogram").read_text())
    vals = pp.resolve(prog, _acqus(ZG))
    p1, d1 = vals["p1"], vals["d1"]
    assert (p1.value, p1.unit, p1.text, p1.source) == (2.1, "µs", "2.1 µs", "acqus")
    assert p1.si == pytest.approx(2.1e-6)
    assert (d1.value, d1.unit, d1.text) == (5.0, "s", "5 s")
    assert vals["ns"].text == "512" and vals["ns"].unit == "count"
    assert vals["td0"].text == "1"
    assert vals["aq"].si == pytest.approx(2048 / (2 * 59523.8095238095))
    assert vals["aq"].text == "17.2 ms" and vals["aq"].source == "derived"
    mc = vals["MCWRK"]
    assert mc.source == "definition" and mc.text == "10 ms" and "0.333333*30m" in mc.note
    assert vals["MCREST"].text == "0 s"
    assert vals["acqt0"].si == pytest.approx(-2.1e-6 * 2 / 3.1416)
    assert vals["ph31"].text == "0 2 2 0 1 3 3 1" and vals["ph31"].note == "8 steps"
    assert vals["ph1"].unit == "× 90°"


def test_timeline_and_describe_shipped_zg():
    prog = pp.parse((ZG / "pulseprogram").read_text())
    vals = pp.resolve(prog, _acqus(ZG))
    tl = pp.timeline(prog, vals)
    assert tl.channels == ["f1", "acq"]
    kinds = [(e.kind, e.label) for e in tl.events]
    assert ("pulse", "p1") in kinds and ("acquire", "go=2") in kinds and ("delay", "d1") in kinds
    p1 = next(e for e in tl.events if e.label == "p1")
    assert (p1.channel, p1.sub, p1.value, p1.width) == ("f1", "ph1", "2.1 µs", pp.W_PULSE)
    d1 = next(e for e in tl.events if e.label == "d1")
    assert d1.value == "5 s" and d1.channel == "" and d1.end == p1.start
    go = next(e for e in tl.events if e.kind == "acquire")
    assert (go.channel, go.sub, go.value) == ("acq", "ph31", "aq 17.2 ms")
    assert go.start == p1.end and go.width == pp.W_ACQ
    ns = next(lp for lp in tl.loops if lp.times == "ns")
    assert ns.start == tl.labels["2"] and ns.end == go.end
    td0 = next(lp for lp in tl.loops if lp.times == "td0")
    assert td0.start == tl.labels["LBLF0"] and td0.end == tl.total
    assert {ns.level, td0.level} == {0, 1} and ns.level == 0    # shortest nearest the rows
    assert ns.text == "× ns = 512" and td0.text == "× td0 = 1"
    assert tl.notes == []
    assert any(e.note == "wr #0" for e in tl.events)
    assert pp.describe(prog, vals) == "zg: d1 – p1 (ph1) – acquire go=2 (ph31); phase cycle of 8 steps"


# ------------------------------------------------------- shipped mp3qdfsz
def test_parse_shipped_mp3qdfsz_four_pulses_one_shaped_nested_loops():
    prog = pp.parse((MQ / "pulseprogram").read_text())
    assert prog.name == "mp3qdfsz"
    assert prog.channels() == ["f1"]
    assert prog.dimension == "2D" and prog.aq_mode == "(F1) States"
    assert prog.prosol == "solids_mqmas"
    assert prog.meta["SUBTYPE"] == "MQMAS" and prog.meta["COMMENT"].startswith("4 pulse")
    assert prog.loopcounters == ["ST1CNT"] and prog.definitions["ST1CNT"] == "trunc(td1 / 2)"
    assert prog.definitions["p2"] == "1s/(cnst31*cnst0)"
    pulses = [s for s in prog.statements if s.kind in ("pulse", "shaped_pulse")]
    assert [(s.kind, s.name, s.power, s.phase) for s in pulses] == [
        ("pulse", "p1", "pl11", "ph1"), ("shaped_pulse", "p2", "sp1", "ph2"),
        ("pulse", "p3", "pl21", "ph3"), ("pulse", "p3", "", "ph4")]
    assert all(s.channel == "f1" for s in pulses)
    delays = [s.name for s in prog.statements if s.kind == "delay" and s.name not in ("MCWRK", "MCREST")]
    assert delays == ["d1", "d0", "d10", "d4"]
    go, = _find(prog, "acquire")
    assert (go.name, go.phase, go.target_label) == ("go=1", "ph31", "1")
    loops = [(s.target_label, s.times_expr) for s in prog.statements if s.kind == "loop"]
    assert loops == [("LBLSTS1", "2"), ("LBLF1", "ST1CNT"), ("LBLAV", "tdav")]
    housekeeping, = [s for s in prog.statements if "ip1" in s.extras]
    assert housekeeping.extras == ["wr #0", "if #0", "zd", "ip1"]
    two, = [s for s in prog.statements if "id10" in s.extras]
    assert two.name == "MCWRK"
    # phase programs: a (12) base, a {0}*24 repetition and a continued ph31
    assert prog.phase_base["ph1"] == 12 and prog.phase_programs["ph1"] == [0, 2, 4, 6, 8, 10]
    assert prog.phase_programs["ph2"] == [0] * 24 + [1] * 24 + [2] * 24 + [3] * 24
    assert prog.phase_programs["ph3"] == [0]
    assert len(prog.phase_programs["ph31"]) == 48
    assert prog.phase_programs["ph31"][24:30] == [2, 0, 2, 0, 2, 0]
    assert prog.phase_cycle_length() == 96
    assert prog.legend["p2"].startswith("=1s/(cnst31*l0)")
    assert prog.legend["spnam1"] == "dfs, file name of sweep"
    assert prog.legend["FnMODE"] == "States"
    assert prog.legend_for("sp1") == "power level for frequency sweep"    # documented as spw1
    assert prog.legend_for("plw11") == prog.legend["pl11"]
    assert prog.legend_for("pl7") == ""
    assert not prog.warnings


def test_resolve_shipped_mp3qdfsz_definitions_and_acqu2s():
    prog = pp.parse((MQ / "pulseprogram").read_text())
    vals = pp.resolve(prog, _acqus(MQ), _acqus(MQ, "acqu2s"))
    p2 = vals["p2"]
    assert p2.source == "definition" and p2.si == pytest.approx(1 / (26000 * 4))
    assert p2.text == "9.615 µs" and "= 1s/(cnst31*cnst0)" in p2.note
    assert "acqus has" not in p2.note                    # the file agrees with the definition
    assert vals["p1"].text == "3.1 µs" and vals["p3"].text == "4 µs"
    assert vals["d0"].text == "1 µs" and vals["d10"].text == "0.1 µs" and vals["d4"].text == "20 µs"
    assert vals["d1"].text == "500 ms"
    assert vals["pl11"].text == "250 W" and vals["pl11"].note == "-24.0 dB"
    assert vals["pl21"].text == "7.5 W" and vals["pl21"].unit == "W"
    assert vals["sp1"].text == "220 W" and "shape dfs" in vals["sp1"].note
    assert vals["spnam1"].text == "dfs" and vals["spnam1"].kind == "string"
    assert vals["cnst31"].text == "26000" and vals["cnst0"].text == "4"
    assert vals["td1"].text == "22" and vals["td1"].source == "acqu2s"
    assert vals["ST1CNT"].text == "11" and vals["ST1CNT"].unit == "count"
    assert vals["ns"].text == "3600" and vals["tdav"].text == "1"
    assert vals["acqt0"].text == "1 µs"                  # 1u*cnst11 with cnst11 = 1
    assert vals["ph1"].unit == "× 30°" and vals["ph1"].note == "6 steps"


def test_timeline_and_describe_shipped_mp3qdfsz():
    prog = pp.parse((MQ / "pulseprogram").read_text())
    vals = pp.resolve(prog, _acqus(MQ), _acqus(MQ, "acqu2s"))
    tl = pp.timeline(prog, vals)
    assert tl.channels == ["f1", "acq"]
    seq = [(e.kind, e.label) for e in tl.events if e.kind != "delay" or not e.label.startswith("MC")]
    assert seq == [("delay", "d1"), ("pulse", "p1"), ("delay", "d0"), ("shaped_pulse", "p2:sp1"),
                   ("delay", "d10"), ("pulse", "p3"), ("delay", "d4"), ("pulse", "p3"),
                   ("acquire", "go=1")]
    shaped = next(e for e in tl.events if e.kind == "shaped_pulse")
    assert shaped.width == pp.W_SHAPED and shaped.value == "9.615 µs" and shaped.sub == "ph2"
    assert shaped.note == "dfs · sp1 220 W"
    p1 = next(e for e in tl.events if e.label == "p1")
    assert p1.note == "pl11 250 W"
    d0 = next(e for e in tl.events if e.label == "d0")
    assert d0.note == "+in0 per increment"               # id0 rides on an MCWRK line
    d10 = next(e for e in tl.events if e.label == "d10")
    assert d10.note == "+in10 per increment"
    # the go loop and the three mc loops overlap without nesting: four levels
    assert sorted(lp.times for lp in tl.loops) == ["2", "ST1CNT", "ns", "tdav"]
    assert sorted(lp.level for lp in tl.loops) == [0, 1, 2, 3]
    outer = max(tl.loops, key=lambda lp: lp.level)
    assert outer.times == "tdav" and outer.start == 0.0 and outer.end == tl.total
    for e in tl.events:
        assert 0 <= e.start and e.end <= tl.total
    assert tl.notes == []
    text = pp.describe(prog, vals)
    assert text.startswith("mp3qdfsz: d1 – p1 (ph1) – d0 – p2:sp1 shaped (ph2) – d10 – p3 (ph3)"
                           " – d4 – p3 (ph4) – acquire go=1 (ph31)")
    assert "loops: ×2 (LBLSTS1), ×ST1CNT (LBLF1), ×tdav (LBLAV)" in text
    assert text.endswith("phase cycle of 96 steps")


# ------------------------------------------------------------- Hahn echo
def test_hahn_echo_definitions_inline_comment_and_rotor_sync():
    prog = pp.parse(HAHN)
    assert prog.name == "hahnecho.nmrfam"
    assert prog.definitions["d6"] == "((1s*l1)/cnst31)-(p1/2)-(p2/2)"
    assert prog.defined_delays == ["del7", "MCWRK", "MCREST"]
    assert prog.legend["del7"] == "calculates approximate d7 value"
    assert prog.legend["cnst2"].startswith("to select between solid echo")
    assert "Author" not in prog.legend                   # prose colons are not a legend
    del7, = _find(prog, "delay", name="del7")
    assert del7.raw.strip().startswith("del7 ;")       # the inline comment is kept in raw
    assert [s.name for s in prog.statements if s.kind == "pulse"] == ["p1", "p2"]
    p2, = _find(prog, "pulse", name="p2")
    assert p2.phase == "ph2" and p2.power == "" and p2.channel == "f1"
    vals = pp.resolve(prog, HAHN_ACQUS)
    assert vals["p2"].text == "6.5 µs" and vals["p2"].source == "definition"
    assert vals["d6"].si == pytest.approx(1 / 35714 - 3.25e-6 / 2 - 6.5e-6 / 2)
    assert vals["d6"].text == "23.13 µs"
    assert vals["del7"].si == pytest.approx(1 / 35714 - 6.5e-6 / 2)
    assert vals["d7"].text == "24.75 µs" and vals["d7"].source == "acqus"
    assert vals["l1"].text == "1" and vals["cnst31"].text == "35714"
    assert vals["acqt0"].text == "0 s"
    tl = pp.timeline(prog, vals)
    order = [e.label for e in tl.events if e.kind in ("pulse", "acquire") or e.label in ("d6", "d7")]
    assert order == ["p1", "d6", "p2", "d7", "go=2"]
    assert pp.describe(prog, vals) == ("hahnecho.nmrfam: del7 – d1 – p1 (ph1) – d6 – p2 (ph2) – d7"
                                       " – acquire go=2 (ph31); phase cycle of 16 steps")


def test_hahn_echo_plain_source_with_ifdef_branches():
    src = HAHN.replace("  (p2 ph2):f1\n", "#ifdef dec\n  (p2 ph2):f1 cpds2:f2\n#else\n  (p2 ph2):f1\n#endif\n")
    prog = pp.parse(src)
    p2s = _find(prog, "pulse", name="p2")
    assert [s.cond for s in p2s] == ["dec", "!dec"]
    dec, = _find(prog, "decouple_on")
    assert dec.cond == "dec" and dec.channel == "f2" and dec.name == "cpds2"


# ---------------------------------------------------- saturation recovery
def test_saturation_recovery_include_goto_print_vd_and_local_loop():
    prog = pp.parse(SATREC)
    assert prog.includes == ["aq_prot.incl"]
    assert prog.dimension == "2D" and prog.aq_mode == "(F1) QF"
    labels = [s.name for s in prog.statements if s.kind == "label"]
    assert labels == ["1", "LBLAV", "2", "LBLF1", "Passaq", "3", "HaltAcqu"]
    cond, plain = _find(prog, "goto")
    assert (cond.name, cond.condition, cond.target_label) == ("if", "aq < 50.1m", "Passaq")
    assert (plain.name, plain.target_label) == ("goto", "HaltAcqu")
    pr, = _find(prog, "misc", name="print")
    assert pr.extras == ["acquisition time exceeds 50m limit!"]
    vd, = _find(prog, "delay", name="vd")
    assert vd.line_no and "recovery delay" not in vd.name
    sat, = [s for s in prog.statements if s.kind == "loop" and s.target_label == "3"]
    assert sat.times_expr == "l20"
    hk, = [s for s in prog.statements if "ivd" in s.extras]
    assert hk.extras == ["wr #0", "if #0", "zd", "ivd"]
    assert [(s.target_label, s.times_expr) for s in prog.statements if s.kind == "loop"] == [
        ("3", "l20"), ("LBLF1", "td1"), ("LBLAV", "tdav")]
    vals = pp.resolve(prog, SATREC_ACQUS, {"TD": 11})
    assert vals["vd"].kind == "list" and vals["vd"].text == "VDLIST"
    assert vals["vd"].note == "delay list T1_standard_array"
    assert vals["l20"].text == "200" and vals["d20"].text == "1 ms" and vals["td1"].text == "11"
    assert vals["acqt0"].text == "-0.25 µs"
    tl = pp.timeline(prog, vals)
    lit = [e for e in tl.events if e.label == "1m"]
    assert len(lit) == 3 and all(e.value == "1 ms" and e.width == pp.W_SHORT for e in lit)
    vde = next(e for e in tl.events if e.label == "vd")
    assert vde.value == "VDLIST" and vde.note == "next list entry per slice"
    l20 = next(lp for lp in tl.loops if lp.times == "l20")
    assert l20.start == tl.labels["3"] and l20.level == 0
    assert any("aq < 50.1m" in n for n in tl.notes)
    assert pp.describe(prog, vals) == (
        "satrect1: d1 – [d20 – p1 (ph4)] ×l20 – vd – p1 (ph1) – acquire go=2 (ph31);"
        " loops: ×td1 (LBLF1), ×tdav (LBLAV); phase cycle of 8 steps")


# ------------------------------------------------------------------ QCPMG
def test_qcpmg_explicit_acquisition_echo_loop_and_decoupling_bar():
    prog = pp.parse(QCPMG)
    assert prog.name == "qcpmg"
    assert [(s.name, s.channel) for s in _find(prog, "set_power")] == [("pl1", "f1"), ("pl12", "f2")]
    adc, = _find(prog, "acquire")
    assert (adc.name, adc.phase, adc.target_label) == ("adc", "ph31", "")
    rc, = [s for s in prog.statements if s.kind == "loop" and s.name.startswith("rcyc")]
    assert (rc.target_label, rc.times_expr) == ("2", "ns")
    echo, = [s for s in prog.statements if s.kind == "loop" and s.name == "lo"]
    assert (echo.target_label, echo.times_expr) == ("3", "l22")
    on, = _find(prog, "decouple_on")
    off, = _find(prog, "decouple_off")
    assert (on.name, on.channel, off.channel) == ("cpds2", "f2", "f2")
    assert [s.extras for s in prog.statements if "eoscnp" in s.extras] == [["eoscnp"]]
    assert "syrec" in [s.name for s in prog.statements if s.kind == "misc"]
    syms = prog.symbols_used()
    assert {"cpdprg2", "pcpd2", "l22", "d3", "d6", "pl12"} <= set(syms)
    vals = pp.resolve(prog, QCPMG_ACQUS)
    assert vals["d3"].si == pytest.approx(200e-6 / 2 - 4e-6 / 2) and vals["d3"].source == "definition"
    assert vals["cpdprg2"].text == "spinal64" and vals["pcpd2"].text == "6.8 µs"
    assert vals["pl12"].text == "20 W" and vals["l22"].text == "64"
    tl = pp.timeline(prog, vals)
    assert tl.channels == ["f1", "f2", "acq"]
    dec = next(e for e in tl.events if e.kind == "decouple")
    acq = next(e for e in tl.events if e.kind == "acquire")
    assert dec.channel == "f2" and dec.label == "cpds2 · spinal64" and dec.value == "pcpd2 6.8 µs"
    assert dec.start < acq.start and dec.end > acq.end       # runs through the train
    p2s = [e for e in tl.events if e.label == "p2"]
    assert len(p2s) == 2 and all(e.channel == "f1" and e.sub == "ph2" for e in p2s)
    l22 = next(lp for lp in tl.loops if lp.times == "l22")
    assert l22.start == tl.labels["3"] and l22.end == p2s[1].end + pp.W_DELAY
    ns = next(lp for lp in tl.loops if lp.times == "ns")
    assert ns.start == tl.labels["2"] and ns.level > l22.level
    text = pp.describe(prog, vals)
    assert "[p2 (ph2) – d6] ×l22" in text
    assert "cpds2 on f2" in text and "acquire adc (ph31)" in text and "do f2" in text
    assert text.endswith("phase cycle of 8 steps")


# ------------------------------------------------------------------ CP/MAS
def test_cpmas_simultaneous_contact_pulses_and_mc_loop():
    prog = pp.parse(CPMAS)
    assert prog.channels() == ["f1", "f2"]
    sim, = _find(prog, "group", name="sim")
    assert [c.channel for c in sim.children] == ["f1", "f2"]
    f1, f2 = sim.children
    assert [(k.kind, k.name, k.power, k.phase) for k in f1.children] == [("pulse", "p15", "pl1", "ph2")]
    assert [(k.kind, k.name, k.power, k.phase) for k in f2.children] == [("shaped_pulse", "p15", "sp0", "ph10")]
    p3, = _find(prog, "pulse", name="p3")
    assert (p3.channel, p3.power, p3.phase) == ("f2", "pl2", "ph1")
    offs = _find(prog, "decouple_off")
    assert [s.line_no for s in offs] == [12, 17]         # `d1 do:f2` and `1m do:f2`
    mc, = [s for s in prog.statements if s.kind == "loop" and s.name == "mc"]
    assert (mc.target_label, mc.times_expr) == ("2", "td0")
    acqus = {"P": [0, 0, 0, 2.5] + [0] * 11 + [2000.0] + [0] * 48, "D": [0, 3.0] + [0] * 62,
             "PLW": [0, 80.0, 90.0] + [0] * 9 + [70.0] + [0] * 51, "SPW": [55.0] + [0] * 63,
             "SPNAM": ["ramp.100"] + [""] * 63, "PCPD": [100, 100, 6.5] + [100] * 7,
             "CPDPRG": ["", "", "spinal64"] + [""] * 6, "NS": 1024, "TD": 2048, "SW_h": 50000.0,
             "TD0": 1, "NUC1": "13C", "NUC2": "1H"}
    vals = pp.resolve(prog, acqus)
    assert vals["p15"].text == "2 ms" and vals["sp0"].text == "55 W" and "ramp.100" in vals["sp0"].note
    tl = pp.timeline(prog, vals)
    assert tl.channels == ["f1", "f2", "acq"]
    contact = [e for e in tl.events if e.label.startswith("p15")]
    assert len(contact) == 2 and contact[0].start == contact[1].start
    assert {e.kind for e in contact} == {"pulse", "shaped_pulse"}
    shaped = next(e for e in contact if e.kind == "shaped_pulse")
    assert shaped.label == "p15:sp0" and shaped.note == "ramp.100 · sp0 55 W"
    nxt = min(e.start for e in tl.events if e.start > contact[0].start)
    assert nxt == contact[0].start + pp.W_SHAPED          # the wider member sets the group width
    dec = next(e for e in tl.events if e.kind == "decouple")
    acq = next(e for e in tl.events if e.kind == "acquire")
    assert dec.start <= acq.start and dec.end >= acq.end
    assert {lp.times for lp in tl.loops} == {"ns", "td0"}
    text = pp.describe(prog, vals)
    assert "p3 (ph1) f2" in text
    assert "[f1: p15 (ph2) | f2: p15:sp0 shaped (ph10)]" in text
    assert "cpds2 on f2 – acquire go=2 (ph31) – do f2" in text


def test_decoupling_during_acquisition_only_and_never_closed():
    src = "1 ze\n2 d1\n  (p1 ph1):f1\n  go=2 ph31 cpd2:f2\n  30m mc #0 to 2 F0(zd)\nexit\nph1=0 2\nph31=0 2\n"
    prog = pp.parse(src)
    tl = pp.timeline(prog, {})
    dec = next(e for e in tl.events if e.kind == "decouple")
    acq = next(e for e in tl.events if e.kind == "acquire")
    assert dec.start == acq.start and dec.end == tl.total
    assert any("never switched off" in n for n in tl.notes)


# ------------------------------------------------------------------- HMQC
def test_hmqc_centered_composite_across_channels():
    prog = pp.parse(HMQC)
    assert prog.name == "pf_hmqc"
    assert prog.definitions["l0"] == "1" and prog.defined_delays == ["zfilter"]
    centre, = _find(prog, "group", name="center")
    assert centre.line_no == 15 and [c.channel for c in centre.children] == ["f2", "f1"]
    assert centre.raw.startswith("(center (p11") and centre.raw.endswith(")")
    f2, f1 = centre.children
    assert [(k.kind, k.name, k.phase) for k in f2.children] == [
        ("pulse", "p11", "ph11"), ("delay", "d0", ""), ("pulse", "p11", "ph12")]
    assert f2.children[1].channel == "f2"
    assert [(k.kind, k.name, k.phase, k.power) for k in f1.children] == [("pulse", "p5", "ph2", "pl1")]
    assert [(s.name, s.channel) for s in _find(prog, "set_power")] == [("pl1", "f1"), ("pl2", "f2")]
    assert prog.phase_programs["ph2"] == [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]
    assert prog.phase_programs["ph31"] == [0, 2, 0, 2, 2, 0, 2, 0]
    acqus = {"P": [0, 7.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.5] + [0] * 52,
             "D": [0, 1.8] + [0] * 62, "PLW": [0, 1.34, 60.5] + [0] * 61,
             "L": [1, 1, 1, 1, 1, 1, 333] + [1] * 9 + [33] + [1] * 15,
             "CNST": [1] * 31 + [33333.0] + [1] * 32, "NS": 4096, "TD": 1024, "SW_h": 200000.0}
    vals = pp.resolve(prog, acqus)
    assert vals["p5"].text == "15 µs" and vals["d0"].text == "15 µs" and vals["d0"].source == "definition"
    assert vals["d6"].si == pytest.approx(333 / 33333 - 7.5e-6 / 2 - 15e-6 / 2 - 2.5e-6)
    assert "zfilter" not in vals                          # defined but unused by this body
    tl = pp.timeline(prog, vals)
    assert tl.channels == ["f1", "f2", "acq"]
    p11 = [e for e in tl.events if e.label == "p11"]
    d0 = next(e for e in tl.events if e.label == "d0")
    p5 = next(e for e in tl.events if e.label == "p5")
    assert len(p11) == 2 and d0.channel == "f2" and p11[0].end == d0.start and d0.end == p11[1].start
    f2_start, f2_end = p11[0].start, p11[1].end
    assert p5.channel == "f1" and p5.start == pytest.approx((f2_start + f2_end) / 2 - pp.W_PULSE / 2)
    after = next(e for e in tl.events if e.label == "d6" and e.start >= f2_end)
    assert after.start == f2_end
    text = pp.describe(prog, vals)
    assert "centre [f2: p11 (ph11) – d0 – p11 (ph12) | f1: p5 (ph2)]" in text


# ------------------------------------------------------- uncompiled source
def test_plain_zg_source_include_and_mc_macro():
    prog = pp.parse(ZG_SOURCE)
    assert prog.name == "zg" and prog.includes == ["Avance.incl"]
    assert prog.legend["d1"] == "relaxation delay; 1-5 * T1"
    st = [(s.kind, s.name) for s in prog.statements]
    assert st == [("label", "1"), ("misc", "ze"), ("label", "2"), ("delay", "30m"), ("delay", "d1"),
                  ("pulse", "p1"), ("acquire", "go=2"), ("delay", "30m"), ("loop", "mc"), ("misc", "exit")]
    mc = prog.statements[-2]
    assert (mc.target_label, mc.times_expr) == ("2", "td0") and mc.extras == ["mc #0 to 2 F0(zd)"]
    assert prog.statements[-3].extras == ["mc #0 to 2 F0(zd)"]
    p1, = _find(prog, "pulse")
    assert (p1.channel, p1.phase) == ("f1", "ph1")
    tl = pp.timeline(prog, pp.resolve(prog, {"P": [0, 2.5], "D": [0, 0.1], "NS": 8, "TD": 1024, "SW_h": 1e5}))
    assert next(e for e in tl.events if e.label == "30m").value == "30 ms"
    assert {lp.times for lp in tl.loops} == {"ns", "td0"}
    assert pp.describe(prog) == "zg: d1 – p1 (ph1) – acquire go=2 (ph31); phase cycle of 8 steps"


def test_mc_macro_with_f1_spec_loops_over_td1():
    prog = pp.parse("1 ze\n2 d1\n  p1 ph1\n  d0\n  p1 ph2\n  go=2 ph31\n  d11 mc #0 to 2 F1PH(ip1, id0)\nexit\n")
    mc, = [s for s in prog.statements if s.kind == "loop"]
    assert (mc.target_label, mc.times_expr) == ("2", "td1")


# ---------------------------------------------------------- phase programs
def test_phase_program_syntax_variants():
    prog = pp.parse("1 ze\n  p1 ph1\n  go=2 ph31\nexit\n"
                    "ph1=(360) 0 90 180 270\n"
                    "ph2 = {0 2}^1^2^3\n"
                    "ph3= {0}*4 {2}*4\n"
                    "ph31=0 2 2 0\n"
                    "     1 3 3 1 ; continued\n"
                    "ph4=0\n")
    assert prog.phase_base["ph1"] == 360 and prog.phase_programs["ph1"] == [0, 90, 180, 270]
    assert prog.phase_programs["ph2"] == [0, 2, 1, 3, 2, 0, 3, 1]
    assert prog.phase_programs["ph3"] == [0, 0, 0, 0, 2, 2, 2, 2]
    assert prog.phase_programs["ph31"] == [0, 2, 2, 0, 1, 3, 3, 1]
    assert prog.phase_programs["ph4"] == [0]
    vals = pp.resolve(prog, {})
    assert vals["ph1"].unit == "× 1°" and vals["ph31"].note == "8 steps"


# ------------------------------------------------------------- odd syntax
def test_gradients_variable_counters_multipliers_and_bare_shaped_pulse():
    src = ("1 ze\n2 d1\n  p16:gp1\n  d16 UNBLKGRAD\n  p1*0.5 ph1\n  p11:sp1 ph2\n  d0*2\n"
           "3 (p2 ph2):f1\n  d20\n  lo to 3 times vc\n  go=2 ph31\n  d11 wr #0 ivc\nexit\n")
    prog = pp.parse(src)
    g, = _find(prog, "pulse", name="p16")
    assert (g.channel, g.power) == ("grad", "gp1")
    assert _find(prog, "delay", name="d16")[0].extras == ["UNBLKGRAD"]
    half, = _find(prog, "pulse", name="p1")
    assert half.duration_expr == "p1*0.5" and half.phase == "ph1"
    sh, = _find(prog, "shaped_pulse")
    assert (sh.name, sh.power, sh.phase, sh.channel) == ("p11", "sp1", "ph2", "f1")
    d0, = _find(prog, "delay", name="d0")
    assert d0.duration_expr == "d0*2"
    lo, = _find(prog, "loop", name="lo")
    assert lo.times_expr == "vc"
    assert not prog.warnings
    vals = pp.resolve(prog, {"P": [0, 4.0] + [0] * 62, "D": [10e-6] + [0] * 63, "GPZ": [0, 50.0]})
    assert vals["gp1"].text == "50 %" and vals["vc"].kind == "list"
    tl = pp.timeline(prog, vals)
    assert "grad" in tl.channels and tl.channels.index("grad") == len(tl.channels) - 2
    assert next(e for e in tl.events if e.label == "p1 ×0.5").value == "2 µs"
    assert next(e for e in tl.events if e.label == "d0 ×2").value == "20 µs"


def test_parser_never_raises_on_garbage_and_reports_unbalanced_parens():
    assert pp.parse("").statements == []
    prog = pp.parse("1 ze\n  (p1 ph1\n  go=2 ph31\nexit\n")
    assert prog.warnings and "unbalanced" in prog.warnings[0]
    prog = pp.parse("1 ze\n  ??? !!! )))\n  go=2 ph31\nexit\n")
    assert any(s.kind == "misc" for s in prog.statements)
    assert _find(prog, "acquire")


def test_unresolved_duration_is_reported_in_the_timeline_notes():
    prog = pp.parse('"d6=1s*l1/cnst31-p1"\n1 ze\n  d1\n  p1 ph1\n  d6\n  go=2 ph31\nexit\nph1=0\nph31=0\n')
    vals = pp.resolve(prog, {"P": [0, 2.0], "D": [0, 1.0], "L": [1, 4]})    # no CNST
    assert vals["d6"].text == "" and "needs cnst31" in vals["d6"].note
    tl = pp.timeline(prog, vals)
    assert next(e for e in tl.events if e.label == "d6").value == ""
    assert any(n.startswith("d6: duration not resolved") and "cnst31" in n for n in tl.notes)


# --------------------------------------------------------------- evaluator
def test_evaluate_arithmetic_units_functions_and_safety():
    env = {"p1": 2.5e-6, "l1": 4, "cnst31": 20000.0, "td1": 22}
    assert pp.evaluate("1s*l1/cnst31", env) == pytest.approx(2e-4)
    assert pp.evaluate("((1s*l1)/cnst31)-(p1/2)", env) == pytest.approx(2e-4 - 1.25e-6)
    assert pp.evaluate("trunc(td1 / 2)", env) == 11
    assert pp.evaluate("-p1*2/3.1416", env) == pytest.approx(-5e-6 / 3.1416)
    assert pp.evaluate("10u", {}) == pytest.approx(1e-5)
    assert pp.evaluate("30m - 30m", {}) == 0
    assert pp.evaluate("0.333333*30m", {}) == pytest.approx(0.00999999)
    assert pp.evaluate("2.5e-3*2", {}) == pytest.approx(5e-3)
    assert pp.evaluate("1s/(cnst31*cnst0)", {"cnst31": 26000, "cnst0": 4}) == pytest.approx(1 / 104000)
    assert pp.evaluate("p1/cnst9", env) is None                  # unknown symbol
    assert pp.evaluate("p1/(l1-4)", env) is None                 # division by zero
    assert pp.evaluate("1s*(", env) is None                      # malformed
    assert pp.evaluate("", env) is None
    assert pp.evaluate("__import__('os').system('x')", env) is None
    assert pp.evaluate("p1 p1", env) is None


def test_value_formatting_and_categories():
    v = pp._make("d1", 3.45e-3, "time", "acqus")
    assert (v.value, v.unit, v.text) == (pytest.approx(3.45), "ms", "3.45 ms")
    assert pp._make("p1", 9.615384615e-6, "time", "acqus").text == "9.615 µs"
    assert pp._make("d1", 3580.0, "time", "acqus").text == "3580 s"
    assert pp._make("pl1", 0.0, "power", "acqus").text == "0 W"
    assert pp._make("pl1", 5.3, "power", "acqus").note == "-7.2 dB"
    assert [pp.category(s) for s in ("p1", "d1", "pl11", "sp1", "spnam1", "l20", "ns", "td1",
                                     "cnst31", "in0", "cpdprg2", "ph31", "gp1", "MCWRK")] == [
        "pulse", "delay", "power", "power", "shape", "loop", "loop", "loop", "constant",
        "increment", "decoupling", "phase", "gradient", "other"]


def test_old_topspin_db_only_power_levels():
    prog = pp.parse("1 ze\n2 d1 pl1:f1\n  p1 ph1\n  go=2 ph31\nexit\nph1=0\nph31=0\n")
    vals = pp.resolve(prog, {"P": [0, 3.0], "D": [0, 2.0], "PL": [120, -3.5]})
    assert vals["pl1"].text == "-3.5 dB" and vals["pl1"].unit == "dB" and vals["pl1"].si is None


# --------------------------------------------------------------------- load
def test_load_shipped_expno_and_channel_nuclei(tmp_path):
    rec = pp.load(ZG)
    assert rec.path == ZG / "pulseprogram" and rec.pulprog == "zg"
    assert rec.program is not None and rec.program.name == "zg"
    assert rec.values["p1"].text == "2.1 µs"
    assert rec.nuclei == {"f1": "11B", "f2": "1H"}
    assert rec.warnings == []
    rec2 = pp.load(MQ / "pdata" / "1")
    assert rec2.expno == MQ and rec2.values["td1"].text == "22"
    assert pp.channel_nuclei({"NUC1": "<27Al>", "NUC2": "1H", "NUC3": "off"}) == {"f1": "27Al", "f2": "1H"}
    assert pp.channel_nuclei(None) == {}
    e = tmp_path / "9"
    e.mkdir()
    (e / "acqus").write_text("##TITLE= Parameter file\n##$PULPROG= <zg>\n##$NUC1= <7Li>\n##END=\n")
    bare = pp.load(e)
    assert bare.program is None and bare.path is None and bare.pulprog == "zg"
    assert bare.warnings == ["no pulse program file in this EXPNO"]
    assert bare.nuclei == {"f1": "7Li"}
