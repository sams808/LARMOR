"""larmor.masrate: the three-source MAS rate (acqus / title / NMRFAM booking
sidecar), the bounded title parse, the majority resolver, the rotor / session
key and the per-session confirmation store the loader applies."""
import inspect
import json
import os
import subprocess
import sys

import pytest

from conftest import BRUKER_1R, require
from larmor import masrate as M
from test_referencing import _expno

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _addenda(expno, rate=None, spinning="true", extra=""):
    """<EXPNO>/experiment_addenda.xml in the 3104 shape (rate in kHz)."""
    rate_tag = ("" if rate is None else
                f"\t<magic_angle_spinning_rate>{rate}</magic_angle_spinning_rate>\n")
    text = ("<experiment_addenda>\n\t<nan_user_id>ssoudani</nan_user_id>\n"
            "\t<session_start_epoch>1777476884</session_start_epoch>\n"
            f"\t<magic_angle_spinning>{spinning}</magic_angle_spinning>\n"
            f"{rate_tag}{extra}"
            "\t<multi_receiver>false</multi_receiver>\n"
            "\t<time_shared>false</time_shared>\n\t<state>solid</state>\n"
            "</experiment_addenda>\n")
    (expno / M.ADDENDA_FILE).write_text(text, encoding="utf-8")


def _ev(acqus, title, booking, **kw):
    return M.MasEvidence(acqus_Hz=acqus, title_Hz=title, booking_Hz=booking, **kw)


@pytest.fixture()
def store(tmp_path, monkeypatch):
    p = tmp_path / "mas_confirmations.jsonl"
    monkeypatch.setenv("LARMOR_MAS_LOG", str(p))
    return p


# ------------------------------------------------------------------ parsing
def test_parse_title_rate_matrix():
    def hz(t):
        return M.parse_title_rate(t).hz

    r = M.parse_title_rate("Rotor SR31648\nMASR 35.714 kHz\nVT ambient")
    assert (r.hz, r.static_zero, r.note, r.repaired) == (35714.0, False, "", False)
    assert r.raw == "MASR 35.714 kHz"
    assert hz("MASR 20kHz") == 20000.0
    assert hz("sample X MAS 12.5 kHz") == 12500.0
    assert hz("MASR = 20 kHz") == 20000.0
    r = M.parse_title_rate("MASR 35714 Hz")
    assert r.hz == 35714.0 and r.repaired is False
    r = M.parse_title_rate("MASR 0 kHz")
    assert r.hz is None and r.static_zero is True
    # the unit typos of the tree: re-read in the other unit and marked
    r = M.parse_title_rate("MASR 35741 kHz")
    assert r.hz == 35741.0 and r.repaired is True
    assert "35741" in r.note and "typo" in r.note
    r = M.parse_title_rate("MASR 35.714 Hz")
    assert r.hz == pytest.approx(35714.0) and r.repaired is True
    # implausible both ways: dropped with a note
    r = M.parse_title_rate("MASR 900 kHz")
    assert r.hz is None and "ignored" in r.note and r.static_zero is False
    # no unit, no prefix, no rate
    assert hz("MAS 20") is None
    assert hz("81Br static") is None
    assert hz("electrostatic sample holder") is None
    assert hz("1 kHz spikelet spacing") is None
    assert hz("") is None and hz(None) is None


def test_read_booking_rate_variants(tmp_path):
    e = tmp_path / "3104"
    e.mkdir()
    _addenda(e, 22)
    assert M.read_booking_rate(e) == (22000.0, True)
    _addenda(e, "35.714")
    assert M.read_booking_rate(e) == (pytest.approx(35714.0), True)
    # the booking form left blank: the flag is informational, never static
    # evidence and never a rate (0 of 214 such files carry one)
    _addenda(e, None, spinning="false")
    assert M.read_booking_rate(e) == (None, False)
    assert M.read_booking_rate(tmp_path / "missing") == (None, None)
    (e / M.ADDENDA_FILE).write_bytes(b"<experiment_addenda><magic_angle")
    assert M.read_booking_rate(e) == (None, None)
    _addenda(e, "")
    assert M.read_booking_rate(e) == (None, True)
    _addenda(e, "n/a")
    assert M.read_booking_rate(e) == (None, True)
    # a Hz value typed into the kHz field is implausible: dropped, flag kept
    _addenda(e, 22000)
    assert M.read_booking_rate(e) == (None, True)


# ----------------------------------------------------------------- resolver
def test_resolve_majority_matrix():
    def r(ev, **kw):
        res = M.resolve(ev, **kw)
        return res.rate_Hz, res.uncertain, res.source

    assert r(_ev(4200, 35714, 35714)) == (35714.0, False, "title+booking")
    assert r(_ev(4200, 20000, 22000)) == (22000.0, True, "highest")
    assert r(_ev(20000, 22000, 20000)) == (20000.0, False, "acqus+booking")
    assert r(_ev(5000, 20000, 10000)) == (20000.0, True, "highest")
    assert r(_ev(None, None, 22000)) == (22000.0, False, "booking")
    assert r(_ev(4200, None, None)) == (4200.0, False, "acqus")
    assert r(_ev(4200, 35714, None)) == (35714.0, True, "highest")
    # all three within 2 %: the (human) title supplies the number
    assert r(_ev(20000, 20300, 20100)) == (20300.0, False, "all")
    assert r(_ev(None, None, None)) == (M.MAS_FALLBACK_HZ, True, "fallback")
    # a repaired title never supplies the number inside an agreeing pair
    rep = _ev(4200, 35741, 35714, title_repaired=True,
              title_note='title says "35741 kHz", read as 35741 Hz (unit typo?)')
    assert r(rep) == (35714.0, False, "title+booking")
    # ... and never settles the rate alone
    assert r(_ev(4200, 35714, None, title_repaired=True)) == (35714.0, True, "highest")
    assert r(_ev(None, 35714, None, title_repaired=True)) == (35714.0, True, "title")
    # a title that was present but unreadable flags whatever else is found
    assert r(_ev(4200, None, None, title_note="ignored")) == (4200.0, True, "acqus")
    # the static matrix bruker's thirteen cases pin
    assert r(_ev(None, None, None), masr_zero=True) == (0.0, True, "static")
    assert r(_ev(None, None, None), masr_zero=True, static_title=True) == (0.0, False, "static")
    assert r(_ev(None, None, None), static_title=True) == (0.0, False, "static")
    assert r(_ev(None, None, None), static_pp=True) == (0.0, True, "static")
    assert r(_ev(14000, None, None), static_pp=True) == (14000.0, True, "acqus")
    assert r(_ev(None, 20000, None), masr_zero=True) == (20000.0, True, "title")
    assert r(_ev(None, 20000, 20000), masr_zero=True) == (20000.0, True, "title+booking")
    # the notes name the outvoted source and the disagreement
    assert "outvoted" in M.resolve(_ev(4200, 35714, 35714)).note
    assert "disagree" in M.resolve(_ev(4200, 20000, 22000)).note


def test_rotor_id_title_then_folder_then_name():
    assert M.rotor_id("Rotor RS2427418\nMASR 20 kHz", "04272026_P5-Bi8-12_SS_ALP") == "RS2427418"
    assert M.rotor_id("11B zg", "01202026_SR31648_Base3Ca_SS_ALP") == "SR31648"
    assert M.rotor_id("Rotor SR31602LS - NaCl", "x") == "SR31602LS"
    assert M.rotor_id("Rotor x", "03232026_P1-Bi0_SS_ALP") == "03232026_P1-Bi0_SS_ALP"
    assert M.rotor_id("Rotor RS", "03232026_P1-Bi0_SS_ALP") == "03232026_P1-Bi0_SS_ALP"
    assert M.rotor_id("19F Hahn echo with 1 rotor period", "Kogarkoite") == "Kogarkoite"
    assert M.rotor_id("", "") == ""


def test_confirmation_key_and_describe(tmp_path):
    root = tmp_path / "2026-05"
    e = _expno(root, "04272026_P5-Bi8-12_SS_ALP", 3104, "31P", 242.792156,
               242.792156, "1777476884", title="Rotor RS2427418\nMASR 20 kHz")
    key = M.confirmation_key(e, "31P", "Rotor RS2427418\nMASR 20 kHz")
    assert key == (str(root), "RS2427418", "31P")
    assert M.describe_key(key) == "2026-05 · rotor RS2427418 · 31P"
    # a 1r file or a pdata folder resolve to the same key
    assert M.confirmation_key(e / "pdata" / "1" / "procs", "<31P>", "Rotor RS2427418") == key
    assert M.confirmation_key(e / "pdata" / "1", "31P", "Rotor RS2427418") == key
    assert M.confirmation_key(e, "", "Rotor RS2427418") is None
    assert M.confirmation_key(tmp_path / "nowhere" / "1", "31P", "") is None
    assert M.confirmation_key("", "31P", "") is None
    # no rotor line: the sample folder stands in and is named as such
    key2 = M.confirmation_key(e, "31P", "no rotor here")
    assert key2[1] == "04272026_P5-Bi8-12_SS_ALP"
    assert M.describe_key(key2) == "2026-05 · folder 04272026_P5-Bi8-12_SS_ALP · 31P"


# -------------------------------------------------------------------- store
def test_store_requires_key_and_exact_evidence(store):
    key = (r"C:\data\2026-05", "RS2427418", "31P")
    ev = _ev(4200, 20000, 22000)
    assert M.lookup(key, ev) is None
    assert M.remember(key, ev, 22000) == store
    assert M.lookup(key, ev)["rate_Hz"] == 22000.0
    # the 05082026 booking of the same rotor: 20 kHz, a different triplet
    assert M.lookup(key, _ev(4200, 20000, 20000)) is None
    assert M.lookup(key, _ev(4200, 20000, None)) is None
    assert M.lookup((key[0], "RS40339", "31P"), ev) is None
    assert M.lookup((key[0], "RS2427418", "27Al"), ev) is None
    assert M.lookup((r"c:\DATA\2026-05", "rs2427418", "31P"), ev)["rate_Hz"] == 22000.0
    assert M.lookup(None, ev) is None
    # the triplet is matched at 1 Hz: a float read of the same values matches
    assert M.lookup(key, _ev(4200.0, 20000.0, 22000.0))["rate_Hz"] == 22000.0
    assert M.lookup(key, _ev(4200.4, 20000.0, 22000.0))["rate_Hz"] == 22000.0
    assert M.lookup(key, _ev(4201, 20000, 22000)) is None


def test_latest_wins_forget_clears_append_only(store):
    key = (r"C:\data\2026-05", "RS2427418", "31P")
    ev = _ev(4200, 20000, 22000)
    M.remember(key, ev, 22000, note="Experiment parameters")
    M.remember(key, ev, 20000)
    assert M.lookup(key, ev)["rate_Hz"] == 20000.0
    M.forget(key, ev)
    assert M.lookup(key, ev) is None
    lines = store.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 3
    first = json.loads(lines[0])
    assert first["rate_Hz"] == 22000.0 and first["action"] == "remember"
    assert first["note"] == "Experiment parameters" and first["larmor"]
    assert json.loads(lines[2])["action"] == "forget"
    # a corrupt line is skipped, a later remember wins again
    with open(store, "a", encoding="utf-8") as f:
        f.write("{not json\n")
    M.remember(key, ev, 22000)
    assert M.lookup(key, ev)["rate_Hz"] == 22000.0


def test_resolve_for_load_provenance_and_confirmation(tmp_path, store):
    from larmor.io import bruker

    root = tmp_path / "2026-05"
    title = "Rotor RS2427418\nMASR 20 kHz"
    e = _expno(root, "04272026_P5-Bi8-12_SS_ALP", 3104, "31P", 242.792156,
               242.792156, "1777476884", title=title)
    _addenda(e, 22)
    meta = bruker._meta_1d({"MASR": 4200, "NUC1": "<31P>"}, title, e)
    r = M.resolve_for_load(meta, e)
    assert r["spin_rate_Hz"] == 22000.0 and r["mas_uncertain"] is True
    assert r["note"] == ""
    b = r["provenance"]
    assert set(b) == {"acqus_Hz", "title_Hz", "title_note", "title_raw",
                      "booking_Hz", "booking_flag", "source", "uncertain",
                      "session", "rotor", "nucleus", "confirmed", "confirmed_note"}
    assert (b["acqus_Hz"], b["title_Hz"], b["booking_Hz"]) == (4200.0, 20000.0, 22000.0)
    assert b["title_raw"] == "MASR 20 kHz" and b["booking_flag"] is True
    assert b["source"] == "highest" and b["uncertain"] is True
    assert (b["session"], b["rotor"], b["nucleus"]) == (str(root), "RS2427418", "31P")
    assert b["confirmed"] is None
    # the store applies: rate, flag, source and a note naming the confirmation
    key = (b["session"], b["rotor"], b["nucleus"])
    M.remember(key, M.MasEvidence.from_dict(b), 20000, note="Experiment parameters")
    r2 = M.resolve_for_load(meta, e)
    assert r2["spin_rate_Hz"] == 20000.0 and r2["mas_uncertain"] is False
    b2 = r2["provenance"]
    assert b2["source"] == "confirmed" and b2["uncertain"] is False
    assert b2["confirmed"] and b2["confirmed"][:2] == "20"
    assert b2["confirmed_note"] == "Experiment parameters"
    assert "confirmed on" in r2["note"] and "rotor RS2427418 · 31P" in r2["note"]
    assert "acqus 4 200 / title 20 000 / booking 22 000" in r2["note"]
    # a certain resolution never consults the store
    _addenda(e, 20)
    meta3 = bruker._meta_1d({"MASR": 4200, "NUC1": "<31P>"}, title, e)
    assert meta3["mas_uncertain"] is False
    M.remember(key, M.MasEvidence.from_dict(meta3["mas_sources"]), 12345)
    r3 = M.resolve_for_load(meta3, e)
    assert r3["spin_rate_Hz"] == 20000.0 and r3["provenance"]["source"] == "title+booking"
    # a meta without the sources block (older reader, QCPMG) still works
    r4 = M.resolve_for_load({"spin_rate_Hz": 5000.0, "mas_uncertain": False,
                             "nucleus": "35Cl"}, "")
    assert r4["spin_rate_Hz"] == 5000.0 and r4["provenance"]["session"] is None


def test_loader_applies_confirmation_on_real_expno(store):
    """2702: acqus 4200 / title 'MASR 20kHz' / booking 22 -- three-way, so
    the loader flags it, and one remembered confirmation clears it."""
    from larmor.loader import load_any

    path = str(require(BRUKER_1R))
    _, _, rec, _, warns = load_any(path)
    assert rec["spin_rate_Hz"] == 22000.0 and rec["mas_uncertain"] is True
    b = rec["provenance"]["mas_rate"]
    assert b["source"] == "highest"
    assert (b["acqus_Hz"], b["title_Hz"], b["booking_Hz"]) == (4200.0, 20000.0, 22000.0)
    assert b["rotor"] == "RS2427418" and b["nucleus"] == "27Al"
    assert b["session"].endswith("2026-05")
    assert sum("booking sidecar says 22000 Hz" in w for w in warns) == 1
    key = (b["session"], b["rotor"], b["nucleus"])
    M.remember(key, M.MasEvidence.from_dict(b), 20000, note="test")
    _, _, rec2, _, warns2 = load_any(path)
    assert rec2["spin_rate_Hz"] == 20000.0 and rec2["mas_uncertain"] is False
    assert rec2["provenance"]["mas_rate"]["source"] == "confirmed"
    assert any("confirmed on" in w for w in warns2)
    # the recipe round-trips the block unchanged
    from larmor.recipe import Recipe
    again = Recipe.from_dict(json.loads(json.dumps(rec2))).to_dict()
    assert again["provenance"]["mas_rate"] == rec2["provenance"]["mas_rate"]


def test_module_is_qt_free_and_exports():
    probe = ("import sys, larmor.masrate; "
             "print(sorted(k for k in sys.modules if k.split('.')[0] in "
             "('PySide6', 'pyqtgraph', 'shiboken6')))")
    out = subprocess.run([sys.executable, "-c", probe], capture_output=True,
                         text=True, cwd=REPO, timeout=120)
    assert out.returncode == 0, out.stderr[-2000:]
    assert out.stdout.strip().splitlines()[-1] == "[]"
    # every function / class defined here and every UPPER constant is exported
    defined = {n for n, o in inspect.getmembers(M)
               if not n.startswith("_") and (
                   (inspect.isfunction(o) or inspect.isclass(o))
                   and getattr(o, "__module__", "") == M.__name__
                   or (n.isupper() and not inspect.ismodule(o)))}
    assert defined == set(M.__all__), defined ^ set(M.__all__)
