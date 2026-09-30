"""The LARMOR-side SR corrections (larmor.referencing.set_override /
override_for / clear_override): the store keys the EXPNO whatever inner
path is given, every set / clear is a record of the referencing log, the
apply decision moves the axis in TopSpin's direction or reports a stale
override -- and the loader applies a live override (ignores a stale one) on
the shipped examples/pCABS2-4/3616, so the desktop, Batch fit, Sequential
fit and the CLI all see the corrected axis. Every store here is a tmp file:
the developer's %LOCALAPPDATA% is never read or written."""
import json
import os
from pathlib import Path

import numpy as np
import pytest

from larmor import loader
from larmor import referencing as R
from larmor.recipe import Recipe

ROOT = Path(__file__).resolve().parents[1]
EXPNO_3616 = ROOT / "examples" / "pCABS2-4" / "3616"


@pytest.fixture()
def store(tmp_path, monkeypatch):
    monkeypatch.setenv("LARMOR_SR_OVERRIDES", str(tmp_path / "store" / "sr_overrides.json"))
    monkeypatch.setenv("LARMOR_REF_LOG", str(tmp_path / "store" / "referencing_log.jsonl"))
    return tmp_path


def _fake_expno(root: Path, sample: str, expno: int) -> Path:
    e = root / sample / str(expno)
    (e / "pdata" / "1").mkdir(parents=True)
    (e / "acqus").write_text("##TITLE= x\n##$NUC1= <27Al>\n##$BF1= 130.3\n##END=\n",
                             encoding="utf-8")
    (e / "pdata" / "1" / "1r").write_bytes(b"\0" * 8)
    return e


def test_store_round_trip_keys_the_expno_and_logs_every_change(store):
    root = store / "2026-05"
    e = _fake_expno(root, "04272026_P1-Bi1-12_SS_ALP", 2702)
    ref = root / "04272026_RS40339_P5-Bi1-12" / "1"
    assert R.overrides_path() == store / "store" / "sr_overrides.json"
    assert R.override_for(e) is None and R.override_for("") is None
    rec = R.set_override(e / "pdata" / "1" / "1r", -135.31, 0.0, reference=str(ref),
                         note="test", nucleus="27Al")
    assert rec["path"] == str(e)                        # keyed by the EXPNO, not the 1r
    probes = [e, e / "pdata" / "1", e / "pdata" / "1" / "1r", str(e) + os.sep]
    if os.name == "nt":
        probes.append(str(e).upper())                   # Windows compares case-blind
    for p in probes:
        got = R.override_for(p)
        assert got and got["new_sr_hz"] == -135.31 and got["old_sr_hz"] == 0.0, p
        assert got["reference"] == str(ref) and got["nucleus"] == "27Al"
    assert R.override_for(root / "other" / "1") is None
    # a second record for the same EXPNO replaces the first
    R.set_override(e, -135.40, 0.0)
    assert R.override_for(e)["new_sr_hz"] == -135.40
    data = json.loads(R.overrides_path().read_text(encoding="utf-8"))
    assert len(data) == 1
    # clearing returns what was dropped; a second clear is a no-op
    assert R.clear_override(e)["new_sr_hz"] == -135.40
    assert R.override_for(e) is None and R.clear_override(e) is None
    assert json.loads(R.overrides_path().read_text(encoding="utf-8")) == {}
    # every set / clear is one record of the referencing log, filed under
    # the EXPNO's session, so previous_audits() sees the whole trail
    recs = [json.loads(ln) for ln in R.log_path().read_text(encoding="utf-8").splitlines()]
    assert [r["action"] for r in recs] == ["override", "override", "override-cleared"]
    assert all(Path(r["session"]) == root for r in recs)
    ent = recs[0]["entries"][0]
    assert (ent["sample"], ent["expno"], ent["nucleus"], ent["procno"]) == \
        ("04272026_P1-Bi1-12_SS_ALP", 2702, "27Al", 1)
    assert ent["old_sr_hz"] == 0.0 and ent["new_sr_hz"] == -135.31 and ent["note"] == "test"
    assert ent["reference"] == str(ref) and Path(ent["path"]) == e
    assert recs[-1]["entries"][0]["new_sr_hz"] == -135.40
    assert len(R.previous_audits(root)) == 3
    # a damaged store reads as empty rather than raising
    R.overrides_path().write_text("{not json", encoding="utf-8")
    assert R.override_for(e) is None


def test_apply_override_moves_the_axis_or_reports_stale():
    ppm = np.linspace(-50.0, 150.0, 5)
    ov = {"old_sr_hz": 0.0, "new_sr_hz": -135.31, "note": "n"}
    out, sr, prov, note = R.apply_override(ov, ppm, 0.2, 156.281744)
    assert np.allclose(out, ppm + 135.31 / 156.281744)      # a smaller SR: higher ppm
    assert sr == -135.31 and prov["source"] == "sr_overrides" and prov["note"] == "n"
    assert prov["axis_shift_ppm"] == pytest.approx(+0.8658, abs=1e-3)
    assert prov["old_sr_hz"] == 0.0 and prov["new_sr_hz"] == -135.31
    assert note.startswith("SR corrected by the referencing audit (was 0.00 Hz)")
    # TopSpin fixed the file since: the recorded old SR no longer matches
    out2, sr2, prov2, note2 = R.apply_override(ov, ppm, -135.30, 156.281744)
    assert out2 is ppm and sr2 == -135.30 and prov2 is None and "stale" in note2
    # right at the tolerance the override still applies
    out3, sr3, prov3, _ = R.apply_override(ov, ppm, 0.5, 156.281744)
    assert prov3 is not None and sr3 == -135.31


def test_loader_applies_a_live_override_and_ignores_a_stale_one(store):
    if not (EXPNO_3616 / "acqus").exists():
        pytest.skip("examples/pCABS2-4/3616 not shipped here")
    ppm0, amp0, rec0, _m, w0 = loader.load_any(EXPNO_3616)
    sr0, sf = float(rec0["sr_hz"]), float(rec0["larmor_frequency_MHz"])
    assert sr0 == pytest.approx(-607.40, abs=0.01)         # (SF - BF1) of the shipped procs
    assert not any("SR corrected" in w or "stale" in w for w in w0)
    assert "referencing" not in rec0["provenance"]
    R.set_override(EXPNO_3616, sr0 + 100.0, sr0, reference="", note="pin", nucleus="27Al")
    ppm1, amp1, rec1, _m, w1 = loader.load_any(EXPNO_3616)
    # the axis moves by exactly -(new - old) / SF, TopSpin's direction
    assert np.allclose(ppm1, ppm0 - 100.0 / sf, atol=1e-9)
    assert np.allclose(ppm1, ppm0 + R.axis_shift_ppm(sr0, sr0 + 100.0, sf), atol=1e-9)
    assert np.array_equal(amp1, amp0)
    assert rec1["sr_hz"] == pytest.approx(sr0 + 100.0)
    prov = rec1["provenance"]["referencing"]
    assert prov["old_sr_hz"] == pytest.approx(sr0) and prov["new_sr_hz"] == pytest.approx(sr0 + 100.0)
    assert prov["source"] == "sr_overrides" and prov["note"] == "pin"
    assert prov["axis_shift_ppm"] == pytest.approx(-100.0 / sf)
    assert any(w.startswith("SR corrected by the referencing audit (was") for w in w1)
    # the 1r path keys the same EXPNO
    ppm1r, *_ = loader.load_any(EXPNO_3616 / "pdata" / "1" / "1r")
    assert np.allclose(ppm1r, ppm1)
    # a recipe saved from the corrected load reopens on the same axis with
    # no "re-referenced" note: stored SR and freshly read SR agree
    r = Recipe.from_dict(rec1)
    p = store / "corrected.recipe.json"
    r.save(p)
    ppm2, _a, rec2, _m, w2 = loader.load_any(p)
    assert np.allclose(ppm2, ppm1) and rec2["sr_hz"] == pytest.approx(sr0 + 100.0)
    assert not any("re-referenced" in w for w in w2)
    # stale: the recorded old SR is not what the file says -> nothing applied
    R.set_override(EXPNO_3616, sr0 + 100.0, sr0 + 50.0)
    ppm3, _a, rec3, _m, w3 = loader.load_any(EXPNO_3616)
    assert np.array_equal(ppm3, ppm0) and rec3["sr_hz"] == pytest.approx(sr0)
    assert "referencing" not in rec3["provenance"]
    assert any("stale" in w for w in w3)
    # the recipe fitted on the corrected axis now reopens on the file's axis
    # (its provenance names the file's SR as the one it started from, so the
    # reload note stays quiet -- provenance._sr_explained's existing rule)
    ppm4, *_ = loader.load_any(p)
    assert np.allclose(ppm4, ppm0)
    R.clear_override(EXPNO_3616)
    ppm5, _a, rec5, _m, w5 = loader.load_any(EXPNO_3616)
    assert np.array_equal(ppm5, ppm0) and not any("stale" in w for w in w5)
