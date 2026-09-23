"""Every tutorial must ship, keep the house style, and name only things that
exist: each bold **Menu > Item** path resolves to a real menu title / action
label of the offscreen MainWindow, each fenced ``larmor <sub>`` line names a
real CLI subcommand (and only real flags), each ``examples/...`` path exists,
each ``module.name(`` in a python fence resolves, external data is declared
honestly, and the Help ▸ Tutorials submenu lists all of them and renders them.

The two tutorials written on the developer's instrument data (4 and 7) also
get a layout check behind conftest.require(): it pins the premises the text
teaches (the acqus-vs-title MAS conflict of the 11B series; the stale
MASR=14000 under WCPMG of the 81Br set) against the files, and skips cleanly
on a machine without them (the real-data banner then names them)."""
import importlib
import inspect
import os
import re
from pathlib import Path

import pytest

from conftest import LAW_CA_11B, MAGLAB_81BR, MAGLAB_35CL, DATA_ROOT, require

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

ROOT = Path(__file__).resolve().parents[1]
TUT_DIR = ROOT / "docs" / "tutorials"

# (file stem, Help ▸ Tutorials title) -- must equal larmor.desktop.app.TUTORIALS
TUTORIALS = [
    ("01-first-fit-27Al-czjzek", "1 · A first fit — ²⁷Al Czjzek"),
    ("02-constraints", "2 · Constraints — fix, bound, link"),
    ("03-figures", "3 · Plotting studio figures"),
    ("04-batch-fitting", "4 · Batch fitting a composition series"),
    ("05-mqmas", "5 · MQMAS — 2D processing & fitting"),
    ("06-error-analysis", "6 · Error analysis — covariance, Monte-Carlo, χ² profile"),
    ("07-static-81Br-wcpmg", "7 · Static wideline — ⁸¹Br WURST-CPMG"),
]
STEMS = [s for s, _ in TUTORIALS]
TOP_MENUS = {"File", "Process", "Decomposition", "View", "Tools", "Plotting", "Help"}


# ------------------------------------------------------------------ helpers
def _text(stem: str) -> str:
    return (TUT_DIR / f"{stem}.md").read_text(encoding="utf-8")


def _norm(label: str) -> str:
    """A Qt label as the user sees it: '&&' -> '&', the mnemonic '&' dropped,
    whitespace collapsed."""
    s = label.replace("&&", "\x00").replace("&", "").replace("\x00", "&")
    return re.sub(r"\s+", " ", s).strip()


def _candidates(label: str) -> set[str]:
    """The strings a tutorial may use for one label: the full label, or the
    label with its trailing '(hint)' explainer removed -- each with and
    without the trailing ellipsis."""
    n = _norm(label)
    out = {n, re.sub(r"\s+\([^()]*\)\s*$", "", n).strip()}
    return out | {c.rstrip("…").strip() for c in out}


def _bold_paths(md: str) -> list[list[str]]:
    """Every **A > B > C** in the text (dialog buttons are bold WITHOUT '>')."""
    return [[seg.strip() for seg in cap.split(" > ")]
            for cap in re.findall(r"\*\*([^*\n]+?)\*\*", md) if " > " in cap]


def _fenced(md: str, lang: str | None = None) -> list[str]:
    """The bodies of ``` fences (optionally only those tagged ``lang``)."""
    out = []
    for m in re.finditer(r"```([\w+-]*)\n(.*?)```", md, flags=re.DOTALL):
        if lang is None or m.group(1) == lang:
            out.append(m.group(2))
    return out


def _menu_labels(win) -> set[str]:
    from PySide6.QtWidgets import QMenu

    labels = set()
    for m in win.menuBar().findChildren(QMenu):
        labels |= _candidates(m.title())
        for a in m.actions():
            if a.text():
                labels |= _candidates(a.text())
    return labels


# ----------------------------------------------------------------- fixtures
@pytest.fixture(scope="module")
def qapp():
    from PySide6.QtWidgets import QApplication

    yield QApplication.instance() or QApplication([])


@pytest.fixture(scope="module")
def win(qapp):
    os.environ["LARMOR_NO_SESSION"] = "1"
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    yield w
    w.close()


# -------------------------------------------------------------------- tests
@pytest.mark.parametrize("stem", STEMS)
def test_tutorial_ships_and_keeps_the_house_style(stem):
    p = TUT_DIR / f"{stem}.md"
    assert p.is_file(), f"docs/tutorials/{stem}.md is missing"
    text = p.read_text(encoding="utf-8")
    assert len(text) > 3000, f"{stem}.md looks too thin to be a tutorial"
    first = text.splitlines()[0]
    m = re.match(r"^# Tutorial (\d+) — ", first)
    assert m, f"{stem}.md must start with '# Tutorial N — …' (got {first!r})"
    assert int(m.group(1)) == STEMS.index(stem) + 1, "tutorial number ≠ position"
    paras = [b.strip() for b in text.split("\n\n") if b.strip()]
    assert paras[1].startswith("*Time:"), "the preamble must open with '*Time:'"
    assert text.count("\n## ") >= 4, "fewer than four numbered sections"
    assert "\t" not in text, "tabs in a tutorial"
    assert "$$" not in text, "display math does not belong in a tutorial"


def test_no_unexpected_tutorial_files():
    assert sorted(p.stem for p in TUT_DIR.glob("*.md")) == STEMS


@pytest.mark.parametrize("stem", STEMS)
def test_tutorial_menu_paths_exist_in_app(stem, win):
    labels = _menu_labels(win)
    bad = []
    for path in _bold_paths(_text(stem)):
        if path[0] not in TOP_MENUS:
            bad.append(" > ".join(path) + "  (not a top-level menu)")
            continue
        for seg in path[1:]:
            if seg.rstrip("…").strip() not in labels:
                bad.append(" > ".join(path) + f"  (no action {seg!r})")
    assert not bad, f"{stem}.md names menu paths the app does not have:\n  " + "\n  ".join(bad)


@pytest.mark.parametrize("stem", STEMS)
def test_tutorial_cli_commands_exist(stem):
    import larmor.cli as cli

    src = inspect.getsource(cli)
    subs = set(re.findall(r'add_parser\("(\w+)"', src))
    flags = set(re.findall(r'"(--[\w-]+)"', src)) | {"-o", "-v"}
    bad = []
    for body in _fenced(_text(stem)):
        for line in body.splitlines():
            m = re.match(r"^\s*larmor\s+([a-z]+)(.*)$", line)
            if not m:
                continue
            if m.group(1) not in subs:
                bad.append(f"{line.strip()!r}: no subcommand {m.group(1)!r}")
            for tok in re.findall(r"(?<![\w-])(-{1,2}[\w-]+)", m.group(2)):
                if tok not in flags and not re.fullmatch(r"-\d+(\.\d+)?", tok):
                    bad.append(f"{line.strip()!r}: unknown flag {tok!r}")
    assert not bad, f"{stem}.md:\n  " + "\n  ".join(bad)


@pytest.mark.parametrize("stem", STEMS)
def test_tutorial_example_paths_exist(stem):
    missing = []
    for tok in re.findall(r"examples/[\w./-]+", _text(stem)):
        tok = tok.rstrip(".,;:")
        if not (ROOT / tok).exists():
            missing.append(tok)
    assert not missing, f"{stem}.md cites examples paths that do not exist: {missing}"


@pytest.mark.parametrize("stem", STEMS)
def test_tutorial_python_symbols_exist(stem):
    bad = []
    for body in _fenced(_text(stem), "python"):
        for mod, name in re.findall(
                r"\b(autofit|chi2map|engine|figures|loader|batchfit|twod|fit|"
                r"qcpmg|qcpmg_fields|staticct|convert|multifit|processing)\.(\w+)\(", body):
            module = importlib.import_module(f"larmor.{mod}")
            if not hasattr(module, name):
                bad.append(f"larmor.{mod}.{name}")
    assert not bad, f"{stem}.md calls symbols that do not exist: {bad}"


@pytest.mark.parametrize("stem", STEMS)
def test_tutorial_declares_external_data(stem):
    text = _text(stem)
    if "<NMR>" in text or "WSU_work" in text:
        preamble = text.split("\n## ", 1)[0]
        assert "not part of the repository" in preamble, \
            f"{stem}.md uses private data without saying so in the preamble"
    assert "C:/Users/" not in text and "C:\\Users\\" not in text, \
        f"{stem}.md carries a username-bearing absolute path"


def test_all_tutorials_in_help_menu(win):
    from PySide6.QtWidgets import QMenu

    from larmor.desktop import app as appmod
    from larmor.desktop.help_dialog import tutorial_path

    assert list(appmod.TUTORIALS) == TUTORIALS
    tut = next(m for m in win.menuBar().findChildren(QMenu)
               if m.title() == "&Tutorials")
    assert [a.text() for a in tut.actions() if a.text()] == [t for _, t in TUTORIALS]
    for stem in STEMS:
        assert tutorial_path(stem) is not None, stem


@pytest.mark.parametrize("stem", STEMS)
def test_tutorial_renders_to_html(stem, qapp):
    from larmor.desktop.mdrender import render_help_html

    html = render_help_html(_text(stem))
    assert "MDMATH" not in html
    assert "<p" in html or "<pre" in html


# ------------------------------------------------------ real-data premises
_NMR = DATA_ROOT / "Desktop/WSU_work/NMR"


def test_tutorial07_dataset_layout():
    from larmor.io import bruker

    expno = require(MAGLAB_81BR / "30")
    exp = bruker.read_expno(expno)
    assert exp.nucleus == "81Br"
    assert "wcpmg" in exp.pulse_program.lower()
    assert float(exp.masr_Hz) == 14000.0
    assert exp.conflicts == []                    # the title carries no kHz
    assert (expno / "fid").exists() and (expno / "pdata" / "1" / "1r").exists()
    d = bruker.read(str(expno))
    assert d.meta.get("mas_uncertain") is True, "the stale MASR must be flagged"
    for k in ("31", "32", "33", "34"):
        assert (MAGLAB_81BR / k / "fid").exists()
    # the two-field 35Cl route: the LAW3Cl0Ca pair and the ten saved datasets
    require(MAGLAB_35CL / "1" / "fid")
    require(_NMR / "NMRFAM/DATA/2026-06_35Cl/06102026_RS40175_LAW3Cl0Ca_SS_ALP/3/fid")
    from larmor.io import spectra
    for k in range(5):
        _, _, m_lo = spectra.read_csv(require(_NMR / f"MagLab/DATA/LAW{k}Ca-3Cl_850_MHz.csv"))
        _, _, m_hi = spectra.read_csv(require(
            _NMR / f"NMRFAM/DATA/2026-06_35Cl/LAW{k}Ca-3Cl_1p1GHz.csv"))
        assert m_lo["larmor_MHz"] == pytest.approx(78.354, rel=1e-3)
        assert m_hi["larmor_MHz"] == pytest.approx(107.811, rel=1e-3)


def test_tutorial04_series_layout():
    from larmor.io import bruker

    require(LAW_CA_11B[0])
    for expno in LAW_CA_11B:
        assert (expno / "pdata" / "1" / "1r").exists(), expno
        exp = bruker.read_expno(expno)
        assert exp.nucleus == "11B"
        assert exp.conflicts and exp.conflicts[0].startswith(
            "MAS rate: acqus says 4200 Hz but the title says 35714 Hz"), exp.conflicts
