"""Unit tests of the command palette's Qt-free matcher (larmor.palette).

Plain strings only: no PySide6 import, no QApplication. The scores pinned
here are the ones the palette dialog relies on for its ranking (word-start
contiguous > interior contiguous > scattered subsequence; ties in menu
order; an empty query is the menu map).
"""
import subprocess
import sys
from pathlib import Path

from larmor.palette import (
    fold, fuzzy_score, rank, strip_hint, strip_mnemonic,
)

REPO = Path(__file__).resolve().parents[1]


def test_strip_mnemonic_and_hint():
    assert strip_mnemonic("&Open…  (spectrum / recipe / 1r / 2rr)") == \
        "Open…  (spectrum / recipe / 1r / 2rr)"
    assert strip_mnemonic("&Integrals && measurements…  (integral, %, FWHM, CoM)") == \
        "Integrals & measurements…  (integral, %, FWHM, CoM)"
    assert strip_mnemonic("Ad&vanced") == "Advanced"
    assert strip_mnemonic("Fit") == "Fit"
    assert strip_hint("Add &line  (pick a model)") == "Add line"
    assert strip_hint("&Apply recipe  (a recent fit → this data)") == "Apply recipe"
    assert strip_hint("Autophase (ACME)") == "Autophase (ACME)"     # single space: a name
    assert strip_hint("Czjzek &width display") == "Czjzek width display"


def test_fold_transliterates_greek_accents_superscripts():
    assert fold("χ² map (parameter pair)…") == "chi2 map (parameter pair)"
    assert "eden" in fold("Literature shift ranges  (Edén 2023)")
    assert fold("Measure Δ (ppm / Hz)").startswith("measure delta")
    assert fold("νrot") == "nurot"
    assert fold("FIT") == "fit"
    assert fold("Read static pattern (C_Q, η)") == "read static pattern (c_q, eta)"


def test_subsequence_required_and_ordered():
    s = fuzzy_score("czj", "Decomposition ▸ Czjzek distribution P(C_Q)…  (what σ stands for)")
    assert isinstance(s, int) and s > 0
    assert fuzzy_score("xyz", "Fit") is None
    assert fuzzy_score("tif", "Fit") is None            # order matters
    assert fuzzy_score("chi2", "Decomposition ▸ χ² map (parameter pair)…") is not None


def test_case_insensitive():
    assert fuzzy_score("FIT", "Decomposition ▸ Fit") == \
        fuzzy_score("fit", "Decomposition ▸ Fit") == \
        fuzzy_score("fit", "DECOMPOSITION ▸ FIT")


def test_word_start_substring_ranks_first_then_shorter():
    # whole word 24 > word-start 21 > interior contiguous 16 == 16 -> shorter first
    assert fuzzy_score("fit", "Fit") == 24
    assert fuzzy_score("fit", "Fits") == 21
    assert fuzzy_score("fit", "Profit") == fuzzy_score("fit", "Refit") == 16
    assert rank("fit", ["Profit", "Fit", "Refit"]) == [1, 2, 0]


def test_whole_word_beats_shorter_word_start_text():
    # measured on the real menus: without the whole-word bonus the shorter
    # 'View ▸ Animate fits' won the tie and 'type fit, Enter' ran the wrong
    # command
    texts = ["View ▸ Animate fits", "Decomposition ▸ Fit  F5",
             "Decomposition ▸ Add fit zone"]
    assert rank("fit", texts) == [1, 2, 0]


def test_consecutive_beats_scattered():
    assert fuzzy_score("fit", "Fit") > fuzzy_score("fit", "Fixture at") > 0


def test_multi_token_all_must_match_any_order():
    assert fuzzy_score("theme dark", "View ▸ Theme ▸ Dark") is not None
    assert fuzzy_score("dark theme", "View ▸ Theme ▸ Dark") is not None
    assert fuzzy_score("theme dark", "View ▸ Theme ▸ Light") is None


def test_empty_query_keeps_menu_order():
    assert rank("", ["b", "a", "c"]) == [0, 1, 2]
    assert rank("   ", ["b", "a", "c"]) == [0, 1, 2]
    assert fuzzy_score("", "anything") == 0


def test_ties_are_stable():
    assert rank("a", ["Alpha", "Apple"]) == [0, 1]
    texts = ["Alpha", "Apple", "Amber", "Fit", "Azure"]
    assert rank("a", texts) == rank("a", texts)


def test_shortcut_token_in_haystack_matches_and_nonmatches_are_dropped():
    texts = ["Decomposition ▸ Simulate  (recompute the model)  F9",
             "Decomposition ▸ Fit  F5"]
    assert rank("f5", texts) == [1]


def test_word_start_pass_falls_back_to_plain_greedy():
    # the word-start pass takes a@5 and dead-ends on 'b'; the plain pass
    # finds a@1, b@3
    assert fuzzy_score("ab", "xa b_a") is not None


def test_realistic_fit_ranking():
    texts = [
        "Decomposition ▸ Advanced ▸ Fit completion threshold…  (Δσ % to stop)",
        "Tools ▸ Batch fit report…  (publication table + plots)",
        "File ▸ Save fit as…  (txt / csv / json / dmfit)  Ctrl+Shift+E",
        "Decomposition ▸ Fit  F5",
        "View ▸ Full spectrum",
    ]
    r = rank("fit", texts)
    assert r[0] == 3                     # shortest of the word-start contiguous hits
    assert 4 not in r                    # 'Full spectrum' has no i after the f
    assert set(r) == {0, 1, 2, 3}


def test_module_is_qt_free():
    """Belt and braces next to tests/test_core_qt_free.py's sweep."""
    code = ("import sys, larmor.palette; print(any(m.split('.')[0] in "
            "('PySide6', 'shiboken6') for m in sys.modules))")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True,
                         text=True, cwd=str(REPO), timeout=120)
    assert out.returncode == 0, out.stderr[-2000:]
    assert out.stdout.strip() == "False"
