"""The "Keep fit parameters?" question: same nucleus keeps dmfit's Yes
default; a different nucleus defaults to No and names both nuclei, so a 27Al
model is never fitted to a 31P spectrum by pressing Enter."""
import pytest

pytest.importorskip("PySide6")


def test_prompt_defaults_by_nucleus():
    from larmor.desktop.mw_files import keep_fit_prompt

    text, yes = keep_fit_prompt("27Al", "27Al")
    assert yes and "Yes = keep the lines" in text
    text, yes = keep_fit_prompt("27Al", "31P")
    assert not yes
    assert "27Al model" in text and "spectrum is 31P" in text
    assert "No = start empty (recommended)" in text
    # unknown nucleus on either side: no basis to refuse, keep the old default
    assert keep_fit_prompt("", "31P")[1] and keep_fit_prompt("27Al", None)[1]
