"""Fuzzy matcher for the desktop command palette (Qt-free).

The palette (larmor.desktop.palette_dialog) lists every menu and toolbar
entry as "Menu ▸ Submenu ▸ Label  (hint)  Shortcut" and filters that list as
the user types. This module holds the text side of that: the mnemonic/hint
conventions of MainWindow._build_menus, a Unicode fold that makes the menus'
Greek letters, accents and superscripts reachable from an ASCII keyboard
("chi2" finds "χ² map", "eden" finds "Edén"), and a small scoring function
that prefers contiguous word-start matches over scattered subsequences.

Kept free of Qt so it can be unit-tested with plain strings; the sweep in
tests/test_core_qt_free.py enforces that.
"""
from __future__ import annotations

import unicodedata
from collections.abc import Sequence

__all__ = [
    "GREEK", "strip_mnemonic", "strip_hint", "fold", "is_word_start",
    "fuzzy_score", "rank",
]

#: ASCII spellings of the Greek letters used in menu labels (χ² map, Δ,
#: νrot, σ, η ...). Keys are the casefolded forms; capitals fold to them.
GREEK: dict[str, str] = {
    "χ": "chi", "δ": "delta", "Δ": "delta", "ν": "nu", "σ": "sigma",
    "η": "eta", "ζ": "zeta", "τ": "tau", "ω": "omega",
}


def strip_mnemonic(text: str) -> str:
    """Remove Qt mnemonics: ``&&`` becomes ``&``, a single ``&`` is dropped."""
    out = []
    i = 0
    while i < len(text):
        ch = text[i]
        if ch == "&":
            if i + 1 < len(text) and text[i + 1] == "&":
                out.append("&")
                i += 2
                continue
            i += 1
            continue
        out.append(ch)
        i += 1
    return "".join(out)


def strip_hint(text: str) -> str:
    """Drop the two-space ``  (hint)`` suffix used by the menu builder, plus
    mnemonics: 'Add &line  (pick a model)' -> 'Add line'. A single-space
    parenthetical such as 'Autophase (ACME)' is part of the name and stays."""
    text = strip_mnemonic(text)
    cut = text.find("  (")
    if cut >= 0:
        text = text[:cut]
    return text.rstrip()


def fold(s: str) -> str:
    """Casefold, strip accents and superscripts (NFKD), transliterate Greek
    letters and drop ellipses, so ASCII typing reaches every menu label."""
    # the ellipsis is removed BEFORE NFKD, which would expand it to "..."
    s = unicodedata.normalize("NFKD", s.replace("…", "").casefold())
    out = []
    for ch in s:
        if unicodedata.combining(ch):
            continue
        out.append(GREEK.get(ch, ch))
    return "".join(out)


def is_word_start(text: str, i: int) -> bool:
    """True when position ``i`` begins a word: index 0, or preceded by a
    non-alphanumeric character (space, ▸, parenthesis, dash, slash ...)."""
    return i == 0 or not text[i - 1].isalnum()


def _subsequence(token: str, hay: str, prefer_word_start: bool) -> list[int] | None:
    """Greedy left-to-right subsequence match; returns the hit positions or
    None. With ``prefer_word_start`` each character takes the nearest
    word-start occurrence when there is one, else the nearest occurrence."""
    hits: list[int] = []
    start = 0
    for ch in token:
        i = -1
        if prefer_word_start:
            j = hay.find(ch, start)
            while j >= 0 and not is_word_start(hay, j):
                j = hay.find(ch, j + 1)
            i = j
        if i < 0:
            i = hay.find(ch, start)
        if i < 0:
            return None
        hits.append(i)
        start = i + 1
    return hits


def _score_token(token: str, hay: str) -> int | None:
    """Score one folded query token against a folded haystack."""
    pos = hay.find(token)
    if pos >= 0:
        score = 10 + 2 * len(token)
        word_start = whole_word = False
        while pos >= 0:
            if is_word_start(hay, pos):
                word_start = True
                end = pos + len(token)
                if end == len(hay) or not hay[end].isalnum():
                    whole_word = True
                    break
            pos = hay.find(token, pos + 1)
        # 'fit' must reach 'Fit' before 'Animate fits' (both word starts;
        # the shorter text would otherwise win the tie)
        return score + (5 if word_start else 0) + (3 if whole_word else 0)
    hits = _subsequence(token, hay, prefer_word_start=True)
    if hits is None:
        hits = _subsequence(token, hay, prefer_word_start=False)
    if hits is None:
        return None
    score = 0
    prev = -2
    for i in hits:
        if is_word_start(hay, i):
            score += 3
        elif i == prev + 1:
            score += 2
        else:
            score += 1
        prev = i
    span = hits[-1] - hits[0] + 1
    score -= (span - len(token)) // 4
    return max(score, 1)


def fuzzy_score(query: str, text: str) -> int | None:
    """Score ``text`` against ``query``; None means no match, a blank query
    scores 0. Both sides are folded; the query is split on whitespace and
    every token must match (in any order), token scores summed. A contiguous
    substring scores 10 + 2*len (+5 when it starts a word, +3 more when it
    is a whole word); otherwise a word-start-preferring greedy subsequence
    scores +3 per word-start hit, +2 per consecutive hit, +1 otherwise,
    minus a small gap penalty."""
    tokens = fold(query).split()
    if not tokens:
        return 0
    hay = fold(text)
    total = 0
    for tok in tokens:
        s = _score_token(tok, hay)
        if s is None:
            return None
        total += s
    return total


def rank(query: str, texts: Sequence[str]) -> list[int]:
    """Indices of the matching ``texts`` ordered by (score desc, length asc,
    original index). A blank query returns every index in the original order,
    so an empty palette query is a complete map of the menus."""
    if not fold(query).split():
        return list(range(len(texts)))
    scored = []
    for i, t in enumerate(texts):
        s = fuzzy_score(query, t)
        if s is not None:
            scored.append((-s, len(t), i))
    scored.sort()
    return [i for _, _, i in scored]
