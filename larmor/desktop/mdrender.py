"""Prettier help rendering for the Markdown manuals.

QTextBrowser can show Markdown, but it can't typeset equations. Here we
(1) render every LaTeX `$...$` / `$$...$$` snippet to a crisp transparent PNG
with matplotlib's mathtext and splice it back in as an inline image, and
(2) supply a stylesheet so the manuals read like a real document rather than
raw Markdown. No web engine, no extra dependencies — matplotlib is already a
core dependency.
"""
from __future__ import annotations

import base64
import io
import re

# body text / heading colours (light theme; the help dialog has a white page)
_INK = "#22303a"
_BRAND = "#1f3a5f"
_math_cache: dict = {}


def _math_datauri(tex: str, fontsize: float, color: str = _INK) -> str:
    key = (tex, fontsize, color)
    if key in _math_cache:
        return _math_cache[key]
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    fig = Figure()
    FigureCanvasAgg(fig)
    fig.patch.set_alpha(0.0)
    fig.text(0.0, 0.0, f"${tex}$", fontsize=fontsize, color=color)
    buf = io.BytesIO()
    try:
        fig.savefig(buf, dpi=150, transparent=True, bbox_inches="tight",
                    pad_inches=0.02)
    except Exception:
        # a malformed expression shouldn't blow up the whole manual
        return ""
    uri = "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")
    _math_cache[key] = uri
    return uri


def _protect_code(md: str):
    """Stash fenced/inline code so `$` inside code is never treated as math."""
    blocks: list[str] = []

    def stash(m):
        blocks.append(m.group(0))
        return f"\x00{len(blocks) - 1}\x00"

    md = re.sub(r"```.*?```", stash, md, flags=re.DOTALL)
    md = re.sub(r"`[^`\n]*`", stash, md)
    return md, blocks


def _restore_code(md: str, blocks: list[str]) -> str:
    return re.sub(r"\x00(\d+)\x00", lambda m: blocks[int(m.group(1))], md)


def render_math(md: str) -> str:
    """Replace $$display$$ and $inline$ math with rendered-image Markdown."""
    md, blocks = _protect_code(md)

    def display(m):
        uri = _math_datauri(m.group(1).strip(), 17)
        return f"\n\n![equation]({uri})\n\n" if uri else m.group(0)

    def inline(m):
        uri = _math_datauri(m.group(1).strip(), 12)
        return f"![eq]({uri})" if uri else m.group(0)

    md = re.sub(r"\$\$(.+?)\$\$", display, md, flags=re.DOTALL)
    md = re.sub(r"(?<!\$)\$(?!\$)(.+?)(?<!\$)\$(?!\$)", inline, md,
                flags=re.DOTALL)
    return _restore_code(md, blocks)


HELP_CSS = f"""
    body {{ color: {_INK}; font-size: 14px; line-height: 155%; }}
    h1 {{ color: {_BRAND}; font-size: 25px; }}
    h2 {{ color: {_BRAND}; font-size: 19px; }}
    h3 {{ color: {_BRAND}; font-size: 15px; }}
    a  {{ color: #0a5a62; text-decoration: none; }}
    code {{ font-family: Consolas, "Courier New", monospace; font-size: 12px;
            background: #eef2ef; color: #16202a; }}
    pre {{ font-family: Consolas, "Courier New", monospace; font-size: 12px;
           background: #f4f6f4; color: #16202a; padding: 8px; }}
    th {{ background: #e7ece8; color: {_BRAND}; padding: 4px 8px;
          text-align: left; }}
    td {{ padding: 4px 8px; }}
    blockquote {{ color: #4a5560; font-style: italic; }}
"""
