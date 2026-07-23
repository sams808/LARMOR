"""Prettier help rendering for the Markdown manuals.

QTextBrowser can show Markdown, but it can't typeset equations, and it drops
raw ``<img>`` tags inside ``setMarkdown`` (so their size can't be controlled).
So we go through HTML instead:

1. render every LaTeX ``$...$`` / ``$$...$$`` snippet to a crisp transparent
   PNG with matplotlib's mathtext (rendered high-res, then scaled *down* in the
   ``<img>`` so it stays sharp and sits at the size of the surrounding text),
2. let Qt turn the Markdown into HTML, then splice the sized ``<img>`` tags in,
3. supply a stylesheet so the manuals read like a real document.

No web engine, no extra dependencies — matplotlib is already a core dependency.
"""
from __future__ import annotations

import base64
import io
import re

# body text / heading colours (light theme; the help dialog has a white page)
_INK = "#22303a"
_BRAND = "#1f3a5f"

# matplotlib renders every snippet at this point-size/DPI, then the <img> is
# scaled down to the target on-screen size below — high render res => crisp.
_RENDER_PT = 16.0
_RENDER_DPI = 220
_FONT_PX = _RENDER_PT * _RENDER_DPI / 72.0     # glyph px in the rendered PNG
_INLINE_PX = 15.0                              # on-screen size of inline math
_DISPLAY_PX = 20.0                             # on-screen size of display math

_math_cache: dict = {}


def _math_png(tex: str, color: str = _INK):
    """Render ``tex`` to a transparent PNG; return (data-uri, width, height)."""
    key = (tex, color)
    if key in _math_cache:
        return _math_cache[key]
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    fig = Figure()
    FigureCanvasAgg(fig)
    fig.patch.set_alpha(0.0)
    fig.text(0.0, 0.0, f"${tex}$", fontsize=_RENDER_PT, color=color)
    buf = io.BytesIO()
    try:
        fig.savefig(buf, dpi=_RENDER_DPI, transparent=True,
                    bbox_inches="tight", pad_inches=0.02)
    except Exception:
        # a malformed expression shouldn't blow up the whole manual
        result = ("", 0, 0)
        _math_cache[key] = result
        return result
    data = buf.getvalue()
    w = int.from_bytes(data[16:20], "big")     # PNG IHDR width  (px)
    h = int.from_bytes(data[20:24], "big")     # PNG IHDR height (px)
    uri = "data:image/png;base64," + base64.b64encode(data).decode("ascii")
    result = (uri, w, h)
    _math_cache[key] = result
    return result


def _img_tag(tex: str, target_px: float, *, inline: bool) -> str:
    uri, w, h = _math_png(tex)
    if not uri:
        return ""
    scale = target_px / _FONT_PX
    dw, dh = max(1, round(w * scale)), max(1, round(h * scale))
    style = ' style="vertical-align:middle"' if inline else ""
    return f'<img src="{uri}" width="{dw}" height="{dh}"{style} />'


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


def render_help_html(md: str) -> str:
    """Markdown (with LaTeX math) -> styled HTML with typeset equations."""
    from PySide6.QtGui import QTextDocument

    md, blocks = _protect_code(md)

    inline_map: dict[str, str] = {}
    display_map: dict[str, str] = {}

    def display(m):
        tag = _img_tag(m.group(1).strip(), _DISPLAY_PX, inline=False)
        if not tag:
            return m.group(0)
        tok = f"MDMATHD{len(display_map)}X"
        display_map[tok] = tag
        return f"\n\n{tok}\n\n"

    def inline(m):
        tag = _img_tag(m.group(1).strip(), _INLINE_PX, inline=True)
        if not tag:
            return m.group(0)
        tok = f"MDMATHI{len(inline_map)}X"
        inline_map[tok] = tag
        return tok

    md = re.sub(r"\$\$(.+?)\$\$", display, md, flags=re.DOTALL)
    md = re.sub(r"(?<!\$)\$(?!\$)(.+?)(?<!\$)\$(?!\$)", inline, md,
                flags=re.DOTALL)
    md = _restore_code(md, blocks)

    doc = QTextDocument()
    doc.setDefaultStyleSheet(HELP_CSS)
    doc.setMarkdown(md, QTextDocument.MarkdownDialectGitHub)
    html = doc.toHtml()

    # centre each display equation (it sits alone in its own paragraph)
    def center(m):
        return (f'<p align="center" style="margin-top:12px;margin-bottom:12px">'
                f'{display_map[m.group(1)]}</p>')

    html = re.sub(r"<p[^>]*>\s*(MDMATHD\d+X)\s*</p>", center, html)
    for tok, tag in display_map.items():          # any that weren't alone
        html = html.replace(tok, tag)
    for tok, tag in inline_map.items():
        html = html.replace(tok, tag)
    return html


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
