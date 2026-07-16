"""LARMOR local web app: dmfit-style interactive fitting in the browser.

Run with `larmor app` (or `python -m larmor.app`), then open the printed URL.
Single-user, local-only by default (binds 127.0.0.1).

The precomputed Czjzek kernel is what makes live parameter adjustment
possible: after the one-time kernel build, every re-simulation costs
milliseconds.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

try:
    from fastapi import FastAPI, HTTPException
    from fastapi.responses import HTMLResponse
    from pydantic import BaseModel
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "the LARMOR app needs the 'app' extra: pip install fastapi uvicorn"
    ) from exc

from larmor import engine
from larmor import fit as fitmod
from larmor.recipe import Recipe

app = FastAPI(title="LARMOR")

_STATIC = Path(__file__).parent / "static"

# per-process spectrum store so /simulate and /fit don't re-read the source
_DATA: dict[str, tuple[np.ndarray, np.ndarray]] = {}


def _norm(path: str | Path) -> str:
    """One canonical key per file, whatever slashes the client sent."""
    return str(Path(path).resolve())


class LoadRequest(BaseModel):
    path: str


class SimulateRequest(BaseModel):
    recipe: dict
    source_path: str


class FitRequest(BaseModel):
    recipe: dict
    source_path: str
    window: tuple[float, float] | None = None


class SaveRequest(BaseModel):
    recipe: dict
    path: str


class FigureTemplateRequest(BaseModel):
    source_path: str
    recipe: dict | None = None


class FigurePreviewRequest(BaseModel):
    spec: dict


class FigureExportRequest(BaseModel):
    spec: dict
    path: str                     # base path, extension added per format
    formats: list[str] = ["png", "svg", "pdf"]
    dpi: int = 600


@app.get("/", response_class=HTMLResponse)
def index() -> str:
    return (_STATIC / "index.html").read_text(encoding="utf-8")


@app.post("/api/load")
def load(req: LoadRequest):
    path = Path(req.path)
    if not path.exists():
        raise HTTPException(404, f"not found: {path}")

    if path.suffix.lower() == ".fxmla":
        from larmor.io import fxmla

        dm = fxmla.read(path)
        if dm.spectrum is None or dm.is_2d:
            raise HTTPException(422, "no 1D experimental data in this fxmla")
        recipe, warnings = fxmla.to_recipe(dm)
        ppm, amp = dm.spectrum.ppm, dm.spectrum.amplitude
        meta = f"dmfit {dm.version} | {dm.comment}"
    else:
        from larmor.io import bruker

        if not bruker.is_expno(path):
            raise HTTPException(422, f"not a .fxmla file or Bruker EXPNO folder: {path}")
        exp = bruker.read_expno(path)
        if exp.processed is None:
            raise HTTPException(422, "EXPNO has no processed pdata/1 to display")
        ppm, amp = exp.processed_ppm, exp.processed.astype(float)
        recipe = Recipe(
            sample=exp.title.splitlines()[0] if exp.title else "",
            source_kind="bruker", source_path=str(path),
            nucleus=exp.nucleus, larmor_frequency_MHz=exp.sfo1_MHz,
            spin_rate_Hz=exp.masr_Hz or 0.0,
        )
        warnings = exp.conflicts
        meta = exp.summary.splitlines()[1]

    order = np.argsort(ppm)
    ppm, amp = np.asarray(ppm)[order], np.asarray(amp)[order]
    _DATA[_norm(path)] = (ppm, amp)
    return {
        "recipe": recipe.to_dict(),
        "warnings": list(warnings),
        "meta": meta,
        "ppm": ppm.tolist(),
        "amp": amp.tolist(),
    }


def _get_data(source_path: str) -> tuple[np.ndarray, np.ndarray]:
    key = _norm(source_path)
    if key not in _DATA:
        load(LoadRequest(path=source_path))
    return _DATA[key]


def _kernel_for(recipe: Recipe, exp_ppm: np.ndarray):
    if engine.needs_kernel(recipe):
        return engine.build_kernel(
            recipe.nucleus, recipe.larmor_frequency_MHz, recipe.spin_rate_Hz)
    return engine.Axis(x_ppm=exp_ppm)


@app.post("/api/simulate")
def simulate(req: SimulateRequest):
    recipe = Recipe.from_dict(req.recipe)
    ppm, _ = _get_data(req.source_path)
    # resolve constraint expressions (e.g. "0.5 * s0.amplitude") so linked
    # parameters draw at their derived values, and bad ones fail here with a
    # clear message rather than mid-fit
    try:
        params = fitmod._make_params(recipe)
        fitmod._apply_params(recipe, params)
    except fitmod.ConstraintError as exc:
        raise HTTPException(422, str(exc))
    kernel = _kernel_for(recipe, ppm)
    per_site = [engine.simulate_site(s, kernel) for s in recipe.sites]
    total = np.sum(per_site, axis=0) if per_site else np.zeros_like(kernel.x_ppm)
    return {
        "x": kernel.x_ppm.tolist(),
        "total": total.tolist(),
        "sites": [y.tolist() for y in per_site],
        "labels": [s.label or s.model for s in recipe.sites],
    }


@app.post("/api/fit")
def run_fit(req: FitRequest):
    recipe = Recipe.from_dict(req.recipe)
    ppm, amp = _get_data(req.source_path)
    kernel = _kernel_for(recipe, ppm)
    try:
        result = fitmod.fit(recipe, ppm, amp, window_ppm=req.window, kernel=kernel)
    except fitmod.ConstraintError as exc:
        raise HTTPException(422, str(exc))
    return {
        "recipe": recipe.to_dict(),
        "rmsd": result.rmsd,
        "report": result.report,
        "frozen": result.frozen_sites or [],
        "at_bounds": result.at_bounds or [],
        "x": result.x_ppm.tolist(),
        "total": result.y_fit.tolist(),
        "sites": [y.tolist() for y in result.per_site],
        "labels": [s.label or s.model for s in recipe.sites],
    }


@app.post("/api/save")
def save(req: SaveRequest):
    target = Path(req.path)
    # data-protection guard: never write into an instrument data folder
    _guard_instrument_dir(target)
    Recipe.from_dict(req.recipe).save(target)
    return {"saved": str(target)}


def _guard_instrument_dir(target: Path) -> None:
    for parent in [target.parent, *target.parent.parents]:
        if (parent / "acqus").exists() or (parent / "fid").exists() or (parent / "ser").exists():
            raise HTTPException(
                403, f"refusing to write inside instrument data folder: {parent}")


@app.post("/api/figure/template")
def figure_template(req: FigureTemplateRequest):
    """Build a sensible starter spec for whatever the source offers."""
    path = Path(req.source_path)
    templates: dict[str, dict] = {}

    is_2d = (path / "acqu2s").exists()
    pdata1 = path / "pdata" / "1"
    if path.suffix.lower() == ".fxmla" or (path.is_dir() and not is_2d):
        traces = [{"path": str(path), "label": "experiment",
                   "color": "black", "linewidth": 0.9}]
        templates["1d"] = {"kind": "1d", "style": "article-wide",
                           "traces": traces}
    if req.recipe and req.recipe.get("sites"):
        # note: the recipe must be SAVED for figure traces to reference it
        templates["1d-fit"] = {
            "kind": "1d", "style": "article-wide",
            "traces": [
                {"path": str(path), "label": "experiment", "color": "black",
                 "linewidth": 0.9},
                {"recipe": "<save recipe first, then put its path here>",
                 "part": "total", "label": "fit", "color": "crimson"},
            ],
        }
    if is_2d:
        templates["2d"] = {
            "kind": "2d", "style": "thesis", "path": str(path),
            "xlabel": "F2 shift (ppm)", "ylabel": "F1 (ppm)",
            "levels": {"mode": "log", "n": 12},
        }
    if (pdata1 / "t1ints.txt").exists():
        templates["satrec"] = {"kind": "series", "mode": "satrec",
                               "style": "article", "path": str(path),
                               "stretched": True}
    if (pdata1 / "redor.txt").exists():
        templates["redor"] = {"kind": "series", "mode": "redor",
                              "style": "article", "path": str(path)}
    return {"templates": templates}


@app.post("/api/figure/preview")
def figure_preview(req: FigurePreviewRequest):
    import base64

    from larmor import figures

    try:
        png = figures.render_png_bytes(req.spec, dpi=110)
    except Exception as exc:
        raise HTTPException(422, f"figure failed: {exc}")
    return {"png_base64": base64.b64encode(png).decode("ascii")}


@app.post("/api/figure/export")
def figure_export(req: FigureExportRequest):
    from larmor import figures

    target = Path(req.path)
    _guard_instrument_dir(target)
    try:
        saved = figures.export(req.spec, target,
                               formats=tuple(req.formats), dpi=req.dpi)
    except Exception as exc:
        raise HTTPException(422, f"figure failed: {exc}")
    return {"saved": saved}


def serve(host: str = "127.0.0.1", port: int = 8642,
          open_browser: bool = False) -> None:  # pragma: no cover
    import uvicorn

    url = f"http://{host}:{port}"
    print(f"LARMOR app: {url}")
    if open_browser:
        import threading
        import webbrowser

        threading.Timer(1.2, lambda: webbrowser.open(url)).start()
    uvicorn.run(app, host=host, port=port, log_level="warning")


if __name__ == "__main__":  # pragma: no cover
    serve()
