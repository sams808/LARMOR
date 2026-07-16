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
    from fastapi.staticfiles import StaticFiles
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
app.mount("/static", StaticFiles(directory=str(_STATIC)), name="static")

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

    if path.suffix.lower() == ".json":
        # a saved LARMOR recipe: load it AND its referenced source data
        recipe = Recipe.load(path)
        if not recipe.source_path or not Path(recipe.source_path).exists():
            raise HTTPException(
                422, f"recipe's source data not found: {recipe.source_path}")
        inner = load(LoadRequest(path=recipe.source_path))
        inner["recipe"] = recipe.to_dict()
        inner["meta"] = f"recipe {path.name} | " + inner["meta"]
        # the client will reference the .json path in later calls
        _DATA[_norm(path)] = _DATA[_norm(recipe.source_path)]
        return inner

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
        x, total, per_site = engine.simulate(recipe, exp_ppm=ppm)
    except (fitmod.ConstraintError, ValueError) as exc:
        raise HTTPException(422, str(exc))
    return {
        "x": x.tolist(),
        "total": total.tolist(),
        "sites": [y.tolist() for y in per_site],
        "labels": [s.label or s.model for s in recipe.sites],
    }


@app.post("/api/fit")
def run_fit(req: FitRequest):
    recipe = Recipe.from_dict(req.recipe)
    ppm, amp = _get_data(req.source_path)
    try:
        result = fitmod.fit(recipe, ppm, amp, window_ppm=req.window)
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


class BrowseRequest(BaseModel):
    path: str = ""


class QuantifyRequest(BaseModel):
    recipe: dict
    window: tuple[float, float] | None = None


class ProcessRequest(BaseModel):
    source_path: str
    ops: list[dict]
    use_raw: bool = False        # start from the raw fid instead of pdata


@app.get("/api/models")
def list_models():
    from larmor import models as model_registry

    return {"models": model_registry.describe_all()}


@app.post("/api/browse")
def browse(req: BrowseRequest):
    """Minimal file browser: folders, EXPNOs, and openable files."""
    import string

    if not req.path:
        drives = [f"{d}:\\" for d in string.ascii_uppercase
                  if Path(f"{d}:\\").exists()]
        return {"path": "", "parent": None, "dirs": drives, "expnos": [],
                "files": []}
    p = Path(req.path)
    if not p.is_dir():
        raise HTTPException(404, f"not a folder: {p}")
    dirs, expnos, files = [], [], []
    try:
        for child in sorted(p.iterdir(), key=lambda c: c.name.lower()):
            try:
                if child.is_dir():
                    if (child / "acqus").exists():
                        expnos.append(child.name)
                    else:
                        dirs.append(child.name)
                elif child.suffix.lower() in (".fxmla", ".fxml", ".json"):
                    files.append(child.name)
            except OSError:
                continue
    except PermissionError:
        raise HTTPException(403, f"no permission to read {p}")
    return {"path": str(p), "parent": str(p.parent) if p.parent != p else None,
            "dirs": dirs, "expnos": expnos, "files": files}


@app.post("/api/quantify")
def run_quantify(req: QuantifyRequest):
    from larmor.quantify import quantify

    recipe = Recipe.from_dict(req.recipe)
    if not recipe.sites:
        raise HTTPException(422, "no sites to quantify")
    try:
        return quantify(recipe, window_ppm=req.window)
    except Exception as exc:
        raise HTTPException(422, f"quantification failed: {exc}")


@app.post("/api/process")
def process(req: ProcessRequest):
    """Apply a processing pipeline; result replaces the working spectrum."""
    from larmor import processing as proc

    key = _norm(req.source_path)
    try:
        if req.use_raw:
            s = proc.from_bruker_fid(req.source_path)
        else:
            ppm, amp = _get_data(req.source_path)
            # recover SFO1 for phase pivoting when available
            sfo1 = 0.0
            from larmor.io import bruker

            if bruker.is_expno(Path(req.source_path)):
                sfo1 = bruker.read_expno(req.source_path).sfo1_MHz
            s = proc.from_processed(ppm, amp, sfo1)
        s = proc.apply(s, req.ops)
    except Exception as exc:
        raise HTTPException(422, f"processing failed: {exc}")
    if s.domain != "freq":
        raise HTTPException(422, "pipeline must end in the frequency domain "
                                 "(add an 'ft' step)")
    x, y = np.asarray(s.x_ppm), s.y.real
    order = np.argsort(x)
    x, y = x[order], y[order]
    _DATA[key] = (x, y)
    return {"ppm": x.tolist(), "amp": y.tolist()}


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
