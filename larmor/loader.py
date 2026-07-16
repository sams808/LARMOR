"""One place that turns any source into (ppm, amplitude, recipe).

Shared by the desktop app, the CLI and the figure studio so they cannot
drift apart. Applies the recipe's stored processing pipeline, which is what
makes a saved fit reproducible end to end: reopening a recipe re-derives the
exact spectrum it was fitted against, from the untouched instrument files.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from larmor.recipe import Recipe


def apply_processing(recipe: Recipe, ppm: np.ndarray, amp: np.ndarray,
                     source_path: str | None = None,
                     ) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Replay recipe.processing. Returns (ppm, amp, notes)."""
    ops = list(recipe.processing or [])
    if not ops:
        return ppm, amp, []
    from larmor import processing as proc
    from larmor.io import bruker

    needs_raw = any(o.get("op") in ("em", "gm", "sine", "traf", "tdeff",
                                    "fcor", "zf", "ft", "lp", "shift_fid",
                                    "swap_echo", "echo_apodize")
                    for o in ops)
    notes = []
    if needs_raw:
        if not (source_path and bruker.is_expno(Path(source_path))):
            raise ValueError(
                "this recipe's processing starts from the raw fid, but the "
                f"source is not a Bruker EXPNO: {source_path}")
        s = proc.from_bruker_fid(source_path)
    else:
        s = proc.from_processed(ppm, amp, recipe.larmor_frequency_MHz or 0.0)
    s = proc.apply(s, ops)
    if s.domain != "freq":
        raise ValueError("recipe processing must end in the frequency domain "
                         "(add an 'ft' step)")
    order = np.argsort(s.x_ppm)
    notes.append(f"replayed {len(ops)} processing step(s) from the recipe")
    return np.asarray(s.x_ppm)[order], s.y.real[order], notes


def load_any(path: str | Path, replay: bool = True):
    """Load a .fxmla, a LARMOR .recipe.json, or a Bruker EXPNO.

    Returns (ppm, amp, recipe_dict, meta, warnings).
    """
    p = Path(path)

    if p.suffix.lower() == ".json":
        recipe = Recipe.load(p)
        if not recipe.source_path or not Path(recipe.source_path).exists():
            raise ValueError(
                f"recipe's source data not found: {recipe.source_path}")
        ppm, amp, _, meta, warnings = load_any(recipe.source_path, replay=False)
        if replay and recipe.processing:
            try:
                ppm, amp, notes = apply_processing(recipe, ppm, amp,
                                                   recipe.source_path)
                warnings = list(warnings) + notes
            except Exception as exc:
                warnings = list(warnings) + [f"processing replay failed: {exc}"]
        return ppm, amp, recipe.to_dict(), f"recipe {p.name} | {meta}", warnings

    if p.suffix.lower() in (".fxmla", ".fxml"):
        from larmor.io import fxmla

        dm = fxmla.read(p)
        if dm.spectrum is None:
            raise ValueError("no experimental data in this fxmla")
        if dm.is_2d:
            raise ValueError(
                "this is a 2D (MQMAS) dmfit file -- open it with the 2D tools "
                "(larmor.twod); the 1D workbench cannot fit it")
        recipe, warnings = fxmla.to_recipe(dm)
        ppm, amp = dm.spectrum.ppm, dm.spectrum.amplitude
        order = np.argsort(ppm)
        return (ppm[order], amp[order], recipe.to_dict(),
                f"dmfit {dm.version} | {dm.comment}", warnings)

    from larmor.io import bruker

    if bruker.is_expno(p):
        exp = bruker.read_expno(p)
        if exp.processed is None:
            raise ValueError(
                "this EXPNO has no processed pdata; process the raw fid "
                "first (Processing panel -> raw fid)")
        ppm, amp = exp.processed_ppm, exp.processed.astype(float)
        order = np.argsort(ppm)
        recipe = Recipe(
            sample=exp.title.splitlines()[0] if exp.title else "",
            source_kind="bruker", source_path=str(p),
            nucleus=exp.nucleus, larmor_frequency_MHz=exp.sfo1_MHz,
            spin_rate_Hz=exp.masr_Hz or 0.0,
        )
        return (ppm[order], amp[order], recipe.to_dict(),
                exp.summary.splitlines()[1], exp.conflicts)

    raise ValueError(f"unrecognized source: {p} (expected .fxmla, "
                     ".recipe.json, or a Bruker EXPNO folder)")
