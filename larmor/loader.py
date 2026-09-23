"""One place that turns any source into (ppm, amplitude, recipe).

Shared by the desktop app, the CLI and the figure studio so they cannot
drift apart. Applies the recipe's stored processing pipeline, which is what
makes a saved fit reproducible end to end: reopening a recipe re-derives the
exact spectrum it was fitted against, from the untouched instrument files.
The sample of a Bruker source is the sample folder without its date / rotor /
operator tokens (or the title's ``Sample …`` line for EXPNO-per-sample
layouts, see ``larmor.io.scan.sample_name``); the title's first line is kept
in the recipe's provenance.
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

    # a chain starting with a time-domain op needs the raw fid; one that first
    # goes through ift (e.g. [hilbert, ift, em, ft], a re-apodized 1r/CSV)
    # replays from the processed arrays -- one rule, shared with the panel
    needs_raw = proc.chain_start_domain(ops) == "time"
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
        src = recipe.source_path
        # a fit made on pdata/2 reopens pdata/2, not pdata/1 (the recipe's
        # acquisition block records the procno; source_path stays the EXPNO
        # because the raw-fid replay and batchfit.align_result key on it)
        try:
            procno = int((recipe.acquisition or {}).get("procno") or 1)
            if procno != 1:
                from larmor.io import bruker as _bruker

                if _bruker.is_expno(src) and (Path(src) / "pdata" / str(procno)).is_dir():
                    src = str(Path(src) / "pdata" / str(procno))
        except (TypeError, ValueError):
            pass
        ppm, amp, fresh, meta, warnings = load_any(src, replay=False)
        # reload checks: the data file's hash, the SR, recorded read-out
        # changes -- status-bar lines through the existing warnings channel
        from larmor import provenance

        stored = recipe.to_dict()
        warnings = (list(warnings) + provenance.verify_source(stored, fresh)
                    + provenance.version_notes(stored))
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

    if p.suffix.lower() in (".csv", ".txt", ".dat"):
        from larmor.io import spectra

        ppm, amp, meta = spectra.read_csv(p)
        from larmor import provenance

        recipe = Recipe(
            sample=meta.get("sample") or p.stem, source_kind="csv",
            source_path=str(p), nucleus=meta.get("nucleus", ""),
            larmor_frequency_MHz=float(meta.get("larmor_MHz", 0.0) or 0.0),
            spin_rate_Hz=float(meta.get("spin_rate_Hz", 0.0) or 0.0),
            source_sha256=provenance.source_sha256(p))
        return ppm, amp, recipe.to_dict(), f"spectrum {p.name}", []

    from larmor.io import varian
    if varian.is_varian(p):
        ppm, amp, meta = varian.read_spectrum(p)
        from larmor import provenance

        recipe = Recipe(
            sample=meta.get("title") or varian._fid_dir(p).name,
            source_kind="varian", source_path=str(varian._fid_dir(p)),
            nucleus=meta.get("nucleus", ""),
            larmor_frequency_MHz=float(meta.get("larmor_MHz", 0.0) or 0.0),
            spin_rate_Hz=0.0, sr_hz=float(meta.get("sr_hz", 0.0) or 0.0),
            source_sha256=provenance.source_sha256(varian._fid_dir(p) / "fid"))
        return (ppm, amp, recipe.to_dict(),
                f"Varian {meta.get('nucleus', '')} (default EM+FT+phase)",
                ["Varian import: a default EM+FT+phase was applied — "
                 "re-process / re-phase as needed"])

    from larmor.io import bruker, scan

    # any Bruker path: a 1r/2rr/fid/ser file, a pdata folder, or an EXPNO
    try:
        ref = bruker.resolve(p)
    except (ValueError, FileNotFoundError):
        ref = None
    if ref is not None:
        data = bruker.read(p)
        if data.ndim == 2:
            what = "arrayed/relaxation" if data.is_pseudo2d else "2D"
            raise ValueError(
                f"this is a {what} dataset — open it with Tools ▸ 2D MQMAS "
                "viewer (or Tools ▸ Relaxation for a series)")
        if data.domain == "time":
            raise ValueError(
                "this is a raw FID — open it with File ▸ Open FID… to process "
                "it (apodize, zero-fill, phase) before the Fourier transform")
        ppm = data.axes[0].values
        amp = data.data
        title = data.meta.get("title", "")
        # a MAS rate confirmed once for this session / rotor / nucleus (and
        # the same three source values) applies here too -- the desktop,
        # Batch fit, Sequential fit and the CLI all load through this branch
        from larmor import masrate

        mas = masrate.resolve_for_load(data.meta, ref.expno)
        # the sample of a Bruker source is the sample folder without its
        # date / rotor / operator tokens (or the title's "Sample …" line for
        # EXPNO-per-sample layouts), never the title's pulse note; the
        # title's first line and the raw folder survive in the provenance
        name = scan.sample_name(Path(ref.expno), title)
        # the flat acquisition record (acqus / procs / title / uxnmr.info /
        # booking sidecar) and the SHA-256 of the exact 1r the fit sees
        from larmor import acquisition, provenance

        block = acquisition.read_block(ref.expno, ref.procno)
        recipe = Recipe(
            sample=name.key,
            source_kind="bruker", source_path=str(ref.expno),
            nucleus=data.nucleus, larmor_frequency_MHz=data.meta["larmor_MHz"],
            spin_rate_Hz=mas["spin_rate_Hz"],
            mas_uncertain=mas["mas_uncertain"],
            sr_hz=data.meta.get("sr_hz", 0.0),
            provenance={"mas_rate": mas["provenance"]},
            source_sha256=provenance.source_sha256(
                Path(ref.expno) / "pdata" / str(ref.procno) / "1r"),
            acquisition=block,
        )
        if (title or "").strip():
            recipe.provenance["title"] = name.title_first
            recipe.provenance["sample_folder"] = name.folder
        warns = list(data.warnings) + ([mas["note"]] if mas["note"] else [])
        return ppm, amp, recipe.to_dict(), data.summary, warns

    raise ValueError(f"unrecognized source: {p} (expected .fxmla, "
                     ".recipe.json, a Bruker 1r/2rr/fid/ser file, or an "
                     "EXPNO / pdata folder)")
