"""LARMOR command line.

    larmor info PATH                     identify + summarize a data source
    larmor import PATH -o recipe.json    convert a dmfit .fxmla to a recipe
    larmor fit recipe.json [-o out.json] [--plot out.png]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


def cmd_info(args: argparse.Namespace) -> int:
    path = Path(args.path)
    if path.suffix.lower() == ".fxmla":
        from larmor.io import fxmla

        dm = fxmla.read(path)
        print(f"dmfit fit file (version {dm.version}, mode {dm.fit_mode!r})")
        if dm.comment:
            print(f"comment: {dm.comment}")
        for dim in dm.dimensions:
            print(f"dimension {dim.label}: {dim.nucleus} at {dim.frequency_MHz:.3f} MHz, "
                  f"{len(dim.lines)} lines")
            for i, line in enumerate(dim.lines):
                pos = line.params.get("pos")
                extras = []
                if "sCZ_CQ" in line.params:
                    extras.append(f"sCZ_CQ={line.params['sCZ_CQ'].value:.0f} kHz")
                if "wid" in line.params:
                    extras.append(f"wid={line.params['wid'].value:.2f} ppm")
                print(f"  [{i}] {line.model_name:<10} "
                      f"pos={pos.value if pos else float('nan'):.2f} ppm  "
                      + "  ".join(extras))
        if dm.spectrum is not None:
            h = dm.spectrum.header
            kind = "2D" if dm.is_2d else "1D"
            print(f"embedded spectrum: {kind}, NP={int(h.get('NP', 0))}"
                  + (f" x NI={int(h['NI'])}" if "NI" in h else ""))
        return 0

    from larmor.io import bruker

    if bruker.is_expno(path):
        exp = bruker.read_expno(path)
        print(exp.summary)
        return 0

    print(f"unrecognized data source: {path}", file=sys.stderr)
    return 1


def cmd_import(args: argparse.Namespace) -> int:
    from larmor.io import fxmla

    dm = fxmla.read(args.path)
    recipe, warnings = fxmla.to_recipe(dm)
    for w in warnings:
        print(f"note: {w}")
    out = Path(args.output or Path(args.path).with_suffix(".recipe.json").name)
    recipe.save(out)
    print(f"recipe written: {out}  ({len(recipe.sites)} sites, "
          f"{recipe.nucleus} at {recipe.larmor_frequency_MHz:.3f} MHz)")
    return 0


def cmd_fit(args: argparse.Namespace) -> int:
    from larmor.recipe import Recipe
    from larmor import fit as fitmod

    recipe = Recipe.load(args.recipe)

    # experimental data comes from the recipe's source reference
    if recipe.source_kind == "fxmla":
        from larmor.io import fxmla

        dm = fxmla.read(recipe.source_path)
        if dm.spectrum is None:
            print("source fxmla holds no experimental data", file=sys.stderr)
            return 1
        exp_ppm, exp_amp = dm.spectrum.ppm, dm.spectrum.amplitude
    elif recipe.source_kind == "bruker":
        from larmor.io import bruker

        exp = bruker.read_expno(recipe.source_path)
        if exp.processed is None:
            print("EXPNO has no processed pdata to fit", file=sys.stderr)
            return 1
        exp_ppm, exp_amp = exp.processed_ppm, exp.processed.astype(float)
    else:
        print(f"unknown source kind {recipe.source_kind!r}", file=sys.stderr)
        return 1

    window = tuple(args.window) if args.window else None
    print("building kernel (one-time, cached per session)...")
    result = fitmod.fit(recipe, exp_ppm, exp_amp, window_ppm=window)

    print(result.report)
    print(f"\nnormalized RMSD: {result.rmsd:.4f}")

    out = Path(args.output or args.recipe)
    recipe.save(out)
    print(f"updated recipe written: {out}")

    if args.plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        hi, lo = recipe.fit_window_ppm
        sel = (exp_ppm >= lo) & (exp_ppm <= hi)
        fig, ax = plt.subplots(figsize=(9, 4.5))
        ax.plot(exp_ppm[sel], exp_amp[sel], "k", lw=0.9, label="experiment")
        ax.plot(result.x_ppm, result.y_fit, "r", lw=1.2, alpha=0.9,
                label=f"LARMOR fit (RMSD {result.rmsd:.3f})")
        for site, ys in zip(recipe.sites, result.per_site):
            ax.plot(result.x_ppm, ys, lw=0.8, alpha=0.6, ls="--", label=site.label)
        ax.set_xlim(hi, lo)
        ax.set_xlabel(f"{recipe.nucleus} shift / ppm")
        ax.legend(fontsize=8)
        ax.set_title(recipe.sample or Path(recipe.source_path).name)
        fig.tight_layout()
        fig.savefig(args.plot, dpi=130)
        print(f"overlay written: {args.plot}")
    return 0


def cmd_app(args: argparse.Namespace) -> int:
    from larmor.app import serve

    serve(host=args.host, port=args.port)
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="larmor", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    p_info = sub.add_parser("info", help="identify and summarize a data source")
    p_info.add_argument("path")
    p_info.set_defaults(func=cmd_info)

    p_imp = sub.add_parser("import", help="convert a dmfit .fxmla to a LARMOR recipe")
    p_imp.add_argument("path")
    p_imp.add_argument("-o", "--output", help="output recipe path")
    p_imp.set_defaults(func=cmd_import)

    p_fit = sub.add_parser("fit", help="refine a recipe against its source data")
    p_fit.add_argument("recipe")
    p_fit.add_argument("-o", "--output", help="output recipe path (default: in place)")
    p_fit.add_argument("--window", nargs=2, type=float, metavar=("HI_PPM", "LO_PPM"))
    p_fit.add_argument("--plot", help="write an overlay PNG")
    p_fit.set_defaults(func=cmd_fit)

    p_app = sub.add_parser("app", help="launch the interactive web app")
    p_app.add_argument("--host", default="127.0.0.1")
    p_app.add_argument("--port", type=int, default=8642)
    p_app.set_defaults(func=cmd_app)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
