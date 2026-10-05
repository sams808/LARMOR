"""Distribution self-check: does THIS installation carry and run everything?

    LARMOR.exe --distcheck [--gui] [--quick]     the frozen app
    python -m larmor.distcheck [--gui]            a development install

``larmor.selftest`` answers one question ("does Fit work here?"). This runs
the whole capability list the way a release must: every module the source
imports anywhere (the manifest, :mod:`larmor.distcheck.manifest`), every
bundled resource, every engine on synthetic data, every export format, the
process pool, the command line, and -- with ``--gui`` -- every window,
menu action and tool of the desktop, offscreen. It writes
``~/LARMOR_distcheck.log`` (``LARMOR_DISTCHECK_LOG`` overrides) and a JSON
twin next to it, and exits 0 only when every stage passed.

Stages are plain callables ``fn(say, ctx)``: ``say(*parts)`` logs a line,
``ctx`` is ``{"out": Path (a scratch folder), "examples": Path | None (the
bundled example data), "quick": bool, "gui": bool, "frozen": bool}``; a
stage passes by returning, fails by raising (the traceback goes to the log).
``core_stages`` holds the environment / imports / resources / exports / pool
/ CLI stages, ``engines`` the Qt-free capabilities (with ``engine_records``,
``engine_specialty`` and ``engine_figures``), ``gui_tools`` the desktop
workflows, ``gui_workflows`` the tool dialogs and ``gui_menus`` the menu
sweep (offscreen). Unless ``--real-settings`` is given, a run never
touches the user's own stores: the preferences go to a scratch .ini
(``LARMOR_SETTINGS_FILE``), %LOCALAPPDATA% and each store's override point
into a sandbox, and ``run`` puts the environment back on the way out.
"""
from __future__ import annotations

import datetime
import json
import os
import shutil
import sys
import tempfile
import time
import traceback
from pathlib import Path

__all__ = ["run", "log_path", "all_stages"]


def log_path() -> Path:
    return Path(os.environ.get("LARMOR_DISTCHECK_LOG")
                or (Path.home() / "LARMOR_distcheck.log"))


class _Log:
    def __init__(self, path: Path):
        self.path = path
        self.f = open(path, "a", encoding="utf-8")

    def __call__(self, *parts):
        line = " ".join(str(p) for p in parts)
        self.f.write(line + "\n")
        self.f.flush()
        if sys.stdout is not None:
            try:
                print(line, flush=True)
            except Exception:                            # noqa: BLE001
                pass

    def close(self):
        self.f.close()


def all_stages(gui: bool) -> list:
    """``[(name, fn)]`` in run order."""
    from larmor.distcheck import core_stages, engines

    out = list(core_stages.STAGES)
    out += list(engines.stages())
    if gui:
        # the workflows first, the dialogs next, the menu sweep last: it
        # triggers every action and leaves the most state behind
        from larmor.distcheck import gui_menus, gui_tools, gui_workflows
        out += list(gui_tools.stages())
        out += list(gui_workflows.stages())
        out += list(gui_menus.stages())
    return out


def _examples_dir() -> Path | None:
    from larmor.distcheck.manifest import resource_root
    p = resource_root() / "examples"
    return p if (p / "pCABS2-4").is_dir() else None


def run(argv=None) -> int:
    """Run the check (``argv`` as on the command line; None = sys.argv);
    returns the exit code. The sandbox variables it sets -- the scratch
    preferences, LOCALAPPDATA, the per-store overrides -- are undone on the
    way out, so an in-process caller (the test suite) gets its environment
    back exactly as it was."""
    env_before = dict(os.environ)
    try:
        return _run(argv)
    finally:
        for k in [k for k in os.environ if k not in env_before]:
            del os.environ[k]
        for k, v in env_before.items():
            if os.environ.get(k) != v:
                os.environ[k] = v


def _run(argv) -> int:
    import argparse

    ap = argparse.ArgumentParser(prog="larmor distcheck", add_help=True,
                                 description=__doc__.split("\n\n")[1])
    ap.add_argument("--gui", action="store_true", help="also sweep the desktop (offscreen)")
    ap.add_argument("--no-gui", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--quick", action="store_true", help="the short form of every stage")
    ap.add_argument("--only", default="", help="comma-separated stage names to run")
    ap.add_argument("--skip", default="", help="comma-separated stage names to skip")
    ap.add_argument("--out", default="", help="keep the scratch outputs in this folder")
    ap.add_argument("--list", action="store_true", help="list the stages and exit")
    ap.add_argument("--real-settings", action="store_true",
                    help="use the user's real preferences and %%LOCALAPPDATA%%\\LARMOR "
                         "stores (default: scratch copies, so nothing of yours is touched)")
    args = ap.parse_args([a for a in (argv if argv is not None else sys.argv[1:])])
    gui = bool(args.gui) and not args.no_gui
    if not args.real_settings:
        # a check run on a user's machine must leave every one of their
        # stores as it found it: the preferences go to a scratch .ini (the
        # GUI stages change themes, recent files, options), and everything
        # LARMOR keeps under %LOCALAPPDATA%\LARMOR -- the crash-recovery
        # session, the kernel cache, aliases, the rename / referencing logs,
        # MAS confirmations, SR overrides -- to a scratch LOCALAPPDATA
        sandbox = Path(tempfile.mkdtemp(prefix="larmor_distcheck_user_"))
        if not os.environ.get("LARMOR_SETTINGS_FILE"):
            os.environ["LARMOR_SETTINGS_FILE"] = str(sandbox / "settings.ini")
        os.environ.setdefault("LARMOR_REAL_LOCALAPPDATA", os.environ.get("LOCALAPPDATA", ""))
        os.environ["LOCALAPPDATA"] = str(sandbox / "LocalAppData")
        (sandbox / "LocalAppData").mkdir(parents=True, exist_ok=True)
        for var, name in (("LARMOR_ALIASES", "aliases.json"),
                          ("LARMOR_RENAME_LOG", "rename_log.jsonl"),
                          ("LARMOR_MAS_LOG", "mas_confirmations.jsonl"),
                          ("LARMOR_REF_LOG", "referencing_log.jsonl"),
                          ("LARMOR_SR_OVERRIDES", "sr_overrides.json")):
            os.environ[var] = str(sandbox / name)
        os.environ.setdefault("LARMOR_NO_SESSION", "1")
    if gui:
        # before any Qt import: no window ever shows, native dialogs become
        # widgets the sweep can dismiss, nothing of the user's session is read
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        os.environ.setdefault("LARMOR_NO_SESSION", "1")
        os.environ.setdefault("LARMOR_NO_KERNEL_WARM", "1")
    stages = all_stages(gui)
    if args.list:
        for name, _fn in stages:
            print(name)
        return 0
    only = {s.strip() for s in args.only.split(",") if s.strip()}
    skip = {s.strip() for s in args.skip.split(",") if s.strip()}
    if only:
        stages = [(n, f) for n, f in stages if n in only]
    if skip:
        stages = [(n, f) for n, f in stages if n not in skip]

    say = _Log(log_path())
    keep = bool(args.out)
    out = Path(args.out) if args.out else Path(tempfile.mkdtemp(prefix="larmor_distcheck_"))
    out.mkdir(parents=True, exist_ok=True)
    ctx = {"out": out, "examples": _examples_dir(), "quick": bool(args.quick),
           "gui": gui, "frozen": bool(getattr(sys, "frozen", False))}
    import larmor
    say(f"=== LARMOR distribution check {datetime.datetime.now().isoformat(timespec='seconds')}"
        f" · larmor {larmor.__version__} · frozen {ctx['frozen']} · gui {gui}"
        f" · quick {ctx['quick']} · {len(stages)} stages · scratch {out}")
    report = {"version": larmor.__version__, "frozen": ctx["frozen"], "gui": gui,
              "started": datetime.datetime.now().isoformat(timespec="seconds"),
              "stages": []}
    failed = []
    t_all = time.perf_counter()
    for name, fn in stages:
        say(f"--- {name}")
        notes: list = []

        def note(*parts, _notes=notes):
            line = " ".join(str(p) for p in parts)
            _notes.append(line)
            say("   ", line)

        t0 = time.perf_counter()
        ok = True
        err = ""
        try:
            fn(note, ctx)
        except BaseException as exc:                     # noqa: BLE001
            ok = False
            err = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
            say(err)
            failed.append(name)
            if isinstance(exc, KeyboardInterrupt):
                break
        dt = time.perf_counter() - t0
        if gui and name.split("-")[0] in ("tool", "dialog", "menu"):
            # every GUI stage starts from an empty application: windows a
            # stage left behind would be restyled by every later theme
            # switch. A worker thread a stage left running is waited for
            # first (and named); one that never ends fails the stage.
            try:
                from larmor.distcheck.gui_common import destroy_all_windows
                n_win, waited, still = destroy_all_windows()
                if waited:
                    note(f"worker thread(s) still running when the stage ended: "
                         f"{', '.join(waited)}" + ("" if still else " (finished after)"))
                if still:
                    ok = False
                    err = err or (f"worker thread(s) still running 30 s after the stage: "
                                  f"{', '.join(still)}")
                    failed.append(name) if name not in failed else None
                if n_win:
                    say(f"    {n_win} window(s) left by {name} destroyed")
            except Exception as exc:                     # noqa: BLE001
                say(f"    window clean-up after {name} failed: {exc!r}")
        say(f"{'ok  ' if ok else 'FAIL'} {name} ({dt:.1f} s)")
        report["stages"].append({"name": name, "ok": ok, "seconds": round(dt, 2),
                                 "notes": notes, "error": err})
    total = time.perf_counter() - t_all
    verdict = (f"=== PASS: {len(stages)} stages in {total:.0f} s" if not failed else
               f"=== FAIL: {len(failed)} of {len(stages)} stages failed in {total:.0f} s: "
               + ", ".join(failed))
    say(verdict)
    say(f"(log: {say.path})")
    report["verdict"] = verdict
    report["failed"] = failed
    try:
        say.path.with_suffix(".json").write_text(json.dumps(report, indent=1),
                                                 encoding="utf-8")
    except OSError:
        pass
    say.close()
    if not keep:
        shutil.rmtree(out, ignore_errors=True)
    return 0 if not failed else 1
