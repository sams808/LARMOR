"""``python -m larmor.distcheck [--gui] [--quick] [--only a,b]``.

Guarded: the pool stage spawns worker processes, and on Windows a spawned
worker re-imports the parent's ``__main__`` -- without the guard every
worker would start a distribution check of its own (and its workers
theirs). ``packaging/launcher.py`` carries the same guard for the exe.
"""
from larmor.distcheck import run

if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    raise SystemExit(run())
