"""Shared CPU-parallel map for embarrassingly-parallel batch work.

Monte-Carlo error trials, chi-square profile scan points, and batch-fit
spectra are all the same shape of problem: many small, fully independent
nonlinear fits whose results don't depend on each other. That is a textbook
fit for OS-level process parallelism -- and specifically NOT a fit for a
GPU: each unit of work is `lmfit`'s serial Levenberg-Marquardt loop running
Python-level control flow over a spectrum of a few thousand points, not a
large batched tensor/matrix operation. Threads don't help either (the GIL
stays held through nearly every line of that loop, since almost none of it
is a single big released-GIL numpy/BLAS call) -- only separate processes
give real parallelism here.

Windows note: process pools use the "spawn" start method, which re-imports
this module (and pickles every argument) in each worker -- so the mapped
function must be a plain, importable, module-level function, and every
argument must be picklable (plain dicts/lists/floats/strings, numpy arrays,
JSON strings -- not a bound method, a lambda, or a live Qt object). See
`larmor/autofit.py`'s `_mc_trial_worker`/`_profile_point_worker` for the
pattern. A **frozen** (PyInstaller) build additionally needs
`multiprocessing.freeze_support()` at its entry point (see
`packaging/launcher.py`) or each worker re-launches the whole app instead
of becoming a worker.
"""
from __future__ import annotations

import os
from concurrent.futures import Future, ProcessPoolExecutor, as_completed

#: below this many items, process-pool startup overhead (a few hundred ms on
#: Windows) would swamp any benefit -- run sequentially instead. Real batch/
#: error-analysis runs are always far above this (tens to hundreds of items);
#: it mainly keeps small/test workloads from paying pool-startup cost for
#: nothing.
MIN_ITEMS_FOR_PROCESSES = 8


def default_worker_count() -> int:
    """Leave one core free for the UI/event loop -- a run that claims every
    last core makes the window itself sluggish while it's in flight."""
    n = os.cpu_count() or 1
    return max(1, n - 1)


def parallel_map(fn, items: list, *, max_workers: int | None = None,
                 should_stop=None, on_result=None, use_processes: bool = True,
                 executor: ProcessPoolExecutor | None = None) -> list:
    """Apply ``fn`` to every item in ``items``; return results in ORIGINAL
    item order (not completion order).

    A result is ``None`` for any item that raised, or that never ran because
    ``should_stop`` cut the run short -- callers decide how to treat a hole
    (a batch-fit spectrum reports its unfit values; an error estimate simply
    drops that trial/point from the statistics).

    ``on_result(index, result)`` fires as each result arrives -- in
    COMPLETION order when running across processes (not item order), so it's
    for live progress ticks ("N of M done"), not for anything that assumes
    submission order.

    ``should_stop()`` is checked after every completion; once true, every
    not-yet-DONE future is cancelled -- BEST EFFORT: ``Future.cancel()`` only
    succeeds for work a worker hasn't picked up yet, and a process pool
    prefetches items into its internal call queue ahead of that, so exactly
    how much gets cancelled depends on timing, not something callers control.
    Already-running/already-dispatched items always finish (no mid-
    computation kill -- the same "finishes its current unit of work"
    contract every other stoppable long-running op in this app already has);
    stopping mainly guarantees "don't keep queuing new work", not "stop
    instantly".

    ``use_processes=False`` runs everything sequentially in the calling
    process, in order -- the default for library-level callers (so existing
    behaviour/timing is unchanged unless a caller opts in), and always used
    for fewer than ``MIN_ITEMS_FOR_PROCESSES`` items regardless of the flag.

    ``executor``: reuse an already-running ``ProcessPoolExecutor`` instead of
    creating (and tearing down) a new one for this call -- pool startup has
    real cost (Windows spawns a fresh interpreter per worker), so a caller
    making MANY of these calls back-to-back (e.g. one profile per parameter,
    for every spectrum in a batch) should create ONE pool up front and pass
    it through every call, closing it only once, at the very end.
    """
    n = len(items)
    results: list = [None] * n
    if n == 0:
        return results

    if not use_processes or (executor is None and n < MIN_ITEMS_FOR_PROCESSES):
        for i, item in enumerate(items):
            if should_stop is not None and should_stop():
                break
            try:
                r = fn(item)
            except Exception:
                r = None
            results[i] = r
            if on_result:
                on_result(i, r)
        return results

    def _run(pool: ProcessPoolExecutor):
        futures: dict[Future, int] = {pool.submit(fn, item): i
                                      for i, item in enumerate(items)}
        stopped = False
        for fut in as_completed(futures):
            i = futures[fut]
            # a future cancelled below (not-yet-started when stop fired)
            # raises CancelledError from .result() -- treated as a hole,
            # same as anything that was never submitted at all
            try:
                r = fut.result()
            except Exception:
                r = None
            results[i] = r
            if on_result:
                on_result(i, r)
            if stopped:
                continue
            if should_stop is not None and should_stop():
                stopped = True
                for f2 in futures:
                    if not f2.done():
                        f2.cancel()

    if executor is not None:
        _run(executor)
    else:
        workers = max(1, min(max_workers or default_worker_count(), n))
        with ProcessPoolExecutor(max_workers=workers) as pool:
            _run(pool)
    return results
