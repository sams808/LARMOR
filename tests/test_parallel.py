"""Regressions for larmor.parallel: the shared process-pool map used by
Monte-Carlo error trials, chi-square profile points, and (soon) batch-fit
spectra. Correctness matters more than speed here -- these tests exist to
prove the parallel path returns EXACTLY what the sequential path returns,
not just that it doesn't crash.

Worker functions must be plain, importable, MODULE-LEVEL functions (picklable
for Windows' spawn start method) -- see larmor/parallel.py's own docstring.
"""
import time
from concurrent.futures import ProcessPoolExecutor

from larmor import parallel


def _square(x):
    return x * x


def _maybe_fail(x):
    if x == 3:
        raise ValueError("boom")
    return x * 2


def _slow_square(x):
    time.sleep(0.15)                   # give should_stop a real window to land
    return x * x


def test_sequential_path_preserves_order_and_values():
    items = list(range(5))
    out = parallel.parallel_map(_square, items, use_processes=False)
    assert out == [0, 1, 4, 9, 16]


def test_small_item_count_stays_sequential_even_with_use_processes_true():
    # below MIN_ITEMS_FOR_PROCESSES: no pool should be needed/created
    items = list(range(3))
    out = parallel.parallel_map(_square, items, use_processes=True)
    assert out == [0, 1, 4]


def test_process_path_matches_sequential_results():
    items = list(range(20))
    seq = parallel.parallel_map(_square, items, use_processes=False)
    par = parallel.parallel_map(_square, items, use_processes=True,
                                max_workers=4)
    assert par == seq                       # SAME order, SAME values


def test_exceptions_become_a_hole_not_a_crash():
    items = list(range(10))
    out = parallel.parallel_map(_maybe_fail, items, use_processes=False)
    assert out[3] is None
    assert out[0] == 0 and out[9] == 18

    out_par = parallel.parallel_map(_maybe_fail, items, use_processes=True,
                                    max_workers=4)
    assert out_par[3] is None
    assert out_par == out


def test_on_result_fires_once_per_item():
    items = list(range(12))
    seen = []
    parallel.parallel_map(_square, items, use_processes=True, max_workers=3,
                          on_result=lambda i, r: seen.append((i, r)))
    assert len(seen) == len(items)
    assert dict(seen) == {i: i * i for i in items}


def test_should_stop_cancels_at_least_some_not_yet_started_items():
    """Cancellation of queued-but-undispatched work is genuinely best-effort
    (ProcessPoolExecutor prefetches into its internal call queue ahead of a
    worker actually picking an item up, and that hand-off race isn't
    something a caller can fully control) -- so this only asserts stopping
    early actually prevents SOME work, and that whatever DID run returned a
    correct (not corrupted) value, rather than pinning an exact cutoff index
    that would just be re-testing ProcessPoolExecutor's internal queue depth."""
    items = list(range(12))
    stop = {"n": 0}

    def should_stop():
        stop["n"] += 1
        return stop["n"] > 1                # stop right after the 1st result

    # one worker + a real per-item delay gives should_stop an actual window
    # to land before every item has already been dispatched/finished
    out = parallel.parallel_map(_slow_square, items, use_processes=True,
                                max_workers=1, should_stop=should_stop)
    ran = [(i, r) for i, r in enumerate(out) if r is not None]
    assert ran                               # at least the first item ran
    assert all(r == i * i for i, r in ran)   # nothing corrupted
    assert any(r is None for r in out)       # stopping early actually cut work


def test_shared_executor_is_reused_across_calls_and_not_closed():
    items = list(range(10))
    pool = ProcessPoolExecutor(max_workers=2)
    try:
        a = parallel.parallel_map(_square, items, executor=pool)
        b = parallel.parallel_map(_square, [1, 2, 3, 4, 5, 6, 7, 8],
                                  executor=pool)
        assert a == [i * i for i in items]
        assert b == [i * i for i in [1, 2, 3, 4, 5, 6, 7, 8]]
        # a THIRD call on the same pool still works -- proof parallel_map
        # never shut it down itself (that's the caller's job)
        c = parallel.parallel_map(_square, [10, 11], executor=pool)
        assert c == [100, 121]
    finally:
        pool.shutdown(wait=True)


def test_default_worker_count_leaves_one_core_free(monkeypatch):
    monkeypatch.setattr(parallel.os, "cpu_count", lambda: 8)
    assert parallel.default_worker_count() == 7
    monkeypatch.setattr(parallel.os, "cpu_count", lambda: 1)
    assert parallel.default_worker_count() == 1     # never below 1
    monkeypatch.setattr(parallel.os, "cpu_count", lambda: None)
    assert parallel.default_worker_count() == 1
