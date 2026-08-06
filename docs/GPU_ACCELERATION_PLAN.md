# GPU acceleration — an opt-in "advanced calculation mode"

Written 2026-08-06 after profiling a real batch fit (7-site 31P, 16 spectra) that
took 15 minutes / ~8000 lmfit evaluations for amplitude+δiso, and got slower again
releasing width. LARMOR must stay installable and runnable on any colleague's
laptop with nothing beyond the standard install — this plan is explicitly for an
**optional, auto-detected, silently-absent-if-unavailable** mode for users who
have a CUDA GPU (the author's Acer Nitro / RTX is the reference machine), not a
new baseline requirement.

## 0 · Relationship to the CPU optimization plan

This is **not** a substitute for the algorithmic fixes in `docs/AUDIT.md` §CPU
optimization (skip the errorbar-retry doubling; restrict the simulation grid to
the fit window; fit each batch spectrum independently instead of one 224+
parameter joint problem). Those land first, are free on every machine, and are a
prerequisite here: GPU acceleration only pays for itself once evaluations are
*batched* (many spectra, or many finite-difference probes, computed together) —
which is exactly the shape the independent-per-spectrum refactor creates. Doing
GPU work before that refactor would mean accelerating a computation shape (one
huge sequential joint problem) that shouldn't exist in the first place.

**Sequencing:** CPU Points 1–3 → profile the actual remaining bottleneck on a
large Czjzek batch + Monte-Carlo run → *then* decide whether this plan is worth
building. Section 6 gives the go/no-go criteria.

## 1 · Where LARMOR's computation is actually GPU-shaped

Profiled by reading (not guessing) `larmor/engine.py`, `larmor/models/*.py`,
`larmor/fit.py`, `larmor/multifit.py`:

| Computation | Shape | GPU-friendly? |
|---|---|---|
| `gauss_lor`, `voigt`, `gl_norm`, `jmultiplet`, `sidebands` (`larmor/models/analytic.py`) | closed-form pointwise formula over ~1000–4000 points | **No** — already microseconds on CPU; a GPU kernel launch costs more than the compute it replaces at this size |
| `quad_ct`/`quad_first`/`quad_csa`/`csa_mas`/`csa_czjzek` (`larmor/models/_singlesite.py`) | one mrsimulator `Simulator.run()` per (rounded) parameter set, LRU-cached | **No** — mrsimulator is a C-extension LARMOR doesn't control; the cache already avoids re-running identical parameter sets |
| **Czjzek kernel reweighting** — `kernel.weights(sigma) @ kernel.K` (`larmor/models/quadrupolar.py::_render_czjzek`, `_render_ext_czjzek`, `_render_amorphous`) | one (grid-points,) row-vector times a fixed (grid-points × spectral-points) matrix, called on **every** evaluation | **Yes, once batched** — this is the one piece of LARMOR's own math that's a genuine linear-algebra op, not a pointwise formula |
| **Monte-Carlo error analysis** (`larmor/autofit.py::monte_carlo_errors`) — `n_trials` (default 200) independent synthetic-noise refits | naturally embarrassingly parallel, identical operation repeated N times | **Yes** — the single best-shaped GPU (or even just multi-core) target in the app |
| `apply_params`/`_pname`/`_key` bookkeeping (`larmor/multifit.py`) | Python dict lookups + string formatting, called on every evaluation | **No** — pure interpreter overhead; a GPU can't help, this is what CPU Point 2 fixes |

**Conclusion:** the target is narrow and specific — batch the Czjzek/ext_czjzek/
Amorphous kernel reweighting step across many simultaneous evaluations (many MC
trials, or many finite-difference Jacobian probes, or many spectra in a batch),
and hand *that* one matmul to the GPU. Nothing else in the current codebase
benefits from a GPU as-is.

## 2 · Library choice

| Option | Pros | Cons | Verdict |
|---|---|---|---|
| **CuPy** | Near drop-in numpy API — `kernel.weights(...) @ kernel.K` becomes the same line on `cupy.ndarray`s; smallest diff against the current code | Needs a CUDA-version-matched wheel (`cupy-cuda11x`/`cupy-cuda12x`); no autodiff | **Primary choice for v1** |
| **PyTorch** | Automatic differentiation — could supply an *exact* Jacobian for the reweighting step instead of finite-differencing it, removing the O(n_params) evaluation cost algorithmically, with GPU as a side effect | Heavier dependency; lmfit expects a plain numpy residual, so a boundary layer is needed either way; bigger paradigm shift | Worth prototyping later for the analytic-Jacobian angle (see §7), not v1 |
| **JAX** | Best-in-class autodiff + JIT + GPU/TPU; mathematically the "right" tool for "many small parametrized simulations, need gradients, want GPU" | Least common in scientific NMR software; steepest adoption/debugging curve for future contributors | Not recommended — optimizes for a problem LARMOR doesn't have (TPUs, huge models); PyTorch gets the autodiff benefit with a gentler curve |
| Numba CUDA (hand-written kernels) | Full control | All of the engineering cost, none of numpy-API convenience | Not recommended — CuPy already gets the matmul onto the GPU with a 1-line change |

**Decision: CuPy for the reweighting matmul. Revisit PyTorch only if/when an
analytic-Jacobian rewrite is separately justified** (§7 — a bigger, CPU-side win
that happens to pair naturally with GPU).

## 3 · Architecture

### 3.1 One isolated backend module, not a forked engine

New module `larmor/gpu.py` (Qt-free, like every other core module):

```python
"""Optional CuPy-backed array backend for the Czjzek-family kernel reweighting
step. Absent by default; every call site falls back to plain numpy identically.
Nothing else in LARMOR imports cupy directly."""

_XP = None          # resolved lazily, never at import time
_DEVICE_NAME = None

def available() -> bool:
    """True once xp() has successfully resolved a working CUDA device."""
    ...

def xp():
    """The array module to use: cupy if a working CUDA device was found and the
    user has GPU mode on, else numpy. Resolved once, cached, never raises —
    any failure (no cupy installed, no device, driver mismatch) falls back."""
    ...

def to_device(a): ...     # numpy -> xp array (no-op if xp() is numpy)
def to_host(a): ...       # xp array -> numpy (no-op if xp() is numpy)
```

`larmor/models/quadrupolar.py` changes **only** at the reweighting line, gated
by a module-level flag the desktop Settings toggle controls — not a parallel
code path per model. Kernel *building* (`engine.build_kernel`, which calls into
mrsimulator) stays on CPU/numpy always — mrsimulator itself is CPU-only, and the
build is a one-time, already-cached cost.

### 3.2 The batching layer is the real work, not the CuPy call

CuPy on a single `(grid,) @ (grid, points)` vector-matrix product buys almost
nothing — the launch overhead argument in §1 applies here too. The actual change
needed is upstream: whichever caller currently does

```python
for each evaluation:
    y = kernel.weights(sigma_i) @ kernel.K       # one at a time
```

needs to instead accumulate a batch of `sigma` values and do

```python
W = xp.stack([kernel_weights(s) for s in sigma_batch])   # (n_batch, grid)
Y = W @ kernel.K                                          # ONE matmul, all at once
```

Two natural batching points, in order of how cleanly they fit today's code:

1. **Monte-Carlo trials** (`autofit.monte_carlo_errors`): today, `n_trials`
   independent `fitmod.fit(...)` calls run sequentially. Restructuring this to
   batch the *reweighting step* across trials means the trials' own optimizers
   still run independently (each trial still needs its own lmfit.minimize with
   its own converging trajectory — trials don't finish in lockstep), so this
   isn't a trivial "stack everything" batch. The practical version: batch
   reweighting **within** the finite-difference Jacobian of a *single* trial's
   fit (same idea as #2 below, applied once per trial) rather than across
   trials. Cross-trial batching would require re-architecting MC as a custom
   batched optimizer loop — worth doing eventually (§7) but is a bigger lift
   than v1 should take on.
2. **The finite-difference Jacobian of one Czjzek fit** (whether a single
   spectrum or, after CPU Point 2, one of the independent per-spectrum batch
   fits): supply `scipy.optimize.least_squares` a custom `jac=` callable that
   perturbs every free parameter **once, all at once**, batches every resulting
   `sigma` value into one `W` matrix, does one `W @ kernel.K`, and assembles the
   Jacobian from the batched result — instead of `least_squares`'s default
   `2-point` scheme calling the residual function once per parameter,
   sequentially. **This is the actual v1 deliverable** — everything else in this
   plan exists to support it.

### 3.3 Precision

Consumer RTX cards are FP32-optimized; FP64 throughput is a fraction of FP32 on
non-datacenter silicon. Plan to run the forward simulation (`kernel.weights(...)
@ kernel.K`) in FP32 on GPU and cast back to FP64 for lmfit's own bookkeeping
(the optimizer's internal linear algebra stays CPU/FP64 either way — only the
model-evaluation step touches the GPU). **Must be validated, not assumed**: a
dedicated test compares FP32-GPU vs FP64-CPU simulated spectra for a
representative Czjzek recipe and asserts the difference is far below the noise
floor before this ships.

## 4 · Staged implementation plan

**Phase 0 — Backend scaffold (no behavior change, safe to ship standalone).**
- `larmor/gpu.py`: `available()`/`xp()`/`to_device()`/`to_host()` as above.
  Lazy import (`import cupy` only inside `xp()`, wrapped in `try/except`, same
  convention already used everywhere in this codebase for optional deps).
- A Settings/Preferences toggle in the desktop app: "Use GPU acceleration
  (experimental)" — off by default, disabled/greyed with an explanatory tooltip
  when `gpu.available()` is False, so the control is honest about whether it can
  do anything on this machine.
- An "About" or Help-menu line reporting detected device name/CUDA version when
  available, for support/debugging.
- Tests: `gpu.available()` returns `False` cleanly on a machine without CUDA
  (this dev's CI/most machines); a GPU-gated test (skipped without a device)
  exercises `xp()`/`to_device()`/`to_host()` round-trips on the author's laptop.

**Phase 1 — Batched Jacobian for Czjzek-family fits.**
- A custom `jac=` callable for `fit.py::fit()`'s `lmfit.minimize(..., method="least_squares")`
  call, used only when (a) GPU mode is on, (b) `gpu.available()`, and (c) the
  recipe's sites are all Czjzek-family (reuses the same allowlist concept as CPU
  Point 3 — deliberately conservative, opt-in per model, not a denylist).
- Batches every free parameter's forward-difference perturbation into one `W`
  matrix, one `W @ kernel.K` on GPU, assembles the Jacobian, returns it to
  `least_squares`.
- Falls back to today's default finite-difference behavior transparently for
  every other case (non-Czjzek models, GPU unavailable, GPU mode off).
- Tests: fitted parameter values and their uncertainties match the
  finite-difference path within a tight numerical tolerance on a synthetic
  Czjzek recipe (this is the correctness gate — a faster wrong answer is worse
  than a slow right one).

**Phase 2 — Wire into Monte-Carlo and batch fitting.**
- `autofit.monte_carlo_errors` and the (by-then-independent, per CPU Point 2)
  per-spectrum batch fit both call `fit.py::fit()` — Phase 1's Jacobian batching
  applies to them automatically, no separate integration work, *provided* CPU
  Point 2 has already landed (independent per-spectrum fits are what makes each
  individual fit's own Jacobian batching meaningful at batch scale).
- Add a `k`/`n`-style progress note ("GPU: batching N-parameter Jacobian") to
  the existing progress-callback text so it's visible when GPU mode is actually
  doing something, for user confidence/debugging.

**Phase 3 — Packaging.**
- `pyproject.toml` optional extra: `pip install larmor[gpu]` pulling in a CuPy
  version range (not pinned to one CUDA version — document how to pick the
  right `cupy-cudaXXx` wheel for a given driver in the README, same as CuPy's
  own install docs).
- The standard PyInstaller-frozen `.exe` distributed to colleagues does **not**
  bundle CuPy — GPU mode starts life as a source/dev-environment feature only.
  Revisit bundling only if there's real demand from someone running the frozen
  build who also has a suitable GPU.

**Phase 4 — Validation & rollout.**
- FP32/FP64 comparison test (§3.3).
- Before/after wall-clock benchmark on the author's actual 27Al Czjzek batch +
  Monte-Carlo workload (the real target workload, not a synthetic microbenchmark)
  documented in this file or a follow-up `docs/GPU_BENCHMARK.md`.
- Ship behind the off-by-default toggle for at least one full batch session
  before considering any change to defaults.

## 5 · Testing strategy

- Every GPU-touching test uses a `pytest.mark.skipif(not gpu.available(), ...)`
  guard — mirrors the existing pattern for other optional/hardware-dependent
  tests in this suite. CI and every colleague's machine skip them silently;
  they only run where they can mean something.
- A **non-GPU-gated** test suite covers the batching/Jacobian-assembly *logic*
  against plain numpy (no CuPy import needed) — the batching math itself (which
  parameters get perturbed, how the Jacobian is assembled from the batched
  result) is backend-agnostic and should be fully covered without hardware.
- Equivalence tests (Phase 1) are the load-bearing ones: GPU-batched Jacobian
  fit results must match today's finite-difference fit results within
  tolerance, always, not just "look reasonable."

## 6 · Go / no-go — when this is actually worth building

Do **not** start Phase 0 until, after CPU Points 1–3 land:
1. A real 27Al Czjzek batch + Monte-Carlo run has been profiled (not guessed),
   and
2. The Czjzek kernel-reweighting step is confirmed to still be the dominant
   remaining cost (as opposed to, say, mrsimulator's own kernel-build time, or
   something CPU Points 1–3 didn't fully address).

If the CPU-side fixes already bring a large Czjzek batch + 200-trial
Monte-Carlo run down to a few minutes, the engineering cost of Phases 0–4 (a
new module, a Settings toggle, a custom Jacobian, packaging, and a whole
GPU-gated test tier to maintain going forward) may not be worth it for a
single-user "nice to have." This plan exists so the decision is well-informed
when that moment comes, not so it's built on spec.

## 7 · Explicitly out of scope for this plan (future ideas, not commitments)

- **Analytic Jacobians for the closed-form models** (`gauss_lor`, `voigt`, …):
  removes the O(n_params) finite-difference cost *algorithmically*, no GPU
  needed, and is arguably higher-value than anything in this document for the
  model family the user's current 31P batch actually uses. A good candidate for
  its own, GPU-independent workplan.
- **PyTorch/autodiff rewrite** of the Czjzek reweighting step, replacing
  hand-batched finite differences with an exact gradient. Natural follow-on to
  Phase 1 if the finite-difference-but-batched approach proves the concept.
- **Cross-trial Monte-Carlo batching** (running many trials' optimizers in
  lockstep so their reweighting steps batch together, not just each trial's own
  internal Jacobian) — bigger optimizer-architecture change, deferred.
- **GPU-accelerated kernel *building*** (the mrsimulator call inside
  `engine.build_kernel`): out of LARMOR's control (mrsimulator is CPU-only
  today) and already a one-time cached cost — not worth chasing unless
  mrsimulator itself grows GPU support upstream.
- **Bundling CuPy into the distributed PyInstaller build** for colleagues —
  deferred until there's demonstrated demand from someone who both runs the
  frozen build and has a suitable GPU.
