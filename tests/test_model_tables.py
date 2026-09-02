"""The hand-maintained model tables must be complete, explicit partitions.

Adding a lineshape means touching six tables outside the registry, and
omission from each degrades behaviour SILENTLY (wrong Jacobian step, slow
full-grid probing, no data-driven seeds, a misapplied width floor). These
tests turn every silent omission into a red test naming the decision the
model's author still has to make.
"""
from larmor import models


def _registry() -> set:
    return set(models.REGISTRY)


def test_fit_jacobian_step_tables_partition_the_registry():
    """Every model is either analytic (scipy's fine diff step) or simulated
    (the coarse SIMULATED_DIFF_STEP) -- explicitly, never by omission."""
    from larmor.fit import _ANALYTIC_MODELS, _SIMULATED_MODELS

    reg = _registry()
    both = _ANALYTIC_MODELS & _SIMULATED_MODELS
    assert not both, f"models claimed as both analytic and simulated: {both}"
    missing = reg - _ANALYTIC_MODELS - _SIMULATED_MODELS
    assert not missing, (
        f"models with no declared Jacobian step: {missing} -- add each to "
        "fit._ANALYTIC_MODELS (closed-form, pointwise) or "
        "fit._SIMULATED_MODELS (grid-based, needs the coarse step)")
    stale = (_ANALYTIC_MODELS | _SIMULATED_MODELS) - reg
    assert not stale, f"tables list unregistered models: {stale}"


def test_grid_restriction_tables_partition_the_registry():
    """Every model either tolerates a window-restricted grid or is audited
    as needing the full experimental axis."""
    from larmor.engine import _GRID_FULL_REQUIRED, _GRID_RESTRICTABLE

    reg = _registry()
    both = _GRID_RESTRICTABLE & _GRID_FULL_REQUIRED
    assert not both, f"models in both grid tables: {both}"
    missing = reg - _GRID_RESTRICTABLE - _GRID_FULL_REQUIRED
    assert not missing, (
        f"models with no grid-restriction audit: {missing} -- add each to "
        "engine._GRID_RESTRICTABLE (pointwise or own-grid) or "
        "engine._GRID_FULL_REQUIRED (derives its simulation from the axis)")
    stale = (_GRID_RESTRICTABLE | _GRID_FULL_REQUIRED) - reg
    assert not stale, f"tables list unregistered models: {stale}"


def test_estimate_seed_tables_partition_the_registry():
    """Every model either declares its breadth parameter for data-driven
    starting values or is explicitly excused."""
    from larmor.estimate import _NO_WIDTH_SEED, _WIDTH_KEY

    reg = _registry()
    keyed = set(_WIDTH_KEY)
    both = keyed & _NO_WIDTH_SEED
    assert not both, f"models both seeded and excused: {both}"
    missing = reg - keyed - _NO_WIDTH_SEED
    assert not missing, (
        f"models with no seeding decision: {missing} -- add a breadth "
        "parameter to estimate._WIDTH_KEY or excuse it in _NO_WIDTH_SEED")
    stale = (keyed | _NO_WIDTH_SEED) - reg
    assert not stale, f"tables list unregistered models: {stale}"
    # and each declared width key must be a real parameter of that model
    for name, (key, _is_cq) in _WIDTH_KEY.items():
        assert key in models.get(name).param_names, (name, key)


def test_peak_fwhm_tables_partition_the_registry():
    """The Eden >=4 ppm amorphous-peak floor applies only where
    shift_fwhm_ppm IS the peak width; every model states its side."""
    from larmor.constraints_util import (_NOT_PEAK_FWHM_MODELS,
                                         _PEAK_FWHM_MODELS)

    reg = _registry()
    peak, non = set(_PEAK_FWHM_MODELS), set(_NOT_PEAK_FWHM_MODELS)
    both = peak & non
    assert not both, f"models on both sides of the width floor: {both}"
    missing = reg - peak - non
    assert not missing, (
        f"models with no width-floor decision: {missing} -- add each to "
        "constraints_util._PEAK_FWHM_MODELS (shift_fwhm_ppm is the peak "
        "FWHM) or _NOT_PEAK_FWHM_MODELS (it is dCS / not a peak width)")
    stale = (peak | non) - reg
    assert not stale, f"tables list unregistered models: {stale}"


def test_param_columns_have_no_dead_entries():
    """PARAM_COLUMNS drives column order (its fallback covers omissions, so
    ordering is the only stake) -- but an entry no model uses is a trap."""
    from larmor.desktop import table

    all_params = {p for name in _registry()
                  for p in models.get(name).param_names}
    dead = [key for key, _label in table.PARAM_COLUMNS
            if key not in all_params]
    assert not dead, f"PARAM_COLUMNS lists parameters no model has: {dead}"
