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


def test_dft_seed_tables_partition_the_registry():
    """The DFT import routes tensor components by ROLE (dft._SEED_KEYS): a
    model states which parameter takes the shielding zeta / eta and C_Q /
    eta_Q, or is declared not seedable. Matching on the shared name 'eta'
    once wrote a quadrupolar eta into csa_mas."""
    from larmor.dft import _NOT_SEEDABLE, _SEED_KEYS

    reg = _registry()
    keyed = set(_SEED_KEYS)
    both = keyed & _NOT_SEEDABLE
    assert not both, f"models both seedable and not: {both}"
    missing = reg - keyed - _NOT_SEEDABLE
    assert not missing, (
        f"models with no DFT seeding decision: {missing} -- add a "
        "(zeta_key, eta_cs_key, cq_key, eta_q_key) entry to dft._SEED_KEYS "
        "or excuse the model in dft._NOT_SEEDABLE")
    stale = (keyed | _NOT_SEEDABLE) - reg
    assert not stale, f"tables list unregistered models: {stale}"
    for name, keys in _SEED_KEYS.items():
        assert len(keys) == 4, name
        pnames = models.get(name).param_names
        for k in keys:
            assert k is None or k in pnames, (name, k)
        zeta, eta_cs, cq, eta_q = keys
        # an eta is only meaningful next to its anisotropy / coupling
        assert (eta_cs is None) or (zeta is not None), name
        assert (eta_q is None) or (cq is not None), name
        # the two etas of one model must be two different parameters
        assert eta_cs is None or eta_q is None or eta_cs != eta_q, name
        if models.get(name).needs_quadrupolar:
            assert cq is not None or zeta is None, (
                f"{name} is quadrupolar but seeds a CSA without a C_Q")


def test_param_columns_have_no_dead_entries():
    """PARAM_COLUMNS drives column order (its fallback covers omissions, so
    ordering is the only stake) -- but an entry no model uses is a trap."""
    from larmor.desktop import table

    all_params = {p for name in _registry()
                  for p in models.get(name).param_names}
    dead = [key for key, _label in table.PARAM_COLUMNS
            if key not in all_params]
    assert not dead, f"PARAM_COLUMNS lists parameters no model has: {dead}"


#: the models whose lb (line_fwhm_ppm) the registry holds at its default
#: unless the user frees it -- dmfit's greyed CzSimple Lb: every distribution
#: model (a Czjzek / Gaussian spread already carries the breadth, so a free
#: lb only ends 'at bounds' in most fits). Every other model that has an lb
#: fits it (a crystalline quad_ct has nothing else to set its width).
DEFAULT_FIXED_LB = frozenset({"czjzek", "czjzek_d", "czjzek_corr",
                              "ext_czjzek", "amorphous"})


def test_default_fixed_partitions_the_models_with_a_line_width():
    """Every model with ``line_fwhm_ppm`` states whether it is fitted by
    default or held at its default (``ParamDef.default_fixed``), and the
    flag never exists without ``vary=False`` behind it."""
    with_lb = {name for name in _registry()
               if "line_fwhm_ppm" in models.get(name).param_names}
    assert DEFAULT_FIXED_LB <= with_lb, DEFAULT_FIXED_LB - with_lb
    held = set()
    for name in _registry():
        for pd in models.get(name).params:
            if pd.default_fixed:
                assert not pd.vary, (name, pd.name)
                assert pd.name == "line_fwhm_ppm", (
                    f"{name}.{pd.name} is default_fixed: add it to this "
                    "partition (only lb is held by default today)")
                held.add(name)
    assert held == DEFAULT_FIXED_LB, (
        f"default_fixed lb differs from DEFAULT_FIXED_LB: "
        f"+{held - DEFAULT_FIXED_LB} -{DEFAULT_FIXED_LB - held}")
    for name in with_lb - DEFAULT_FIXED_LB:
        pd = next(p for p in models.get(name).params if p.name == "line_fwhm_ppm")
        assert pd.vary, f"{name}.line_fwhm_ppm is fixed without default_fixed"
    # a fresh site of a held-lb model comes out pinned, at the default
    for name in DEFAULT_FIXED_LB:
        p = models.get(name).defaults()["line_fwhm_ppm"]
        assert p.vary is False
        assert p.value == next(pd.default for pd in models.get(name).params
                               if pd.name == "line_fwhm_ppm")
    # the JSON dump carries the flag for the web UI
    dumped = {m["name"]: m for m in models.describe_all()}
    assert {n for n, m in dumped.items()
            if any(p["default_fixed"] for p in m["params"])} == DEFAULT_FIXED_LB
