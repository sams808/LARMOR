"""Lineshape model registry. Importing this package registers all models."""
from larmor.models.base import (  # noqa: F401
    Model, ParamDef, SimContext, REGISTRY, describe_all, get, register,
)
from larmor.models import analytic, quadrupolar, csa  # noqa: F401  (register)
