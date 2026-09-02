"""Lineshape model registry. Importing this package registers all models."""
from larmor.models.base import (
    Model, ParamDef, SimContext, REGISTRY, describe_all, get, register,
)
# importing a model module registers its models; a new module must be added
# here or it never registers (docs/development-notes.md section 6)
from larmor.models import analytic, quadrupolar, csa, external

__all__ = [
    "Model", "ParamDef", "SimContext", "REGISTRY", "describe_all", "get",
    "register", "analytic", "quadrupolar", "csa", "external",
]
