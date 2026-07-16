"""Model registry: every lineshape model LARMOR can fit, self-describing.

A model declares its parameters (name, short key for lmfit, default, unit,
human description) and how to render itself on a ppm axis. The app UI, the
fit engine, and the constraint translator all read this registry -- adding a
new model in Python makes it appear everywhere, including the "add site"
buttons in the app.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np


@dataclass(frozen=True)
class ParamDef:
    name: str          # recipe/JSON name, e.g. "isotropic_chemical_shift_ppm"
    key: str           # short lmfit-safe key, e.g. "pos"
    default: float
    unit: str = ""
    description: str = ""
    min: float | None = None
    max: float | None = None
    vary: bool = True


@dataclass(frozen=True)
class SimContext:
    """Everything a model needs to know about the experiment."""

    nucleus: str
    larmor_MHz: float
    spin_rate_Hz: float
    x_ppm: np.ndarray           # ascending axis the site must be rendered on


@dataclass(frozen=True)
class Model:
    name: str                       # registry key, e.g. "czjzek"
    label: str                      # human name for the UI
    description: str
    params: tuple[ParamDef, ...]
    render: Callable[[dict, SimContext], np.ndarray]  # values-by-name -> y
    needs_quadrupolar: bool = False  # UI hint: only offer for I > 1/2

    @property
    def param_names(self) -> tuple[str, ...]:
        return tuple(p.name for p in self.params)

    def key_of(self, param_name: str) -> str:
        for p in self.params:
            if p.name == param_name:
                return p.key
        raise KeyError(f"{self.name} has no parameter {param_name!r}")

    def defaults(self) -> dict:
        """Fresh recipe-style params dict for a new site of this model."""
        from larmor.recipe import Param

        return {p.name: Param(p.default, vary=p.vary, min=p.min, max=p.max)
                for p in self.params}


REGISTRY: dict[str, Model] = {}


def register(model: Model) -> Model:
    REGISTRY[model.name] = model
    return model


def get(name: str) -> Model:
    try:
        return REGISTRY[name]
    except KeyError:
        raise ValueError(
            f"unknown site model {name!r} (available: {sorted(REGISTRY)})") from None


def describe_all() -> list[dict]:
    """JSON-friendly registry dump for the app UI."""
    return [
        {
            "name": m.name,
            "label": m.label,
            "description": m.description,
            "needs_quadrupolar": m.needs_quadrupolar,
            "params": [
                {"name": p.name, "key": p.key, "default": p.default,
                 "unit": p.unit, "description": p.description,
                 "min": p.min, "max": p.max, "vary": p.vary}
                for p in m.params
            ],
        }
        for m in REGISTRY.values()
    ]
