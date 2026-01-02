"""
Lightweight enzyme kinetics utilities for bayesian_enzyme_v4.py

This module provides a minimal, self‑contained implementation that matches the
API expected by bayesian_enzyme_v4.py:

- EnzymeKinetics class with an .integrate(duration_days, membrane_turnover)
  method returning (NAA_molar, Cho_molar)
- compute_protection_factor(xi, beta_xi)
- coherence_modulation(coh, gamma)
- ENZYME placeholder (for potential future use)

Implementation notes
- Designed to work both with plain Python/NumPy floats and Aesara tensors used by
  PyMC. Operations switch to aesara.tensor (at) when any input is a tensor.
- The model here is a simple, dimensionally consistent surrogate capturing the
  intended qualitative dependencies:
  • Higher protection and coherence increase NAA and mildly decrease Cho
  • Higher membrane turnover increases Cho (membrane synthesis marker)
  • Viral damage reduces effective capacity
- Baseline concentrations are set so that, after division by creatine = 8 mM in
  the caller, the healthy condition roughly matches clinical data (NAA/Cr ≈ 1.105,
  Cho/Cr ≈ 0.225).

This is NOT a biochemical simulator; it is an algebraic surrogate intended to be
compatible with PyMC’s computation graph.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np

try:
    # Prefer PyTensor (PyMC v5+); fall back to Aesara if unavailable
    import pytensor.tensor as at
    from pytensor.tensor.var import TensorVariable
except Exception:
    try:
        import aesara.tensor as at
        from aesara.tensor.var import TensorVariable
    except Exception:  # pragma: no cover — allows use without PyMC/pytensor/aesara installed
        at = None
        TensorVariable = tuple()  # dummy


# Baseline molar concentrations chosen to match observed metabolite ratios after
# dividing by creatine (8 mM) in the caller.
NAA_BASE = 8.84e-3  # ~8.84 mM → NAA/Cr ≈ 1.105 in healthy
CHO_BASE = 1.80e-3  # ~1.80 mM → Cho/Cr ≈ 0.225 in healthy

XI_REF = 0.80e-9  # 0.8 nm reference correlation length


def _is_tensor(x) -> bool:
    """Return True if x is an Aesara tensor-like object."""
    if at is None:
        return False
    return isinstance(x, TensorVariable) or getattr(x, "owner", None) is not None


def _pow(x, y):
    if _is_tensor(x) or _is_tensor(y):
        return at.power(x, y)
    return np.power(x, y)


def _mul(a, b):
    if _is_tensor(a) or _is_tensor(b):
        return a * b  # at handles operator overload
    return a * b


def _div(a, b):
    if _is_tensor(a) or _is_tensor(b):
        return a / b
    return a / b


def _add(a, b):
    if _is_tensor(a) or _is_tensor(b):
        return a + b
    return a + b


def _as_tensor(x):
    if _is_tensor(x):
        return x
    if at is not None:
        try:
            return at.as_tensor_variable(x)
        except Exception:
            return x
    return x


def compute_protection_factor(xi, beta_xi: float, xi_ref: float = XI_REF):
    """
    Protection factor from correlation length xi.
    Pi_xi = (xi_ref / xi) ** beta_xi
    Works with floats or Aesara tensors.
    """
    xi = _as_tensor(xi)
    beta_xi = _as_tensor(beta_xi)
    return _pow(_div(xi_ref, xi), beta_xi)


def coherence_modulation(coh, gamma: float):
    """
    Coherence modulation factor.
    eta_coh = coh ** gamma
    """
    coh = _as_tensor(coh)
    gamma = _as_tensor(gamma)
    return _pow(coh, gamma)


@dataclass
class EnzymeKinetics:
    Pi_xi: float
    eta_coh: float
    viral_damage_factor: float = 1.0  # 1.0 = no damage

    def integrate(self, duration_days: float, membrane_turnover: float) -> Tuple[object, object]:
        """
        Algebraic surrogate for metabolite steady-state over the specified duration.
        Returns (NAA_molar, Cho_molar). Compatible with Aesara if inputs are tensors.
        """
        Pi = _as_tensor(self.Pi_xi)
        eta = _as_tensor(self.eta_coh)
        vd = _as_tensor(self.viral_damage_factor)
        mt = _as_tensor(membrane_turnover)

        # Effective capacity — geometric blend; viral damage scales linearly
        eff = _mul(_mul(_pow(Pi, 0.6), _pow(eta, 0.4)), vd)

        # NAA increases with effective capacity
        NAA = _mul(NAA_BASE, eff)

        # Cho increases with membrane turnover and decreases weakly with eff
        Cho = _mul(CHO_BASE, _mul(mt, _pow(_div(1.0, eff), 0.2)))

        return NAA, Cho


# Placeholder — exposed for API completeness; not used in v4 code
class ENZYME:
    NAA = "NAA"
    CHO = "CHO"
