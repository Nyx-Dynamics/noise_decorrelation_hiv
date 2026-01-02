"""
Final calibrated model (root wrapper)
=====================================

This module re-exports the canonical implementation from
quantum.final_calibrated_model so it can be imported or executed from
the repository root as:

  python -m final_calibrated_model

It provides the same public API:
- coherence_to_NAA_optimized(coherence, xi, sigma_r)
- run_full_model(condition)
- validate_model()
- plot_xi_dependence()
"""
from __future__ import annotations

from quantum.final_calibrated_model import (
    coherence_to_NAA_optimized,
    run_full_model,
    validate_model,
    plot_xi_dependence,
)

__all__ = [
    "coherence_to_NAA_optimized",
    "run_full_model",
    "validate_model",
    "plot_xi_dependence",
]


if __name__ == "__main__":
    # Run validation and a quick diagnostic plot to match the behavior of the
    # canonical module when executed directly.
    validate_model()
    try:
        plot_xi_dependence()
    except Exception:
        # Plotting is optional (e.g., headless environment)
        pass
