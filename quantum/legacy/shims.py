"""
Stable import shims for cylindrical-aware metrics used by SSE Option B.

Preference order:
1) script/Tegmarks_Cat_source_code_3/metrics.py
2) script/tegmarks_cat_source_code/code/metrics.py

Exposes:
- calculate_probability_density
- calculate_coherence
- calculate_dispersion_metrics

Notes:
- These functions incorporate cylindrical volume elements (R*dr*dz*2π) as used in legacy code.
"""

# Try preferred source first
_METRICS_SOURCE = None
try:  # Preferred cohesive implementation
    from Legacy.Tegmarks_Cat_source_code_3.metrics import (
        calculate_probability_density,
        calculate_coherence,
        calculate_dispersion_metrics,
    )  # type: ignore
    _METRICS_SOURCE = "script/Tegmarks_Cat_source_code_3/metrics.py"
except Exception:  # Fallback to alternate legacy path
    try:
        from Legacy.tegmarks_cat_source_code.code.metrics import (
            calculate_probability_density,
            calculate_coherence,
            calculate_dispersion_metrics,
        )  # type: ignore
        _METRICS_SOURCE = "script/tegmarks_cat_source_code/code/metrics.py"
    except Exception as e:
        raise ImportError(
            "Unable to import cylindrical metrics from preferred or fallback locations. "
            "Ensure legacy metrics modules are present."
        ) from e


def get_metrics_source() -> str:
    """Return the actual metrics module path that was imported via this shim."""
    return _METRICS_SOURCE or "unknown"
