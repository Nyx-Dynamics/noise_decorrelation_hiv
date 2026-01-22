#!/usr/bin/env python3
"""
Check Bayesian stack installation: imports and versions.

Usage:
  PYTHONPATH=. python Extra/check_bayes_env.py

Exit codes:
  0 = all good, imports succeeded
  1 = at least one import failed
"""
import sys

pkgs = [
    ("numpy", "numpy"),
    ("scipy", "scipy"),
    ("matplotlib", "matplotlib"),
    ("pymc", "pymc"),
    ("pytensor", "pytensor"),
    ("arviz", "arviz"),
    ("pandas", "pandas"),
    ("xarray", "xarray"),
    ("netcdf4", "netCDF4"),
]

print("Bayesian stack environment check")
print("-" * 40)

ok = True
versions = {}
for name, mod in pkgs:
    try:
        m = __import__(mod)
        ver = getattr(m, "__version__", getattr(m, "__VERSION__", "(unknown)"))
        versions[name] = ver
        print(f"{name:>10}: {ver}")
    except Exception as e:
        ok = False
        print(f"{name:>10}: IMPORT FAILED — {e}")

# Quick compatibility hint for SciPy/ArviZ
sci = versions.get("scipy")
avz = versions.get("arviz")
if sci and avz:
    try:
        from packaging.version import Version
        sci_v = Version(sci)
        arv_v = Version(avz)
        # Known issue: older ArviZ releases used scipy.signal.gaussian removed in SciPy>=1.11
        if sci_v >= Version("1.11.0"):
            # Recommend ArviZ >= 0.16 as a safe floor; 0.17.1 current works well
            if arv_v < Version("0.16.0"):
                print("\nWarning: SciPy >= 1.11 detected with ArviZ < 0.16.\n         Consider upgrading ArviZ: pip install --upgrade arviz")
    except Exception:
        pass

print("-" * 40)

# Explicit NumPy 2.x check (PyTensor expects NumPy 1.x config API)
try:
    import numpy as _np
    _major = int(_np.__version__.split(".")[0])
    if _major >= 2:
        ok = False
        print(
            "Detected NumPy %s (>=2.0). PyTensor (PyMC backend) may fail with AttributeError on numpy.__config__.get_info.\n"
            "Fix: pip install \"numpy<2.0\" in your venv or run: make install (uses requirements.txt pins)."
            % _np.__version__
        )
except Exception:
    pass

print("Result:", "OK" if ok else "FAILED")

sys.exit(0 if ok else 1)
