#!/usr/bin/env python3
"""
Shared utilities for the absolute and ratio pipelines.

Functions provided:
- parse_time_slice: interpret user-provided week filters for Valcour.
- parse_regions: normalize regions argument to a normalized set.
- list_curated_files: discover input files from a curated directory using
  optional manifest(s); falls back to pattern-based discovery if manifest
  is not a simple text/csv list.
- write_run_manifest: write a manifest.json alongside outputs capturing
  CLI args, resolved inputs, filters, and (optionally) git hash.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Dict, Any, Optional, Set
import json
import re


def parse_time_slice(spec: str) -> Optional[List[int]]:
    """Parse a time-slice specification for Valcour weeks.

    Supported forms:
    - "all": returns None (no filtering)
    - "baseline": [0]
    - "acute_window": range 0..12 inclusive
    - comma list: "0,4,12,24"
    - range w/ step: "0-24:4" -> [0,4,8,12,16,20,24]
    """
    s = (spec or "").strip().lower()
    if s in ("", "all"):
        return None
    if s == "baseline":
        return [0]
    if s == "acute_window":
        return list(range(0, 13))

    # range with step
    m = re.fullmatch(r"\s*(\d+)\s*-\s*(\d+)\s*:\s*(\d+)\s*", s)
    if m:
        start, end, step = map(int, m.groups())
        if step <= 0:
            step = 1
        vals = list(range(start, end + 1, step))
        return vals

    # comma separated
    parts = [p.strip() for p in s.split(",") if p.strip()]
    try:
        vals = sorted({int(p) for p in parts})
        return vals
    except ValueError:
        # Unknown spec; treat as no filter
        return None


def parse_regions(spec: str | Iterable[str] | None) -> Optional[Set[str]]:
    """Normalize region selection to an upper-case set, or None for all.
    Recognizes common Valcour region codes.
    """
    if spec is None:
        return None
    if isinstance(spec, str):
        s = spec.strip()
        if s.lower() in ("all", ""):
            return None
        parts = [p.strip().upper() for p in s.split(",") if p.strip()]
    else:
        parts = [str(p).strip().upper() for p in spec if str(p).strip()]

    # Standardize to short codes where possible
    norm = []
    for r in parts:
        x = r.upper()
        # Accept both long and short forms
        if x in {"BG", "BASAL GANGLIA", "BASAL_GANGLIA"}:
            norm.append("BG")
        elif x in {"FWM", "FRONTAL WHITE MATTER", "FRONTAL_WHITE_MATTER"}:
            norm.append("FWM")
        elif x in {"PCC", "POSTERIOR CINGULATE CORTEX", "POSTERIOR_CINGULATE_CORTEX"}:
            norm.append("PCC")
        elif x in {"PFC", "PREFRONTAL CORTEX", "PREFRONTAL_CORTEX"}:
            norm.append("PFC")
        else:
            norm.append(x)
    return set(norm)


def _read_manifest_lines(path: Path) -> Optional[List[str]]:
    """Try to read a simple text or CSV manifest listing filenames.
    Returns a list of non-empty lines if recognized; otherwise None.
    """
    try:
        if path.suffix.lower() in {".txt", ".csv", ".lst", ".list", ".md"}:
            data = path.read_text(encoding="utf-8", errors="ignore")
            lines = [ln.strip() for ln in data.splitlines() if ln.strip() and not ln.strip().startswith("#")]
            return lines or None
    except Exception:
        pass
    return None


def list_curated_files(
    directory: Path,
    patterns: Iterable[str] = ("*.csv", "*.xlsx"),
    manifest_names: Iterable[str] = ("manifest.txt", "manifest.csv", "Manifest.txt", "Manifest.csv")
) -> List[Path]:
    """Return curated files to load.

    - If a recognizable text manifest exists (txt/csv), only files listed in it are returned.
    - If only a docx manifest exists, it is recorded but ignored for discovery (cannot parse here).
    - If no manifest, return files that match provided patterns.
    """
    directory = Path(directory)
    if not directory.exists():
        return []

    # Look for simple manifests first
    for name in manifest_names:
        mpath = directory / name
        if mpath.exists():
            lines = _read_manifest_lines(mpath)
            if lines:
                resolved: List[Path] = []
                for ln in lines:
                    p = (directory / ln).resolve()
                    if p.exists() and p.is_file():
                        resolved.append(p)
                return resolved

    # Fallback: glob by patterns
    files: List[Path] = []
    for pat in patterns:
        files.extend(directory.glob(pat))
    # Filter out temporary or manifest files
    files = [f for f in files if f.is_file() and not f.name.lower().endswith("manifest.docx")]
    return sorted({f.resolve() for f in files})


def write_run_manifest(results_dir: Path, info: Dict[str, Any]) -> Path:
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    out_path = results_dir / "manifest.json"
    try:
        with out_path.open("w", encoding="utf-8") as f:
            json.dump(info, f, indent=2, sort_keys=True, default=str)
    except Exception:
        # As a last resort, write a minimal version
        try:
            with out_path.open("w", encoding="utf-8") as f:
                json.dump({k: str(v) for k, v in info.items()}, f, indent=2)
        except Exception:
            pass
    return out_path
