"""
Build a figures_manifest.json that maps each manuscript figure to its source run
(directory + run_id) and the model/ratio/era it represents.

Heuristics:
- Figures live under figures/figures/*.png
- We infer model/ratio/era and suffix from the filename, then locate the most
  recent run under results/<model>/<ratio>/<era>/ unless --run-id is provided.
- Known filename patterns:
  - *_3_1_1_both_v3_6.png  → model=bayesian_v3_6, ratio=3_1_1, era=both
  - *_3_1_1_pre_v3_6.png   → model=bayesian_v3_6, ratio=3_1_1, era=pre_modern
  - *_3_1_1_post_v3_6.png  → model=bayesian_v3_6, ratio=3_1_1, era=post_modern
  - ode_phase_ratios_3_1_1_both.png → model=enzyme_v3 (ode), ratio=3_1_1, era=both

You can override run selection with --force-run-id mapping per (model,ratio,era).
"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional, Dict, List, Tuple


@dataclass
class FigureEntry:
    filename: str
    model: str
    ratio: str
    era: str
    run_id: Optional[str]
    run_dir: Optional[str]


def parse_args():
    p = argparse.ArgumentParser(description="Build figures_manifest.json for manuscript figures")
    p.add_argument("--figdir", default="figures/figures", help="Directory with manuscript figures")
    p.add_argument("--out", default="figures_manifest.json", help="Output JSON path")
    # Optional explicit run id overrides in form model:ratio:era=RUNID
    p.add_argument("--force-run-id", action="append", default=[], help="Override run id mapping, e.g. bayesian_v3_6:3_1_1:both=20251123T223826Z_a8c794cc")
    return p.parse_args()


def latest_run_dir(base: Path) -> Optional[Path]:
    if not base.exists():
        return None
    runs = [p for p in base.iterdir() if p.is_dir()]
    if not runs:
        return None
    runs.sort()
    return runs[-1]


def infer_meta(filename: str) -> Optional[Tuple[str, str, str]]:
    # Return (model, ratio, era)
    name = filename
    if name.endswith("_3_1_1_both_v3_6.png"):
        return ("bayesian_v3_6", "3_1_1", "both")
    if name.endswith("_3_1_1_pre_v3_6.png"):
        return ("bayesian_v3_6", "3_1_1", "pre_modern")
    if name.endswith("_3_1_1_post_v3_6.png"):
        return ("bayesian_v3_6", "3_1_1", "post_modern")
    if name == "ode_phase_ratios_3_1_1_both.png" or name == "ode_phase_ratios_3_1_1_both_ode.png":
        return ("enzyme_v3", "3_1_1", "both")
    # Try broader patterns
    if "v3_6" in name and "3_1_1" in name:
        if "_pre_v3_6" in name:
            return ("bayesian_v3_6", "3_1_1", "pre_modern")
        if "_post_v3_6" in name:
            return ("bayesian_v3_6", "3_1_1", "post_modern")
        if "_both_v3_6" in name:
            return ("bayesian_v3_6", "3_1_1", "both")
    return None


def build_manifest(figdir: Path, overrides: Dict[Tuple[str, str, str], str]) -> List[FigureEntry]:
    entries: List[FigureEntry] = []
    for f in sorted(figdir.glob("*.png")):
        meta = infer_meta(f.name)
        if not meta:
            # Non-data/composite figure; include with composite tag and N/A fields
            entries.append(FigureEntry(filename=str(f), model="composite", ratio="n/a", era="n/a", run_id=None, run_dir=None))
            continue
        model, ratio, era = meta
        # Determine results base
        if model == "enzyme_v3":
            base = Path("results/enzyme_v3") / ratio / era
        else:
            base = Path("results") / model / ratio / era
        run_id = overrides.get((model, ratio, era))
        run_dir = None
        if run_id:
            rd = base / run_id
            if rd.exists():
                run_dir = str(rd)
            else:
                # fallback to latest
                latest = latest_run_dir(base)
                run_dir = str(latest) if latest else None
                run_id = latest.name if latest else None
        else:
            latest = latest_run_dir(base)
            run_dir = str(latest) if latest else None
            run_id = latest.name if latest else None
        entries.append(FigureEntry(filename=str(f), model=model, ratio=ratio, era=era, run_id=run_id, run_dir=run_dir))
    return entries


def main():
    args = parse_args()
    figdir = Path(args.figdir)
    figdir.mkdir(parents=True, exist_ok=True)

    overrides: Dict[Tuple[str, str, str], str] = {}
    for item in args.force_run_id:
        # format: model:ratio:era=RUNID
        try:
            lhs, runid = item.split("=", 1)
            model, ratio, era = lhs.split(":", 2)
            overrides[(model, ratio, era)] = runid
        except Exception:
            print(f"WARN: could not parse --force-run-id '{item}' (expected model:ratio:era=RUNID)")

    entries = build_manifest(figdir, overrides)
    out = Path(args.out)
    out.write_text(json.dumps([asdict(e) for e in entries], indent=2))
    print(f"Wrote {out} with {len(entries)} entries")


if __name__ == "__main__":
    main()
