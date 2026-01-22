#!/usr/bin/env python3
"""
File organization script for Zenodo package
Copies and organizes files from project and local sources to Zenodo structure
"""

import os
import shutil
from pathlib import Path
import json
import argparse
from datetime import datetime

# For generating figures manifest
try:
    from quantum.build_figures_manifest import build_manifest
except Exception:
    build_manifest = None

EXCLUDED_DIRS = {".venv", "venv", "__pycache__", ".mypy_cache", ".pytest_cache", ".ipynb_checkpoints"}
EXCLUDED_NAMES = {".DS_Store"}


def newer(src: Path, dst: Path) -> bool:
    try:
        return src.stat().st_mtime > dst.stat().st_mtime
    except FileNotFoundError:
        return True


def discover_and_copy(sources, target_dir: Path, dry_run: bool, force: bool, max_size_mb: int):
    """Discover scattered files by keywords/patterns and copy into Zenodo structure.

    Overwrite policy: overwrite if newer or if force=True.
    Size policy: skip files larger than max_size_mb (if > 0).
    """
    target_dir.mkdir(parents=True, exist_ok=True)

    # Ensure subdirs exist
    subdirs = [
        "01_manuscript",
        "02_data_raw",
        "03_data_processed",
        "04_model_outputs",
        "05_figures",
        "06_code",
        "07_cross_validation",
        "08_sensitivity",
        "09_supplementary",
        "10_quantum_hand",
        "11_local_development",
    ]
    for sd in subdirs:
        (target_dir / sd).mkdir(parents=True, exist_ok=True)

    # Define routing rules by keywords
    routes = [
        # Manuscript PDFs
        ({"noise_v2_8_main.pdf", "noise_v2_8_supplemental.pdf"}, None, "01_manuscript"),
        # v3_6 outputs / traces / predictions
        ({"v3_6", "posterior_predictive"}, None, "04_model_outputs"),
        ({"trace_v3_6"}, {".nc"}, "04_model_outputs"),
        # v4 related
        ({"v4"}, None, "04_model_outputs"),
        # v2 processed data
        ({"v2", "sse"}, None, "03_data_processed"),
        ({"v2", "array"}, None, "03_data_processed"),
        ({"tegmark"}, None, "03_data_processed"),
        # Figures
        ({"forest_hdi", "posterior_", "ppc_bar", "trace_plots", "group_constraints"}, {".png", ".pdf"}, "05_figures"),
        # Cross-validation and Valcour
        ({"cv_", "valcour"}, None, "07_cross_validation"),
        # Sensitivity
        ({"nopseudo", "sensitivity"}, None, "08_sensitivity"),
        # Supplementary tables
        ({"TableS", "Supplementary_Figure_"}, None, "09_supplementary"),
        # Code snippets of interest
        ({"bayesian_v3_6", "bayesian_enzyme_v4", "enzyme_kinetics", "final_calibrated_model_v2"}, {".py"}, "06_code"),
        ({"Makefile", "Makefile.txt", "README.md", "venv_info.txt", "run_info.json"}, None, "06_code"),
    ]

    copied = {sd: [] for sd in subdirs}
    skipped = []

    def match_route(fname_lower: str, suffix: str):
        for keywords, exts, dest in routes:
            if any(kw in fname_lower for kw in keywords):
                if exts is None or suffix in exts:
                    return dest
        return None

    for src_root in sources:
        src_root = Path(src_root).expanduser().resolve()
        if not src_root.exists():
            continue
        for dirpath, dirnames, filenames in os.walk(src_root):
            # prune excluded dirs in-place
            dirnames[:] = [d for d in dirnames if d not in EXCLUDED_DIRS]
            for name in filenames:
                if name in EXCLUDED_NAMES:
                    continue
                src_path = Path(dirpath) / name
                try:
                    if max_size_mb > 0 and src_path.stat().st_size > max_size_mb * 1024 * 1024:
                        skipped.append(str(src_path))
                        continue
                except FileNotFoundError:
                    continue
                lower = name.lower()
                dest_sub = match_route(lower, src_path.suffix.lower())
                if not dest_sub:
                    continue
                rel = name
                dst_path = target_dir / dest_sub / rel
                # If collision, keep original name but place; could add disambiguation if needed
                action = None
                if dst_path.exists():
                    if force or newer(src_path, dst_path):
                        action = "update"
                    else:
                        action = "skip"
                else:
                    action = "copy"

                if action == "skip":
                    continue
                if dry_run:
                    copied[dest_sub].append({"source": str(src_path), "dest": str(dst_path), "action": action})
                else:
                    dst_path.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(src_path, dst_path)
                    copied[dest_sub].append({"source": str(src_path), "dest": str(dst_path), "action": action})

    report = {
        "timestamp": datetime.now().isoformat(),
        "target": str(target_dir),
        "copied_counts": {k: len(v) for k, v in copied.items()},
        "copied": copied,
        "skipped_large": skipped,
        "max_size_mb": max_size_mb,
        "force": force,
        "dry_run": dry_run,
        "sources": [str(Path(s)) for s in sources],
    }

    report_path = target_dir / ("file_copy_report_dry_run.json" if dry_run else "file_copy_report.json")
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)

    print("\n" + "=" * 60)
    print("COPY SUMMARY")
    print(f"Target: {target_dir}")
    print(f"Total copied (by category): {report['copied_counts']}")
    if skipped:
        print(f"Skipped (size>{max_size_mb}MB): {len(skipped)}")
    print(f"Report saved to: {report_path}")

    return report


def create_supplementary_readme(target_dir: Path):
    """Create README for files from Quantum Hand project."""

    readme = """# Quantum Hand Project Data

## Overview
This directory contains data from the "Quantum Hand" project, which represents
the early v2 model development phase of the noise-mediated neuroprotection study.

## Contents

### Model Development
- Early iterations of the Bayesian model (v2.0-v2.5)
- Parameter optimization studies
- Convergence diagnostics from development phase

### Validation Studies
- Initial cross-validation attempts
- Parameter sensitivity explorations
- Model comparison studies

## Integration with Main Analysis
These files provide historical context for model development and demonstrate
the evolution of the analytical approach from initial hypothesis to final
validated model (v3.6).

## Note
Some files may contain intermediate results that were superseded by the
final analysis. They are included for completeness and reproducibility.
"""

    quantum_dir = target_dir / "10_quantum_hand"
    quantum_dir.mkdir(exist_ok=True)
    (quantum_dir / "README.md").write_text(readme)
    print("✓ Created Quantum Hand project directory")


def create_local_data_readme(target_dir: Path):
    """Create README for local GitHub data."""

    readme = """# Local Development Data

## Overview
This directory should contain data from:
~/Documents/GitHub/noise_decorrelation_HIV

## Expected Contents

### Code Repository
- Version control history
- Development notebooks
- Testing scripts
- Local configuration files

### Additional Analyses
- Exploratory data analyses
- Alternative model specifications
- Visualization experiments

## Instructions
To add local data:
1. Navigate to ~/Documents/GitHub/noise_decorrelation_HIV
2. Copy relevant files to this directory
3. Update the file manifest
4. Generate checksums

## Note
Ensure no sensitive or private data is included before upload.
"""

    local_dir = target_dir / "11_local_development"
    local_dir.mkdir(exist_ok=True)
    (local_dir / "README.md").write_text(readme)
    print("✓ Created local development directory")


def parse_args():
    parser = argparse.ArgumentParser(description="Organize files into Zenodo package")
    parser.add_argument("--target", type=Path, default=Path("/Users/acdmbpmax/Desktop/zenodo_package"), help="Target Zenodo package directory")
    parser.add_argument("--sources", type=Path, nargs="*", default=None, help="Source roots to scan (defaults to project root and Results_v2)")
    parser.add_argument("--dry-run", action="store_true", help="Don’t copy, only report planned actions")
    parser.add_argument("--force", action="store_true", help="Force overwrite even if destination is newer")
    parser.add_argument("--max-size-mb", type=int, default=800, help="Skip files larger than this size (MB). 0 disables size check")
    # One-step pinning for figures manifest: model:ratio:era=RUNID (can be passed multiple times)
    parser.add_argument("--force-run-id", action="append", default=[], help="Override run id mapping for figures manifest, e.g. bayesian_v3_6:3_1_1:both=20251125T160846Z_ce3b0657; can be repeated")
    return parser.parse_args()


if __name__ == "__main__":
    print("=" * 60)
    print("ORGANIZING FILES FOR ZENODO PACKAGE")
    print("=" * 60)

    args = parse_args()

    # Default sources: this repo root and the user’s Documents repo plus Results_v2
    repo_root = Path(__file__).resolve().parent
    default_sources = [
        repo_root,
        Path("/Users/acdmbpmax/Documents/GitHub/noise_decorrelation_HIV"),
        Path("/Users/acdmbpmax/Documents/GitHub/noise_decorrelation_HIV/Results_v2"),
    ]
    sources = args.sources if args.sources else default_sources

    # Discover and copy
    report = discover_and_copy(sources, args.target, args.dry_run, args.force, args.max_size_mb)

    # Create additional directories/readmes
    create_supplementary_readme(args.target)
    create_local_data_readme(args.target)

    # Generate figures manifest inside 05_figures
    figures_dir = args.target / "05_figures"
    try:
        if build_manifest is None:
            print("WARN: Could not import quantum.build_figures_manifest; skipping figures manifest generation.")
        else:
            print("Generating figures_manifest.json in 05_figures ...")
            # Parse --force-run-id list into overrides dict: (model, ratio, era) -> run_id
            overrides = {}
            for item in args.force_run_id:
                try:
                    lhs, runid = item.split("=", 1)
                    model, ratio, era = lhs.split(":", 2)
                    overrides[(model, ratio, era)] = runid
                except Exception:
                    print(f"WARN: could not parse --force-run-id '{item}' (expected model:ratio:era=RUNID)")
            entries = build_manifest(figures_dir, overrides=overrides)
            manifest_path = figures_dir / "figures_manifest.json"
            manifest_path.write_text(json.dumps([e.__dict__ for e in entries], indent=2))
            print(f"✓ Wrote {manifest_path} with {len(entries)} entries")
    except Exception as e:
        print(f"WARN: Failed to generate figures manifest: {e}")

    print("\n" + "=" * 60)
    print("FILE ORGANIZATION COMPLETE")
    print("=" * 60)
    print("\nRemaining manual steps:")
    print("1. Verify that v2 SSE/array/Tegmark and v3_6 traces/results are present")
    print("2. Add any missing Quantum Hand or local development data as needed")
    print("3. Run zenodo.py --out /Users/acdmbpmax/Desktop/zenodo_package --checksums [--zip]")
    print("4. Review file_copy_report.json for any misses/exclusions")