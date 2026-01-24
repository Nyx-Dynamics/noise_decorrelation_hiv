#!/usr/bin/env python3
"""
Data Collector for HIV Noise Neuroprotection Study
Gathers scattered files from multiple locations into Zenodo package
"""

import os
import shutil
from pathlib import Path
from datetime import datetime


class DataCollector:
    def __init__(self):
        self.zenodo_dir = Path("~/PycharmProjects/noise_decorrelation_HIV/zenodo_package").expanduser()
        self.found_files = []
        self.missing_files = []

        # Define search locations
        self.search_paths = [
            Path("~/PycharmProjects/noise_decorrelation_HIV").expanduser(),
            Path("~/Documents/GitHub/noise_decorrelation_HIV").expanduser(),
            Path("~/Desktop").expanduser(),
            Path("~/Downloads").expanduser(),
            Path("~/Library/Mobile Documents/com~apple~CloudDocs").expanduser(),  # iCloud
            Path("/Volumes/").expanduser(),  # External drives
        ]

        print("🔍 DATA COLLECTOR FOR ZENODO PACKAGE")
        print("=" * 50)
        print(f"Target: {self.zenodo_dir}")
        print()

    def find_manuscript_files(self):
        """Locate manuscript PDFs and LaTeX files."""
        print("📄 SEARCHING FOR MANUSCRIPT FILES...")
        print("-" * 40)

        manuscript_files = [
            "noise_v2_8_main.pdf",
            "noise_v2_8_supplemental.pdf",
            "noise_v2_7_main.pdf",  # Earlier versions
            "noise_v2_7_supplemental.pdf",
            "noise_v2_1_main.tex",
            "noise_v2_1_supplemental.tex"
        ]

        for filename in manuscript_files:
            self.search_file(filename, "01_manuscript")

    def find_data_files(self):
        """Locate all CSV data files."""
        print("\n📊 SEARCHING FOR DATA FILES...")
        print("-" * 40)

        # Raw data tables
        raw_data = ["TableS1*.csv", "TableS2*.csv", "TableS3*.csv"]
        for pattern in raw_data:
            self.search_pattern(pattern, "02_data_raw")

        # Processed data
        processed_patterns = [
            "group_inputs*.csv",
            "group_likelihood*.csv",
            "summary_*.csv",
            "results_v3_6*.csv",
            "posterior_predictive*.csv",
            "valcour*.csv"
        ]
        for pattern in processed_patterns:
            self.search_pattern(pattern, "03_data_processed")

    def find_code_files(self):
        """Locate Python scripts and notebooks."""
        print("\n💻 SEARCHING FOR CODE FILES...")
        print("-" * 40)

        code_files = [
            "bayesian_v3_6_corrected_local.py",
            "bayesian_enzyme_v4.py",
            "enzyme_kinetics.py",
            "final_calibrated_model_v2.py",
            "Makefile",
            "requirements.txt",
            "*.ipynb"
        ]

        for pattern in code_files:
            self.search_pattern(pattern, "06_code")

    def find_figures(self):
        """Locate all figure files."""
        print("\n🎨 SEARCHING FOR FIGURES...")
        print("-" * 40)

        figure_patterns = [
            "forest_hdi*.png",
            "posterior*.png",
            "trace_plots*.png",
            "ppc_bar*.png",
            "group_constraints*.png",
            "cv_*.png",  # Cross-validation → 07
            "*nopseudo*.png"  # Sensitivity → 08
        ]

        for pattern in figure_patterns:
            if "cv_" in pattern:
                target = "07_cross_validation"
            elif "nopseudo" in pattern:
                target = "08_sensitivity"
            else:
                target = "05_figures"

            self.search_pattern(pattern, target)

    def search_file(self, filename, target_dir):
        """Search for a specific file."""
        for search_path in self.search_paths:
            if not search_path.exists():
                continue

            # Search recursively
            for file_path in search_path.rglob(filename):
                if file_path.is_file():
                    print(f"  ✓ Found: {file_path}")
                    self.found_files.append((file_path, target_dir))
                    return

        print(f"  ✗ Not found: {filename}")
        self.missing_files.append(filename)

    def search_pattern(self, pattern, target_dir):
        """Search for files matching a pattern."""
        found_count = 0

        for search_path in self.search_paths:
            if not search_path.exists():
                continue

            # Handle wildcards
            if "*" in pattern:
                for file_path in search_path.glob(pattern):
                    if file_path.is_file():
                        if found_count == 0:
                            print(f"  ✓ Found {pattern}:")
                        print(f"     {file_path}")
                        self.found_files.append((file_path, target_dir))
                        found_count += 1

        if found_count == 0:
            print(f"  ✗ No files matching: {pattern}")

    def copy_to_zenodo(self, dry_run=True):
        """Copy found files to Zenodo package."""
        if dry_run:
            print("\n📋 DRY RUN - Would copy these files:")
            print("-" * 40)
        else:
            print("\n📦 COPYING FILES TO ZENODO PACKAGE...")
            print("-" * 40)

        for source_path, target_dir in self.found_files:
            target_path = self.zenodo_dir / target_dir / source_path.name

            if dry_run:
                print(f"  {source_path.name} → {target_dir}/")
            else:
                target_path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source_path, target_path)
                print(f"  ✓ Copied {source_path.name} → {target_dir}/")

        if dry_run:
            print(f"\nTotal files to copy: {len(self.found_files)}")
            print("Run with --copy to actually copy files")

    def generate_report(self):
        """Generate a collection report."""
        report_path = Path("zenodo_collection_report.txt")

        with open(report_path, 'w') as f:
            f.write("ZENODO DATA COLLECTION REPORT\n")
            f.write(f"Generated: {datetime.now()}\n")
            f.write("=" * 50 + "\n\n")

            f.write(f"FOUND FILES ({len(self.found_files)}):\n")
            for source, target in self.found_files:
                f.write(f"  {source} → {target}/\n")

            f.write(f"\nMISSING FILES ({len(self.missing_files)}):\n")
            for filename in self.missing_files:
                f.write(f"  - {filename}\n")

            f.write("\nSEARCH LOCATIONS:\n")
            for path in self.search_paths:
                status = "✓" if path.exists() else "✗"
                f.write(f"  {status} {path}\n")

        print(f"\n📄 Report saved to: {report_path}")


def main():
    import sys

    collector = DataCollector()

    # Search for all file types
    collector.find_manuscript_files()
    collector.find_data_files()
    collector.find_code_files()
    collector.find_figures()

    # Check if --copy flag is present
    if "--copy" in sys.argv:
        collector.copy_to_zenodo(dry_run=False)
    else:
        collector.copy_to_zenodo(dry_run=True)

    # Generate report
    collector.generate_report()

    print("\n" + "=" * 50)
    print("✨ COLLECTION COMPLETE!")
    print("=" * 50)

    if len(collector.missing_files) > 0:
        print(f"\n⚠️  Some files not found. Check report for details.")

    if "--copy" not in sys.argv:
        print("\n💡 TIP: Run with --copy flag to actually copy files:")
        print("   python3 collect_data.py --copy")


if __name__ == "__main__":
    main()