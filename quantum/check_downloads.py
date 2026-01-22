#!/usr/bin/env python3
"""
Check what HIV MRS related files are in Downloads
Helps identify what you have before organizing
"""

from pathlib import Path

DOWNLOADS = Path.home() / "Downloads"

# Keywords to look for in filenames
KEYWORDS = [
    'NAA', 'MRS', 'HIV', 'CHANG', 'SAILASUTA', 'VALCOUR', 'YOUNG',
    'DAHMANI', 'Table', 'posterior', 'summary', 'MASTER', 'CRITICAL'
]


def main():
    print("=" * 60)
    print("Checking Downloads for HIV MRS Data Files")
    print("=" * 60)
    print(f"Location: {DOWNLOADS}\n")

    if not DOWNLOADS.exists():
        print(f"❌ Downloads folder not found!")
        return

    # Get all files in Downloads
    all_files = list(DOWNLOADS.glob('*'))

    # Filter for relevant files
    relevant = []
    for f in all_files:
        if f.is_file():
            name = f.name.upper()
            if any(keyword.upper() in name for keyword in KEYWORDS):
                relevant.append(f)

    if not relevant:
        print("❌ No MRS/HIV related files found in Downloads\n")
        print("Expected files like:")
        print("  - CHANG_2002_EXTRACTED.csv")
        print("  - VALCOUR_2015_INDIVIDUAL_PATIENTS.csv")
        print("  - MASTER_HIV_MRS_DATABASE_v2.csv")
        print("  - etc.")
        return

    print(f"✓ Found {len(relevant)} relevant files:\n")

    # Group by type
    csvs = [f for f in relevant if f.suffix == '.csv']
    excels = [f for f in relevant if f.suffix in ['.xlsx', '.xls']]
    docs = [f for f in relevant if f.suffix in ['.docx', '.DOCX']]
    others = [f for f in relevant if f not in csvs + excels + docs]

    if csvs:
        print("CSV FILES:")
        for f in sorted(csvs):
            size_kb = f.stat().st_size / 1024
            print(f"  ✓ {f.name} ({size_kb:.1f} KB)")

    if excels:
        print("\nEXCEL FILES:")
        for f in sorted(excels):
            size_kb = f.stat().st_size / 1024
            print(f"  ✓ {f.name} ({size_kb:.1f} KB)")

    if docs:
        print("\nDOCX/TABLE FILES:")
        for f in sorted(docs):
            size_kb = f.stat().st_size / 1024
            print(f"  ✓ {f.name} ({size_kb:.1f} KB)")

    if others:
        print("\nOTHER FILES:")
        for f in sorted(others):
            print(f"  • {f.name}")

    print("\n" + "=" * 60)
    print("Ready to organize these files!")
    print("Run: python organize_from_downloads.py")


if __name__ == "__main__":
    main()