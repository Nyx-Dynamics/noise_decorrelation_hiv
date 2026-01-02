#!/usr/bin/env python3
"""
Bayesian Model v3.6 - EXPANDED WITH RATIO DATA
Including Valcour (absolute) + Young + Sailasuta (ratios converted to absolute)

Updated: November 16, 2025
Investigator: AC (Nyx Dynamics LLC)
"""

import numpy as np
import pandas as pd
import pymc as pm
import arviz as az
import matplotlib.pyplot as plt
from scipy import stats
import argparse
from datetime import datetime
import json
import platform
from pathlib import Path
import warnings
import sys

warnings.filterwarnings('ignore')

print("\n" + "╔" + "=" * 78 + "╗")
print("║" + " " * 15 + "BAYESIAN MODEL v3.6 - WITH RATIO DATA" + " " * 22 + "║")
print("║" + " " * 78 + "║")
print("║" + " Valcour + Chang (absolute → ratios) + Young + Sailasuta (ratios)" + " " * 7 + "║")
print("╚" + "=" * 78 + "╝\n")

# ============================================================================
# FILE PATH CONFIGURATION
# ============================================================================

# Auto-detect if we're running locally or on Claude's system
script_dir = Path(__file__).resolve().parent

# Try to find the data directory
possible_data_dirs = [
    script_dir / "data" / "curated",  # If script is in quantum/quantum/
    script_dir.parent / "data" / "curated",  # If script is in quantum/
    script_dir.parent.parent / "data" / "curated",  # One more level up
    Path("/mnt/user-data/uploads"),  # Claude's system
]

data_dir = None
for possible_dir in possible_data_dirs:
    if possible_dir.exists():
        data_dir = possible_dir
        print(f"✅ Found data directory: {data_dir}")
        break

if data_dir is None:
    print("❌ Could not find data directory. Please set manually.")
    print("Looking for:")
    for pd in possible_data_dirs:
        print(f"  - {pd}")
    sys.exit(1)

# Define file paths
VALCOUR_DIR = data_dir / "absolute"
RATIO_DIR = data_dir / "ratio"

VALCOUR_FILES = [
    VALCOUR_DIR / "valcour_2015_week_0.xlsx",
    VALCOUR_DIR / "valcour_2015_week_4.xlsx",
    VALCOUR_DIR / "valcour_2015_week_12.xlsx",
    VALCOUR_DIR / "valcour_2015_week_24.xlsx",
]

YOUNG_FILE = RATIO_DIR / "YOUNG_2014_CROSS_SECTIONAL_DATA.csv"
SAILASUTA_FILE = RATIO_DIR / "Sailasuta_2012.csv"
CHANG_FILE = data_dir / "absolute" / "CHANG_2002_EXTRACTED.csv"
EXTRA_GROUP_FILE = RATIO_DIR / "all_extracted_hiv_neurodata.csv"
# New combined BG meta-analysis file (Winston et al. + Dahmani supplement)
COMBINED_BG_FILE = RATIO_DIR / "bg_combined_with_winston_dahmani.csv"

print(f"\n📂 Data paths configured:")
print(f"   Absolute data: {VALCOUR_DIR}")
print(f"   Ratio data: {RATIO_DIR}")

# ============================================================================
# REFERENCE CREATINE VALUES (used for abs↔ratio conversions when needed)
# ============================================================================
CR_REFERENCE = {
    'BG': 8.0,  # Basal Ganglia
    'FWM': 6.8,  # Frontal White Matter
    'PGM': 7.8,  # Posterior Gray Matter / Parietal
    'FGM': 7.8,  # Frontal Gray Matter
    'AC': 7.8,  # Anterior Cingulate
    'OGM': 7.5  # Occipital Gray Matter
}


def convert_ratio_to_absolute(ratio_value, region):
    """Convert metabolite/Cr ratio to absolute concentration (mM)"""
    cr_ref = CR_REFERENCE.get(region, 7.5)
    return ratio_value * cr_ref

def convert_absolute_to_ratio(abs_value, region):
    """Convert absolute concentration (mM) to metabolite/Cr ratio using reference Cr."""
    cr_ref = CR_REFERENCE.get(region, 7.5)
    if pd.isna(abs_value):
        return np.nan
    return abs_value / cr_ref


# ============================================================================
# LOAD DATA
# ============================================================================

print("\n" + "=" * 80)
print("LOADING DATA")
print("=" * 80)

# ----------------------------------------------------------------------------
# CLI ARGUMENTS
# ----------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Bayesian v3.6 (ratio scale) with optional Valcour ablation, tagging, NetCDF saving, and density/extended plots, and timestamped outputs")
parser.add_argument("--exclude-valcour", action="store_true",
                    help="Exclude Valcour individual acute observations from the ACUTE pool (Young/Sailasuta kept).")
parser.add_argument("--tag", type=str, default="",
                    help="Optional tag to suffix output filenames for side-by-side comparisons (e.g., with_valcour, no_valcour).")
parser.add_argument("--run-label", type=str, default="",
                    help="Optional label for this run; if provided, used for the run folder name instead of timestamp.")
parser.add_argument("--no-timestamp", action="store_true",
                    help="Disable timestamp suffixing of artifact filenames and run folder nesting.")
parser.add_argument("--no-save-trace", action="store_true",
                    help="Disable saving full posterior trace to NetCDF (enabled by default).")
parser.add_argument("--save-posterior-csv", action="store_true",
                    help="Additionally save a CSV of flattened posterior samples for key parameters.")
parser.add_argument("--plot-densities", action="store_true",
                    help="Generate ArviZ posterior KDE plots for key parameters and save to figures directory.")
parser.add_argument("--plots-extended", action="store_true",
                    help="Generate extended figures: PPC comparison, trace plots, forest plots, Valcour degradation ribbon (if available).")
parser.add_argument("--plots-extra", action="store_true",
                    help="Alias for --plots-extended.")
parser.add_argument("--compare-trace", type=str, default="",
                    help="Optional path to a NetCDF trace file to overlay densities against (e.g., the ablation run).")
parser.add_argument("--validate-valcour", action="store_true",
                    help="Run a Valcour-manuscript validation pass (week and control comparisons, p-values, Fig.2-style plots) and exit.")
parser.add_argument("--fig2", action="store_true",
                    help="Force generation of Fig.2-style plots (used with --validate-valcour; implied when --validate-valcour is set).")
parser.add_argument("--week-pvals", action="store_true",
                    help="Print detailed week-level statistical test results during validation.")
parser.add_argument("--valcour-aux", type=str, default="off", choices=["off","time","time+vl","plasma","csf","both"],
                    help=(
                        "Run Valcour-only auxiliary model(s) on weekly means (fully decoupled):\n"
                        "  off     : disabled (default)\n"
                        "  time    : week effect only (no VL)\n"
                        "  time+vl : legacy alias for plasma arm with VL\n"
                        "  plasma  : plasma VL arm (logpVL preferred; fallback pVL→log10)\n"
                        "  csf     : CSF VL arm (logcVL preferred; fallback cVL→log10)\n"
                        "  both    : run plasma and CSF arms sequentially"
                    ))
parser.add_argument("--regions", type=str, default="all", choices=["BG","all"],
                    help=(
                        "Scope of curated group constraints: \n"
                        "  BG  : restrict group-mean constraints to Basal Ganglia only (BG).\n"
                        "  all : include group-mean constraints from all regions (default)."
                    ))
parser.add_argument("--exclude-region", type=str, default="",
                    help=(
                        "Exclude a specific region from curated group constraints (case-insensitive). "
                        "Examples: --exclude-region BG, --exclude-region FWM, --exclude-region FGM, --exclude-region AC."
                    ))
parser.add_argument("--no-acute-pseudo", action="store_true",
                    help=(
                        "Do not replicate acute Young/Sailasuta BG means into pseudo-observations. "
                        "Instead, use proper group-mean Normal likelihoods on the acute phase with reported SEs."
                    ))
parser.add_argument("--exclude-source", type=str, default="",
                    help=(
                        "Exclude curated group constraints from a specific source/study name (case-insensitive substring). "
                        "Examples: --exclude-source Young_2014, --exclude-source Sailasuta_2012, --exclude-source Mohamed, --exclude-source Boban, --exclude-source Dahmani, --exclude-source Chang."
                    ))
args, unknown = parser.parse_known_args()

if args.exclude_valcour:
    print("\n🧪 Configuration: EXCLUDING Valcour acute individuals from ACUTE pool (for ablation)")
else:
    print("\n🧪 Configuration: Including Valcour acute individuals in ACUTE pool (default)")

output_tag = ("_" + args.tag.strip()) if args.tag and len(args.tag.strip()) > 0 else ""
if output_tag:
    print(f"   Output files will include tag: '{output_tag}'")
print(f"   Save NetCDF trace: {'NO' if args.no_save_trace else 'YES'}")
if args.plot_densities:
    print("   Will generate posterior density plots (ArviZ KDE)")
if args.compare_trace:
    print(f"   Will overlay densities against: {args.compare_trace}")
if args.plots_extended or args.plots_extra:
    print("   Will generate extended plots (PPC, trace, forest, Valcour ribbon if available)")
if args.valcour_aux and args.valcour_aux != "off":
    print(f"   Valcour auxiliary analysis enabled: {args.valcour_aux} (fully decoupled from primary model)")
print(f"   Group constraints scope: {'BG-only' if args.regions=='BG' else 'multi-region (all)'}")
if args.exclude_region.strip():
    print(f"   Excluding curated region from constraints: {args.exclude_region.strip()}")
if args.no_acute_pseudo:
    print("   Acute pseudo-observations: DISABLED (use group-mean likelihoods for acute Young/Sailasuta)")
if args.exclude_source.strip():
    print(f"   Excluding curated source from constraints: {args.exclude_source.strip()}")

# Timestamp/run labeling
run_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
use_timestamp = not args.no_timestamp
run_name = args.run_label.strip() if args.run_label.strip() else (f"run_{run_timestamp}" if use_timestamp else (args.tag.strip() if args.tag.strip() else "run"))
print(f"   Run identifier: {run_name}")

# Prepare results directories early if running validation/figures without full inference
base_results_dir = script_dir / "results_v3_6"
runs_dir = base_results_dir / "runs"
run_dir = runs_dir / run_name
figures_dir = run_dir / "figures"
if args.validate_valcour or args.fig2:
    base_results_dir.mkdir(exist_ok=True)
    runs_dir.mkdir(exist_ok=True)
    run_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)

# ------------------------------
# Load Valcour data (absolute concentrations)
# ------------------------------
print("\n📁 Loading Valcour 2015 data (absolute concentrations)...")

valcour_dfs = []
for week_file in VALCOUR_FILES:
    if week_file.exists():
        try:
            df = pd.read_excel(str(week_file))
            week_name = week_file.stem.replace('valcour_2015_', '')
            df['Week'] = week_name
            valcour_dfs.append(df)
            print(f"   ✅ Loaded {week_name}: {len(df)} patients")
        except Exception as e:
            print(f"   ⚠ Error reading {week_file.name}: {e}")
    else:
        print(f"   ⚠ File not found: {week_file.name}")

if valcour_dfs:
    valcour_df = pd.concat(valcour_dfs, ignore_index=True)
    print(f"\n✅ Total Valcour records: {len(valcour_df)}")
else:
    print("   ⚠ No Valcour data loaded")
    valcour_df = pd.DataFrame()

# ------------------------------
# Load Young 2014 data (ratios)
# ------------------------------
print("\n📁 Loading Young 2014 data (ratios, keep as ratios)...")

if YOUNG_FILE.exists():
    try:
        young_df = pd.read_csv(str(YOUNG_FILE))

        # Parse the ratio data - extract NAA/Cr and Cho/Cr as-is
        young_processed = []

        for _, row in young_df.iterrows():
            if row['Metabolite'] == 'NAA/Cr':
                young_processed.append({
                    'Study': 'Young_2014',
                    'Phase': 'Acute' if row['Phase'] == 'Primary' else row['Phase'],
                    'Region': row['Region'],
                    'n': row['n'],
                    'NAA_ratio': row['Ratio_Median'],
                    'NAA_SE_ratio': row['SE'],
                    'Source': 'Ratio'
                })
            elif row['Metabolite'] == 'Cho/Cr':
                # Find matching NAA entry to update
                for proc in young_processed:
                    if (proc['Region'] == row['Region'] and
                            proc['Phase'] == ('Acute' if row['Phase'] == 'Primary' else row['Phase']) and
                            proc['n'] == row['n']):
                        proc['Cho_ratio'] = row['Ratio_Median']
                        proc['Cho_SE_ratio'] = row['SE']
                        break

        young_clean = pd.DataFrame(young_processed)
        print(f"✅ Loaded Young 2014: {len(young_clean)} group means")
        print(f"   Phases: {young_clean['Phase'].value_counts().to_dict()}")
        print(f"   Regions: {young_clean['Region'].value_counts().to_dict()}")

    except Exception as e:
        print(f"⚠ Error loading Young 2014: {e}")
        young_clean = pd.DataFrame()
else:
    print(f"⚠ File not found: {YOUNG_FILE}")
    young_clean = pd.DataFrame()

# ------------------------------
# Load Sailasuta 2012 data (ratios)
# ------------------------------
print("\n📁 Loading Sailasuta 2012 data (ratios, keep as ratios)...")

if SAILASUTA_FILE.exists():
    try:
        sailasuta_df = pd.read_csv(str(SAILASUTA_FILE), skiprows=1)

        # Process Sailasuta data (means reported as ratios)
        sailasuta_processed = []

        for phase in sailasuta_df['Group'].unique():
            for region in sailasuta_df['Brain_Region'].unique():
                phase_region_data = sailasuta_df[
                    (sailasuta_df['Group'] == phase) &
                    (sailasuta_df['Brain_Region'] == region)
                    ]

                naa_row = phase_region_data[phase_region_data['Metabolite'] == 'NAA']
                cho_row = phase_region_data[phase_region_data['Metabolite'] == 'tCho']

                if not naa_row.empty:
                    naa_value = naa_row['Value'].values[0]  # already a ratio
                    naa_sd = naa_row['SD'].values[0]
                    n = int(naa_row['n'].values[0])

                    naa_ratio = naa_value
                    naa_se_ratio = naa_sd / np.sqrt(n)

                    entry = {
                        'Study': 'Sailasuta_2012',
                        'Phase': phase,
                        'Region': region,
                        'n': n,
                        'NAA_ratio': naa_ratio,
                        'NAA_SE_ratio': naa_se_ratio,
                        'Source': 'Ratio'
                    }

                    if not cho_row.empty:
                        cho_value = cho_row['Value'].values[0]
                        cho_sd = cho_row['SD'].values[0]
                        entry['Cho_ratio'] = cho_value
                        entry['Cho_SE_ratio'] = cho_sd / np.sqrt(n)

                    sailasuta_processed.append(entry)

        sailasuta_clean = pd.DataFrame(sailasuta_processed)
        print(f"✅ Loaded Sailasuta 2012: {len(sailasuta_clean)} group means")
        print(f"   Phases: {sailasuta_clean['Phase'].value_counts().to_dict()}")
        print(f"   Regions: {sailasuta_clean['Region'].value_counts().to_dict()}")

    except Exception as e:
        print(f"⚠ Error loading Sailasuta 2012: {e}")
        sailasuta_clean = pd.DataFrame()
else:
    print(f"⚠ File not found: {SAILASUTA_FILE}")
    sailasuta_clean = pd.DataFrame()

# ------------------------------
# Load Chang 2002 reference data (absolute)
# ------------------------------
print("\n📁 Loading Chang 2002 reference data (absolute → ratios)...")

if CHANG_FILE.exists():
    try:
        chang_df = pd.read_csv(str(CHANG_FILE))
        print(f"✅ Loaded Chang 2002: {len(chang_df)} measurements")
    except Exception as e:
        print(f"⚠ Error loading Chang 2002: {e}")
        chang_df = pd.DataFrame()
else:
    print(f"⚠ File not found: {CHANG_FILE}")
    chang_df = pd.DataFrame()

# ------------------------------
# Load additional curated group-level BG data (ratios)
# Robust to schema variants in all_extracted_hiv_neurodata.csv
# ------------------------------
print("\n📁 Loading additional curated group-level data (ratios, optional)...")
extra_group_df = pd.DataFrame()
extra_group_registry = []

def _parse_se_from_token(token: str) -> float:
    try:
        # Accept 'SE=0.03', 'SD=0.09', '0.03', 0.03
        if pd.isna(token):
            return np.nan
        if isinstance(token, (int, float)):
            return float(token)
        s = str(token).strip()
        if '=' in s:
            s = s.split('=', 1)[1]
        # remove non-numeric chars except . and e-+
        s = ''.join(ch for ch in s if (ch.isdigit() or ch in ['.','e','E','-','+']))
        return float(s) if s not in ('', None) else np.nan
    except Exception:
        return np.nan

def _canon_phase(ph: str) -> str:
    if pd.isna(ph):
        return ''
    s = str(ph).strip()
    # remove parenthetical qualifiers, keep first token
    if '(' in s:
        s = s.split('(', 1)[0].strip()
    s = s.title()
    # map common variants
    if s.startswith('Acute') or s.startswith('Early') or s.startswith('Primary'):
        return 'Acute'
    if s.startswith('Chronic'):
        return 'Chronic'
    if s.startswith('Control') or s in ['Ctrl','Controls']:
        return 'Control'
    return s

def _canon_metabolite(m: str) -> str:
    if pd.isna(m):
        return ''
    s = str(m).replace(' ', '').upper()
    # Normalize several variants
    if s in ['NAA/CR', 'NAA_CR', 'NAA\CR', 'NAA']:
        return 'NAA/Cr'
    if s in ['CHO/CR', 'CHO_CR', 'TCHO/CR', 'TCHO', 'CHO']:
        return 'Cho/Cr'
    return ''  # unsupported metrics ignored

def _canon_region(r: str) -> str:
    """Canonicalize region names for group-level ingestion.
    Option A: allow multi-region constraints (not just BG).
    Returns short upper-case tokens such as BG, FWM, PGM, AC, OGM, PCC, etc.
    Unknown regions return their upper-case token so they can be listed/audited,
    but non-NAA/Cr metrics will still be excluded by _canon_metabolite.
    """
    if pd.isna(r):
        return ''
    s = str(r).strip()
    su = s.upper()
    # Common normalizations
    if su in ['BASAL GANGLIA', 'BASALGANGLIA']:
        return 'BG'
    if su in ['FRONTAL WHITE MATTER', 'FRONTALWM', 'F WM', 'FRONTAL WM']:
        return 'FWM'
    if su in ['POSTERIOR GRAY MATTER', 'PARIETAL GM', 'PARIETAL', 'PARIETAL CORTEX GM', 'PGM', 'PCG', 'PCC', 'POSTERIOR CINGULATE CORTEX']:
        return 'PGM'
    if su in ['ANTERIOR CINGULATE', 'AC', 'ACC']:
        return 'AC'
    if su in ['OCCIPITAL GRAY MATTER', 'OGM', 'OCCIPITAL GM']:
        return 'OGM'
    if su in ['FRONTAL CORTEX', 'FC']:
        return 'FC'
    if su in ['FRONTAL SUBCORTICAL WHITE MATTER', 'FSWM', 'FSWM (RIGHT)', 'FSWM RIGHT']:
        return 'FSWM'
    if su.startswith('BG'):
        return 'BG'
    # Default to upper-case short token (strip spaces and punctuation)
    return ''.join(ch for ch in su if ch.isalnum())

def _try_load_extra_csv(path: Path, source_name: str = None) -> pd.DataFrame:
    try:
        df_raw = pd.read_csv(str(path))
    except Exception as e:
        print(f"   ⚠ Error loading {path.name}: {e}")
        return pd.DataFrame()

    # Harmonize column names
    cols = {c: c for c in df_raw.columns}
    rename_map = {}
    for c in df_raw.columns:
        cu = c.strip().lower()
        if cu in ['phase', 'hiv phase', 'group']:
            rename_map[c] = 'Phase'
        elif cu in ['region', 'brain region']:
            rename_map[c] = 'Region'
        elif cu in ['metabolite', 'analyte']:
            rename_map[c] = 'Metabolite'
        elif cu in ['mean', 'value', 'ratio', 'mean_ratio']:
            rename_map[c] = 'Mean'
        elif cu in ['se', 'sem', 'se/sd', 'se_sd', 'se/ sd', 'se/ sd ', 'se/SD', 'se/sd', 'se\sd', 'se or sd', 'se/sd/n']:
            rename_map[c] = 'SE_SD'
        elif cu in ['se/sd']:
            rename_map[c] = 'SE_SD'
        elif cu in ['n', 'N']:
            rename_map[c] = 'n'
        elif cu in ['study', 'source']:
            rename_map[c] = 'Study'
        elif cu in ['notes']:
            rename_map[c] = 'Notes'
    df = df_raw.rename(columns=rename_map).copy()

    # Create unified fields
    df['Study'] = df.get('Study', pd.Series(['']*len(df))).astype(str)
    df['Phase'] = df.get('Phase', pd.Series(['']*len(df))).apply(_canon_phase)
    df['Region'] = df.get('Region', pd.Series(['']*len(df))).apply(_canon_region)
    df['Metabolite'] = df.get('Metabolite', pd.Series(['']*len(df))).apply(_canon_metabolite)
    # Mean
    df['Mean'] = pd.to_numeric(df.get('Mean', pd.Series([np.nan]*len(df))), errors='coerce')
    # n
    df['n'] = pd.to_numeric(df.get('n', pd.Series([np.nan]*len(df))), errors='coerce')
    # SE or SD
    # Support both 'SE' column and the combined 'SE/SD' column from the uploaded file
    if 'SE' in df.columns:
        df['SE_raw'] = df['SE']
    elif 'SE_SD' in df.columns:
        df['SE_raw'] = df['SE_SD']
    else:
        df['SE_raw'] = np.nan

    se_vals = []
    included_flags = []
    reasons = []
    for i, row in df.iterrows():
        se = _parse_se_from_token(row.get('SE_raw', np.nan))
        reason = ''
        inc = False
        # If token indicated SD (e.g., original had 'SD=0.09'), compute SE using n
        token = row.get('SE_raw', '')
        token_s = str(token) if not pd.isna(token) else ''
        if (pd.isna(se) or se <= 0) and ('SD' in token_s.upper()):
            # Try to parse numeric again (already done); if still nan, cannot compute
            sd = _parse_se_from_token(token_s)
            n = row.get('n', np.nan)
            if np.isfinite(sd) and np.isfinite(n) and n > 1:
                se = float(sd) / float(np.sqrt(n))
                reason = 'SE computed from SD/sqrt(n)'
        # Basic inclusion criteria
        if row['Region'] == 'BG' and row['Metabolite'] in ['NAA/Cr', 'Cho/Cr'] and np.isfinite(row['Mean']):
            # Phase must be one of the modeled phases
            phase_ok = row['Phase'] in ['Acute', 'Chronic', 'Control']
            if phase_ok:
                if not np.isfinite(se) or se <= 0:
                    reason = reason + ('; ' if reason else '') + 'Missing/invalid SE'
                    inc = False
                else:
                    inc = True
            else:
                reason = 'Unsupported phase'
        else:
            reason = 'Unsupported region/metabolite or missing mean'
        se_vals.append(se if np.isfinite(se) else np.nan)
        included_flags.append(bool(inc))
        extra_group_registry.append({
            'Study': row.get('Study',''),
            'Phase': row.get('Phase',''),
            'Region': row.get('Region',''),
            'Metabolite': row.get('Metabolite',''),
            'Mean': row.get('Mean', np.nan),
            'SE_raw': row.get('SE_raw', np.nan),
            'SE_final': se if np.isfinite(se) else np.nan,
            'n': row.get('n', np.nan),
            'Included': inc,
            'Reason': reason,
            'SourceFile': (source_name if source_name else path.name),
            'DuplicateKey': '',
            'DedupReason': ''
        })
    df['SE'] = se_vals
    df['Included'] = included_flags
    # Keep only rows intended for ingestion
    df_in = df[df['Included']].copy()
    if not df_in.empty:
        print(f"✅ Loaded extra group-level data: {len(df_in)} usable BG rows (from {path.name})")
    else:
        print(f"   ℹ No usable BG rows found in {path.name} (see registry for details)")
    # Return only standardized, usable columns
    # Track provenance
    df_in['SourceFile'] = (source_name if source_name else path.name)
    return df_in[['Study','Phase','Region','Metabolite','Mean','SE','n','SourceFile']]

loaded_sources = []

# Load new combined BG meta-analysis first (preferred source for overlaps)
combined_df = pd.DataFrame()
if COMBINED_BG_FILE.exists():
    combined_df = _try_load_extra_csv(COMBINED_BG_FILE, source_name=COMBINED_BG_FILE.name)
    if not combined_df.empty:
        loaded_sources.append('combined')
        try:
            inc_summary = combined_df.groupby(['Phase','Metabolite']).size().to_dict()
            print(f"   ➕ Usable BG constraints (combined file): {inc_summary}")
        except Exception:
            pass
else:
    print(f"   ℹ Optional combined file not found: {COMBINED_BG_FILE.name}")

# Load previously curated extras (legacy/all_extracted)
legacy_df = pd.DataFrame()
if EXTRA_GROUP_FILE.exists():
    legacy_df = _try_load_extra_csv(EXTRA_GROUP_FILE, source_name=EXTRA_GROUP_FILE.name)
elif (RATIO_DIR / 'bg_only_hiv_neurodata.csv').exists():
    alt_path = RATIO_DIR / 'bg_only_hiv_neurodata.csv'
    print(f"   ℹ Trying alternate file: {alt_path.name}")
    legacy_df = _try_load_extra_csv(alt_path, source_name=alt_path.name)
else:
    print(f"   ℹ Optional file not found: {EXTRA_GROUP_FILE.name}")

# Merge with de-duplication (prefer combined file; do not duplicate Dahmani)
def _canonical_study_name(s: str) -> str:
    if pd.isna(s):
        return ''
    s0 = str(s).strip().lower()
    # Normalize common variants for Dahmani to a single token to avoid duplicates
    if 'dahmani' in s0:
        return 'dahmani_2021'
    return s0

def _dedup_extras(df_list):
    merged = pd.concat([d for d in df_list if isinstance(d, pd.DataFrame) and not d.empty], ignore_index=True)
    if merged.empty:
        return merged
    keep_rows = []
    best_by_key = {}
    # Helper for priority
    def priority(row):
        study_norm = _canonical_study_name(row['Study'])
        # Prefer combined file entries
        src_priority = 1 if str(row.get('SourceFile','')).startswith('bg_combined') else 0
        # Valid SE and larger n are better; smaller SE also desirable
        se_ok = 1 if (np.isfinite(row['SE']) and row['SE'] > 0) else 0
        n_val = int(row['n']) if pd.notna(row['n']) and row['n'] == row['n'] else 0
        return (src_priority, se_ok, n_val, -float(row['SE']) if (np.isfinite(row['SE']) and row['SE']>0) else -1e9)
    for idx, row in merged.iterrows():
        key = (
            _canonical_study_name(row['Study']),
            row['Phase'], row['Region'], row['Metabolite']
        )
        pr = priority(row)
        if key not in best_by_key:
            best_by_key[key] = (pr, idx)
        else:
            if pr > best_by_key[key][0]:
                best_by_key[key] = (pr, idx)
    # Mark registry dedup reasons
    idx_keep = set(i for (_, i) in best_by_key.values())
    for i, reg in enumerate(extra_group_registry):
        # Find matching row in merged by key
        key = (
            _canonical_study_name(reg.get('Study','')),
            reg.get('Phase',''), reg.get('Region',''), reg.get('Metabolite','')
        )
        # Determine if this registry entry corresponds to a kept row
        # We don't have row indices mapping directly; use inclusion + source and key to guess
        is_kept = False
        # If it's not included initially, leave as is
        if reg.get('Included', False):
            # Check if any kept row matches this key and source
            for _, idx in best_by_key.values():
                mr = merged.iloc[idx]
                k2 = (
                    _canonical_study_name(mr['Study']), mr['Phase'], mr['Region'], mr['Metabolite']
                )
                if k2 == key and reg.get('SourceFile','') == mr.get('SourceFile',''):
                    is_kept = True
                    break
            if not is_kept:
                reg['Included'] = False
                reg['DedupReason'] = 'Dropped duplicate in favor of preferred source/SE/n'
                reg['DuplicateKey'] = '|'.join(key)
        else:
            # leave non-included entries
            pass
    kept_rows = merged.iloc[list(idx_keep)].copy()
    return kept_rows

extra_group_df = _dedup_extras([combined_df, legacy_df])
if not extra_group_df.empty:
    # Optional: exclude curated source by user request
    if args.exclude_source.strip():
        key = args.exclude_source.strip().lower()
        try:
            before = len(extra_group_df)
            extra_group_df = extra_group_df[~extra_group_df['Study'].astype(str).str.lower().str.contains(key, na=False)].copy()
            removed = before - len(extra_group_df)
            if removed > 0:
                print(f"   🔎 Excluded {removed} rows from curated constraints matching source '~{key}~'.")
        except Exception:
            pass
    try:
        inc_summary_all = extra_group_df.groupby(['Phase','Metabolite']).size().to_dict()
        print(f"   ➕ Final usable BG group constraints after de-dup: {inc_summary_all}")
    except Exception:
        pass

# ---------------------------------------------------------------------------
# Option A: augment group constraints with multi-region entries from Young 2014
# and Sailasuta 2012 for Chronic/Control (exclude BG here to avoid duplicating
# the BG entries already handled elsewhere). These will be treated as additional
# independent Normal constraints tied to phase-level means.
# ---------------------------------------------------------------------------
try:
    extra_known_rows = []
    # Helper to add a row
    def _add_known(study, phase, region, metabolite, mean, se, n, source):
        if not (np.isfinite(mean) and np.isfinite(se) and se > 0):
            return
        extra_known_rows.append({
            'Study': study,
            'Phase': _canon_phase(phase),
            'Region': _canon_region(region),
            'Metabolite': _canon_metabolite(metabolite),
            'Mean': float(mean),
            'SE': float(se),
            'n': int(n) if pd.notna(n) else np.nan,
            'SourceFile': source
        })

    # Young 2014: use all regions except BG; Chronic & Control only; NAA/Cr and Cho/Cr
    if 'young_clean' in locals() and isinstance(young_clean, pd.DataFrame) and not young_clean.empty:
        for _, r in young_clean.iterrows():
            region = str(r.get('Region',''))
            if _canon_region(region) == 'BG':
                continue  # skip BG here to avoid double counting
            phase = r.get('Phase','')
            if phase not in ['Chronic','Control']:
                continue
            n = r.get('n', np.nan)
            # NAA/Cr
            if pd.notna(r.get('NAA_ratio', np.nan)) and pd.notna(r.get('NAA_SE_ratio', np.nan)):
                _add_known('Young_2014', phase, region, 'NAA/Cr', r['NAA_ratio'], r['NAA_SE_ratio'], n, 'YOUNG_2014_CROSS_SECTIONAL_DATA.csv')
            # Cho/Cr (optional)
            if 'Cho_ratio' in r and pd.notna(r.get('Cho_ratio', np.nan)) and pd.notna(r.get('Cho_SE_ratio', np.nan)):
                _add_known('Young_2014', phase, region, 'Cho/Cr', r['Cho_ratio'], r['Cho_SE_ratio'], n, 'YOUNG_2014_CROSS_SECTIONAL_DATA.csv')

    # Sailasuta 2012: include non-BG regions (e.g., OGM) for Chronic/Control
    if 'sailasuta_clean' in locals() and isinstance(sailasuta_clean, pd.DataFrame) and not sailasuta_clean.empty:
        for _, r in sailasuta_clean.iterrows():
            region = str(r.get('Region',''))
            if _canon_region(region) == 'BG':
                continue
            phase = r.get('Phase','')
            if phase not in ['Chronic','Control']:
                continue
            n = r.get('n', np.nan)
            if pd.notna(r.get('NAA_ratio', np.nan)) and pd.notna(r.get('NAA_SE_ratio', np.nan)):
                _add_known('Sailasuta_2012', phase, region, 'NAA/Cr', r['NAA_ratio'], r['NAA_SE_ratio'], n, 'Sailasuta_2012.csv')
            if 'Cho_ratio' in r and pd.notna(r.get('Cho_ratio', np.nan)) and pd.notna(r.get('Cho_SE_ratio', np.nan)):
                _add_known('Sailasuta_2012', phase, region, 'Cho/Cr', r['Cho_ratio'], r['Cho_SE_ratio'], n, 'Sailasuta_2012.csv')

    if len(extra_known_rows) > 0:
        extra_known_df = pd.DataFrame(extra_known_rows)
        # Concat and rely on same de-dup policy against combined/legacy if overlaps exist
        extra_group_df = pd.concat([extra_group_df, extra_known_df], ignore_index=True) if not extra_group_df.empty else extra_known_df
        # Recompute counts
        extra_group_counts = {
            'NAA/Cr': {'Acute': 0, 'Chronic': 0, 'Control': 0},
            'Cho/Cr': {'Acute': 0, 'Chronic': 0, 'Control': 0}
        }
        for _, r in extra_group_df.iterrows():
            m = str(r['Metabolite']).strip()
            p = str(r['Phase']).strip().title()
            if m in extra_group_counts and p in extra_group_counts[m]:
                extra_group_counts[m][p] += 1
        # Save inclusion table
        try:
            base_results_dir.mkdir(exist_ok=True)
            (base_results_dir / 'runs').mkdir(exist_ok=True)
            run_dir.mkdir(exist_ok=True)
        except Exception:
            pass
        incl_path = run_dir / f"group_inputs_included{output_tag}.csv"
        extra_group_df[['Study','Phase','Region','Metabolite','Mean','SE','n','SourceFile']].to_csv(str(incl_path), index=False)
        # Print quick summary
        try:
            inc_summary_all = extra_group_df.groupby(['Phase','Metabolite']).size().to_dict()
            print(f"   ➕ Final usable group constraints (multi-region): {inc_summary_all}")
            # One-line summary: counts and sum of group n by Phase×Metabolite
            try:
                # Ensure numeric n
                df_n = extra_group_df.copy()
                df_n['n'] = pd.to_numeric(df_n['n'], errors='coerce')
                sum_n = df_n.groupby(['Phase','Metabolite'])['n'].sum(min_count=1).to_dict()
                def _cnt(ph, met):
                    return inc_summary_all.get((ph, met), 0)
                def _sum(ph, met):
                    v = sum_n.get((ph, met), np.nan)
                    return int(v) if np.isfinite(v) else 0
                line = (
                    "   ➕ Included group constraints — "
                    f"NAA/Cr: Acute={_cnt('Acute','NAA/Cr')} (sum n={_sum('Acute','NAA/Cr')}), "
                    f"Chronic={_cnt('Chronic','NAA/Cr')} (sum n={_sum('Chronic','NAA/Cr')}), "
                    f"Control={_cnt('Control','NAA/Cr')} (sum n={_sum('Control','NAA/Cr')}); "
                    f"Cho/Cr: Acute={_cnt('Acute','Cho/Cr')} (sum n={_sum('Acute','Cho/Cr')}), "
                    f"Chronic={_cnt('Chronic','Cho/Cr')} (sum n={_sum('Chronic','Cho/Cr')}), "
                    f"Control={_cnt('Control','Cho/Cr')} (sum n={_sum('Control','Cho/Cr')})"
                )
                print(line)
            except Exception:
                pass
        except Exception:
            pass
except Exception as e:
    print(f"   ⚠ Failed to augment multi-region group constraints: {e}")

# Precompute convenience artifacts for logging/audit and later summaries
extra_group_counts = {
    'NAA/Cr': {'Acute': 0, 'Chronic': 0, 'Control': 0},
    'Cho/Cr': {'Acute': 0, 'Chronic': 0, 'Control': 0}
}
if not extra_group_df.empty:
    try:
        for _, r in extra_group_df.iterrows():
            m = str(r['Metabolite']).strip()
            p = str(r['Phase']).strip().title()
            if m in extra_group_counts and p in extra_group_counts[m]:
                extra_group_counts[m][p] += 1
    except Exception:
        pass

# Apply region exclusion and scope filter (BG-only vs all), then save/print inclusion table
try:
    if not extra_group_df.empty:
        # Exclude region if requested
        if args.exclude_region.strip():
            def _canon_for_cmp(r: str) -> str:
                ru = str(r).upper().strip()
                # Map common aliases
                if ru in ["CG", "CINGULATE", "CINGULATE GYRUS", "CINGULATE CORTEX", "ACC", "AC"]:
                    return "AC"
                return ru
            excl = _canon_for_cmp(args.exclude_region)
            extra_group_df = extra_group_df[extra_group_df['Region'].apply(lambda x: _canon_for_cmp(x) != excl)].copy()
        # Filter by requested regions scope
        if args.regions == 'BG':
            extra_group_df = extra_group_df[extra_group_df['Region'] == 'BG'].copy()
        inclusion_tbl = extra_group_df.copy()
        inclusion_tbl = inclusion_tbl[['Study','Phase','Region','Metabolite','Mean','SE','n','SourceFile']]
        # Print a compact table preview (first 8 rows)
        preview = inclusion_tbl.head(8).to_string(index=False)
        print("\n🧾 Group constraints included (" + ("BG-only" if args.regions=='BG' else "multi-region") + ", ratios):\n" + preview + ("\n   …" if len(inclusion_tbl) > 8 else ""))
        # Ensure run_dir exists for early saves as well
        try:
            base_results_dir.mkdir(exist_ok=True)
            (base_results_dir / 'runs').mkdir(exist_ok=True)
            run_dir.mkdir(exist_ok=True)
        except Exception:
            pass
        incl_path = run_dir / f"group_inputs_included{output_tag}.csv"
        inclusion_tbl.to_csv(str(incl_path), index=False)
except Exception:
    pass

# ============================================================================
# PREPARE DATA FOR BASAL GANGLIA (PRIMARY FOCUS)
# ============================================================================

print("\n" + "=" * 80)
print("PREPARING BASAL GANGLIA DATA (working on NAA/Cr and Cho/Cr ratios)")
print("=" * 80)

# ------------------------------
# Valcour: Individual acute patients (BASAL GANGLIA ONLY)
# ------------------------------
if not valcour_df.empty:
    print(f"\n📋 Valcour columns available: {len(valcour_df.columns)} total")

    # Specifically look for Basal Ganglia (BG) columns
    # Valcour format: BGNAA, BGCho, BGCr, etc.
    naa_col = None
    cho_col = None

    # Priority order: BG > other regions
    bg_naa_variants = ['BGNAA', 'BG_NAA', 'NAA_BG']
    bg_cho_variants = ['BGCho', 'BG_Cho', 'Cho_BG']
    bg_cr_variants = ['BGCr', 'BG_Cr', 'Cr_BG']

    # First try exact BG matches
    for col in valcour_df.columns:
        if col in bg_naa_variants or col.upper() in bg_naa_variants:
            naa_col = col
        if col in bg_cho_variants or col.upper() in bg_cho_variants:
            cho_col = col
        if col in bg_cr_variants or str(col).upper() in [s.upper() for s in bg_cr_variants]:
            cr_col = col

    # If not found, search for BG prefix
    if not ('cr_col' in locals() and cr_col) or not naa_col or not cho_col:
        for col in valcour_df.columns:
            col_upper = str(col).upper()
            if 'BG' in col_upper and 'NAA' in col_upper and not naa_col:
                naa_col = col
            if 'BG' in col_upper and 'CHO' in col_upper and not cho_col:
                cho_col = col
            if 'BG' in col_upper and ('CR' in col_upper or 'CREATINE' in col_upper) and not ('cr_col' in locals() and cr_col):
                cr_col = col

    if naa_col and cho_col:
        print(f"   ✅ Found Basal Ganglia columns: NAA='{naa_col}', Cho='{cho_col}'" + (f", Cr='{cr_col}'" if 'cr_col' in locals() and cr_col else " (Cr not found; will use reference)"))

        # Use weeks 0, 4, 12, and 24 — all considered acute timeframe
        # Filter for HIV+ patients only (HIV column = 1)
        acute_weeks = ['week_0', 'week_4', 'week_12', 'week_24']
        week_filter = valcour_df['Week'].isin(acute_weeks)
        hiv_filter = valcour_df['HIV'] == 1  # HIV+ patients only
        naa_filter = valcour_df[naa_col].notna()
        cho_filter = valcour_df[cho_col].notna()

        valcour_bg_acute = valcour_df[week_filter & hiv_filter & naa_filter & cho_filter].copy()

        # Log per-week counts to show inclusion beyond week_0
        try:
            counts_per_week = valcour_bg_acute['Week'].value_counts().to_dict()
            print(f"   Included Valcour weeks (HIV+ with BG NAA & Cho present): {counts_per_week}")
        except Exception:
            pass

        # Compute ratios per individual; prefer measured Cr, fallback to reference BG Cr
        if 'cr_col' in locals() and cr_col in valcour_bg_acute.columns and valcour_bg_acute[cr_col].notna().any():
            cr_vals = valcour_bg_acute[cr_col].values
            naa_ratio_vals = valcour_bg_acute[naa_col].values / cr_vals
            cho_ratio_vals = valcour_bg_acute[cho_col].values / cr_vals
        else:
            cr_ref = CR_REFERENCE.get('BG', 8.0)
            naa_ratio_vals = valcour_bg_acute[naa_col].values / cr_ref
            cho_ratio_vals = valcour_bg_acute[cho_col].values / cr_ref

        naa_ratio_obs_acute_valcour = naa_ratio_vals
        cho_ratio_obs_acute_valcour = cho_ratio_vals
        n_acute_valcour = len(naa_ratio_obs_acute_valcour)

        print(f"\n📊 VALCOUR ACUTE BG RATIOS (HIV+ weeks 0/4/12/24, n={n_acute_valcour}):")
        print(f"   NAA/Cr: {naa_ratio_obs_acute_valcour.mean():.3f} ± {naa_ratio_obs_acute_valcour.std():.3f}")
        print(f"   Cho/Cr: {cho_ratio_obs_acute_valcour.mean():.3f} ± {cho_ratio_obs_acute_valcour.std():.3f}")

        # ------------------------------------------------------------------
        # Build Valcour week-level aggregates to study degradation over 24w
        # NOTE: This is for auxiliary analysis only. The PRIMARY model does
        # NOT use VL or these aggregates; it remains VL-free by design.
        # ------------------------------------------------------------------
        # Helper: deterministic VL selection and cleaning per modality
        def _select_vl_column(df_modality: pd.DataFrame, modality: str):
            """Return (col_name, is_log10) for requested modality 'plasma' or 'csf'."""
            cols = [str(c) for c in df_modality.columns]
            if modality == 'plasma':
                if 'logpVL' in cols: return 'logpVL', True
                if 'pVL' in cols:    return 'pVL', False
                # legacy names fallback (none expected in full sheet)
                if 'VL_log10' in cols: return 'VL_log10', True
                if 'VL' in cols:       return 'VL', False
            if modality == 'csf':
                if 'logcVL' in cols: return 'logcVL', True
                if 'cVL' in cols:    return 'cVL', False
            return None, None

        def _weekly_agg_for_arm(df_arm: pd.DataFrame, modality: str):
            """Compute week-level aggregates for NAA/Cr and chosen VL arm.
            Returns dict with arrays and audit rows.
            """
            use_col, is_log10 = _select_vl_column(df_arm, modality)
            if use_col is None:
                print(f"   ℹ No {modality.upper()} VL column detected; this arm will run time-only if requested.")
            # Map week labels to numeric weeks
            week_map = {'week_0': 0, 'week_4': 4, 'week_12': 12, 'week_24': 24}
            week_labels, week_nums, week_mean_naa, week_se_naa, week_n = [], [], [], [], []
            week_mean_vl_log10 = []
            audit_rows = []
            for wl, wn in week_map.items():
                wk = df_arm[df_arm['Week'] == wl]
                if len(wk) == 0:
                    continue
                # NAA ratio for this week
                if 'cr_col' in locals() and cr_col in wk.columns and wk[cr_col].notna().any():
                    cr_vals_wk = wk[cr_col].replace(0, np.nan).values
                    naa_r = wk[naa_col].values / cr_vals_wk
                else:
                    naa_r = wk[naa_col].values / CR_REFERENCE.get('BG', 8.0)
                naa_r = naa_r[~pd.isna(naa_r)]
                if len(naa_r) == 0:
                    continue
                n_i = len(naa_r)
                mu_i = float(np.mean(naa_r))
                se_i = float(np.std(naa_r, ddof=1) / np.sqrt(max(1, n_i))) if n_i > 1 else 0.05

                # VL handling for this week
                vl_mu = np.nan
                n_raw = n_clean = n_clipped = 0
                if use_col is not None and use_col in wk.columns:
                    vl_vals_raw = pd.to_numeric(wk[use_col], errors='coerce').replace([np.inf, -np.inf], np.nan)
                    vl_vals_raw = vl_vals_raw.dropna().values
                    n_raw = vl_vals_raw.size
                    if n_raw > 0:
                        if is_log10:
                            vl_vals = np.asarray(vl_vals_raw, dtype=float)
                        else:
                            with np.errstate(divide='ignore'):
                                vl_vals = np.log10(np.clip(vl_vals_raw.astype(float), 1.0, None))
                        # clip to plausible range [1, 8]
                        before = vl_vals.copy()
                        vl_vals = np.clip(vl_vals, 1.0, 8.0)
                        n_clipped = int(np.sum(before != vl_vals))
                        # drop non-finite just in case
                        vl_vals = vl_vals[np.isfinite(vl_vals)]
                        n_clean = vl_vals.size
                        if n_clean > 0:
                            vl_mu = float(np.mean(vl_vals))
                week_labels.append(wl)
                week_nums.append(wn)
                week_mean_naa.append(mu_i)
                week_se_naa.append(se_i)
                week_n.append(n_i)
                week_mean_vl_log10.append(vl_mu)
                audit_rows.append({
                    'week_label': wl,
                    'vl_column': use_col if use_col else '',
                    'modality': modality,
                    'n_raw': int(n_raw),
                    'n_clean': int(n_clean),
                    'n_clipped': int(n_clipped),
                    'mean_log10': float(vl_mu) if np.isfinite(vl_mu) else np.nan
                })

            return {
                'labels': week_labels,
                'nums': np.array(week_nums, dtype=float),
                'naa_mean': np.array(week_mean_naa, dtype=float),
                'naa_se': np.array(week_se_naa, dtype=float),
                'n': np.array(week_n, dtype=int),
                'vl_mean_log10': np.array(week_mean_vl_log10, dtype=float),
                'audit': audit_rows,
                'chosen_col': use_col,
                'is_log10': bool(is_log10) if use_col is not None else False,
                'modality': modality
            }

        # Legacy heuristic aggregation removed. We will compute per-arm aggregates below.

        # Legacy single-set weekly arrays for backward-compatibility (used by time/time+vl alias)
        _legacy = _weekly_agg_for_arm(valcour_bg_acute, modality='plasma')
        valcour_week_labels = _legacy['labels']
        valcour_week_nums = _legacy['nums']
        valcour_week_mean_naa = _legacy['naa_mean']
        valcour_week_se_naa = _legacy['naa_se']
        valcour_week_n = _legacy['n']
        valcour_week_mean_vl_log10 = _legacy['vl_mean_log10']

        # Save a small CSV snapshot for validation utilities (legacy plasma default)
        if len(valcour_week_nums) > 0:
            try:
                week_df_tmp = pd.DataFrame({
                    'week_label': valcour_week_labels,
                    'week_num': valcour_week_nums,
                    'n': valcour_week_n,
                    'NAA_over_Cr_mean': valcour_week_mean_naa,
                    'SE': valcour_week_se_naa,
                    'VL_log10_mean': valcour_week_mean_vl_log10
                })
                snap_path = run_dir / f"valcour_bg_week_summary{output_tag}.csv"
                base_results_dir.mkdir(exist_ok=True)
                (base_results_dir / "runs").mkdir(exist_ok=True)
                run_dir.mkdir(exist_ok=True)
                week_df_tmp.to_csv(str(snap_path), index=False)
            except Exception:
                pass
        else:
            print("   ⚠ No week-level aggregates could be computed for Valcour.")
    else:
        print("\n⚠ Could not identify Basal Ganglia NAA/Cho columns in Valcour data")
        print(
            f"   Available region columns: {[c for c in valcour_df.columns if 'NAA' in c.upper() or 'CHO' in c.upper()]}")
        naa_ratio_obs_acute_valcour = np.array([])
        cho_ratio_obs_acute_valcour = np.array([])
        n_acute_valcour = 0
        # Empty placeholders for degradation analysis
        valcour_week_labels = []
        valcour_week_nums = np.array([])
        valcour_week_mean_naa = np.array([])
        valcour_week_se_naa = np.array([])
        valcour_week_n = np.array([])
        valcour_week_mean_vl_log10 = np.array([])
else:
    naa_ratio_obs_acute_valcour = np.array([])
    cho_ratio_obs_acute_valcour = np.array([])
    n_acute_valcour = 0
    print("\n⚠ No Valcour acute data available")
    # Empty placeholders for degradation analysis
    valcour_week_labels = []
    valcour_week_nums = np.array([])
    valcour_week_mean_naa = np.array([])
    valcour_week_se_naa = np.array([])
    valcour_week_n = np.array([])
    valcour_week_mean_vl_log10 = np.array([])

# ----------------------------------------------------------------------------
# Validation helpers (Valcour manuscript replication & Fig.2 plots)
# ----------------------------------------------------------------------------

def _detect_id_column(df: pd.DataFrame):
    for c in ['SubjectID', 'ID', 'ParticipantID', 'Participant', 'PID', 'StudyID']:
        if c in df.columns:
            return c
    return None

def _find_region_columns(df: pd.DataFrame, region: str, metabolite: str):
    # Build candidate lists
    r = region.upper()
    m = metabolite.upper()
    # Handle PGM/PCG naming
    region_aliases = [r]
    if r == 'PGM':
        region_aliases += ['PCG', 'PCC', 'PCG']
    if r == 'FGM':
        region_aliases += ['FG']
    candidates = []
    for col in df.columns:
        cu = str(col).upper()
        if any(ra in cu for ra in region_aliases) and m in cu:
            candidates.append(col)
    # Cr columns per region
    cr_candidates = []
    for col in df.columns:
        cu = str(col).upper()
        if any(ra in cu for ra in region_aliases) and ('CR' in cu or 'CREAT' in cu):
            cr_candidates.append(col)
    naa_col = candidates[0] if candidates else None
    cr_col = cr_candidates[0] if cr_candidates else None
    return naa_col, cr_col

def _compute_ratio_series(df: pd.DataFrame, value_col: str, cr_col: str, region: str):
    if value_col is None or value_col not in df.columns:
        return pd.Series([np.nan]*len(df))
    if cr_col and cr_col in df.columns and df[cr_col].notna().any():
        cr_vals = pd.to_numeric(df[cr_col], errors='coerce').replace(0, np.nan)
        val = pd.to_numeric(df[value_col], errors='coerce')
        return val / cr_vals
    else:
        ref = CR_REFERENCE.get(region, 7.5)
        val = pd.to_numeric(df[value_col], errors='coerce')
        return val / ref

def _week_group_stats(s: pd.Series):
    s = pd.to_numeric(s, errors='coerce').dropna()
    n = int(s.shape[0])
    if n == 0:
        return n, np.nan, np.nan
    mu = float(s.mean())
    se = float(s.std(ddof=1) / np.sqrt(n)) if n > 1 else np.nan
    return n, mu, se

def _collect_controls(region: str, metabolite: str):
    # Prefer Young Control, then Sailasuta Control, then Chang Control
    means = []
    ses = []
    ns = []
    # Young
    try:
        y = young_clean[(young_clean['Region'] == region) & (young_clean['Phase'] == 'Control')]
        if not y.empty and f"{'NAA' if metabolite=='NAA' else 'Cho'}_ratio" in y.columns:
            means.append(float(y[f"{'NAA' if metabolite=='NAA' else 'Cho'}_ratio"].values[0]))
            ses.append(float(y[f"{'NAA' if metabolite=='NAA' else 'Cho'}_SE_ratio"].values[0]))
            ns.append(int(y['n'].values[0]))
    except Exception:
        pass
    # Sailasuta
    try:
        s = sailasuta_clean[(sailasuta_clean['Region'] == region) & (sailasuta_clean['Phase'] == 'Control')]
        if not s.empty and f"{'NAA' if metabolite=='NAA' else 'Cho'}_ratio" in s.columns:
            means.append(float(s[f"{'NAA' if metabolite=='NAA' else 'Cho'}_ratio"].values[0]))
            se_col = f"{'NAA' if metabolite=='NAA' else 'Cho'}_SE_ratio"
            if se_col in s.columns:
                ses.append(float(s[se_col].values[0]))
                ns.append(int(s['n'].values[0]))
    except Exception:
        pass
    # Chang (absolute -> ratio)
    try:
        if not chang_df.empty:
            r = region
            meta = 'NAA' if metabolite=='NAA' else 'Cho'
            dfc = chang_df[(chang_df['Phase']=='Control') & (chang_df['Region']==r) & (chang_df['Metabolite']==meta)]
            if not dfc.empty:
                ref = CR_REFERENCE.get(region, 7.5)
                means.append(float(dfc['Mean'].values[0]/ref))
                if 'SE' in dfc.columns and not pd.isna(dfc['SE'].values[0]):
                    ses.append(float(dfc['SE'].values[0]/ref))
                    ns.append(int(dfc['n'].values[0]) if 'n' in dfc.columns else 25)
    except Exception:
        pass
    # Combine by inverse-variance weighting when SE available; else simple mean
    if len(means)==0:
        return np.nan, np.nan, 0
    if len(ses)>0 and all([se>0 for se in ses]):
        w = np.array([1.0/(se**2) for se in ses])
        mu = float(np.sum(w*np.array(means))/np.sum(w))
        se = float(np.sqrt(1.0/np.sum(w)))
        n_eff = int(np.round(np.sum(ns))) if ns else 0
        return mu, se, n_eff
    else:
        mu = float(np.mean(means))
        return mu, np.nan, sum(ns) if ns else 0

def run_valcour_validation(valcour_all: pd.DataFrame):
    if valcour_all.empty:
        print("⚠ Validation skipped: no Valcour data loaded.")
        return
    print("\n════════════════════════════════════════════════════════════════════════════════")
    print("VALCOUR MANUSCRIPT VALIDATION (weeks 0/4/12/24; ratios)")
    print("════════════════════════════════════════════════════════════════════════════════")
    id_col = _detect_id_column(valcour_all)
    if id_col:
        print(f"   Detected subject ID column: {id_col}")
    else:
        print("   ⚠ No subject ID column detected; will use unpaired tests.")
    acute_weeks = ['week_0','week_4','week_12','week_24']
    hiv_mask = (valcour_all['HIV']==1) if 'HIV' in valcour_all.columns else pd.Series([True]*len(valcour_all))
    df = valcour_all[valcour_all['Week'].isin(acute_weeks) & hiv_mask].copy()
    results_rows = []
    # Panels: (a) FGM NAA, (b) FWM NAA, (c) PGM/PCG NAA, (d) BG Cho
    panels = [
        ('FGM','NAA','a'),
        ('FWM','NAA','b'),
        ('PGM','NAA','c'),
        ('BG','Cho','d')
    ]
    fig2_data_rows = []
    for region, metabolite, panel in panels:
        val_col, cr_col_r = _find_region_columns(df, region, metabolite)
        if val_col is None:
            print(f"   ⚠ {region} {metabolite} column not found; skipping panel {panel}.")
            continue
        df[f'{region}_{metabolite}_ratio'] = _compute_ratio_series(df, val_col, cr_col_r, region)
        # Week stats
        week_stats = {}
        for wl in acute_weeks:
            sub = df[df['Week']==wl][f'{region}_{metabolite}_ratio']
            n, mu, se = _week_group_stats(sub)
            week_stats[wl] = (n, mu, se)
            fig2_data_rows.append({'panel':panel,'region':region,'metabolite':metabolite,'week':wl,'n':n,'mean':mu,'se':se})
        # Controls
        ctrl_mu, ctrl_se, ctrl_n = _collect_controls(region, metabolite)
        fig2_data_rows.append({'panel':panel,'region':region,'metabolite':metabolite,'week':'Control','n':ctrl_n,'mean':ctrl_mu,'se':ctrl_se})
        # 0 vs 24
        w0 = df[df['Week']=='week_0'][f'{region}_{metabolite}_ratio']
        w24 = df[df['Week']=='week_24'][f'{region}_{metabolite}_ratio']
        p_w0_w24 = np.nan
        try:
            if id_col and id_col in df.columns:
                merged = df[df['Week'].isin(['week_0','week_24'])][[id_col,'Week',f'{region}_{metabolite}_ratio']].dropna()
                piv = merged.pivot_table(index=id_col, columns='Week', values=f'{region}_{metabolite}_ratio')
                paired = piv.dropna()
                if paired.shape[0] > 2:
                    t, p_w0_w24 = stats.ttest_rel(paired['week_0'], paired['week_24'])
                else:
                    t, p_w0_w24 = stats.ttest_ind(w0.dropna(), w24.dropna(), equal_var=False)
            else:
                t, p_w0_w24 = stats.ttest_ind(w0.dropna(), w24.dropna(), equal_var=False)
        except Exception:
            p_w0_w24 = np.nan
        # 24 vs Control (Welch using SEs if available)
        n24, mu24, se24 = week_stats['week_24'] if 'week_24' in week_stats else (0,np.nan,np.nan)
        p_24_ctrl = np.nan
        try:
            if np.isfinite(mu24) and np.isfinite(ctrl_mu):
                if np.isfinite(se24) and np.isfinite(ctrl_se) and se24>0 and ctrl_se>0:
                    se_diff = np.sqrt(se24**2 + ctrl_se**2)
                    z = (mu24 - ctrl_mu) / se_diff
                    p_24_ctrl = 2*(1 - stats.norm.cdf(abs(z)))
                else:
                    # fallback to sample Welch test if we have raw samples; otherwise NaN
                    p_24_ctrl = np.nan
        except Exception:
            p_24_ctrl = np.nan
        results_rows.append({
            'panel': panel,
            'region': region,
            'metabolite': metabolite,
            'p_week0_vs_week24': p_w0_w24,
            'week24_mean': mu24,
            'control_mean': ctrl_mu,
            'p_week24_vs_control': p_24_ctrl,
            'n_week0': week_stats.get('week_0',(0,np.nan,np.nan))[0],
            'n_week24': n24,
            'n_control': ctrl_n
        })
        if args.week_pvals:
            print(f"   [{panel}] {region} {metabolite}: p(0 vs 24) = {p_w0_w24:.3g}; p(24 vs Ctrl) = {p_24_ctrl:.3g}")
    # Save validation CSVs and Fig.2 data
    try:
        res_df = pd.DataFrame(results_rows)
        res_path = run_dir / f"valcour_replication_stats{output_tag}.csv"
        res_df.to_csv(str(res_path), index=False)
        print(f"✅ Saved Valcour replication stats: {res_path}")
        fig2_df = pd.DataFrame(fig2_data_rows)
        fig2_data_path = run_dir / f"valcour_fig2_boxplot_data{output_tag}.csv"
        fig2_df.to_csv(str(fig2_data_path), index=False)
    except Exception as e:
        print(f"⚠ Could not save validation CSVs: {e}")
    # Create Fig.2-like plots if requested
    if args.fig2 or args.validate_valcour:
        try:
            fig, axes = plt.subplots(2,2, figsize=(9,6))
            panel_map = {('FGM','NAA'):(0,0),'FWM_NAA':(0,1),'PGM_NAA':(1,0),'BG_Cho':(1,1)}
            def panel_ax(region, metabolite):
                key = (region, metabolite)
                if key==( 'FGM','NAA'):
                    return axes[0,0]
                if key==( 'FWM','NAA'):
                    return axes[0,1]
                if key==( 'PGM','NAA'):
                    return axes[1,0]
                return axes[1,1]
            # Build per-panel box data from fig2_df (synthesize control samples when only mean/SE available)
            for region, metabolite, title in [('FGM','NAA','a) FGM NAA'),('FWM','NAA','b) FWM NAA'),('PGM','NAA','c) PGM/PCG NAA'),('BG','Cho','d) BG Choline')]:
                ax = panel_ax(region, metabolite)
                sub = fig2_df[(fig2_df['region']==region) & (fig2_df['metabolite']==metabolite)]
                groups = ['week_0','week_12','week_24','Control']
                box_data = []
                for g in groups:
                    if g=='Control' and np.isfinite(sub[sub['week']=='Control']['mean']).any():
                        row = sub[sub['week']=='Control'].iloc[0]
                        n = int(row['n']) if pd.notna(row['n']) and int(row['n'])>1 else 20
                        mu = float(row['mean'])
                        se = float(row['se']) if pd.notna(row['se']) and row['se']>0 else 0.05
                        np.random.seed(123)
                        box_data.append(np.random.normal(mu, se*np.sqrt(n), n))
                    else:
                        vals = df[df['Week']==g][f'{region}_{metabolite}_ratio'] if g!='Control' else pd.Series(dtype=float)
                        box_data.append(pd.to_numeric(vals, errors='coerce').dropna().values)
                ax.boxplot(box_data, labels=['W0','W12','W24','Co'])
                ax.set_title(title)
                ax.set_ylabel(f"{metabolite}/Cr" if metabolite!='Cho' else "Cho/Cr")
            plt.tight_layout()
            fig2_path = figures_dir / f"fig2_valcour_like{output_tag}.png"
            plt.savefig(str(fig2_path), dpi=220)
            plt.close()
            print(f"✅ Saved Fig.2-style plot: {fig2_path}")
        except Exception as e:
            print(f"⚠ Could not generate Fig.2-style plots: {e}")

# If user requested validation pass, run it now and exit (model remains unchanged)
if args.validate_valcour:
    run_valcour_validation(valcour_df if 'valcour_df' in locals() else pd.DataFrame())
    # Save minimal run metadata for this validation-only run
    try:
        run_info = {
            "timestamp": run_timestamp,
            "run_name": run_name,
            "tag": args.tag,
            "validate_only": True,
            "data_counts": {
                "valcour_rows": int(len(valcour_df)) if 'valcour_df' in locals() else 0
            }
        }
        base_results_dir.mkdir(exist_ok=True)
        runs_dir.mkdir(exist_ok=True)
        run_dir.mkdir(exist_ok=True)
        with open(run_dir/"run_info.json","w") as f:
            json.dump(run_info, f, indent=2)
        print(f"🧾 Saved run metadata: {run_dir/ 'run_info.json'}")
    except Exception:
        pass
    print("\n✅ Validation pass complete. Exiting before Bayesian inference as requested (--validate-valcour).\n")
    sys.exit(0)

# ------------------------------
# Young: Group means for BG
# ------------------------------
if not young_clean.empty:
    young_bg = young_clean[young_clean['Region'] == 'BG'].copy()

    young_bg_acute = young_bg[young_bg['Phase'] == 'Acute']
    young_bg_chronic = young_bg[young_bg['Phase'] == 'Chronic']
    young_bg_control = young_bg[young_bg['Phase'] == 'Control']

    print(f"\n📊 YOUNG 2014 BG RATIOS:")
    if not young_bg_acute.empty:
        print(f"   Acute: n={young_bg_acute['n'].values[0]}, NAA/Cr={young_bg_acute['NAA_ratio'].values[0]:.3f}")
    if not young_bg_chronic.empty:
        print(f"   Chronic: n={young_bg_chronic['n'].values[0]}, NAA/Cr={young_bg_chronic['NAA_ratio'].values[0]:.3f}")
    if not young_bg_control.empty:
        print(f"   Control: n={young_bg_control['n'].values[0]}, NAA/Cr={young_bg_control['NAA_ratio'].values[0]:.3f}")
else:
    young_bg_acute = pd.DataFrame()
    young_bg_chronic = pd.DataFrame()
    young_bg_control = pd.DataFrame()

# ------------------------------
# Sailasuta: Group means for BG
# ------------------------------
if not sailasuta_clean.empty:
    sailasuta_bg = sailasuta_clean[sailasuta_clean['Region'] == 'BG'].copy()

    sailasuta_bg_acute = sailasuta_bg[sailasuta_bg['Phase'] == 'Acute']
    sailasuta_bg_chronic = sailasuta_bg[sailasuta_bg['Phase'] == 'Chronic']
    sailasuta_bg_control = sailasuta_bg[sailasuta_bg['Phase'] == 'Control']

    print(f"\n📊 SAILASUTA 2012 BG RATIOS:")
    if not sailasuta_bg_acute.empty:
        print(f"   Acute: n={sailasuta_bg_acute['n'].values[0]}, NAA/Cr={sailasuta_bg_acute['NAA_ratio'].values[0]:.3f}")
    if not sailasuta_bg_chronic.empty:
        print(f"   Chronic: n={sailasuta_bg_chronic['n'].values[0]}, NAA/Cr={sailasuta_bg_chronic['NAA_ratio'].values[0]:.3f}")
    if not sailasuta_bg_control.empty:
        print(f"   Control: n={sailasuta_bg_control['n'].values[0]}, NAA/Cr={sailasuta_bg_control['NAA_ratio'].values[0]:.3f}")
else:
    sailasuta_bg_acute = pd.DataFrame()
    sailasuta_bg_chronic = pd.DataFrame()
    sailasuta_bg_control = pd.DataFrame()

# ------------------------------
# Combine data for modeling
# ------------------------------

# ACUTE: Combine Valcour individual ratios + Young/Sailasuta group means (ratios)
# Optionally exclude Valcour via CLI flag
if args.exclude_valcour:
    acute_naa_ratio_list = []
    acute_cho_ratio_list = []
    n_acute_valcour_effective = 0
else:
    acute_naa_ratio_list = [naa_ratio_obs_acute_valcour] if n_acute_valcour > 0 else []
    acute_cho_ratio_list = [cho_ratio_obs_acute_valcour] if n_acute_valcour > 0 else []
    n_acute_valcour_effective = n_acute_valcour

# Add Young acute as weighted observations (repeat by sample size) unless disabled
if not args.no_acute_pseudo and not young_bg_acute.empty:
    n_young = int(young_bg_acute['n'].values[0])
    naa_young = float(young_bg_acute['NAA_ratio'].values[0])
    se_young = float(young_bg_acute['NAA_SE_ratio'].values[0])
    # Create pseudo-observations with appropriate spread
    np.random.seed(42)
    acute_naa_ratio_list.append(np.random.normal(naa_young, se_young * np.sqrt(n_young), n_young))

    if 'Cho_ratio' in young_bg_acute.columns and not pd.isna(young_bg_acute['Cho_ratio'].values[0]):
        cho_young = float(young_bg_acute['Cho_ratio'].values[0])
        cho_se_young = float(young_bg_acute['Cho_SE_ratio'].values[0])
        acute_cho_ratio_list.append(np.random.normal(cho_young, cho_se_young * np.sqrt(n_young), n_young))

# Add Sailasuta acute unless disabled
if not args.no_acute_pseudo and not sailasuta_bg_acute.empty:
    n_sail = int(sailasuta_bg_acute['n'].values[0])
    naa_sail = float(sailasuta_bg_acute['NAA_ratio'].values[0])
    se_sail = float(sailasuta_bg_acute['NAA_SE_ratio'].values[0])
    np.random.seed(43)
    acute_naa_ratio_list.append(np.random.normal(naa_sail, se_sail * np.sqrt(n_sail), n_sail))

    if 'Cho_ratio' in sailasuta_bg_acute.columns and not pd.isna(sailasuta_bg_acute['Cho_ratio'].values[0]):
        cho_sail = float(sailasuta_bg_acute['Cho_ratio'].values[0])
        cho_se_sail = float(sailasuta_bg_acute['Cho_SE_ratio'].values[0])
        acute_cho_ratio_list.append(np.random.normal(cho_sail, cho_se_sail * np.sqrt(n_sail), n_sail))

naa_ratio_obs_acute = np.concatenate(acute_naa_ratio_list) if acute_naa_ratio_list else np.array([])
cho_ratio_obs_acute = np.concatenate(acute_cho_ratio_list) if acute_cho_ratio_list else np.array([])

# CHRONIC: Combine group means
chronic_naa_ratio_list = []
chronic_cho_ratio_list = []

if not young_bg_chronic.empty:
    chronic_naa_ratio_list.append(young_bg_chronic['NAA_ratio'].values)
    if 'Cho_ratio' in young_bg_chronic.columns and not pd.isna(young_bg_chronic['Cho_ratio'].values[0]):
        chronic_cho_ratio_list.append(young_bg_chronic['Cho_ratio'].values)

if not sailasuta_bg_chronic.empty:
    chronic_naa_ratio_list.append(sailasuta_bg_chronic['NAA_ratio'].values)
    if 'Cho_ratio' in sailasuta_bg_chronic.columns and not pd.isna(sailasuta_bg_chronic['Cho_ratio'].values[0]):
        chronic_cho_ratio_list.append(sailasuta_bg_chronic['Cho_ratio'].values)

naa_ratio_obs_chronic = np.concatenate(chronic_naa_ratio_list) if chronic_naa_ratio_list else np.array([1.10])
cho_ratio_obs_chronic = np.concatenate(chronic_cho_ratio_list) if chronic_cho_ratio_list else np.array([0.27])

# CONTROL: Use literature reference (Chang absolute) converted to ratios
if not chang_df.empty:
    chang_bg_control = chang_df[(chang_df['Phase'] == 'Control') & (chang_df['Region'] == 'BG')]
    if not chang_bg_control.empty:
        cr_ref_bg = CR_REFERENCE.get('BG', 8.0)
        chang_naa_abs = chang_bg_control[chang_bg_control['Metabolite'] == 'NAA']['Mean'].values
        chang_cho_abs = chang_bg_control[chang_bg_control['Metabolite'] == 'Cho']['Mean'].values
        naa_ratio_obs_control = chang_naa_abs / cr_ref_bg if len(chang_naa_abs) > 0 else np.array([9.55 / cr_ref_bg])
        cho_ratio_obs_control = chang_cho_abs / cr_ref_bg if len(chang_cho_abs) > 0 else np.array([2.18 / cr_ref_bg])
    else:
        cr_ref_bg = CR_REFERENCE.get('BG', 8.0)
        naa_ratio_obs_control = np.array([9.55 / cr_ref_bg])
        cho_ratio_obs_control = np.array([2.18 / cr_ref_bg])
else:
    cr_ref_bg = CR_REFERENCE.get('BG', 8.0)
    naa_ratio_obs_control = np.array([9.55 / cr_ref_bg])
    cho_ratio_obs_control = np.array([2.18 / cr_ref_bg])

print("\n" + "=" * 80)
print("FINAL COMBINED DATA FOR MODEL (all on ratio scale)")
print("=" * 80)

if len(naa_ratio_obs_acute) > 0:
    print(f"\n✅ ACUTE: n={len(naa_ratio_obs_acute)}")
    print(f"   NAA/Cr: {naa_ratio_obs_acute.mean():.3f} ± {naa_ratio_obs_acute.std():.3f}")
    if len(cho_ratio_obs_acute) > 0:
        print(f"   Cho/Cr: {cho_ratio_obs_acute.mean():.3f} ± {cho_ratio_obs_acute.std():.3f}")
else:
    print(f"\n⚠ ACUTE: n=0 (no data available)")

print(f"\n✅ CHRONIC: n={len(naa_ratio_obs_chronic)}")
print(f"   NAA/Cr: {naa_ratio_obs_chronic.mean():.3f} ± {naa_ratio_obs_chronic.std():.3f}")
if len(cho_ratio_obs_chronic) > 0:
    print(f"   Cho/Cr: {cho_ratio_obs_chronic.mean():.3f} ± {cho_ratio_obs_chronic.std():.3f}")

print(f"\n✅ CONTROL: n={len(naa_ratio_obs_control)}")
print(f"   NAA/Cr: {naa_ratio_obs_control.mean():.3f}")
print(f"   Cho/Cr: {cho_ratio_obs_control.mean():.3f}")

# Abort if no acute individual data unless acute group-mean mode is enabled
if len(naa_ratio_obs_acute) == 0:
    if args.no_acute_pseudo:
        print("\n⚠ ACUTE: No individual acute observations. Proceeding with acute group-mean constraints only (no_acute_pseudo).")
    else:
        print("\n❌ ERROR: No acute data available for modeling!")
        print("   Consider enabling --no-acute-pseudo to use acute group-mean constraints, or include Valcour individuals.")
        sys.exit(1)

# ============================================================================
# MODEL PARAMETERS
# ============================================================================

print("\n" + "=" * 80)
print("MODEL PARAMETERS")
print("=" * 80)

# Microtubule parameters
L_MT = 2000e-9  # Microtubule length (m)
k_B = 1.38e-23  # Boltzmann constant
T = 310  # Temperature (K)
hbar = 1.055e-34  # Reduced Planck constant

# Enzyme kinetics parameters
V_max_base = 100.0  # Base enzyme velocity (nmol/min/mg)
K_m = 50.0  # Michaelis constant (μM)
S_0 = 100.0  # Substrate concentration (μM)

print("\n🔬 Physical Constants:")
print(f"   Microtubule length: {L_MT * 1e9:.0f} nm")
print(f"   Temperature: {T} K")
print(f"   Enzyme V_max baseline: {V_max_base} nmol/min/mg")

# ============================================================================
# BAYESIAN MODEL v3.6
# ============================================================================

print("\n" + "=" * 80)
print("BUILDING BAYESIAN MODEL v3.6 (ratio scale)")
print("=" * 80)

with pm.Model() as model:
    # ========================================================================
    # PRIORS: Noise correlation lengths
    # ========================================================================

    # Use LogNormal priors for positivity and better geometry
    ξ_acute    = pm.LogNormal('ξ_acute',    mu=np.log(0.6), sigma=0.15)
    ξ_chronic  = pm.LogNormal('ξ_chronic',  mu=np.log(0.8), sigma=0.15)
    ξ_control  = pm.LogNormal('ξ_control',  mu=np.log(0.5), sigma=0.15)

    # ========================================================================
    # PROTECTION MECHANISM
    # ========================================================================

    # Unconstrained Normal prior (ξ positivity handled above)
    β_ξ = pm.Normal('β_ξ', mu=-2.0, sigma=0.6)


    def quantum_protection_factor(ξ, L_MT, β_ξ):
        """Quantum coherence protection: Γ_eff = exp(β_ξ * ξ/L_MT)"""
        return pm.math.exp(β_ξ * (ξ * 1e-9) / L_MT)


    Γ_acute = quantum_protection_factor(ξ_acute, L_MT, β_ξ)
    Γ_chronic = quantum_protection_factor(ξ_chronic, L_MT, β_ξ)
    Γ_control = quantum_protection_factor(ξ_control, L_MT, β_ξ)

    # ========================================================================
    # METABOLIC MODEL
    # ========================================================================

    # Baseline NAA/Cr (control): LogNormal centered near literature control
    r_NAA_ratio_baseline = pm.LogNormal('r_NAA_ratio_baseline', mu=np.log(1.10), sigma=0.10)

    # Phase-specific modulation (positive multipliers)
    α_acute = pm.LogNormal('α_acute',   mu=np.log(1.00), sigma=0.15)
    α_chronic = pm.LogNormal('α_chronic', mu=np.log(0.95), sigma=0.12)

    # Predicted NAA/Cr ratios
    NAA_ratio_control_mean = r_NAA_ratio_baseline
    NAA_ratio_acute_mean = r_NAA_ratio_baseline * α_acute * Γ_acute
    NAA_ratio_chronic_mean = r_NAA_ratio_baseline * α_chronic * Γ_chronic

    # Store for posterior analysis
    NAA_control_pred = pm.Deterministic('NAA_ratio_control_mean', NAA_ratio_control_mean)
    NAA_acute_pred = pm.Deterministic('NAA_ratio_acute_mean', NAA_ratio_acute_mean)
    NAA_chronic_pred = pm.Deterministic('NAA_ratio_chronic_mean', NAA_ratio_chronic_mean)

    # ========================================================================
    # CHOLINE MODEL
    # ========================================================================

    Cho_ratio_baseline = pm.TruncatedNormal('Cho_ratio_baseline', mu=0.27, sigma=0.05,
                                            lower=0.1, upper=0.5)

    # Inflammation increases Cho
    Cho_ratio_control_mean = Cho_ratio_baseline
    Cho_ratio_acute_mean = Cho_ratio_baseline * 1.15  # Acute inflammation
    Cho_ratio_chronic_mean = Cho_ratio_baseline * 1.10  # Chronic inflammation

    # ========================================================================
    # LIKELIHOOD
    # ========================================================================

    # Individual-level acute observations
    σ_naa_acute = pm.HalfNormal('σ_naa_acute', sigma=0.20)
    σ_cho_acute = pm.HalfNormal('σ_cho_acute', sigma=0.1)

    if len(naa_ratio_obs_acute) > 0:
        # Robust Student-t likelihood for acute individual data
        ν_naa = pm.Exponential('ν_naa', lam=1/10)
        NAA_acute_obs = pm.StudentT('NAA_ratio_acute_obs',
                                    nu=ν_naa + 2.0,  # ensure df > 2
                                    mu=NAA_ratio_acute_mean,
                                    sigma=σ_naa_acute,
                                    observed=naa_ratio_obs_acute)

    if len(cho_ratio_obs_acute) > 0:
        Cho_acute_obs = pm.Normal('Cho_ratio_acute_obs',
                                  mu=Cho_ratio_acute_mean,
                                  sigma=σ_cho_acute,
                                  observed=cho_ratio_obs_acute)

    # Group-level chronic observations
    if len(naa_ratio_obs_chronic) > 0:
        NAA_chronic_obs = pm.Normal('NAA_ratio_chronic_obs',
                                    mu=NAA_ratio_chronic_mean,
                                    sigma=0.15,
                                    observed=naa_ratio_obs_chronic)

    if len(cho_ratio_obs_chronic) > 0:
        Cho_chronic_obs = pm.Normal('Cho_ratio_chronic_obs',
                                    mu=Cho_ratio_chronic_mean,
                                    sigma=0.1,
                                    observed=cho_ratio_obs_chronic)

    # Control observations
    NAA_control_obs = pm.Normal('NAA_ratio_control_obs',
                                mu=NAA_ratio_control_mean,
                                sigma=0.1,
                                observed=naa_ratio_obs_control)

    Cho_control_obs = pm.Normal('Cho_ratio_control_obs',
                                mu=Cho_ratio_control_mean,
                                sigma=0.08,
                                observed=cho_ratio_obs_control)

    # ------------------------------------------------------------------------
    # Additional group-mean constraints from curated CSV (BG only, ratios)
    # These are independent Normal terms on reported group means with SEs.
    # Does not create pseudo-observations and helps stabilize chronic/control.
    # ------------------------------------------------------------------------
    if 'extra_group_df' in locals() and isinstance(extra_group_df, pd.DataFrame) and not extra_group_df.empty:
        try:
            extra_any = extra_group_df[(extra_group_df['SE'].notna()) & (extra_group_df['Mean'].notna())]
            for _, row in extra_any.iterrows():
                study = str(row['Study']).strip().replace(' ', '_')
                phase = str(row['Phase']).strip().title()
                metabolite = str(row['Metabolite']).strip()
                mean_val = float(row['Mean'])
                se_val = float(row['SE']) if float(row['SE']) > 0 else None
                if not se_val:
                    continue
                if metabolite.upper() in ['NAA/CR', 'NAA_CR', 'NAA']:
                    if phase == 'Acute':
                        mu_node = NAA_ratio_acute_mean
                    elif phase == 'Chronic':
                        mu_node = NAA_ratio_chronic_mean
                    else:
                        mu_node = NAA_ratio_control_mean
                    pm.Normal(f"NAA_{study}_{phase}_mean_obs", mu=mu_node, sigma=se_val, observed=mean_val)
                elif metabolite.upper() in ['CHO/CR', 'CHO_CR', 'CHO', 'TCHO/CR', 'TCHO']:
                    if phase == 'Acute':
                        mu_node = Cho_ratio_acute_mean
                    elif phase == 'Chronic':
                        mu_node = Cho_ratio_chronic_mean
                    else:
                        mu_node = Cho_ratio_control_mean
                    pm.Normal(f"Cho_{study}_{phase}_mean_obs", mu=mu_node, sigma=se_val, observed=mean_val)
        except Exception as _e:
            print(f"   ⚠ Skipped adding some extra group means due to parsing error: {_e}")

    # If acute pseudo-obs are disabled, add BG acute group-mean constraints directly for Young/Sailasuta
    if args.no_acute_pseudo:
        try:
            if not young_bg_acute.empty:
                _mu = float(young_bg_acute['NAA_ratio'].values[0])
                _se = float(young_bg_acute['NAA_SE_ratio'].values[0])
                pm.Normal("NAA_Young_BG_Acute_mean_obs", mu=NAA_ratio_acute_mean, sigma=_se, observed=_mu)
                if 'Cho_ratio' in young_bg_acute.columns and not pd.isna(young_bg_acute['Cho_ratio'].values[0]):
                    _mu_c = float(young_bg_acute['Cho_ratio'].values[0])
                    _se_c = float(young_bg_acute['Cho_SE_ratio'].values[0])
                    pm.Normal("Cho_Young_BG_Acute_mean_obs", mu=Cho_ratio_acute_mean, sigma=_se_c, observed=_mu_c)
            if not sailasuta_bg_acute.empty:
                _mu = float(sailasuta_bg_acute['NAA_ratio'].values[0])
                _se = float(sailasuta_bg_acute['NAA_SE_ratio'].values[0])
                pm.Normal("NAA_Sailasuta_BG_Acute_mean_obs", mu=NAA_ratio_acute_mean, sigma=_se, observed=_mu)
                if 'Cho_ratio' in sailasuta_bg_acute.columns and not pd.isna(sailasuta_bg_acute['Cho_ratio'].values[0]):
                    _mu_c = float(sailasuta_bg_acute['Cho_ratio'].values[0])
                    _se_c = float(sailasuta_bg_acute['Cho_SE_ratio'].values[0])
                    pm.Normal("Cho_Sailasuta_BG_Acute_mean_obs", mu=Cho_ratio_acute_mean, sigma=_se_c, observed=_mu_c)
        except Exception as _e:
            print(f"   ⚠ Failed to add acute group-mean constraints (BG): {_e}")

    # ========================================================================
    # NOTE: Valcour 0–24w auxiliary analysis has been moved OUTSIDE the
    # primary model and fully decoupled. The primary BG model is VL-free.
    # ========================================================================

    # ========================================================================
    # HYPOTHESIS TEST
    # ========================================================================

    Δξ = pm.Deterministic('Δξ', ξ_chronic - ξ_acute)

print("\n✅ Model built successfully!")
print(f"   Free parameters: {len(model.free_RVs)}")
# Extend observations summary with group-mean constraint counts
try:
    gm_naa = extra_group_counts.get('NAA/Cr', {}) if 'extra_group_counts' in locals() else {}
    gm_cho = extra_group_counts.get('Cho/Cr', {}) if 'extra_group_counts' in locals() else {}
    obs_line = (
        f"   Observations: {len(naa_ratio_obs_acute)} acute + {len(naa_ratio_obs_chronic)} chronic + {len(naa_ratio_obs_control)} control"
    )
    if gm_naa or gm_cho:
        obs_line += (
            " | Group means — NAA/Cr: "
            f"Acute={gm_naa.get('Acute',0)}, Chronic={gm_naa.get('Chronic',0)}, Control={gm_naa.get('Control',0)}; "
            "Cho/Cr: "
            f"Acute={gm_cho.get('Acute',0)}, Chronic={gm_cho.get('Chronic',0)}, Control={gm_cho.get('Control',0)}"
        )
    print(obs_line)
except Exception:
    print(f"   Observations: {len(naa_ratio_obs_acute)} acute + {len(naa_ratio_obs_chronic)} chronic + {len(naa_ratio_obs_control)} control")

# ============================================================================
# INFERENCE
# ============================================================================

print("\n" + "=" * 80)
print("RUNNING MCMC INFERENCE")
print("=" * 80)

with model:
    print("\n⏳ Sampling (this may take 5-10 minutes)...")
    trace = pm.sample(
        draws=2500,
        tune=4500,
        chains=4,
        cores=4,
        target_accept=0.995,
        return_inferencedata=True,
        idata_kwargs={"log_likelihood": True},
        random_seed=42
    )

    print("\n✅ Sampling complete!")

    # Quick check: ensure log_likelihood is present for WAIC/LOO
    try:
        has_llk = hasattr(trace, 'log_likelihood') and trace.log_likelihood is not None
        llk_vars = list(trace.log_likelihood.data_vars) if has_llk else []
        print(f"\n🔎 WAIC/LOO readiness: log_likelihood present = {has_llk}; observed terms = {len(llk_vars)}")
    except Exception:
        pass

    print("\n⏳ Generating posterior predictive samples...")
    ppc = pm.sample_posterior_predictive(trace, random_seed=42)
    print("✅ Posterior predictive complete!")

# ============================================================================
# RESULTS
# ============================================================================

print("\n" + "=" * 80)
print("BAYESIAN INFERENCE RESULTS (ratio scale)")
print("=" * 80)

posterior = trace.posterior

# Key parameters
ξ_acute_samples = posterior['ξ_acute'].values.flatten()
ξ_chronic_samples = posterior['ξ_chronic'].values.flatten()
Δξ_samples = posterior['Δξ'].values.flatten()
β_ξ_samples = posterior['β_ξ'].values.flatten()

print("\n🔬 CORRELATION LENGTH PARAMETERS:")
print(f"\n   ξ_acute = {ξ_acute_samples.mean():.3f} ± {ξ_acute_samples.std():.3f} nm")
print(f"   95% HDI: [{np.percentile(ξ_acute_samples, 2.5):.3f}, {np.percentile(ξ_acute_samples, 97.5):.3f}]")

print(f"\n   ξ_chronic = {ξ_chronic_samples.mean():.3f} ± {ξ_chronic_samples.std():.3f} nm")
print(f"   95% HDI: [{np.percentile(ξ_chronic_samples, 2.5):.3f}, {np.percentile(ξ_chronic_samples, 97.5):.3f}]")

print(f"\n   Δξ = ξ_chronic - ξ_acute = {Δξ_samples.mean():.3f} ± {Δξ_samples.std():.3f} nm")
print(f"   95% HDI: [{np.percentile(Δξ_samples, 2.5):.3f}, {np.percentile(Δξ_samples, 97.5):.3f}]")

# Hypothesis test
P_acute_shorter = (Δξ_samples > 0).sum() / len(Δξ_samples)
print(f"\n✅ P(ξ_acute < ξ_chronic) = {P_acute_shorter:.4f}")

if P_acute_shorter > 0.999:
    print("   *** HYPOTHESIS STRONGLY SUPPORTED ***")
elif P_acute_shorter > 0.99:
    print("   *** HYPOTHESIS SUPPORTED ***")
elif P_acute_shorter > 0.95:
    print("   ** HYPOTHESIS LIKELY **")

print(f"\n🔬 PROTECTION FACTOR EXPONENT:")
print(f"   β_ξ = {β_ξ_samples.mean():.2f} ± {β_ξ_samples.std():.2f}")

# Predicted NAA/Cr
naa_acute_pred = posterior['NAA_ratio_acute_mean'].values.flatten()
naa_chronic_pred = posterior['NAA_ratio_chronic_mean'].values.flatten()
naa_control_pred = posterior['NAA_ratio_control_mean'].values.flatten()

print(f"\n📊 PREDICTED vs OBSERVED NAA/Cr:")
print(f"\n   Acute: {naa_acute_pred.mean():.3f} (pred) vs {naa_ratio_obs_acute.mean():.3f} (obs)")
print(f"   Error: {abs(naa_acute_pred.mean() - naa_ratio_obs_acute.mean()) / max(1e-6, naa_ratio_obs_acute.mean()) * 100:.1f}%")

print(f"\n   Chronic: {naa_chronic_pred.mean():.3f} (pred) vs {naa_ratio_obs_chronic.mean():.3f} (obs)")
print(f"   Error: {abs(naa_chronic_pred.mean() - naa_ratio_obs_chronic.mean()) / max(1e-6, naa_ratio_obs_chronic.mean()) * 100:.1f}%")

print(f"\n   Control: {naa_control_pred.mean():.3f} (pred) vs {naa_ratio_obs_control.mean():.3f} (obs)")
print(f"   Error: {abs(naa_control_pred.mean() - naa_ratio_obs_control.mean()) / max(1e-6, naa_ratio_obs_control.mean()) * 100:.1f}%")

# ---------------------------------------------------------------------------
# Posterior predictive by bins (acute tails diagnostic)
# ---------------------------------------------------------------------------
try:
    if 'NAA_ratio_acute_obs' in ppc and len(naa_ratio_obs_acute) > 0:
        # Align posterior predictive draws with observed order
        # ppc['NAA_ratio_acute_obs'] has shape (draws, chains, n_obs) or (draws, n_obs)
        ppc_arr = ppc['NAA_ratio_acute_obs']
        arr = ppc_arr
        # Flatten draws/chains dimensions
        arr = np.asarray(arr)
        if arr.ndim == 3:
            # (chains, draws, n) or (draws, chains, n) depending on PyMC version
            # Try to move n to last and flatten others
            if arr.shape[0] < 10 and arr.shape[1] > 10:
                # assume (chains, draws, n)
                arr = arr.reshape(arr.shape[0]*arr.shape[1], arr.shape[2])
            else:
                # assume (draws, chains, n)
                arr = arr.reshape(arr.shape[0]*arr.shape[1], arr.shape[2])
        elif arr.ndim == 2:
            # (draws, n)
            pass
        else:
            # unexpected shape — skip
            arr = None
        if arr is not None:
            obs = np.asarray(naa_ratio_obs_acute).astype(float)
            # Define quartile bins on observed acute NAA/Cr
            try:
                qs = np.quantile(obs, [0.25, 0.5, 0.75])
            except Exception:
                qs = [np.median(obs)]*3
            bins = [-np.inf, qs[0], qs[1], qs[2], np.inf]
            labels = ['Q1 (low)','Q2','Q3','Q4 (high)']
            bin_idx = np.digitize(obs, bins) - 1
            rows = []
            for b in range(4):
                idx = np.where(bin_idx == b)[0]
                if idx.size == 0:
                    continue
                obs_mean = float(np.mean(obs[idx]))
                # Predicted mean per-observation (average over draws)
                pred_means = np.mean(arr[:, idx], axis=0)
                pred_mean = float(np.mean(pred_means))
                mae = float(np.mean(np.abs(pred_means - obs[idx])))
                rows.append({
                    'bin': labels[b],
                    'n': int(idx.size),
                    'obs_mean': obs_mean,
                    'pred_mean': pred_mean,
                    'MAE': mae
                })
            if rows:
                ppc_bins_df = pd.DataFrame(rows)
                ppc_bins_path = run_dir / f"ppc_bins{output_tag}.csv"
                # Ensure directories exist (in case called before save block)
                try:
                    base_results_dir.mkdir(exist_ok=True)
                    (base_results_dir / 'runs').mkdir(exist_ok=True)
                    run_dir.mkdir(exist_ok=True)
                    figures_dir.mkdir(exist_ok=True)
                except Exception:
                    pass
                ppc_bins_df.to_csv(str(ppc_bins_path), index=False)
                # Simple bar plot
                try:
                    x = np.arange(len(ppc_bins_df))
                    width = 0.35
                    fig, ax = plt.subplots(figsize=(6,3.2))
                    ax.bar(x - width/2, ppc_bins_df['pred_mean'], width, label='Pred')
                    ax.bar(x + width/2, ppc_bins_df['obs_mean'], width, label='Obs')
                    ax.set_xticks(x)
                    ax.set_xticklabels(ppc_bins_df['bin'])
                    ax.set_ylabel('NAA/Cr (acute)')
                    ax.set_title('Posterior predictive by observed NAA/Cr quartiles')
                    ax.legend()
                    plt.tight_layout()
                    outp = figures_dir / f"ppc_bins{output_tag}.png"
                    plt.savefig(str(outp), dpi=220)
                    plt.close()
                except Exception:
                    pass
                print(f"   Saved PPC by bins: {ppc_bins_path} (figure in figures/)")
except Exception as _e:
    print(f"⚠ PPC-by-bins diagnostic failed: {_e}")

# ---------------------------------------------------------------------------
# Post-fit audit: how well the model fits the included group constraints
# ---------------------------------------------------------------------------
try:
    group_fit_rows = []
    # Posterior means and HDIs for phase-level predictions
    def _post_summary(arr):
        arr = np.asarray(arr).flatten()
        h = az.hdi(arr, hdi_prob=0.95)
        return float(np.mean(arr)), float(h[0]), float(h[1])

    naa_mu_acute, naa_hdi_lo_acute, naa_hdi_hi_acute = _post_summary(posterior['NAA_ratio_acute_mean'].values)
    naa_mu_chronic, naa_hdi_lo_chronic, naa_hdi_hi_chronic = _post_summary(posterior['NAA_ratio_chronic_mean'].values)
    naa_mu_control, naa_hdi_lo_control, naa_hdi_hi_control = _post_summary(posterior['NAA_ratio_control_mean'].values)

    # Cho summaries if present in posterior deterministics
    cho_mu_acute = cho_hdi_lo_acute = cho_hdi_hi_acute = np.nan
    cho_mu_chronic = cho_hdi_lo_chronic = cho_hdi_hi_chronic = np.nan
    cho_mu_control = cho_hdi_lo_control = cho_hdi_hi_control = np.nan
    try:
        if 'Cho_ratio_acute_mean' in posterior and 'Cho_ratio_chronic_mean' in posterior and 'Cho_ratio_control_mean' in posterior:
            cho_mu_acute, cho_hdi_lo_acute, cho_hdi_hi_acute = _post_summary(posterior['Cho_ratio_acute_mean'].values)
            cho_mu_chronic, cho_hdi_lo_chronic, cho_hdi_hi_chronic = _post_summary(posterior['Cho_ratio_chronic_mean'].values)
            cho_mu_control, cho_hdi_lo_control, cho_hdi_hi_control = _post_summary(posterior['Cho_ratio_control_mean'].values)
    except Exception:
        pass

    if 'extra_group_df' in locals() and isinstance(extra_group_df, pd.DataFrame) and not extra_group_df.empty:
        for _, row in extra_group_df.iterrows():
            study = str(row['Study'])
            phase = str(row['Phase']).title()
            metab = str(row['Metabolite'])
            mean_obs = float(row['Mean'])
            se_obs = float(row['SE']) if pd.notna(row['SE']) else np.nan
            n_obs = int(row['n']) if pd.notna(row['n']) else np.nan
            if metab == 'NAA/Cr':
                if phase == 'Acute':
                    mu, lo, hi = naa_mu_acute, naa_hdi_lo_acute, naa_hdi_hi_acute
                elif phase == 'Chronic':
                    mu, lo, hi = naa_mu_chronic, naa_hdi_lo_chronic, naa_hdi_hi_chronic
                else:
                    mu, lo, hi = naa_mu_control, naa_hdi_lo_control, naa_hdi_hi_control
            elif metab == 'Cho/Cr':
                if phase == 'Acute':
                    mu, lo, hi = cho_mu_acute, cho_hdi_lo_acute, cho_hdi_hi_acute
                elif phase == 'Chronic':
                    mu, lo, hi = cho_mu_chronic, cho_hdi_lo_chronic, cho_hdi_hi_chronic
                else:
                    mu, lo, hi = cho_mu_control, cho_hdi_lo_control, cho_hdi_hi_control
            else:
                continue
            residual = mean_obs - mu if np.isfinite(mu) else np.nan
            z = residual / se_obs if (np.isfinite(residual) and np.isfinite(se_obs) and se_obs > 0) else np.nan
            group_fit_rows.append({
                'Study': study,
                'Phase': phase,
                'Metabolite': metab,
                'Mean_obs': mean_obs,
                'SE': se_obs,
                'n': n_obs,
                'mu_phase_post_mean': mu,
                'mu_phase_HDI_low': lo,
                'mu_phase_HDI_high': hi,
                'residual': residual,
                'z_score': z
            })

    if group_fit_rows:
        fit_df = pd.DataFrame(group_fit_rows)
        fit_path = run_dir / f"group_likelihood_fit{output_tag}.csv"
        try:
            fit_df.to_csv(str(fit_path), index=False)
            print(f"\n🧪 Group constraint fit summary saved: {fit_path}")
            # Print a quick summary by phase/metabolite
            try:
                by_pm = fit_df.groupby(['Phase','Metabolite']).size().to_dict()
                print(f"   Used in likelihood (counts by Phase×Metabolite): {by_pm}")
                worst = fit_df.reindex(fit_df['z_score'].abs().sort_values(ascending=False).index).head(3)
                print("   Top |z| residuals:\n" + worst[['Study','Phase','Metabolite','Mean_obs','mu_phase_post_mean','SE','z_score']].to_string(index=False))
            except Exception:
                pass

            # ------------------------------------------------------------------
            # Visualization: overlay group means (+SE) vs posterior phase means
            # Saves: figures/group_constraints_fit<tag>.png
            # ------------------------------------------------------------------
            try:
                # Ensure figures directory exists
                try:
                    base_results_dir.mkdir(exist_ok=True)
                    runs_dir.mkdir(exist_ok=True)
                    run_dir.mkdir(exist_ok=True)
                    figures_dir.mkdir(exist_ok=True)
                except Exception:
                    pass

                # Helper to plot a single metabolite panel
                def _plot_panel(ax, metabolite_name: str,
                                 post_means: dict,
                                 post_hdi: dict,
                                 data_subset: pd.DataFrame):
                    phases = ['Acute','Chronic','Control']
                    x = np.arange(len(phases))
                    # Posterior mean and HDI per phase
                    mu = [post_means.get(p, np.nan) for p in phases]
                    lo = [post_hdi.get(p, (np.nan, np.nan))[0] for p in phases]
                    hi = [post_hdi.get(p, (np.nan, np.nan))[1] for p in phases]
                    # Plot posterior bands (HDI) and means
                    for i, p in enumerate(phases):
                        if np.isfinite(mu[i]) and np.isfinite(lo[i]) and np.isfinite(hi[i]):
                            ax.fill_between([i-0.15, i+0.15], [lo[i], lo[i]], [hi[i], hi[i]],
                                            color='C2', alpha=0.2, linewidth=0)
                            ax.errorbar([i], [mu[i]],
                                        yerr=[[mu[i]-lo[i]],[hi[i]-mu[i]]],
                                        fmt='o', color='C2', capsize=3, label='Posterior mean ±95% HDI' if i==0 else None)
                    # Overlay group means with SE (jittered)
                    if not data_subset.empty:
                        for i, p in enumerate(phases):
                            subp = data_subset[data_subset['Phase']==p]
                            if subp.empty:
                                continue
                            # Jitter positions
                            n_points = len(subp)
                            if n_points == 1:
                                xs = [x[i]]
                            else:
                                jitter = np.linspace(-0.12, 0.12, n_points)
                                xs = x[i] + jitter
                            ys = subp['Mean_obs'].values.astype(float)
                            ses = subp['SE'].values.astype(float)
                            ax.errorbar(xs, ys, yerr=ses, fmt='o', color='C0', alpha=0.8, capsize=3,
                                        label='Group mean ± SE' if i==0 else None)
                    ax.set_xticks(x)
                    ax.set_xticklabels(['Acute','Chronic','Control'])
                    ax.set_ylabel(metabolite_name)
                    ax.grid(True, axis='y', alpha=0.2)

                # Build posterior summaries dictionaries
                post_means_naa = {
                    'Acute': naa_mu_acute,
                    'Chronic': naa_mu_chronic,
                    'Control': naa_mu_control,
                }
                post_hdi_naa = {
                    'Acute': (naa_hdi_lo_acute, naa_hdi_hi_acute),
                    'Chronic': (naa_hdi_lo_chronic, naa_hdi_hi_chronic),
                    'Control': (naa_hdi_lo_control, naa_hdi_hi_control),
                }
                have_cho_post = np.isfinite(cho_mu_acute) and np.isfinite(cho_mu_chronic) and np.isfinite(cho_mu_control)
                if have_cho_post:
                    post_means_cho = {'Acute': cho_mu_acute, 'Chronic': cho_mu_chronic, 'Control': cho_mu_control}
                    post_hdi_cho = {
                        'Acute': (cho_hdi_lo_acute, cho_hdi_hi_acute),
                        'Chronic': (cho_hdi_lo_chronic, cho_hdi_hi_chronic),
                        'Control': (cho_hdi_lo_control, cho_hdi_hi_control),
                    }

                # Prepare data subsets for plotting
                naa_df = fit_df[fit_df['Metabolite'] == 'NAA/Cr'].copy()
                cho_df = fit_df[fit_df['Metabolite'] == 'Cho/Cr'].copy()

                # Decide figure layout based on availability
                if not cho_df.empty and have_cho_post:
                    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
                    _plot_panel(axes[0], 'NAA/Cr', post_means_naa, post_hdi_naa, naa_df)
                    axes[0].set_title('Group constraints vs posterior — NAA/Cr (BG)')
                    _plot_panel(axes[1], 'Cho/Cr', post_means_cho, post_hdi_cho, cho_df)
                    axes[1].set_title('Group constraints vs posterior — Cho/Cr (BG)')
                    # Build a shared legend
                    handles, labels = [], []
                    for ax in axes:
                        for h, l in zip(*ax.get_legend_handles_labels()):
                            if l not in labels:
                                handles.append(h); labels.append(l)
                    if handles:
                        fig.legend(handles, labels, loc='upper center', ncol=2, frameon=False)
                    plt.tight_layout(rect=[0, 0, 1, 0.92])
                else:
                    fig, ax = plt.subplots(1, 1, figsize=(5.5, 4))
                    _plot_panel(ax, 'NAA/Cr', post_means_naa, post_hdi_naa, naa_df)
                    ax.set_title('Group constraints vs posterior — NAA/Cr (BG)')
                    # Legend
                    handles, labels = ax.get_legend_handles_labels()
                    if handles:
                        ax.legend(handles, labels, loc='upper center', ncol=2, frameon=False)
                    plt.tight_layout()

                out_path = figures_dir / f"group_constraints_fit{output_tag}.png"
                plt.savefig(str(out_path), dpi=220)
                plt.close()
                print(f"   Saved: {out_path}")
            except Exception as _e:
                print(f"⚠ Could not generate group constraints fit figure: {_e}")
        except Exception as _e:
            print(f"⚠ Could not save group constraint fit CSV: {_e}")
except Exception as _e:
    print(f"⚠ Post-fit group constraint audit failed: {_e}")

# ============================================================================
# CONVERGENCE DIAGNOSTICS
# ============================================================================

print("\n" + "=" * 80)
print("CONVERGENCE DIAGNOSTICS")
print("=" * 80)

# R-hat (should be < 1.01)
rhat = az.rhat(trace)
# Extract max R-hat across all variables
rhat_values = []
for var in rhat.data_vars:
    rhat_var = rhat[var].values
    if np.isscalar(rhat_var):
        rhat_values.append(float(rhat_var))
    else:
        rhat_values.append(float(np.max(rhat_var)))
rhat_max = max(rhat_values)

print(f"\n   Max R-hat: {rhat_max:.4f}")
if rhat_max < 1.01:
    print("   ✅ Excellent convergence (R-hat < 1.01)")
elif rhat_max < 1.05:
    print("   ✅ Good convergence (R-hat < 1.05)")
else:
    print("   ⚠ May need more samples (R-hat > 1.05)")

# Effective sample size
ess = az.ess(trace)
# Extract min ESS across all variables
ess_values = []
for var in ess.data_vars:
    ess_var = ess[var].values
    if np.isscalar(ess_var):
        ess_values.append(float(ess_var))
    else:
        ess_values.append(float(np.min(ess_var)))
ess_min = min(ess_values)

print(f"\n   Min ESS: {ess_min:.0f}")
if ess_min > 400:
    print("   ✅ Adequate effective sample size")
else:
    print("   ⚠ Low ESS, consider more samples")

# ============================================================================
# SAVE RESULTS
# ============================================================================

print("\n" + "=" * 80)
print("SAVING RESULTS")
print("=" * 80)

# Create output directories (base and run-specific)
base_results_dir = script_dir / "results_v3_6"
base_results_dir.mkdir(exist_ok=True)

runs_dir = base_results_dir / "runs"
runs_dir.mkdir(exist_ok=True)

run_dir = runs_dir / run_name
run_dir.mkdir(exist_ok=True)
figures_dir = run_dir / "figures"
figures_dir.mkdir(exist_ok=True)

# Save summary statistics (primary: run folder; also copy to base with timestamp if enabled)
try:
    # Save registry of extra group-level inputs (if any were parsed)
    if 'extra_group_registry' in locals() and isinstance(extra_group_registry, list) and len(extra_group_registry) > 0:
        reg_df = pd.DataFrame(extra_group_registry)
        reg_path = run_dir / f'group_inputs_registry{output_tag}.csv'
        reg_df.to_csv(str(reg_path), index=False)
        print(f"✅ Saved group inputs registry: {reg_path}")
except Exception as e:
    print(f"⚠ Could not save group inputs registry: {e}")

summary_path = run_dir / f'summary{output_tag}.csv'
summary_df = az.summary(trace, hdi_prob=0.95)
summary_df.to_csv(str(summary_path))
print(f"✅ Saved summary: {summary_path}")
if use_timestamp:
    summary_base_ts = base_results_dir / f"summary{output_tag}_{run_timestamp}.csv"
    try:
        summary_df.to_csv(str(summary_base_ts))
        print(f"   ↳ Also saved (timestamped): {summary_base_ts}")
    except Exception:
        pass

# Create results summary
results_summary = {
    'Parameter': ['ξ_acute', 'ξ_chronic', 'Δξ', 'β_ξ', 'P(ξ_acute < ξ_chronic)'],
    'Mean': [
        ξ_acute_samples.mean(),
        ξ_chronic_samples.mean(),
        Δξ_samples.mean(),
        β_ξ_samples.mean(),
        P_acute_shorter
    ],
    'SD': [
        ξ_acute_samples.std(),
        ξ_chronic_samples.std(),
        Δξ_samples.std(),
        β_ξ_samples.std(),
        np.nan
    ],
    'HDI_2.5%': [
        np.percentile(ξ_acute_samples, 2.5),
        np.percentile(ξ_chronic_samples, 2.5),
        np.percentile(Δξ_samples, 2.5),
        np.percentile(β_ξ_samples, 2.5),
        np.nan
    ],
    'HDI_97.5%': [
        np.percentile(ξ_acute_samples, 97.5),
        np.percentile(ξ_chronic_samples, 97.5),
        np.percentile(Δξ_samples, 97.5),
        np.percentile(β_ξ_samples, 97.5),
        np.nan
    ]
}

results_df = pd.DataFrame(results_summary)
results_path = run_dir / f'results_v3_6_ratio_scale{output_tag}.csv'
results_df.to_csv(str(results_path), index=False)
print(f"✅ Saved results: {results_path}")
if use_timestamp:
    results_base_ts = base_results_dir / f"results_v3_6_ratio_scale{output_tag}_{run_timestamp}.csv"
    try:
        results_df.to_csv(str(results_base_ts), index=False)
        print(f"   ↳ Also saved (timestamped): {results_base_ts}")
    except Exception:
        pass

# Save posterior predictive comparison
ppc_summary = {
    'condition': ['Acute', 'Chronic', 'Control'],
    'NAA/Cr_pred': [naa_acute_pred.mean(), naa_chronic_pred.mean(), naa_control_pred.mean()],
    'NAA/Cr_obs': [naa_ratio_obs_acute.mean(), naa_ratio_obs_chronic.mean(), naa_ratio_obs_control.mean()],
    'n_obs': [len(naa_ratio_obs_acute), len(naa_ratio_obs_chronic), len(naa_ratio_obs_control)]
}
ppc_df = pd.DataFrame(ppc_summary)
ppc_path = run_dir / f'posterior_predictive_comparison_ratio{output_tag}.csv'
ppc_df.to_csv(str(ppc_path), index=False)
print(f"✅ Saved posterior predictive: {ppc_path}")
if use_timestamp:
    ppc_base_ts = base_results_dir / f"posterior_predictive_comparison_ratio{output_tag}_{run_timestamp}.csv"
    try:
        ppc_df.to_csv(str(ppc_base_ts), index=False)
        print(f"   ↳ Also saved (timestamped): {ppc_base_ts}")
    except Exception:
        pass

# ==========================================================================
# FULLY DECOUPLED VALCOUR AUXILIARY MODEL (optional; does NOT affect primary)
# ==========================================================================
if args.valcour_aux and args.valcour_aux != 'off':
    try:
        if 'valcour_bg_acute' in locals() and not valcour_bg_acute.empty:
            print("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
            print(f"VALCOUR AUXILIARY ANALYSIS (mode={args.valcour_aux}) — fully decoupled from primary model")
            print("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")

            def _run_aux_for_arm(arm_modality: str, arm_tag: str):
                agg = _weekly_agg_for_arm(valcour_bg_acute, modality=arm_modality)
                labels = agg['labels']; nums = agg['nums']; naa_mu = agg['naa_mean']; naa_se = agg['naa_se']; vl_mu = agg['vl_mean_log10']
                # Save audits
                try:
                    audit_df = pd.DataFrame(agg['audit'])
                    audit_path = run_dir / f"valcour_week_vl_summary_valcour_{arm_tag}{output_tag}.csv"
                    audit_df.to_csv(str(audit_path), index=False)
                except Exception:
                    pass
                if nums.size == 0:
                    print(f"ℹ {arm_modality.upper()} arm skipped: no week-level aggregates available.")
                    return
                # Save week summary for this arm
                try:
                    week_df = pd.DataFrame({
                        'week_label': labels,
                        'week_num': nums,
                        'n': agg['n'],
                        'NAA_over_Cr_mean': naa_mu,
                        'SE': naa_se,
                        'VL_log10_mean': vl_mu
                    })
                    week_path = run_dir / f"valcour_bg_week_summary_valcour_{arm_tag}{output_tag}.csv"
                    week_df.to_csv(str(week_path), index=False)
                except Exception:
                    pass
                # Prepare arrays for modeling
                week_scaled = (nums / 24.0).astype(float)
                se_vals = naa_se.astype(float) if naa_se.size > 0 else np.full_like(week_scaled, 0.05, dtype=float)
                include_vl = np.isfinite(vl_mu).any()
                if include_vl:
                    vl_vals = vl_mu.astype(float)
                    if not np.isfinite(vl_vals).any() or np.nanmax(vl_vals) <= 0:
                        include_vl = False
                with pm.Model() as aux_model:
                    pm_week = pm.Data(f'Aux_{arm_tag}_week_scaled', week_scaled)
                    pm_se = pm.Data(f'Aux_{arm_tag}_week_SE', se_vals)
                    if include_vl:
                        pm_vl = pm.Data(f'Aux_{arm_tag}_VL_log10', vl_vals)
                        vl_ref = float(np.nanmean(vl_vals)) if np.isfinite(vl_vals).any() else 4.0
                    # Priors
                    mu_base = pm.LogNormal(f'μ_valcour_base_{arm_tag}', mu=np.log(1.05), sigma=0.15)
                    eta_t = pm.Normal(f'η_t_per24w_aux_{arm_tag}', mu=0.0, sigma=0.2)
                    if include_vl:
                        eta_vl = pm.Normal(f'η_vl_per_log10_aux_{arm_tag}', mu=0.0, sigma=0.2)
                    # Mean structure
                    if include_vl:
                        delta = pm.math.exp(eta_t * pm_week + eta_vl * (pm_vl - vl_ref))
                    else:
                        delta = pm.math.exp(eta_t * pm_week)
                    pm.Deterministic(f'Aux_Valcour_delta_week_{arm_tag}', delta)
                    pred = mu_base * delta
                    pm.Deterministic(f'Aux_Valcour_pred_week_{arm_tag}', pred)
                    # Noise and likelihood
                    sigma_week = pm.HalfNormal(f'σ_valcour_week_aux_{arm_tag}', sigma=0.05)
                    sigma_obs = pm.math.sqrt(pm_se**2 + sigma_week**2)
                    _ = pm.Normal(f'Aux_Valcour_weekly_NAA_ratio_obs_{arm_tag}', mu=pred, sigma=sigma_obs,
                                  observed=naa_mu.astype(float))
                    aux_trace = pm.sample(draws=1000, tune=2000, target_accept=0.98, chains=4, cores=4,
                                          return_inferencedata=True, random_seed=42)
                # Save outputs
                aux_trace_path = run_dir / f"aux_valcour_trace_valcour_{arm_tag}{output_tag}.nc"
                try:
                    aux_trace.to_netcdf(str(aux_trace_path))
                    print(f"✅ Saved Valcour {arm_modality} auxiliary trace: {aux_trace_path}")
                except Exception as e:
                    print(f"⚠ Could not save Valcour {arm_modality} auxiliary trace: {e}")
                aux_summary = az.summary(aux_trace, hdi_prob=0.95)
                aux_sum_path = run_dir / f"aux_valcour_summary_valcour_{arm_tag}{output_tag}.csv"
                aux_summary.to_csv(str(aux_sum_path))
                print(f"✅ Saved Valcour {arm_modality} auxiliary summary: {aux_sum_path}")
                # Ribbon figure
                try:
                    pred_arr = aux_trace.posterior[f'Aux_Valcour_pred_week_{arm_tag}'].values
                    pred_flat = pred_arr.reshape(-1, pred_arr.shape[-1])
                    mean_line = np.mean(pred_flat, axis=0)
                    hdi = az.hdi(pred_flat, hdi_prob=0.95)
                    lo = hdi[:, 0]; hi = hdi[:, 1]
                    weeks = nums
                    fig, ax = plt.subplots(figsize=(6.5, 3.6))
                    ax.plot(weeks, mean_line, color="#4C78A8", lw=2, label="Predicted NAA/Cr (aux)")
                    ax.fill_between(weeks, lo, hi, color="#4C78A8", alpha=0.2, label="95% HDI")
                    ax.errorbar(weeks, naa_mu, yerr=se_vals, fmt='o', color="#F58518", label="Observed mean ± SE", zorder=3)
                    ax.set_xlabel("Week"); ax.set_ylabel("NAA/Cr (BG)")
                    mode_label = f"{arm_modality.upper()} VL"
                    ax.set_title(f"Valcour auxiliary ({mode_label}){(' ' + output_tag) if output_tag else ''}")
                    ax.legend(); plt.tight_layout()
                    aux_fig = figures_dir / f"aux_valcour_ribbon_valcour_{arm_tag}{output_tag}.png"
                    plt.savefig(str(aux_fig), dpi=220); plt.close()
                    print(f"✅ Saved: {aux_fig}")
                except Exception as e:
                    print(f"⚠ Could not create Valcour {arm_modality} ribbon plot: {e}")

            # Route modes
            mode = args.valcour_aux
            if mode == 'time':
                # Use legacy arrays without VL
                if 'valcour_week_nums' in locals() and isinstance(valcour_week_nums, np.ndarray) and valcour_week_nums.size > 0:
                    _run_aux_for_arm('plasma', 'plasma')  # will include VL if present; acceptable for legacy
                else:
                    print("ℹ Valcour auxiliary skipped: no week-level aggregates available.")
            elif mode in ('time+vl', 'plasma'):
                _run_aux_for_arm('plasma', 'plasma')
            elif mode == 'csf':
                _run_aux_for_arm('csf', 'csf')
            elif mode == 'both':
                _run_aux_for_arm('plasma', 'plasma')
                _run_aux_for_arm('csf', 'csf')
            else:
                print("ℹ Unknown auxiliary mode; skipping.")
        else:
            print("ℹ Valcour auxiliary skipped: no week-level aggregates available.")
    except Exception as e:
        print(f"⚠ Valcour auxiliary analysis failed: {e}")

# ---------------------------------------------------------------------------
# Save full posterior trace to NetCDF and optional CSV of samples
# ---------------------------------------------------------------------------
try:
    if not args.no_save_trace:
        netcdf_path = run_dir / f"trace{output_tag}.nc"
        az.to_netcdf(trace, str(netcdf_path))
        try:
            size_mb = netcdf_path.stat().st_size / (1024 * 1024)
            print(f"✅ Saved full posterior trace (NetCDF): {netcdf_path} ({size_mb:.1f} MB)")
        except Exception:
            print(f"✅ Saved full posterior trace (NetCDF): {netcdf_path}")

        if args.save_posterior_csv:
            print("📝 Saving flattened posterior samples CSV for key parameters…")
            key_vars = [
                'β_ξ', 'ξ_acute', 'ξ_chronic', 'ξ_control', 'Δξ',
                'NAA_ratio_acute_mean', 'NAA_ratio_chronic_mean', 'NAA_ratio_control_mean'
            ]
            rows = {}
            for var in key_vars:
                if var in trace.posterior:
                    vals = trace.posterior[var].values  # (chain, draw, ...)
                    flat = vals.reshape(-1, vals.shape[-1]) if vals.ndim == 3 else vals.reshape(-1)
                    # If last dim >1 (e.g., vector), save separate columns per index
                    if flat.ndim == 1:
                        rows[var] = flat
                    else:
                        for i in range(flat.shape[1]):
                            rows[f"{var}[{i}]"] = flat[:, i]
            poster_csv_df = pd.DataFrame(rows)
            posterior_csv_path = run_dir / f"posterior_samples{output_tag}.csv"
            poster_csv_df.to_csv(str(posterior_csv_path), index=False)
            print(f"✅ Saved posterior samples CSV: {posterior_csv_path}")
    else:
        print("ℹ Skipped saving NetCDF trace (--no-save-trace set)")
except Exception as e:
    print(f"⚠ Could not save NetCDF/CSV posterior outputs: {e}")

# ---------------------------------------------------------------------------
# Posterior density plots (KDE) and optional overlay against a comparison trace
# ---------------------------------------------------------------------------
try:
    if args.plot_densities:
        print("\n🖼 Generating posterior density (KDE) plots…")
        var_names = [
            'β_ξ', 'ξ_acute', 'ξ_chronic', 'ξ_control', 'Δξ',
            'NAA_ratio_acute_mean', 'NAA_ratio_chronic_mean', 'NAA_ratio_control_mean'
        ]

        # Multi-panel overview for available vars
        vars_available = [v for v in var_names if v in trace.posterior]
        if len(vars_available) > 0:
            az.style.use("arviz-darkgrid")
            ax = az.plot_posterior(trace, var_names=vars_available, kind='kde', point_estimate='mean', hdi_prob=0.95)
            plt.suptitle(f"Posterior (KDE) — ratio scale{(' ' + output_tag) if output_tag else ''}")
            overview_path = figures_dir / f"posterior_overview_kde{output_tag}.png"
            plt.tight_layout()
            plt.savefig(str(overview_path), dpi=200)
            plt.close()
            print(f"✅ Saved: {overview_path}")

        # Individual parameter plots
        for v in var_names:
            if v in trace.posterior:
                az.style.use("arviz-darkgrid")
                ax = az.plot_posterior(trace, var_names=[v], kind='kde', point_estimate='mean', hdi_prob=0.95)
                plt.title(f"{v} — KDE{(' ' + output_tag) if output_tag else ''}")
                outp = figures_dir / f"posterior_{v}_kde{output_tag}.png"
                plt.tight_layout()
                plt.savefig(str(outp), dpi=220)
                plt.close()
                print(f"   Saved: {outp}")

        # Overlay against comparison NetCDF if provided
        if args.compare_trace:
            try:
                compare_path = Path(args.compare_trace)
                if not compare_path.exists():
                    print(f"⚠ Comparison trace not found: {compare_path}")
                else:
                    other_idata = az.from_netcdf(str(compare_path))
                    curr_label = args.tag if args.tag else 'current'
                    comp_label = compare_path.stem
                    for v in var_names:
                        if (v in trace.posterior) and (v in other_idata.posterior):
                            az.style.use("arviz-darkgrid")
                            _ = az.plot_density([trace, other_idata], var_names=[v], data_labels=[curr_label, comp_label], shade=0.1)
                            plt.title(f"Overlay KDE: {v} — {curr_label} vs {comp_label}")
                            outp = figures_dir / f"overlay_{v}{output_tag}.png"
                            plt.tight_layout()
                            plt.savefig(str(outp), dpi=220)
                            plt.close()
                            print(f"   Saved overlay: {outp}")
            except Exception as e:
                print(f"⚠ Could not generate overlay density plots: {e}")
except Exception as e:
    print(f"⚠ Could not generate posterior density plots: {e}")

# ---------------------------------------------------------------------------
# Extended plots: PPC bars/violins, trace and forest plots, Valcour ribbon
# ---------------------------------------------------------------------------
try:
    if args.plots_extended or args.plots_extra or args.plot_densities:
        # 1) PPC bar with error annotations
        try:
            fig, ax = plt.subplots(figsize=(6, 4))
            conds = ppc_df['condition'].tolist()
            pred = ppc_df['NAA/Cr_pred'].values
            obs = ppc_df['NAA/Cr_obs'].values
            x = np.arange(len(conds))
            w = 0.35
            ax.bar(x - w/2, obs, width=w, color="#4C78A8", alpha=0.8, label="Observed")
            ax.bar(x + w/2, pred, width=w, color="#F58518", alpha=0.8, label="Predicted")
            for i, (p, o) in enumerate(zip(pred, obs)):
                if o != 0:
                    pe = abs(p - o) / abs(o) * 100.0
                    ax.text(i + w/2, max(p, o) + 0.02, f"{pe:.1f}%", ha='center', va='bottom', fontsize=9)
            ax.set_xticks(x)
            ax.set_xticklabels(conds)
            ax.set_ylabel("NAA/Cr")
            ax.set_title(f"PPC — Predicted vs Observed (ratio){(' ' + output_tag) if output_tag else ''}")
            ax.legend()
            plt.tight_layout()
            ppc_fig_path = figures_dir / f"ppc_bar_ratio{output_tag}.png"
            plt.savefig(str(ppc_fig_path), dpi=220)
            plt.close()
            print(f"✅ Saved: {ppc_fig_path}")
        except Exception as e:
            print(f"⚠ Could not create PPC bar chart: {e}")

        # 2) ArviZ trace plot for key vars
        try:
            key_vars = ['β_ξ', 'ξ_acute', 'ξ_chronic', 'ξ_control']
            avail = [v for v in key_vars if v in trace.posterior]
            if avail:
                az.style.use("arviz-darkgrid")
                ax = az.plot_trace(trace, var_names=avail, compact=True)
                plt.suptitle(f"Trace plots{(' ' + output_tag) if output_tag else ''}")
                trace_fig_path = figures_dir / f"trace_plots{output_tag}.png"
                plt.tight_layout()
                plt.savefig(str(trace_fig_path), dpi=200)
                plt.close()
                print(f"✅ Saved: {trace_fig_path}")
        except Exception as e:
            print(f"⚠ Could not create trace plots: {e}")

        # 3) Forest plot (HDIs)
        try:
            forest_vars = ['β_ξ', 'ξ_acute', 'ξ_chronic', 'ξ_control', 'Δξ']
            availf = [v for v in forest_vars if v in trace.posterior]
            if availf:
                az.style.use("arviz-darkgrid")
                ax = az.plot_forest(trace, var_names=availf, combined=True, hdi_prob=0.95)
                plt.title(f"Forest (95% HDI){(' ' + output_tag) if output_tag else ''}")
                forest_path = figures_dir / f"forest_hdi{output_tag}.png"
                plt.tight_layout()
                plt.savefig(str(forest_path), dpi=200)
                plt.close()
                print(f"✅ Saved: {forest_path}")
        except Exception as e:
            print(f"⚠ Could not create forest plot: {e}")

        # 4) Valcour degradation ribbon plot
        try:
            if 'Valcour_pred_week' in trace.posterior and 'valcour_week_nums' in locals() and isinstance(valcour_week_nums, np.ndarray) and valcour_week_nums.size > 0:
                pred_arr = trace.posterior['Valcour_pred_week'].values  # c,d,w
                pred_flat = pred_arr.reshape(-1, pred_arr.shape[-1])
                # Compute mean and HDI across samples
                mean_line = np.mean(pred_flat, axis=0)
                hdi = az.hdi(pred_flat, hdi_prob=0.95)
                lo = hdi[:, 0]
                hi = hdi[:, 1]
                weeks = valcour_week_nums
                fig, ax = plt.subplots(figsize=(6.5, 3.6))
                ax.plot(weeks, mean_line, color="#4C78A8", lw=2, label="Predicted NAA/Cr")
                ax.fill_between(weeks, lo, hi, color="#4C78A8", alpha=0.2, label="95% HDI")
                # overlay observed weekly means
                if isinstance(valcour_week_mean_naa, np.ndarray) and valcour_week_mean_naa.size == weeks.size:
                    ax.scatter(weeks, valcour_week_mean_naa, color="#F58518", label="Observed mean", zorder=3)
                ax.set_xlabel("Week")
                ax.set_ylabel("NAA/Cr (BG)")
                ax.set_title(f"Valcour 0–24w degradation (BG){(' ' + output_tag) if output_tag else ''}")
                ax.legend()
                plt.tight_layout()
                valcour_fig = figures_dir / f"valcour_degradation_ribbon{output_tag}.png"
                plt.savefig(str(valcour_fig), dpi=220)
                plt.close()
                print(f"✅ Saved: {valcour_fig}")
        except Exception as e:
            print(f"⚠ Could not create Valcour ribbon plot: {e}")
except Exception as e:
    print(f"⚠ Extended plotting failed: {e}")
# Save Valcour 0–24w degradation analysis outputs (if available)
# ---------------------------------------------------------------------------
try:
    if 'valcour_week_nums' in locals() and isinstance(valcour_week_nums, np.ndarray) and valcour_week_nums.size > 0 \
            and 'Valcour_delta_week' in trace.posterior and 'Valcour_pred_week' in trace.posterior:
        print("\n📝 Saving Valcour 0–24w degradation analysis results…")

        # Extract posterior arrays: shape (chain, draw, week)
        delta_arr = trace.posterior['Valcour_delta_week'].values  # c,d,w
        pred_arr = trace.posterior['Valcour_pred_week'].values    # c,d,w

        # Combine chain and draw
        delta_flat = delta_arr.reshape(-1, delta_arr.shape[-1])
        pred_flat = pred_arr.reshape(-1, pred_arr.shape[-1])

        # Compute means and 95% HDIs per week
        def hdi_95(x):
            hdi = az.hdi(x, hdi_prob=0.95)
            return float(hdi[0]), float(hdi[1])

        rows = []
        for i in range(delta_flat.shape[1]):
            d_mean = float(np.mean(delta_flat[:, i]))
            p_mean = float(np.mean(pred_flat[:, i]))
            d_lo, d_hi = hdi_95(delta_flat[:, i])
            p_lo, p_hi = hdi_95(pred_flat[:, i])
            rows.append({
                'week_label': valcour_week_labels[i] if i < len(valcour_week_labels) else f'week_idx_{i}',
                'week_num': float(valcour_week_nums[i]) if i < len(valcour_week_nums) else float('nan'),
                'VL_mean_log10': float(valcour_week_mean_vl_log10[i]) if i < len(valcour_week_mean_vl_log10) else float('nan'),
                'n_week': int(valcour_week_n[i]) if i < len(valcour_week_n) else 0,
                'delta_mean': d_mean,
                'delta_hdi2.5%': d_lo,
                'delta_hdi97.5%': d_hi,
                'pred_NAA/Cr_mean': p_mean,
                'pred_hdi2.5%': p_lo,
                'pred_hdi97.5%': p_hi
            })

        valcour_deg_df = pd.DataFrame(rows)
        valcour_deg_path = run_dir / f'valcour_degradation_0_24w{output_tag}.csv'
        valcour_deg_df.to_csv(str(valcour_deg_path), index=False)
        print(f"✅ Saved Valcour degradation summary: {valcour_deg_path}")
        if use_timestamp:
            try:
                valcour_deg_base_ts = base_results_dir / f"valcour_degradation_0_24w{output_tag}_{run_timestamp}.csv"
                valcour_deg_df.to_csv(str(valcour_deg_base_ts), index=False)
                print(f"   ↳ Also saved (timestamped): {valcour_deg_base_ts}")
            except Exception:
                pass
    else:
        print("ℹ Valcour degradation outputs not generated (no week aggregates or posterior deterministics missing).")
except Exception as e:
    print(f"⚠ Could not save Valcour degradation outputs: {e}")

# ---------------------------------------------------------------------------
# Save run metadata
# ---------------------------------------------------------------------------
try:
    run_info = {
        "timestamp": run_timestamp,
        "run_name": run_name,
        "tag": args.tag,
        "exclude_valcour": bool(args.exclude_valcour),
        "versions": {
            "python": sys.version.split(" (", 1)[0],
            "platform": platform.platform(),
            "pymc": getattr(pm, "__version__", "unknown"),
            "arviz": getattr(az, "__version__", "unknown"),
            "numpy": getattr(np, "__version__", "unknown"),
            "pandas": getattr(pd, "__version__", "unknown"),
        },
        "data_counts": {
            "n_acute_total": int(len(naa_ratio_obs_acute)),
            "n_valcour": int(n_acute_valcour) if 'n_acute_valcour' in locals() else 0,
            "n_young_acute": int(young_bg_acute['n'].values[0]) if 'young_bg_acute' in locals() and not young_bg_acute.empty else 0,
            "n_sailasuta_acute": int(sailasuta_bg_acute['n'].values[0]) if 'sailasuta_bg_acute' in locals() and not sailasuta_bg_acute.empty else 0,
            "n_chronic_points": int(len(naa_ratio_obs_chronic)),
            "n_control_points": int(len(naa_ratio_obs_control)),
        },
        "paths": {
            "base_results_dir": str(base_results_dir),
            "run_dir": str(run_dir),
            "figures_dir": str(figures_dir),
        },
        "cli": vars(args),
    }
    run_info_path = run_dir / "run_info.json"
    with open(run_info_path, "w") as f:
        json.dump(run_info, f, indent=2)
    print(f"🧾 Saved run metadata: {run_info_path}")
except Exception as e:
    print(f"⚠ Could not save run metadata: {e}")

print("\n" + "=" * 80)
print("✅ ANALYSIS COMPLETE USING RATIO DATA (Valcour/Chang converted, Young/Sailasuta native)")
print("=" * 80)

print(f"\n🎯 FINAL SAMPLE SIZES:")
print(f"   Total acute observations: n={len(naa_ratio_obs_acute)}")
if (not args.exclude_valcour) and n_acute_valcour > 0:
    print(f"   - Valcour individuals: {n_acute_valcour}")
if not young_bg_acute.empty:
    print(f"   - Young 2014: {int(young_bg_acute['n'].values[0])}")
if not sailasuta_bg_acute.empty:
    print(f"   - Sailasuta 2012: {int(sailasuta_bg_acute['n'].values[0])}")

print(f"\n🎯 KEY FINDING:")
print(f"   P(ξ_acute < ξ_chronic) = {P_acute_shorter:.4f}")
print(f"   Δξ = {Δξ_samples.mean():.3f} ± {Δξ_samples.std():.3f} nm")
print(f"   β_ξ = {β_ξ_samples.mean():.2f} ± {β_ξ_samples.std():.2f}")

print("\n" + "=" * 80)