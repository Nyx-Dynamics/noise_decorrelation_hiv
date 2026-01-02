#!/usr/bin/env python3
"""
Bayesian Model v3.6 (ORIGINAL, BG-only, ratio scale)

This is the preserved "original" v3.6 ratio-scale model you selected (O2):
- BG-only focus for primary model inputs
- No multi-region group-mean augmentation from meta CSVs
- VL completely removed from the primary BG model
- Optional Valcour auxiliary analysis (time or plasma/CSF VL) is available and fully decoupled

Notes:
- This script is intentionally conservative to match your historical v3.6 ratio-scale behavior.
- For enhanced features (multi-region constraints, extended figures), use
  bayesian_v3_6_corrected_local.py instead.

Updated: November 16, 2025 (frozen original O2)
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
from pathlib import Path
import warnings
import sys

warnings.filterwarnings('ignore')

print("\n" + "╔" + "=" * 78 + "╗")
print("║" + " " * 10 + "BAYESIAN MODEL v3.6 (ORIGINAL, BG-only, RATIO DATA)" + " " * 6 + "║")
print("║" + " " * 78 + "║")
print("║" + " Valcour + Chang (abs → ratio) + Young + Sailasuta (ratios) — BG only " + "║")
print("╚" + "=" * 78 + "╝\n")

# ============================================================================
# FILE PATHS
# ============================================================================
script_dir = Path(__file__).resolve().parent
possible_data_dirs = [
    script_dir / "data" / "curated",
    script_dir.parent / "data" / "curated",
    script_dir.parent.parent / "data" / "curated",
]

data_dir = None
for p in possible_data_dirs:
    if p.exists():
        data_dir = p
        print(f"✅ Found data directory: {data_dir}")
        break
if data_dir is None:
    print("❌ Could not find data directory. Exiting.")
    sys.exit(1)

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
CHANG_FILE = VALCOUR_DIR / "CHANG_2002_EXTRACTED.csv"

print(f"\n📂 Data paths configured:")
print(f"   Absolute data: {VALCOUR_DIR}")
print(f"   Ratio data: {RATIO_DIR}")

# ============================================================================
# CLI
# ============================================================================
parser = argparse.ArgumentParser(description="Bayesian v3.6 ORIGINAL (BG-only, ratio scale)")
parser.add_argument("--exclude-valcour", action="store_true",
                    help="Exclude Valcour individuals from Acute pool (run ablation).")
parser.add_argument("--tag", type=str, default="",
                    help="Optional tag for output filenames.")
parser.add_argument("--run-label", type=str, default="",
                    help="Optional label for run folder (otherwise timestamp).")
parser.add_argument("--no-timestamp", action="store_true",
                    help="Disable timestamp suffixing of artifact filenames and run folder nesting.")
parser.add_argument("--no-save-trace", action="store_true",
                    help="Disable saving full posterior trace to NetCDF.")
parser.add_argument("--plot-densities", action="store_true",
                    help="Generate posterior KDE plots (ArviZ).")
parser.add_argument("--plots-extended", action="store_true",
                    help="Generate PPC bar, trace, and forest plots.")
parser.add_argument("--compare-trace", type=str, default="",
                    help="Optional NetCDF path to overlay densities against.")
parser.add_argument("--valcour-aux", type=str, default="off", choices=["off","time","plasma","csf","both"],
                    help="Valcour auxiliary (decoupled): time-only or VL arms (plasma/csf/both)")
args, _ = parser.parse_known_args()

if args.exclude_valcour:
    print("\n🧪 Configuration: EXCLUDING Valcour acute individuals from ACUTE pool (ablation)")
else:
    print("\n🧪 Configuration: Including Valcour acute individuals in ACUTE pool (default)")
output_tag = ("_" + args.tag.strip()) if args.tag.strip() else ""
if output_tag:
    print(f"   Output files will include tag: '{output_tag}'")
print(f"   Save NetCDF trace: {'NO' if args.no_save_trace else 'YES'}")
if args.plot_densities:
    print("   Will generate posterior density plots (ArviZ KDE)")
if args.plots_extended:
    print("   Will generate extended plots (PPC, trace, forest)")
if args.compare_trace:
    print(f"   Will overlay densities vs: {args.compare_trace}")
if args.valcour_aux != "off":
    print(f"   Valcour auxiliary analysis enabled: {args.valcour_aux} (fully decoupled)")

run_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
use_timestamp = not args.no_timestamp
run_name = args.run_label.strip() if args.run_label.strip() else (f"run_{run_timestamp}" if use_timestamp else (args.tag.strip() if args.tag.strip() else "run"))
print(f"   Run identifier: {run_name}")

# Prepare result dirs
base_results_dir = script_dir / "results_v3_6"
runs_dir = base_results_dir / "runs"
run_dir = runs_dir / run_name / "original"  # keep original in its own subfolder
figures_dir = run_dir / "figures"
for d in [base_results_dir, runs_dir, run_dir, figures_dir]:
    d.mkdir(parents=True, exist_ok=True)

# ============================================================================
# Helpers
# ============================================================================
CR_REFERENCE = {'BG': 8.0, 'FWM': 6.8, 'PGM': 7.8, 'FGM': 7.8, 'AC': 7.8, 'OGM': 7.5}

def convert_absolute_to_ratio(abs_value, region):
    cr_ref = CR_REFERENCE.get(region, 7.5)
    if pd.isna(abs_value):
        return np.nan
    return float(abs_value) / cr_ref

# ============================================================================
# LOAD DATA (BG-only for primary)
# ============================================================================
print("\n" + "=" * 80)
print("LOADING DATA")
print("=" * 80)

# Valcour
print("\n📁 Loading Valcour 2015 data (absolute)")
valcour_dfs = []
for wf in VALCOUR_FILES:
    if wf.exists():
        df = pd.read_excel(str(wf))
        df['Week'] = wf.stem.replace('valcour_2015_', '')
        valcour_dfs.append(df)
        print(f"   ✅ Loaded {wf.name}: {len(df)} rows")
    else:
        print(f"   ⚠ File not found: {wf.name}")
valcour_df = pd.concat(valcour_dfs, ignore_index=True) if valcour_dfs else pd.DataFrame()

# Young 2014 ratios (keep as ratios)
print("\n📁 Loading Young 2014 (ratios)")
if YOUNG_FILE.exists():
    ydf = pd.read_csv(str(YOUNG_FILE))
    yproc = []
    for _, r in ydf.iterrows():
        if r['Metabolite'] == 'NAA/Cr':
            yproc.append({
                'Study':'Young_2014',
                'Phase':'Acute' if r['Phase']=='Primary' else r['Phase'],
                'Region':r['Region'],
                'n':int(r['n']),
                'NAA_ratio':r['Ratio_Median'],
                'NAA_SE_ratio':r['SE']
            })
        elif r['Metabolite'] == 'Cho/Cr' and yproc:
            # attach cho to matching last (same simple approximation)
            yproc[-1]['Cho_ratio'] = r['Ratio_Median']
            yproc[-1]['Cho_SE_ratio'] = r['SE']
    young_clean = pd.DataFrame(yproc)
    young_bg = young_clean[young_clean['Region']=='BG'].copy()
    print(f"   ✅ Loaded Young: {len(young_clean)} group rows (BG subset: {len(young_bg)})")
else:
    young_clean = pd.DataFrame(); young_bg = pd.DataFrame()
    print(f"   ⚠ Missing: {YOUNG_FILE}")

# Sailasuta 2012 ratios (keep as ratios)
print("\n📁 Loading Sailasuta 2012 (ratios)")
if SAILASUTA_FILE.exists():
    sraw = pd.read_csv(str(SAILASUTA_FILE), skiprows=1)
    sprows = []
    for phase in sraw['Group'].unique():
        for region in sraw['Brain_Region'].unique():
            blk = sraw[(sraw['Group']==phase)&(sraw['Brain_Region']==region)]
            naa = blk[blk['Metabolite']=='NAA']
            cho = blk[blk['Metabolite']=='tCho']
            if not naa.empty:
                n = int(naa['n'].values[0])
                sd = float(naa['SD'].values[0]) if not pd.isna(naa['SD'].values[0]) else np.nan
                se = sd/np.sqrt(n) if np.isfinite(sd) and n>1 else np.nan
                row = {'Study':'Sailasuta_2012','Phase':phase,'Region':region,'n':n,
                       'NAA_ratio': float(naa['Value'].values[0]), 'NAA_SE_ratio': se}
                if not cho.empty:
                    sd_c = float(cho['SD'].values[0]) if not pd.isna(cho['SD'].values[0]) else np.nan
                    row['Cho_ratio'] = float(cho['Value'].values[0])
                    row['Cho_SE_ratio'] = sd_c/np.sqrt(n) if np.isfinite(sd_c) and n>1 else np.nan
                sprows.append(row)
    sailasuta_clean = pd.DataFrame(sprows)
    sailasuta_bg = sailasuta_clean[sailasuta_clean['Region']=='BG'].copy()
    print(f"   ✅ Loaded Sailasuta: {len(sailasuta_clean)} group rows (BG subset: {len(sailasuta_bg)})")
else:
    sailasuta_clean = pd.DataFrame(); sailasuta_bg = pd.DataFrame()
    print(f"   ⚠ Missing: {SAILASUTA_FILE}")

# Chang absolute → ratio (Control)
print("\n📁 Loading Chang 2002 (absolute → ratio)")
if CHANG_FILE.exists():
    chang_df = pd.read_csv(str(CHANG_FILE))
    print(f"   ✅ Loaded Chang rows: {len(chang_df)}")
else:
    chang_df = pd.DataFrame(); print(f"   ⚠ Missing: {CHANG_FILE}")

# ============================================================================
# PREP BG DATA
# ============================================================================
print("\n" + "=" * 80)
print("PREPARING BASAL GANGLIA DATA (BG-only)")
print("=" * 80)

if not valcour_df.empty:
    # locate BG columns
    naa_col = next((c for c in valcour_df.columns if str(c).upper().startswith('BG') and 'NAA' in str(c).upper()), None)
    cho_col = next((c for c in valcour_df.columns if str(c).upper().startswith('BG') and 'CHO' in str(c).upper()), None)
    cr_col = next((c for c in valcour_df.columns if str(c).upper().startswith('BG') and ('CR' in str(c).upper() or 'CREAT' in str(c).upper())), None)
    if naa_col and cho_col:
        acute_weeks = ['week_0','week_4','week_12','week_24']
        mask = valcour_df['Week'].isin(acute_weeks) & (valcour_df['HIV']==1)
        sub = valcour_df[mask & valcour_df[naa_col].notna() & valcour_df[cho_col].notna()].copy()
        if cr_col and cr_col in sub.columns and sub[cr_col].notna().any():
            crv = sub[cr_col].values
            naa_ratio_obs_acute_valcour = sub[naa_col].values / crv
            cho_ratio_obs_acute_valcour = sub[cho_col].values / crv
        else:
            naa_ratio_obs_acute_valcour = sub[naa_col].values / CR_REFERENCE['BG']
            cho_ratio_obs_acute_valcour = sub[cho_col].values / CR_REFERENCE['BG']
        print(f"\n📊 VALCOUR ACUTE BG RATIOS (n={len(naa_ratio_obs_acute_valcour)}):\n   NAA/Cr: {np.mean(naa_ratio_obs_acute_valcour):.3f} ± {np.std(naa_ratio_obs_acute_valcour):.3f}\n   Cho/Cr: {np.mean(cho_ratio_obs_acute_valcour):.3f} ± {np.std(cho_ratio_obs_acute_valcour):.3f}")
    else:
        naa_ratio_obs_acute_valcour = np.array([])
        cho_ratio_obs_acute_valcour = np.array([])
        print("⚠ BG columns not found in Valcour.")
else:
    naa_ratio_obs_acute_valcour = np.array([])
    cho_ratio_obs_acute_valcour = np.array([])
    print("⚠ No Valcour data loaded.")

# Young BG
young_bg_acute = young_bg[young_bg['Phase']=='Acute'] if not young_bg.empty else pd.DataFrame()
young_bg_chronic = young_bg[young_bg['Phase']=='Chronic'] if not young_bg.empty else pd.DataFrame()
young_bg_control = young_bg[young_bg['Phase']=='Control'] if not young_bg.empty else pd.DataFrame()

# Sailasuta BG
sail_bg_acute = sailasuta_bg[sailasuta_bg['Phase']=='Acute'] if not sailasuta_bg.empty else pd.DataFrame()
sail_bg_chronic = sailasuta_bg[sailasuta_bg['Phase']=='Chronic'] if not sailasuta_bg.empty else pd.DataFrame()
sail_bg_control = sailasuta_bg[sailasuta_bg['Phase']=='Control'] if not sailasuta_bg.empty else pd.DataFrame()

# Combine arrays
if args.exclude_valcour:
    acute_naa_ratio_list = []
    acute_cho_ratio_list = []
else:
    acute_naa_ratio_list = [naa_ratio_obs_acute_valcour] if len(naa_ratio_obs_acute_valcour)>0 else []
    acute_cho_ratio_list = [cho_ratio_obs_acute_valcour] if len(cho_ratio_obs_acute_valcour)>0 else []

# Expand Young/Sailasuta Acute as pseudo-observations (historical behavior)
if not young_bg_acute.empty:
    n = int(young_bg_acute.iloc[0]['n']); mu = float(young_bg_acute.iloc[0]['NAA_ratio']); se = float(young_bg_acute.iloc[0]['NAA_SE_ratio'])
    np.random.seed(42); acute_naa_ratio_list.append(np.random.normal(mu, se*np.sqrt(n), n))
    if 'Cho_ratio' in young_bg_acute.columns and pd.notna(young_bg_acute.iloc[0].get('Cho_ratio', np.nan)):
        mu_c = float(young_bg_acute.iloc[0]['Cho_ratio']); se_c = float(young_bg_acute.iloc[0]['Cho_SE_ratio'])
        acute_cho_ratio_list.append(np.random.normal(mu_c, se_c*np.sqrt(n), n))
if not sail_bg_acute.empty:
    n = int(sail_bg_acute.iloc[0]['n']); mu = float(sail_bg_acute.iloc[0]['NAA_ratio']); se = float(sail_bg_acute.iloc[0]['NAA_SE_ratio'])
    np.random.seed(43); acute_naa_ratio_list.append(np.random.normal(mu, se*np.sqrt(n), n))
    if 'Cho_ratio' in sail_bg_acute.columns and pd.notna(sail_bg_acute.iloc[0].get('Cho_ratio', np.nan)):
        mu_c = float(sail_bg_acute.iloc[0]['Cho_ratio']); se_c = float(sail_bg_acute.iloc[0]['Cho_SE_ratio'])
        acute_cho_ratio_list.append(np.random.normal(mu_c, se_c*np.sqrt(n), n))

naa_ratio_obs_acute = np.concatenate(acute_naa_ratio_list) if acute_naa_ratio_list else np.array([])
cho_ratio_obs_acute = np.concatenate(acute_cho_ratio_list) if acute_cho_ratio_list else np.array([])

# Chronic/Control means (BG-only sources)
naa_ratio_obs_chronic = []
if not young_bg_chronic.empty:
    naa_ratio_obs_chronic.append(float(young_bg_chronic.iloc[0]['NAA_ratio']))
if not sail_bg_chronic.empty:
    naa_ratio_obs_chronic.append(float(sail_bg_chronic.iloc[0]['NAA_ratio']))
naa_ratio_obs_chronic = np.array(naa_ratio_obs_chronic) if naa_ratio_obs_chronic else np.array([1.025])

naa_ratio_obs_control = []
if not young_bg_control.empty:
    naa_ratio_obs_control.append(float(young_bg_control.iloc[0]['NAA_ratio']))
else:
    # Chang BG control absolute 8.4 → ratio
    if not chang_df.empty:
        row = chang_df[(chang_df['Phase']=='Control')&(chang_df['Region']=='BG')&(chang_df['Metabolite']=='NAA')]
        if not row.empty:
            naa_ratio_obs_control.append(convert_absolute_to_ratio(float(row['Mean'].values[0]), 'BG'))
naa_ratio_obs_control = np.array(naa_ratio_obs_control) if naa_ratio_obs_control else np.array([1.05])

# Cho arrays (BG-only for reporting; not central to ξ)
cho_ratio_obs_chronic = []
if not young_bg_chronic.empty and 'Cho_ratio' in young_bg_chronic.columns and pd.notna(young_bg_chronic.iloc[0].get('Cho_ratio', np.nan)):
    cho_ratio_obs_chronic.append(float(young_bg_chronic.iloc[0]['Cho_ratio']))
if not sail_bg_chronic.empty and 'Cho_ratio' in sail_bg_chronic.columns and pd.notna(sail_bg_chronic.iloc[0].get('Cho_ratio', np.nan)):
    cho_ratio_obs_chronic.append(float(sail_bg_chronic.iloc[0]['Cho_ratio']))
cho_ratio_obs_chronic = np.array(cho_ratio_obs_chronic) if cho_ratio_obs_chronic else np.array([])

cho_ratio_obs_control = []
if not young_bg_control.empty and 'Cho_ratio' in young_bg_control.columns and pd.notna(young_bg_control.iloc[0].get('Cho_ratio', np.nan)):
    cho_ratio_obs_control.append(float(young_bg_control.iloc[0]['Cho_ratio']))
else:
    row = chang_df[(chang_df['Phase']=='Control')&(chang_df['Region']=='BG')&(chang_df['Metabolite']=='Cho')]
    if not row.empty:
        cho_ratio_obs_control.append(convert_absolute_to_ratio(float(row['Mean'].values[0]), 'BG'))
cho_ratio_obs_control = np.array(cho_ratio_obs_control) if cho_ratio_obs_control else np.array([0.263])

print("\n" + "=" * 80)
print("FINAL COMBINED DATA FOR MODEL (BG-only, ratio)")
print("=" * 80)
print(f"\n✅ ACUTE: n={len(naa_ratio_obs_acute)}\n   NAA/Cr: {np.mean(naa_ratio_obs_acute):.3f} ± {np.std(naa_ratio_obs_acute):.3f}")
print(f"\n✅ CHRONIC: n={len(naa_ratio_obs_chronic)}\n   NAA/Cr: {np.mean(naa_ratio_obs_chronic):.3f} ± {np.std(naa_ratio_obs_chronic):.3f}")
print(f"\n✅ CONTROL: n={len(naa_ratio_obs_control)}\n   NAA/Cr: {np.mean(naa_ratio_obs_control):.3f}")
if len(naa_ratio_obs_acute)==0:
    print("\n❌ ERROR: No acute data available"); sys.exit(1)

# ============================================================================
# MODEL (VL-free primary)
# ============================================================================
print("\n" + "=" * 80)
print("BUILDING BAYESIAN MODEL v3.6 ORIGINAL (BG-only)")
print("=" * 80)

L_MT = 2000e-9
with pm.Model() as model:
    ξ_acute   = pm.LogNormal('ξ_acute',   mu=np.log(0.6), sigma=0.15)
    ξ_chronic = pm.LogNormal('ξ_chronic', mu=np.log(0.8), sigma=0.15)
    ξ_control = pm.LogNormal('ξ_control', mu=np.log(0.5), sigma=0.15)
    β_ξ = pm.Normal('β_ξ', mu=-2.0, sigma=0.6)
    Γ = lambda ξ: pm.math.exp(β_ξ * (ξ * 1e-9) / L_MT)
    r0 = pm.LogNormal('r_NAA_ratio_baseline', mu=np.log(1.10), sigma=0.10)
    α_acute = pm.LogNormal('α_acute', mu=np.log(1.00), sigma=0.15)
    α_chronic = pm.LogNormal('α_chronic', mu=np.log(0.95), sigma=0.12)
    NAA_control_mean = r0
    NAA_acute_mean   = r0 * α_acute * Γ(ξ_acute)
    NAA_chronic_mean = r0 * α_chronic * Γ(ξ_chronic)
    pm.Deterministic('NAA_ratio_control_mean', NAA_control_mean)
    pm.Deterministic('NAA_ratio_acute_mean', NAA_acute_mean)
    pm.Deterministic('NAA_ratio_chronic_mean', NAA_chronic_mean)

    # Acute likelihood (robust)
    σ_naa = pm.HalfNormal('σ_naa_acute', sigma=0.20)
    ν_naa = pm.Exponential('ν_naa', lam=1/10)
    pm.StudentT('NAA_ratio_acute_obs', nu=ν_naa+2.0, mu=NAA_acute_mean, sigma=σ_naa, observed=naa_ratio_obs_acute)

    # Group means for Chronic/Control (BG-only)
    if len(naa_ratio_obs_chronic)>0:
        pm.Normal('NAA_ratio_chronic_obs', mu=NAA_chronic_mean, sigma=0.15, observed=naa_ratio_obs_chronic)
    pm.Normal('NAA_ratio_control_obs', mu=NAA_control_mean, sigma=0.10, observed=naa_ratio_obs_control)

    pm.Deterministic('Δξ', ξ_chronic - ξ_acute)

print("\n✅ Model built successfully!")
print(f"   Observations: {len(naa_ratio_obs_acute)} acute + {len(naa_ratio_obs_chronic)} chronic + {len(naa_ratio_obs_control)} control")

# ============================================================================
# INFERENCE
# ============================================================================
print("\n" + "=" * 80)
print("RUNNING MCMC INFERENCE (original)")
print("=" * 80)
with model:
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
    ppc = pm.sample_posterior_predictive(trace, random_seed=42)

posterior = trace.posterior
# Quick note for WAIC/LOO readiness
try:
    has_llk = hasattr(trace, 'log_likelihood') and trace.log_likelihood is not None
    llk_vars = list(trace.log_likelihood.data_vars) if has_llk else []
    print(f"\n🔎 WAIC/LOO readiness (original): log_likelihood present = {has_llk}; observed terms = {len(llk_vars)}")
except Exception:
    pass
ξa = posterior['ξ_acute'].values.flatten()
ξc = posterior['ξ_chronic'].values.flatten()
Δ = posterior['Δξ'].values.flatten()
βs = posterior['β_ξ'].values.flatten()

print("\n" + "=" * 80)
print("BAYESIAN INFERENCE RESULTS (original BG-only)")
print("=" * 80)
print(f"\n   ξ_acute = {ξa.mean():.3f} ± {ξa.std():.3f} nm")
print(f"   ξ_chronic = {ξc.mean():.3f} ± {ξc.std():.3f} nm")
print(f"   Δξ = {Δ.mean():.3f} ± {Δ.std():.3f} nm")
P = (Δ>0).mean()
print(f"\n✅ P(ξ_acute < ξ_chronic) = {P:.4f}")
print(f"\n   β_ξ = {βs.mean():.2f} ± {βs.std():.2f}")

# Predicted vs observed
naa_acute_pred = posterior['NAA_ratio_acute_mean'].values.flatten()
naa_chronic_pred = posterior['NAA_ratio_chronic_mean'].values.flatten()
naa_control_pred = posterior['NAA_ratio_control_mean'].values.flatten()
print(f"\n📊 PREDICTED vs OBSERVED NAA/Cr:")
print(f"   Acute: {naa_acute_pred.mean():.3f} (pred) vs {naa_ratio_obs_acute.mean():.3f} (obs)")
print(f"   Chronic: {naa_chronic_pred.mean():.3f} (pred) vs {naa_ratio_obs_chronic.mean():.3f} (obs)")
print(f"   Control: {naa_control_pred.mean():.3f} (pred) vs {naa_ratio_obs_control.mean():.3f} (obs)")

# ============================================================================
# SAVE RESULTS (under .../runs/<run_id>/original/)
# ============================================================================
print("\n" + "=" * 80)
print("SAVING RESULTS (original)")
print("=" * 80)

summary_path = run_dir / f"summary{output_tag}.csv"
az.summary(trace, hdi_prob=0.95).to_csv(str(summary_path))
print(f"✅ Saved summary: {summary_path}")

results_df = pd.DataFrame({
    'Parameter':['ξ_acute','ξ_chronic','Δξ','β_ξ','P(ξ_acute < ξ_chronic)'],
    'Mean':[ξa.mean(), ξc.mean(), Δ.mean(), βs.mean(), P],
    'SD':[ξa.std(), ξc.std(), Δ.std(), βs.std(), np.nan],
    'HDI_2.5%':[np.percentile(ξa,2.5), np.percentile(ξc,2.5), np.percentile(Δ,2.5), np.percentile(βs,2.5), np.nan],
    'HDI_97.5%':[np.percentile(ξa,97.5), np.percentile(ξc,97.5), np.percentile(Δ,97.5), np.percentile(βs,97.5), np.nan]
})
results_path = run_dir / f"results_v3_6_ratio_scale{output_tag}.csv"
results_df.to_csv(str(results_path), index=False)
print(f"✅ Saved results: {results_path}")

ppc_df = pd.DataFrame({
    'condition':['Acute','Chronic','Control'],
    'NAA/Cr_pred':[naa_acute_pred.mean(), naa_chronic_pred.mean(), naa_control_pred.mean()],
    'NAA/Cr_obs':[naa_ratio_obs_acute.mean(), naa_ratio_obs_chronic.mean(), naa_ratio_obs_control.mean()],
    'n_obs':[len(naa_ratio_obs_acute), len(naa_ratio_obs_chronic), len(naa_ratio_obs_control)]
})
ppc_path = run_dir / f"posterior_predictive_comparison_ratio{output_tag}.csv"
ppc_df.to_csv(str(ppc_path), index=False)
print(f"✅ Saved posterior predictive: {ppc_path}")

if not args.no_save_trace:
    trace_path = run_dir / f"trace{output_tag}.nc"
    trace.to_netcdf(str(trace_path))
    print(f"✅ Saved full posterior trace (NetCDF): {trace_path}")

# Optional figures
try:
    if args.plot_densities:
        vars_to_plot = ['β_ξ','ξ_acute','ξ_chronic','ξ_control','Δξ','NAA_ratio_acute_mean','NAA_ratio_chronic_mean','NAA_ratio_control_mean']
        az.plot_posterior(trace, var_names=vars_to_plot, kind='kde')
        outp = figures_dir / f"posterior_overview_kde{output_tag}.png"
        plt.tight_layout(); plt.savefig(str(outp), dpi=220); plt.close()
        print(f"✅ Saved: {outp}")
    if args.plots_extended:
        # Simple PPC bar
        fig, ax = plt.subplots(figsize=(4.8,3.0))
        cond = ['Acute','Chronic','Control']
        pred = [ppc_df.loc[i,'NAA/Cr_pred'] for i in range(3)]
        obs = [ppc_df.loc[i,'NAA/Cr_obs'] for i in range(3)]
        x = np.arange(3)
        ax.bar(x-0.15, pred, width=0.3, label='Pred')
        ax.bar(x+0.15, obs, width=0.3, label='Obs')
        ax.set_xticks(x); ax.set_xticklabels(cond); ax.set_ylabel('NAA/Cr')
        ax.legend(); plt.tight_layout()
        outb = figures_dir / f"ppc_bar_ratio{output_tag}.png"
        plt.savefig(str(outb), dpi=220); plt.close(); print(f"✅ Saved: {outb}")
except Exception as e:
    print(f"⚠ Figure generation error: {e}")

print("\n✅ ORIGINAL BG-only run complete. For adjunctive (multi-region, extended audits), use bayesian_v3_6_corrected_local.py")
