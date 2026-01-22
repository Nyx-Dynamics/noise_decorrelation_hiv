"""
Enhanced Bayesian Inference v2.2 (FINAL)
=========================================

FINAL FIX: Remove homeostatic floor entirely
- Let natural quantum coupling + astrocyte compensation work
- NAA_base increased to 1.20 for better healthy baseline
- Trust the biological mechanisms without artificial constraints

Your v2.1 results showed:
- Astrocyte compensation works perfectly (1.183)
- Floor was interfering with natural dynamics
- Model wants chronic NAA = 0.76 × 1.183 = 0.90, not forced 0.939
"""

from __future__ import annotations

import os
import argparse
import numpy as np
import pandas as pd

import pymc as pm
import pytensor.tensor as pt
import arviz as az

# --------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------
class CONST:
    xi_baseline = 0.8e-9
    sigma_r_regular = 0.38e-9
    NAA_baseline = 1.105
    Cho_baseline = 0.225


# --------------------------------------------------------------------------------------
# Data
# --------------------------------------------------------------------------------------
CONDITIONS = ("healthy", "acute_HIV", "chronic_HIV")

OBS_NAA = {"healthy": 1.105, "acute_HIV": 1.135, "chronic_HIV": 1.005}
OBS_CHO = {"healthy": 0.225, "acute_HIV": 0.245, "chronic_HIV": 0.235}

COHERENCE_BASE = {"healthy": 0.85, "acute_HIV": 0.84, "chronic_HIV": 0.73}
SIGMA_R = {
    "healthy": CONST.sigma_r_regular,
    "acute_HIV": CONST.sigma_r_regular * 1.05,
    "chronic_HIV": CONST.sigma_r_regular * 1.4,
}
MEMBRANE_RATES = {
    "healthy": (1.0, 1.0),
    "acute_HIV": (2.0, 2.5),
    "chronic_HIV": (1.2, 1.1),
}


# --------------------------------------------------------------------------------------
# Forward Models (FINAL - No Floor)
# --------------------------------------------------------------------------------------

def coherence_from_xi_nonlinear(xi, coherence_base, xi_floor, xi_ceiling):
    """Nonlinear ξ → coherence with floor."""
    xi_normalized = (xi - xi_floor) / (xi_ceiling - xi_floor)
    xi_normalized = pt.clip(xi_normalized, 0.0, 1.0)
    
    C_floor = 0.65
    coherence = C_floor + (coherence_base - C_floor) * (1 - xi_normalized) ** 2
    return coherence


def forward_NAA_final(coherence_base, xi, sigma_r, condition,
                     coh_exp, xi_exp, deloc_exp, NAA_base, astrocyte_comp,
                     xi_floor, xi_ceiling):
    """
    FINAL NAA model - NO artificial floor.
    
    Let quantum coupling + astrocyte compensation determine NAA naturally.
    """
    # 1. Nonlinear coherence
    coherence_effective = coherence_from_xi_nonlinear(
        xi, coherence_base, xi_floor, xi_ceiling
    )
    
    # 2. Quantum coupling
    coherence_term = (coherence_effective / 0.85) ** coh_exp
    xi_protection = (CONST.xi_baseline / xi) ** xi_exp
    deloc_term = (sigma_r / CONST.sigma_r_regular) ** deloc_exp
    
    NAA_quantum = NAA_base * coherence_term * xi_protection * deloc_term
    
    # 3. Astrocyte compensation (chronic only)
    if condition == 'chronic_HIV':
        NAA_total = NAA_quantum * astrocyte_comp
    else:
        NAA_total = NAA_quantum
    
    # 4. NO FLOOR - trust the mechanisms
    return NAA_total


def forward_Cho(damage_rate, repair_rate, k_turnover):
    """Choline dynamics."""
    turnover_factor = (damage_rate + repair_rate) - 1.0
    return CONST.Cho_baseline * (1.0 + k_turnover * turnover_factor)


# --------------------------------------------------------------------------------------
# Bayesian Model (FINAL)
# --------------------------------------------------------------------------------------

def build_model_final():
    """Final model - no artificial constraints."""
    
    naa_obs = np.array([OBS_NAA[c] for c in CONDITIONS])
    cho_obs = np.array([OBS_CHO[c] for c in CONDITIONS])
    
    with pm.Model() as model:
        # ξ parameters
        xi_floor = pm.TruncatedNormal("xi_floor", mu=0.35e-9, sigma=0.1e-9, 
                                      lower=0.2e-9, upper=0.5e-9)
        xi_ceiling = pm.TruncatedNormal("xi_ceiling", mu=0.8e-9, sigma=0.1e-9,
                                        lower=0.6e-9, upper=1.0e-9)
        
        xi_healthy = pm.TruncatedNormal("xi_healthy", mu=0.75e-9, sigma=0.1e-9,
                                        lower=xi_floor, upper=xi_ceiling)
        xi_acute = pm.TruncatedNormal("xi_acute", mu=0.4e-9, sigma=0.1e-9,
                                      lower=xi_floor, upper=0.6e-9)
        xi_chronic = pm.TruncatedNormal("xi_chronic", mu=0.8e-9, sigma=0.1e-9,
                                        lower=0.5e-9, upper=xi_ceiling)
        
        # Coupling exponents
        coh_exp = pm.TruncatedNormal("coh_exp", mu=2.5, sigma=0.5, lower=0.5, upper=5.0)
        xi_exp = pm.TruncatedNormal("xi_exp", mu=0.3, sigma=0.2, lower=0.0, upper=1.5)
        deloc_exp = pm.TruncatedNormal("deloc_exp", mu=0.2, sigma=0.1, lower=0.0, upper=1.0)
        
        # FINAL: Even higher NAA_base to reach healthy levels naturally
        NAA_base = pm.TruncatedNormal("NAA_base", mu=1.20, sigma=0.10,
                                     lower=1.0, upper=1.40)
        k_turnover = pm.TruncatedNormal("k_turnover", mu=0.02, sigma=0.01,
                                       lower=0.0, upper=0.1)
        
        # Astrocyte compensation
        astrocyte_comp = pm.TruncatedNormal("astrocyte_comp", mu=1.18, sigma=0.05,
                                           lower=1.05, upper=1.30)
        
        # Observation noise
        sigma_NAA = pm.HalfNormal("sigma_NAA", sigma=0.05)
        sigma_Cho = pm.HalfNormal("sigma_Cho", sigma=0.01)
        
        # Forward model
        naa_preds = []
        cho_preds = []
        
        for cond in CONDITIONS:
            if cond == "healthy":
                xi_cond = xi_healthy
            elif cond == "acute_HIV":
                xi_cond = xi_acute
            else:
                xi_cond = xi_chronic
            
            # FINAL: No floor
            naa_pred = forward_NAA_final(
                coherence_base=COHERENCE_BASE[cond],
                xi=xi_cond,
                sigma_r=SIGMA_R[cond],
                condition=cond,
                coh_exp=coh_exp,
                xi_exp=xi_exp,
                deloc_exp=deloc_exp,
                NAA_base=NAA_base,
                astrocyte_comp=astrocyte_comp,
                xi_floor=xi_floor,
                xi_ceiling=xi_ceiling
            )
            naa_preds.append(naa_pred)
            
            dmg, rep = MEMBRANE_RATES[cond]
            cho_pred = forward_Cho(dmg, rep, k_turnover)
            cho_preds.append(cho_pred)
        
        naa_preds = pt.stack(naa_preds)
        cho_preds = pt.stack(cho_preds)
        
        pm.Normal("NAA_obs", mu=naa_preds, sigma=sigma_NAA, observed=naa_obs)
        pm.Normal("Cho_obs", mu=cho_preds, sigma=sigma_Cho, observed=cho_obs)
        
        # Derived quantities
        pm.Deterministic("delta_xi", xi_chronic - xi_acute)
        pm.Deterministic("xi_healthy_nm", xi_healthy * 1e9)
        pm.Deterministic("xi_acute_nm", xi_acute * 1e9)
        pm.Deterministic("xi_chronic_nm", xi_chronic * 1e9)
        pm.Deterministic("xi_floor_nm", xi_floor * 1e9)
        pm.Deterministic("xi_ceiling_nm", xi_ceiling * 1e9)
        
        coh_healthy = coherence_from_xi_nonlinear(xi_healthy, 0.85, xi_floor, xi_ceiling)
        coh_acute = coherence_from_xi_nonlinear(xi_acute, 0.84, xi_floor, xi_ceiling)
        coh_chronic = coherence_from_xi_nonlinear(xi_chronic, 0.73, xi_floor, xi_ceiling)
        
        pm.Deterministic("coh_healthy_eff", coh_healthy)
        pm.Deterministic("coh_acute_eff", coh_acute)
        pm.Deterministic("coh_chronic_eff", coh_chronic)
    
    return model


# --------------------------------------------------------------------------------------
# Run Inference
# --------------------------------------------------------------------------------------

def run_inference_final(draws=3000, tune=1500, chains=4, target_accept=0.92, seed=42):
    """Run FINAL inference without artificial floors."""
    
    os.makedirs("results/bayesian_v2_final", exist_ok=True)
    
    print("="*80)
    print(" ENHANCED BAYESIAN INFERENCE v2.2 (FINAL)")
    print("="*80)
    print("\nFINAL VERSION:")
    print("  ✓ No artificial floor - trust natural dynamics")
    print("  ✓ NAA_base prior: μ=1.20 (higher baseline)")
    print("  ✓ Astrocyte compensation: natural regulation")
    print("  ✓ Nonlinear ξ-coherence: biological floor at 0.65")
    print()
    
    model = build_model_final()
    
    with model:
        idata = pm.sample(
            draws=draws, tune=tune, chains=chains,
            target_accept=target_accept, random_seed=seed,
            return_inferencedata=True, progressbar=True,
        )
    
    trace_path = "results/bayesian_v2_final/trace_final.nc"
    idata.to_netcdf(trace_path)
    
    var_names = [
        "coh_exp", "xi_exp", "deloc_exp", "NAA_base", "k_turnover",
        "astrocyte_comp", "sigma_NAA", "sigma_Cho",
        "xi_healthy_nm", "xi_acute_nm", "xi_chronic_nm",
        "xi_floor_nm", "xi_ceiling_nm", "delta_xi",
        "coh_healthy_eff", "coh_acute_eff", "coh_chronic_eff"
    ]
    
    summary = az.summary(idata, var_names=var_names, round_to=4)
    summary_path = "results/bayesian_v2_final/summary_final.csv"
    summary.to_csv(summary_path)
    
    print("\n" + "="*80)
    print(" POSTERIOR SUMMARY")
    print("="*80)
    print(summary)
    
    # P(ξ_acute < ξ_chronic)
    xi_acute_vals = idata.posterior["xi_acute_nm"].values.reshape(-1)
    xi_chronic_vals = idata.posterior["xi_chronic_nm"].values.reshape(-1)
    p_order = float(np.mean(xi_acute_vals < xi_chronic_vals))
    
    # Posterior predictive
    post = idata.posterior
    med = {k: float(np.median(post[k].values)) 
           for k in ["coh_exp", "xi_exp", "deloc_exp", "NAA_base", 
                     "k_turnover", "astrocyte_comp", "xi_floor", "xi_ceiling"]}
    
    preds = []
    for cond in CONDITIONS:
        if cond == "healthy":
            xi_med = float(np.median(post["xi_healthy_nm"].values)) * 1e-9
        elif cond == "acute_HIV":
            xi_med = float(np.median(post["xi_acute_nm"].values)) * 1e-9
        else:
            xi_med = float(np.median(post["xi_chronic_nm"].values)) * 1e-9
        
        # Coherence
        xi_norm = (xi_med - med["xi_floor"]) / (med["xi_ceiling"] - med["xi_floor"])
        xi_norm = np.clip(xi_norm, 0.0, 1.0)
        coh_eff = 0.65 + (COHERENCE_BASE[cond] - 0.65) * (1 - xi_norm) ** 2
        
        # Quantum
        naa_quantum = (
            med["NAA_base"]
            * (coh_eff / 0.85) ** med["coh_exp"]
            * (CONST.xi_baseline / xi_med) ** med["xi_exp"]
            * (SIGMA_R[cond] / CONST.sigma_r_regular) ** med["deloc_exp"]
        )
        
        # Compensation (no floor)
        if cond == "chronic_HIV":
            naa = naa_quantum * med["astrocyte_comp"]
        else:
            naa = naa_quantum
        
        # Choline
        dmg, rep = MEMBRANE_RATES[cond]
        cho = CONST.Cho_baseline * (1.0 + med["k_turnover"] * ((dmg + rep) - 1.0))
        
        preds.append({
            "condition": cond,
            "NAA_pred": naa,
            "Cho_pred": cho,
            "NAA_obs": OBS_NAA[cond],
            "Cho_obs": OBS_CHO[cond],
            "error_NAA_%": 100 * (naa - OBS_NAA[cond]) / OBS_NAA[cond],
            "error_Cho_%": 100 * (cho - OBS_CHO[cond]) / OBS_CHO[cond]
        })
    
    preds_df = pd.DataFrame(preds)
    preds_path = "results/bayesian_v2_final/posterior_predictive_final.csv"
    preds_df.to_csv(preds_path, index=False)
    
    print("\n" + "="*80)
    print(" POSTERIOR PREDICTIVE CHECK")
    print("="*80)
    print(preds_df.to_string(index=False))
    
    print("\n" + "="*80)
    print(" KEY RESULTS")
    print("="*80)
    print(f"P(ξ_acute < ξ_chronic) = {p_order:.4f}")
    print(f"\nPosterior Medians:")
    print(f"  NAA_base: {med['NAA_base']:.3f}")
    print(f"  Astrocyte comp: {med['astrocyte_comp']:.3f}")
    print(f"  ξ floor: {med['xi_floor']*1e9:.2f} nm")
    print(f"  ξ ceiling: {med['xi_ceiling']*1e9:.2f} nm")
    
    print(f"\nNAA Errors:")
    for i, cond in enumerate(CONDITIONS):
        print(f"  {cond:12s}: {preds_df.loc[i, 'error_NAA_%']:+.1f}%")
    
    # Save summary
    with open("results/bayesian_v2_final/results_summary.txt", "w") as f:
        f.write("="*80 + "\n")
        f.write(" FINAL MODEL v2.2 - RESULTS\n")
        f.write("="*80 + "\n\n")
        f.write(f"P(ξ_acute < ξ_chronic) = {p_order:.4f}\n\n")
        f.write("Posterior Medians:\n")
        for k, v in med.items():
            if "xi" in k and k != "xi_exp":
                f.write(f"  {k}: {v*1e9:.3f} nm\n")
            else:
                f.write(f"  {k}: {v:.3f}\n")
        f.write("\n" + preds_df.to_string(index=False) + "\n")
    
    print(f"\nResults saved to: results/bayesian_v2_final/")
    
    return {"idata": idata, "summary": summary, "predictions": preds_df,
            "p_xi_order": p_order, "medians": med}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--draws", type=int, default=3000)
    parser.add_argument("--tune", type=int, default=1500)
    parser.add_argument("--chains", type=int, default=4)
    parser.add_argument("--target-accept", type=float, default=0.92)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    
    run_inference_final(args.draws, args.tune, args.chains, args.target_accept, args.seed)


if __name__ == "__main__":
    main()
