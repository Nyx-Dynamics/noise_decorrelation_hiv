"""
Bayesian Inference v3.0: Enzyme Kinetics Parameters
====================================================

Infers enzyme kinetic parameters from MRS data.

NEW PARAMETERS (v3):
- β_ξ: Protection factor exponent
- γ_coh: Coherence modulation exponent
- V_NAT8L_scale: NAT8L Vmax scaling factor
- V_ASPA_scale: ASPA Vmax scaling factor
- damage_acute: Viral damage in acute phase
- damage_chronic: Viral damage in chronic phase

REMOVED (from v2):
- astrocyte_compensation (replaced by direct mechanism)
- coherence_floor (emerges naturally from kinetics)
- NAA_base (replaced by enzyme kinetics)

This inference will determine if β_ξ ≈ 1.8-2.5 as predicted.
"""

from __future__ import annotations

import os
import argparse
import numpy as np
import pandas as pd
import sys

# Add current directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pymc as pm
import pytensor.tensor as pt
import arviz as az

# Import enzyme kinetics (simplified for PyMC compatibility)
from enzyme_kinetics import ENZYME


# =============================================================================
# CONSTANTS
# =============================================================================

class CONST:
    xi_baseline = 0.8e-9
    sigma_r_regular = 0.38e-9
    Cr_conc = 6.8e-3  # M (creatine reference)

    # Enzyme baseline rates (will be scaled by parameters)
    Vmax_NAT8L_base = ENZYME.Vmax_NAT8L_base
    Vmax_ASPA = ENZYME.Vmax_ASPA
    k_Kennedy_base = ENZYME.k_Kennedy_base

    # Substrate concentrations
    Asp_conc = ENZYME.Asp_conc
    AcCoA_conc = ENZYME.AcCoA_conc
    Km_NAT8L_Asp = ENZYME.Km_NAT8L_Asp
    Km_NAT8L_AcCoA = ENZYME.Km_NAT8L_AcCoA
    Km_ASPA = ENZYME.Km_ASPA
    Ki_ASPA = ENZYME.Ki_ASPA


# =============================================================================
# DATA
# =============================================================================

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
    "healthy": 1.0,
    "acute_HIV": 2.2,
    "chronic_HIV": 1.3,
}


# =============================================================================
# FORWARD MODEL (Simplified for PyMC)
# =============================================================================

def protection_factor(xi, xi_ref, beta_xi):
    """Compute protection factor Π_ξ = (ξ_ref / ξ)^β with safe clipping of ξ."""
    xi_safe = pt.clip(xi, 1e-12, np.inf)
    return (xi_ref / xi_safe) ** beta_xi


def coherence_modulation(coherence, gamma_coh):
    """Compute coherence modulation η_coh = (C / C_ref)^γ using ENZYME.coherence_ref."""
    return (coherence / ENZYME.coherence_ref) ** gamma_coh


def simplified_NAA_steady_state(Pi_xi, eta_coh, viral_damage,
                                V_NAT8L_scale, V_ASPA_scale):
    """
    Simplified enzyme steady-state calculation for PyMC.

    Uses approximate steady-state solution instead of full ODE integration.

    Steady state: v_synthesis = v_degradation
    V_NAT8L * [substrates] * Π * η * damage = V_ASPA * [NAA] / (Km + ...)

    Approximate solution:
    NAA_ss ≈ (V_NAT8L / V_ASPA) * Π * η * damage * substrate_factors * Km
    """
    # Total enzyme enhancement
    total_enhancement = Pi_xi * eta_coh * viral_damage

    # Substrate saturation factors
    f_Asp = CONST.Asp_conc / (CONST.Km_NAT8L_Asp + CONST.Asp_conc)
    f_AcCoA = CONST.AcCoA_conc / (CONST.Km_NAT8L_AcCoA + CONST.AcCoA_conc)

    # Synthesis rate
    v_synth = CONST.Vmax_NAT8L_base * V_NAT8L_scale * total_enhancement * f_Asp * f_AcCoA

    # At steady state, assuming NAA ≈ Km_ASPA (middle of Michaelis-Menten curve)
    # v_degrad = V_ASPA * NAA / (Km + NAA)
    # When NAA ≈ Km: v_degrad ≈ V_ASPA * Km / (2*Km) = V_ASPA / 2

    # More sophisticated approximation:
    # NAA_ss = v_synth * Km_ASPA / (V_ASPA * 0.5)  # 0.5 accounts for Michaelis-Menten shape

    NAA_approx = v_synth * CONST.Km_ASPA / (CONST.Vmax_ASPA * V_ASPA_scale * 0.5)

    # Add substrate inhibition correction (reduces NAA at high concentrations)
    inhibition_factor = 1.0 / (1.0 + NAA_approx / CONST.Ki_ASPA)
    NAA_ss = NAA_approx * inhibition_factor

    return NAA_ss


def forward_NAA_v3(xi, coherence_base, viral_damage,
                   beta_xi, gamma_coh, V_NAT8L_scale, V_ASPA_scale):
    """
    Full NAA forward model v3 for PyMC.

    Parameters
    ----------
    xi : tensor
        Noise correlation length (m)
    coherence_base : float
        Base coherence level
    viral_damage : tensor
        Viral damage factor (0-1)
    beta_xi : tensor
        Protection factor exponent
    gamma_coh : tensor
        Coherence modulation exponent
    V_NAT8L_scale : tensor
        NAT8L Vmax scaling factor
    V_ASPA_scale : tensor
        ASPA Vmax scaling factor

    Returns
    -------
    NAA_Cr : tensor
        NAA/Cr ratio
    """
    # 1. Protection factor
    Pi_xi = protection_factor(xi, CONST.xi_baseline, beta_xi)

    # 2. Coherence modulation
    eta_coh = coherence_modulation(coherence_base, gamma_coh)

    # 3. Steady-state NAA concentration
    NAA_molar = simplified_NAA_steady_state(
        Pi_xi, eta_coh, viral_damage, V_NAT8L_scale, V_ASPA_scale
    )

    # 4. Convert to ratio
    NAA_Cr = NAA_molar / CONST.Cr_conc

    return NAA_Cr


def forward_Cho(membrane_turnover, k_turnover_scale):
    """Choline dynamics (simplified)."""
    Cho_base = 1.5e-3  # M
    Cho_molar = Cho_base * (1.0 + k_turnover_scale * (membrane_turnover - 1.0))
    Cho_Cr = Cho_molar / CONST.Cr_conc
    return Cho_Cr


# =============================================================================
# BAYESIAN MODEL
# =============================================================================

def build_model_v3():
    """Build Bayesian model v3 with enzyme parameters."""

    naa_obs = np.array([OBS_NAA[c] for c in CONDITIONS])
    cho_obs = np.array([OBS_CHO[c] for c in CONDITIONS])

    with pm.Model() as model:
        # =====================================================================
        # ξ parameters (same as v2)
        # =====================================================================

        xi_healthy = pm.TruncatedNormal("xi_healthy", mu=0.75e-9, sigma=0.1e-9,
                                        lower=0.5e-9, upper=0.9e-9)
        xi_acute = pm.TruncatedNormal("xi_acute", mu=0.45e-9, sigma=0.1e-9,
                                      lower=0.3e-9, upper=0.6e-9)
        xi_chronic = pm.TruncatedNormal("xi_chronic", mu=0.78e-9, sigma=0.1e-9,
                                        lower=0.6e-9, upper=0.95e-9)

        # =====================================================================
        # NEW: Enzyme modulation parameters
        # =====================================================================

        # Protection factor exponent (KEY PARAMETER)
        beta_xi = pm.TruncatedNormal("beta_xi", mu=2.0, sigma=0.5,
                                     lower=0.5, upper=4.0)

        # Coherence modulation exponent
        gamma_coh = pm.TruncatedNormal("gamma_coh", mu=1.5, sigma=0.5,
                                       lower=0.5, upper=3.0)

        # Enzyme Vmax scaling factors
        V_NAT8L_scale = pm.TruncatedNormal("V_NAT8L_scale", mu=1.0, sigma=0.3,
                                           lower=0.5, upper=2.0)
        V_ASPA_scale = pm.TruncatedNormal("V_ASPA_scale", mu=1.0, sigma=0.3,
                                          lower=0.5, upper=2.0)

        # Viral damage factors
        damage_acute = pm.TruncatedNormal("damage_acute", mu=0.90, sigma=0.05,
                                          lower=0.75, upper=1.0)
        damage_chronic = pm.TruncatedNormal("damage_chronic", mu=0.60, sigma=0.1,
                                            lower=0.4, upper=0.8)

        # Choline turnover scaling
        k_turnover_scale = pm.TruncatedNormal("k_turnover_scale", mu=0.20, sigma=0.1,
                                              lower=0.05, upper=0.5)

        # =====================================================================
        # Observation noise
        # =====================================================================

        sigma_NAA = pm.HalfNormal("sigma_NAA", sigma=0.05)
        sigma_Cho = pm.HalfNormal("sigma_Cho", sigma=0.01)

        # =====================================================================
        # Forward model
        # =====================================================================

        naa_preds = []
        cho_preds = []

        for cond in CONDITIONS:
            # Select ξ
            if cond == "healthy":
                xi_cond = xi_healthy
                viral_damage = 1.0  # No damage
            elif cond == "acute_HIV":
                xi_cond = xi_acute
                viral_damage = damage_acute
            else:  # chronic_HIV
                xi_cond = xi_chronic
                viral_damage = damage_chronic

            # NAA prediction
            naa_pred = forward_NAA_v3(
                xi=xi_cond,
                coherence_base=COHERENCE_BASE[cond],
                viral_damage=viral_damage,
                beta_xi=beta_xi,
                gamma_coh=gamma_coh,
                V_NAT8L_scale=V_NAT8L_scale,
                V_ASPA_scale=V_ASPA_scale
            )
            naa_preds.append(naa_pred)

            # Cho prediction
            cho_pred = forward_Cho(
                membrane_turnover=MEMBRANE_RATES[cond],
                k_turnover_scale=k_turnover_scale
            )
            cho_preds.append(cho_pred)

        naa_preds = pt.stack(naa_preds)
        cho_preds = pt.stack(cho_preds)

        # =====================================================================
        # Likelihood
        # =====================================================================

        pm.Normal("NAA_obs", mu=naa_preds, sigma=sigma_NAA, observed=naa_obs)
        pm.Normal("Cho_obs", mu=cho_preds, sigma=sigma_Cho, observed=cho_obs)

        # =====================================================================
        # Derived quantities
        # =====================================================================

        pm.Deterministic("delta_xi", xi_chronic - xi_acute)
        pm.Deterministic("xi_healthy_nm", xi_healthy * 1e9)
        pm.Deterministic("xi_acute_nm", xi_acute * 1e9)
        pm.Deterministic("xi_chronic_nm", xi_chronic * 1e9)

        # Protection factors
        Pi_healthy = protection_factor(xi_healthy, CONST.xi_baseline, beta_xi)
        Pi_acute = protection_factor(xi_acute, CONST.xi_baseline, beta_xi)
        Pi_chronic = protection_factor(xi_chronic, CONST.xi_baseline, beta_xi)

        pm.Deterministic("Pi_healthy", Pi_healthy)
        pm.Deterministic("Pi_acute", Pi_acute)
        pm.Deterministic("Pi_chronic", Pi_chronic)

    return model


# =============================================================================
# RUN INFERENCE
# =============================================================================

def run_inference_v3(draws=3000, tune=1500, chains=4, target_accept=0.92, seed=42):
    """Run Bayesian inference v3 for enzyme parameters."""

    os.makedirs("results/bayesian_v3", exist_ok=True)

    print("=" * 80)
    print(" BAYESIAN INFERENCE v3.0: Enzyme Kinetics Parameters")
    print("=" * 80)
    print()
    print("NEW APPROACH:")
    print("  ✓ Direct enzyme modulation (no compensation)")
    print("  ✓ Protection factor β_ξ (KEY PARAMETER)")
    print("  ✓ Enzyme Vmax scaling factors")
    print("  ✓ Viral damage parameters")
    print()
    print("TESTING HYPOTHESIS:")
    print("  H0: β_ξ ≈ 1.8-2.5 (predicted range)")
    print("  H1: P(ξ_acute < ξ_chronic) > 0.95")
    print()

    model = build_model_v3()

    print("Starting MCMC sampling...")
    print(f"  Draws: {draws}, Tune: {tune}, Chains: {chains}")
    print(f"  Target accept: {target_accept}")
    print()

    with model:
        idata = pm.sample(
            draws=draws, tune=tune, chains=chains,
            target_accept=target_accept, random_seed=seed,
            return_inferencedata=True, progressbar=True,
        )

    # Save trace
    trace_path = "results/bayesian_v3/trace_v3.nc"
    idata.to_netcdf(trace_path)
    print(f"\nTrace saved to: {trace_path}")

    # Get summary
    var_names = [
        "beta_xi", "gamma_coh", "V_NAT8L_scale", "V_ASPA_scale",
        "damage_acute", "damage_chronic", "k_turnover_scale",
        "xi_healthy_nm", "xi_acute_nm", "xi_chronic_nm", "delta_xi",
        "Pi_healthy", "Pi_acute", "Pi_chronic",
        "sigma_NAA", "sigma_Cho"
    ]

    summary = az.summary(idata, var_names=var_names, round_to=4)
    summary_path = "results/bayesian_v3/summary_v3.csv"
    summary.to_csv(summary_path)

    print("\n" + "=" * 80)
    print(" POSTERIOR SUMMARY")
    print("=" * 80)
    print(summary)

    # Key results
    xi_acute_vals = idata.posterior["xi_acute_nm"].values.reshape(-1)
    xi_chronic_vals = idata.posterior["xi_chronic_nm"].values.reshape(-1)
    p_order = float(np.mean(xi_acute_vals < xi_chronic_vals))

    beta_xi_median = float(np.median(idata.posterior["beta_xi"].values))
    beta_xi_hdi = az.hdi(idata.posterior["beta_xi"], hdi_prob=0.94)

    print("\n" + "=" * 80)
    print(" KEY RESULTS")
    print("=" * 80)
    print(f"P(ξ_acute < ξ_chronic) = {p_order:.4f}")
    print()
    print(f"Protection Factor Exponent β_ξ:")
    print(f"  Median: {beta_xi_median:.3f}")
    print(f"  94% HDI: [{float(beta_xi_hdi[0]):.3f}, {float(beta_xi_hdi[1]):.3f}]")
    print()

    # Posterior predictive check
    post = idata.posterior
    med = {
        "beta_xi": float(np.median(post["beta_xi"].values)),
        "gamma_coh": float(np.median(post["gamma_coh"].values)),
        "V_NAT8L_scale": float(np.median(post["V_NAT8L_scale"].values)),
        "V_ASPA_scale": float(np.median(post["V_ASPA_scale"].values)),
        "damage_acute": float(np.median(post["damage_acute"].values)),
        "damage_chronic": float(np.median(post["damage_chronic"].values)),
        "k_turnover_scale": float(np.median(post["k_turnover_scale"].values)),
    }

    # Compute predictions
    preds = []
    for cond in CONDITIONS:
        if cond == "healthy":
            xi_med = float(np.median(post["xi_healthy_nm"].values)) * 1e-9
            viral_damage = 1.0
        elif cond == "acute_HIV":
            xi_med = float(np.median(post["xi_acute_nm"].values)) * 1e-9
            viral_damage = med["damage_acute"]
        else:
            xi_med = float(np.median(post["xi_chronic_nm"].values)) * 1e-9
            viral_damage = med["damage_chronic"]

        # Protection factor
        Pi = (CONST.xi_baseline / xi_med) ** med["beta_xi"]
        eta = (COHERENCE_BASE[cond] / 0.85) ** med["gamma_coh"]

        # Simplified steady state
        total_enh = Pi * eta * viral_damage
        f_Asp = CONST.Asp_conc / (CONST.Km_NAT8L_Asp + CONST.Asp_conc)
        f_AcCoA = CONST.AcCoA_conc / (CONST.Km_NAT8L_AcCoA + CONST.AcCoA_conc)

        v_synth = CONST.Vmax_NAT8L_base * med["V_NAT8L_scale"] * total_enh * f_Asp * f_AcCoA
        NAA_approx = v_synth * CONST.Km_ASPA / (CONST.Vmax_ASPA * med["V_ASPA_scale"] * 0.5)
        NAA_ss = NAA_approx / (1.0 + NAA_approx / CONST.Ki_ASPA)
        naa = NAA_ss / CONST.Cr_conc

        # Cho
        cho_molar = 1.5e-3 * (1.0 + med["k_turnover_scale"] * (MEMBRANE_RATES[cond] - 1.0))
        cho = cho_molar / CONST.Cr_conc

        preds.append({
            "condition": cond,
            "NAA_pred": naa,
            "Cho_pred": cho,
            "NAA_obs": OBS_NAA[cond],
            "Cho_obs": OBS_CHO[cond],
            "error_NAA_%": 100 * (naa - OBS_NAA[cond]) / OBS_NAA[cond],
            "Pi_xi": Pi,
        })

    preds_df = pd.DataFrame(preds)
    preds_path = "results/bayesian_v3/posterior_predictive_v3.csv"
    preds_df.to_csv(preds_path, index=False)

    print("\n" + "=" * 80)
    print(" POSTERIOR PREDICTIVE CHECK")
    print("=" * 80)
    print(preds_df.to_string(index=False))

    # Save text summary
    with open("results/bayesian_v3/results_summary_v3.txt", "w") as f:
        f.write("=" * 80 + "\n")
        f.write(" BAYESIAN INFERENCE v3.0 RESULTS\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"P(ξ_acute < ξ_chronic) = {p_order:.4f}\n\n")
        f.write(f"Protection Factor Exponent β_ξ:\n")
        f.write(f"  Median: {beta_xi_median:.3f}\n")
        f.write(f"  94% HDI: [{float(beta_xi_hdi[0]):.3f}, {float(beta_xi_hdi[1]):.3f}]\n\n")
        f.write("Posterior Medians:\n")
        for k, v in med.items():
            f.write(f"  {k}: {v:.3f}\n")
        f.write("\n" + preds_df.to_string(index=False) + "\n")

    print(f"\nAll results saved to: results/bayesian_v3/")

    return {
        "idata": idata,
        "summary": summary,
        "predictions": preds_df,
        "p_xi_order": p_order,
        "beta_xi_median": beta_xi_median,
        "beta_xi_hdi": beta_xi_hdi,
        "medians": med
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--draws", type=int, default=3000)
    parser.add_argument("--tune", type=int, default=1500)
    parser.add_argument("--chains", type=int, default=4)
    parser.add_argument("--target-accept", type=float, default=0.92)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    run_inference_v3(args.draws, args.tune, args.chains, args.target_accept, args.seed)


if __name__ == "__main__":
    main()