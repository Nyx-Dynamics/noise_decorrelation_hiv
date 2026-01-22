"""
Bayesian parameter inference (PyMC v4) to optimize coupling exponents while preserving mechanism.

- Infers exponents for NAA coupling: coherence exponent, xi exponent, delocalization exponent
- Infers Cho turnover coupling coefficient
- Includes latent ξ per condition (healthy, acute, chronic) with informative priors
- Fits to Sailasuta et al. (2012) target NAA/Cr and Cho/Cr values
- Reports P(xi_acute < xi_chronic)

Usage:
    python -m quantum.bayesian_optimization --draws 2000 --tune 1000 --chains 4 --target-accept 0.9

Outputs:
    results/bayesian/trace.nc                     (inference data)
    results/bayesian/summary.csv                  (posterior summary)
    results/bayesian/posterior_predictive.csv     (PPD means and 95% CI)
    Console: P(xi_acute < xi_chronic)
"""

from __future__ import annotations

import argparse
import os

import arviz as az
import numpy as np
import pandas as pd
# PyMC v4
import pymc as pm
import pytensor.tensor as pt

from .coupling_functions import CONST

# --------------------------------------------------------------------------------------
# Data: Sailasuta et al. (2012) target means for BG/FC averages (approx)
# --------------------------------------------------------------------------------------
CONDITIONS = ("healthy", "acute_HIV", "chronic_HIV")

OBS_NAA = {
    "healthy": 1.105,
    "acute_HIV": 1.135,
    "chronic_HIV": 1.005,
}

OBS_CHO = {
    "healthy": 0.225,
    "acute_HIV": 0.245,
    "chronic_HIV": 0.235,
}

# Fixed inputs per condition (based on final_calibrated_model.py)
COHERENCE = {
    "healthy": 0.85,
    "acute_HIV": 0.84,
    "chronic_HIV": 0.73,
}

SIGMA_R = {
    "healthy": CONST.sigma_r_regular,
    "acute_HIV": CONST.sigma_r_regular * 1.05,
    "chronic_HIV": CONST.sigma_r_regular * 1.4,
}

# Membrane turnover proxies used by choline_dynamics in the existing model
MEMBRANE_RATES = {
    "healthy": (1.0, 1.0),  # damage, repair
    "acute_HIV": (2.0, 2.5),
    "chronic_HIV": (1.2, 1.1),
}


# --------------------------------------------------------------------------------------
# Forward models
# --------------------------------------------------------------------------------------

def forward_NAA(coherence: pt.TensorVariable,
                xi: pt.TensorVariable,
                sigma_r: pt.TensorVariable,
                coh_exp: pt.TensorVariable,
                xi_exp: pt.TensorVariable,
                deloc_exp: pt.TensorVariable,
                NAA_base: pt.TensorVariable) -> pt.TensorVariable:
    """
    Mechanism-preserving NAA/Cr forward function that mirrors the calibrated form
    but with learnable exponents and base.

    NAA = NAA_base * (coherence/0.85)^coh_exp * (xi_baseline/xi)^xi_exp * (sigma_r/sigma_r0)^deloc_exp
    """
    sigma_r0 = CONST.sigma_r_regular
    xi0 = CONST.xi_baseline
    return (
        NAA_base
        * (coherence / 0.85) ** coh_exp
        * (xi0 / xi) ** xi_exp
        * (sigma_r / sigma_r0) ** deloc_exp
    )


def forward_Cho(membrane_damage: pt.TensorVariable,
                membrane_repair: pt.TensorVariable,
                cho_baseline: pt.TensorVariable,
                k_turnover: pt.TensorVariable) -> pt.TensorVariable:
    """
    Cho/Cr forward function consistent with choline_dynamics but with learnable
    turnover sensitivity coefficient k_turnover (0.1 in the deterministic model).

    Cho = Cho_baseline * (1 + k_turnover * ( (damage + repair) - 1 ))
    """
    net_turnover = membrane_damage + membrane_repair
    return cho_baseline * (1.0 + k_turnover * (net_turnover - 1.0))


# --------------------------------------------------------------------------------------
# PyMC model
# --------------------------------------------------------------------------------------

def build_model():
    with pm.Model() as model:
        # Priors for NAA coupling exponents
        coh_exp = pm.Normal("coh_exp", mu=3.0, sigma=0.5)
        xi_exp = pm.Normal("xi_exp", mu=0.5, sigma=0.2)
        deloc_exp = pm.Normal("deloc_exp", mu=0.15, sigma=0.1)
        NAA_base = pm.Normal("NAA_base", mu=1.10, sigma=0.1)

        # Priors for ξ per condition (in meters), weakly informative but centered per hypothesis
        xi_healthy = pm.TruncatedNormal("xi_healthy", mu=0.8e-9, sigma=0.05e-9, lower=0.3e-9, upper=1.2e-9)
        xi_acute = pm.TruncatedNormal("xi_acute", mu=0.4e-9, sigma=0.06e-9, lower=0.2e-9, upper=1.0e-9)
        xi_chronic = pm.TruncatedNormal("xi_chronic", mu=0.8e-9, sigma=0.06e-9, lower=0.3e-9, upper=1.2e-9)

        # Cho coupling coefficient (turnover sensitivity); positive
        k_turnover = pm.TruncatedNormal("k_turnover", mu=0.10, sigma=0.05, lower=0.0, upper=0.5)

        # Observation noise (shared for simplicity; could be per-modality)
        sigma_NAA = pm.HalfNormal("sigma_NAA", sigma=0.05)
        sigma_Cho = pm.HalfNormal("sigma_Cho", sigma=0.02)

        # Deterministic predictions per condition
        xi_map = {
            "healthy": xi_healthy,
            "acute_HIV": xi_acute,
            "chronic_HIV": xi_chronic,
        }

        # NAA predictions
        naa_preds = []
        cho_preds = []
        naa_obs = []
        cho_obs = []

        for cond in CONDITIONS:
            coh = COHERENCE[cond]
            sig = SIGMA_R[cond]
            xi_c = xi_map[cond]

            # NAA
            naa_pred = forward_NAA(
                coherence=pt.as_tensor_variable(coh),
                xi=xi_c,
                sigma_r=pt.as_tensor_variable(sig),
                coh_exp=coh_exp,
                xi_exp=xi_exp,
                deloc_exp=deloc_exp,
                NAA_base=NAA_base,
            )
            naa_preds.append(naa_pred)
            naa_obs.append(OBS_NAA[cond])

            # Cho using membrane rates proxy
            dmg, rep = MEMBRANE_RATES[cond]
            cho_pred = forward_Cho(
                membrane_damage=pt.as_tensor_variable(dmg),
                membrane_repair=pt.as_tensor_variable(rep),
                cho_baseline=pt.as_tensor_variable(CONST.Cho_baseline),
                k_turnover=k_turnover,
            )
            cho_preds.append(cho_pred)
            cho_obs.append(OBS_CHO[cond])

        # Stack
        naa_preds = pt.stack(naa_preds)
        cho_preds = pt.stack(cho_preds)
        naa_obs = np.array(naa_obs, dtype=float)
        cho_obs = np.array(cho_obs, dtype=float)

        # Likelihoods (vectorized over 3 conditions)
        pm.Normal("NAA_obs", mu=naa_preds, sigma=sigma_NAA, observed=naa_obs)
        pm.Normal("Cho_obs", mu=cho_preds, sigma=sigma_Cho, observed=cho_obs)

        # Mechanistic metric
        pm.Deterministic("delta_xi", xi_chronic - xi_acute)
        pm.Deterministic("xi_healthy_nm", xi_healthy * 1e9)
        pm.Deterministic("xi_acute_nm", xi_acute * 1e9)
        pm.Deterministic("xi_chronic_nm", xi_chronic * 1e9)

    return model


# --------------------------------------------------------------------------------------
# Run inference
# --------------------------------------------------------------------------------------

def run_inference(draws: int = 2000, tune: int = 1000, chains: int = 4, target_accept: float = 0.9, seed: int | None = 42):
    os.makedirs("results/bayesian", exist_ok=True)

    model = build_model()
    with model:
        idata = pm.sample(
            draws=draws,
            tune=tune,
            chains=chains,
            target_accept=target_accept,
            random_seed=seed,
            return_inferencedata=True,
            progressbar=True,
        )

    # Save trace
    trace_path = "results/bayesian/trace.nc"
    idata.to_netcdf(trace_path)

    # Summary
    summary = az.summary(idata, var_names=["coh_exp", "xi_exp", "deloc_exp", "NAA_base", "k_turnover", "sigma_NAA", "sigma_Cho", "xi_healthy_nm", "xi_acute_nm", "xi_chronic_nm", "delta_xi"], round_to=4)
    summary_path = "results/bayesian/summary.csv"
    summary.to_csv(summary_path)

    # Compute probability P(xi_acute < xi_chronic)
    xi_acute_vals = idata.posterior["xi_acute_nm"].values.reshape(-1)
    xi_chronic_vals = idata.posterior["xi_chronic_nm"].values.reshape(-1)
    p_order = float(np.mean(xi_acute_vals < xi_chronic_vals))

    # Posterior predictive on means for reporting (deterministic re-evaluation not needed; we can use posterior means of preds)
    # Here we re-compute point predictions using posterior medians for readability
    post = idata.posterior
    med = {k: float(np.median(post[k].values)) for k in ["coh_exp", "xi_exp", "deloc_exp", "NAA_base", "k_turnover"]}

    # Point predictions using medians
    preds = []
    for cond in CONDITIONS:
        xi_med = float(np.median(post[f"xi_{'acute' if cond=='acute_HIV' else ('chronic' if cond=='chronic_HIV' else 'healthy')}_nm"].values)) * 1e-9
        naa = (
            med["NAA_base"]
            * (COHERENCE[cond] / 0.85) ** med["coh_exp"]
            * (CONST.xi_baseline / xi_med) ** med["xi_exp"]
            * (SIGMA_R[cond] / CONST.sigma_r_regular) ** med["deloc_exp"]
        )
        dmg, rep = MEMBRANE_RATES[cond]
        cho = CONST.Cho_baseline * (1.0 + med["k_turnover"] * ((dmg + rep) - 1.0))
        preds.append({"condition": cond, "NAA_pred": naa, "Cho_pred": cho, "NAA_obs": OBS_NAA[cond], "Cho_obs": OBS_CHO[cond]})

    preds_df = pd.DataFrame(preds)
    preds_path = "results/bayesian/posterior_predictive.csv"
    preds_df.to_csv(preds_path, index=False)

    print("\nBayesian inference complete.")
    print(f"Trace saved to: {trace_path}")
    print(f"Summary saved to: {summary_path}")
    print(f"Posterior predictive saved to: {preds_path}")
    print(f"P(xi_acute < xi_chronic) = {p_order:.3f}")

    # Also return values for programmatic use
    return {
        "idata": idata,
        "summary_path": summary_path,
        "p_xi_acute_lt_chronic": p_order,
        "preds": preds_df,
    }


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Bayesian inference of coupling parameters with PyMC v4")
    parser.add_argument("--draws", type=int, default=2000)
    parser.add_argument("--tune", type=int, default=1000)
    parser.add_argument("--chains", type=int, default=4)
    parser.add_argument("--target-accept", type=float, default=0.9)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    run_inference(draws=args.draws, tune=args.tune, chains=args.chains, target_accept=args.target_accept, seed=args.seed)


if __name__ == "__main__":
    main()
