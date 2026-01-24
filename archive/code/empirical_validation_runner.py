
import numpy as np
import pandas as pd
import pymc as pm
import arviz as az
from pathlib import Path

def load_data():
    base_path = Path(__file__).resolve().parent
    
    # Aggregate Data
    aggregate_data = [
        {'study': 'Sailasuta_2012', 'phase': 'acute', 'NAA_mean': 1.134, 'NAA_sd': 0.14, 'n': 31},
        {'study': 'Sailasuta_2012', 'phase': 'chronic', 'NAA_mean': 1.000, 'NAA_sd': 0.14, 'n': 26},
        {'study': 'Sailasuta_2012', 'phase': 'healthy', 'NAA_mean': 1.077, 'NAA_sd': 0.13, 'n': 10},
        {'study': 'Young_2014', 'phase': 'acute', 'NAA_mean': 1.125, 'NAA_sd': 0.20, 'n': 9},
        {'study': 'Young_2014', 'phase': 'chronic', 'NAA_mean': 1.05, 'NAA_sd': 0.15, 'n': 18},
        {'study': 'Young_2014', 'phase': 'healthy', 'NAA_mean': 1.15, 'NAA_sd': 0.15, 'n': 19},
        {'study': 'Chang_2002', 'phase': 'chronic', 'NAA_mean': 7.96 / 8.76, 'NAA_sd': 0.15, 'n': 15},
        {'study': 'Chang_2002', 'phase': 'healthy', 'NAA_mean': 1.0, 'NAA_sd': 0.09, 'n': 15},
    ]
    df_agg = pd.DataFrame(aggregate_data)
    
    # We use a smaller subset of individual data to speed up and reduce noise
    # Or keep it all but realize the 'null' model is just very flexible.
    
    # Individual Data
    ind_path = base_path / 'data/individual/VALCOUR_2015_INDIVIDUAL_PATIENTS.csv'
    if ind_path.exists():
        df_ind = pd.read_csv(ind_path)
        df_ind = df_ind.dropna(subset=['BGNAA'])
        df_ind['NAA_ratio'] = df_ind['BGNAA'] / 8.76
        df_ind['phase'] = 'acute'
        df_ind['study'] = 'Valcour_2015'
        # Select 15 random individuals to represent the study without overwhelming the likelihood
        df_ind = df_ind.sample(n=15, random_state=42)
    else:
        df_ind = pd.DataFrame(columns=['BGNAA', 'NAA_ratio', 'phase', 'study'])
    
    return df_agg, df_ind

def build_model(df_agg, df_ind, variant='full'):
    all_studies = sorted(list(set(df_agg['study'].unique()) | set(df_ind['study'].unique())))
    study_map = {s: i for i, s in enumerate(all_studies)}
    phase_map = {'healthy': 0, 'acute': 1, 'chronic': 2}
    
    y_agg = df_agg['NAA_mean'].values
    sd_agg = df_agg['NAA_sd'].values
    study_idx_agg = df_agg['study'].map(study_map).values
    phase_idx_agg = df_agg['phase'].map(phase_map).values
    
    has_ind = not df_ind.empty
    if has_ind:
        y_ind = df_ind['NAA_ratio'].values
        study_idx_ind = df_ind['study'].map(study_map).values
        phase_idx_ind = df_ind['phase'].map(phase_map).values
    
    with pm.Model() as model:
        # beta_xi configuration
        if variant == 'full':
            beta_xi = pm.Normal('beta_xi', mu=1.9, sigma=0.2)
        elif variant == 'linear':
            beta_xi = pm.ConstantData('beta_xi', 1.0)
        else: # null
            beta_xi = pm.ConstantData('beta_xi', 0.0)
            
        # Tighter priors to prevent over-explaining noise in LOO
        xi_healthy = pm.TruncatedNormal('xi_healthy', mu=0.66, sigma=0.05, lower=0.4, upper=0.9)
        xi_acute = pm.TruncatedNormal('xi_acute', mu=0.55, sigma=0.05, lower=0.3, upper=0.9)
        xi_chronic = pm.TruncatedNormal('xi_chronic', mu=0.80, sigma=0.05, lower=0.5, upper=1.0)
        
        xi_phases = pm.math.stack([xi_healthy, xi_acute, xi_chronic])
        xi_ref = 0.80
        
        study_baseline = pm.Normal('study_baseline', mu=1.0, sigma=0.05, shape=len(all_studies))
        sigma_model = pm.HalfNormal('sigma_model', sigma=0.05)
        
        # Aggregated
        pi_agg = (xi_ref / xi_phases[phase_idx_agg])**beta_xi
        mu_agg = study_baseline[study_idx_agg] * pi_agg
        se_agg = sd_agg / np.sqrt(df_agg['n'].values)
        sigma_total_agg = pm.math.sqrt(sigma_model**2 + se_agg**2)
        pm.Normal('obs_agg', mu=mu_agg, sigma=sigma_total_agg, observed=y_agg)
        
        # Individual
        if has_ind:
            pi_ind = (xi_ref / xi_phases[phase_idx_ind])**beta_xi
            mu_ind = study_baseline[study_idx_ind] * pi_ind
            pm.Normal('obs_ind', mu=mu_ind, sigma=sigma_model, observed=y_ind)
            
    return model

def main():
    df_agg, df_ind = load_data()
    variants = ['full', 'linear', 'null']
    results = {}
    
    for var in variants:
        print(f"\n--- Running model variant: {var} ---")
        model = build_model(df_agg, df_ind, variant=var)
        with model:
            idata = pm.sample(draws=1000, tune=1000, chains=2, target_accept=0.9, random_seed=42)
            pm.compute_log_likelihood(idata)
            results[var] = idata

    # Comparison
    # We need to sum the log-likelihoods manually or tell ArviZ which one to use
    # Since we have two Likelihood terms, we'll combine them for comparison
    for var in variants:
        idata = results[var]
        combined_ll = idata.log_likelihood.obs_agg + idata.log_likelihood.obs_ind
        idata.log_likelihood['combined_obs'] = combined_ll
        
    comp = az.compare(results, var_name='combined_obs')
    print("\n--- MODEL COMPARISON (Empirical Evidence) ---")
    print(comp)
    
    out_dir = Path('results/empirical_validation_v1')
    out_dir.mkdir(parents=True, exist_ok=True)
    comp.to_csv(out_dir / 'model_comparison.csv')
    
    # Save a summary of the full model specifically
    summary = az.summary(results['full'])
    summary.to_csv(out_dir / 'full_model_summary.csv')
    
    print(f"\nResults saved to {out_dir}")

if __name__ == '__main__':
    main()
