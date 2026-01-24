import arviz as az
idata = az.from_netcdf('results/bayesian_v3_6/3_1_1/both/20251123T223826Z_a8c794cc/trace_v3_6.nc')
print('max R-hat:', az.rhat(idata).to_array().max().values)
print('min ESS bulk:', az.ess(idata, method='bulk').to_array().min().values)
print('min ESS tail:', az.ess(idata, method='tail').to_array().min().values)