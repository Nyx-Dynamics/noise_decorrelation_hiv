import importlib

print('Testing enzyme_kinetics...')
ek = importlib.import_module('enzyme_kinetics')
print('compute_protection_factor:', ek.compute_protection_factor(0.8e-9, 0.5))
E = ek.EnzymeKinetics(Pi_xi=1.0, eta_coh=1.0)
print('integrate:', E.integrate(30, membrane_turnover=1.2))

print('Testing final_calibrated_model wrapper...')
fm = importlib.import_module('final_calibrated_model')
print('run_full_model healthy:', fm.run_full_model('healthy')['NAA_Cr'])

print('Testing final_calibrated_model_v2...')
fm2 = importlib.import_module('final_calibrated_model_v2')
print('run_full_model_v2 healthy:', fm2.run_full_model_v2('healthy')['NAA_Cr'])
