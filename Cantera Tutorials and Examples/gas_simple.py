import cantera as ct

gas = ct.Solution('gri30.yaml')
major = gas['CH4', 'O2', 'CO2','H2O','N2']
cp_major = major.partial_molar_cp
wdot_major = major.net_production_rates

gas()
