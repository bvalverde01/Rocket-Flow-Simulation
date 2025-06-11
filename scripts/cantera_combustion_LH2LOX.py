import cantera as ct

gas = ct.Solution('gri30.yaml')

fuel = 'H2:2'
ox = 'O2:1'

gas.set_equivalence_ratio(phi = 1.0, fuel=fuel, oxidizer=ox)
print(gas.stoich_air_fuel_ratio('H2:2', 'O2:1'))