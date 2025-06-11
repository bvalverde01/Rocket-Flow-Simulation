import cantera as ct

# Load the mechanism
gas = ct.Solution("nDodecane_Reitz.yaml")

# Inputs
OF = 2.27                        # oxidizer-to-fuel mass ratio
T0 = 298.15                      # K
P0 = 7e6                         # Pa
fuel = 'C12H26'
oxidizer = 'O2'

# Get molecular weights
MW_fuel = gas.molecular_weights[gas.species_index(fuel)]
MW_ox = gas.molecular_weights[gas.species_index(oxidizer)]

# Calculate mole ratio from mass ratio
# OF = (n_ox * MW_ox) / (n_fuel * MW_fuel)
# let n_fuel = 1 mol
n_fuel = 1.0
n_ox = (OF * MW_fuel) / MW_ox

# Set the gas composition
gas.TPX = T0, P0, {fuel: n_fuel, oxidizer: n_ox}

# Equilibrate at constant pressure and enthalpy (adiabatic)
gas.equilibrate('HP')

# Output
print(f"Adiabatic flame temperature: {gas.T:.2f} K")
print(f"Pressure: {gas.P:.2f} Pa")
print(f"Density: {gas.density:.3f} kg/m³")
print(f"gamma: {gas.cp_mass/gas.cv_mass:.4f}")
print(f"Molecular weight: {gas.mean_molecular_weight:.3f} kg/kmol")
print(f"Speed of sound: {gas.sound_speed:.2f} m/s")

# Print top species by mole fraction
print("\nMajor species:")
for sp, x in sorted(zip(gas.species_names, gas.X), key=lambda p: -p[1]):
    if x > 5e-5:
        print(f"{sp:>10s}: {x:.5f}")
