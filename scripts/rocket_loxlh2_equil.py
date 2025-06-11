from rocketcea.cea_obj import CEA_Obj

# Define propellants
ispObj = CEA_Obj(oxName='LOX', fuelName='RP-1')

# Define test parameters
Pc = 70.0  # chamber pressure in bar
MR = 2.27  # mixture ratio
eps = 15   # expansion ratio

# Get mole fractions
molWtD, moleFracD = ispObj.get_SpeciesMoleFractions(Pc=Pc, MR=MR, eps=eps, frozen=1, frozenAtThroat=0, min_fraction=5e-06)

# Print mole fractions at chamber (index 1)
print("Major Species at Chamber (Mole Fractions > 1%):")
for species, fracs in moleFracD.items():
    if fracs[1] > 0.01:
        print(f"{species:>8}: {fracs[1]:.6f}")
