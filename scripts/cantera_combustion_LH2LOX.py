import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

gas = ct.Solution('gri30.yaml')
n = 200
fuel = 'H2:2'
ox = 'O2:1'
P = 7e6
Patm = 101325
T0 = np.zeros(n)
phis = np.zeros(n)
gamma = np.zeros(n)
Rspec = np.zeros(n)
OFact = np.zeros(n)
ve = np.zeros(n)
Isp = np.zeros(n)
g0 = 9.807
F = 8e6

for i, phi in enumerate((np.linspace(0.6, 5, n))):
    gas.set_equivalence_ratio(phi = phi, fuel=fuel, oxidizer=ox)
    gas.TP = 298.15, P
    gas.equilibrate('HP')
    cp = gas.cp_mass
    cv = gas.cv_mass
    gamma_i = cp / cv
    gamma[i] = gamma_i
    T0_i = gas.T
    T0[i] = T0_i
    phis[i] = phi
    R_i = cp - cv
    Rspec[i] = R_i
    exp_term = (Patm / P)**((gamma_i - 1.0) / gamma_i)
    ve[i] = np.sqrt( (2 * gamma_i) / (gamma_i - 1.0) * R_i * T0_i * (1.0 - exp_term) )
    Isp[i] = ve[i]/g0
    

OFsto = gas.stoich_air_fuel_ratio('H2:2', 'O2:1')
OFact = OFsto/phis
mdot = F/ve
At = (((mdot * np.sqrt(T0))/P)*np.sqrt(Rspec/gamma)*((gamma+1)/2)**((gamma+1)/(2*(gamma-1))))
print('Mass flow rate')
print(mdot)
print('Throat Area')
print(At)

fig1, ax = plt.subplots()
ax.grid()
h1 = ax.plot(OFact, Isp)
ax.set(xlabel='OF Ratio', ylabel='Isp [s]')
plt.show()

fig2, ax2 = plt.subplots()
ax2.grid()
h1 = ax2.plot(At, mdot)
ax2.set(xlabel='Throat Area', ylabel='Mass Flow Rate')
plt.show()

