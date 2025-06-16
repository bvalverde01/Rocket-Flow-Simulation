import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from atmosphere import grid_profile
from scipy.optimize import root_scalar

P0_array = np.arange(1e6, 7e6, 0.5e6)
phi = 1.6
fuel = 'H2:2'
ox = 'O2:1'
gas = ct.Solution('gri30.yaml')
n = 200
natm = 85
gamma = np.zeros(len(P0_array))
nepsilon = n
epsilons = np.linspace(10, 50, nepsilon)
atm = grid_profile(n=natm, h_min=0, h_max=84000)
Patm = atm['P']
Me = np.zeros((len(P0_array), nepsilon))
Pe = np.zeros((len(P0_array), nepsilon))
deltaP = np.zeros((len(P0_array), nepsilon, natm))
mean_deltaP = np.zeros((len(P0_array), nepsilon))

def area_ratio(Me, gamma):
    term = (2/(gamma+1)) * (1 + (gamma-1)/2 * Me**2)
    exponent = (gamma+1)/(2*(gamma-1))
    return (1/Me) * term**exponent

def solve_Me_from_eps(epsilon, gamma):
    # Residual function: f(Me) = area_ratio(Me, gamma) - epsilon
    def residual(Me):
        return area_ratio(Me, gamma) - epsilon

    # Solve using Brent's method with bracketed guess (Mach > 1)
    sol = root_scalar(residual, bracket=[1.01, 50.0], method='brentq', xtol=1e-6, rtol=1e-8)

    if not sol.converged:
        raise RuntimeError(f"Mach solve failed for ε={epsilon:.2f}, γ={gamma:.3f}")
    
    return sol.root

for i, P0 in enumerate(P0_array):
    gas.set_equivalence_ratio(phi = phi, fuel=fuel, oxidizer=ox)
    gas.TP = 288.15, P0
    gas.equilibrate('HP')
    cp = gas.cp_mass
    cv = gas.cv_mass
    gamma_i = cp / cv
    gamma[i] = gamma_i
    for j, epsilon in enumerate(epsilons):
        Me_ij = solve_Me_from_eps(epsilon, gamma_i)
        Me[i,j] = Me_ij
        Pe_ij= P0*(1 + ((gamma_i-1)/2)*(Me_ij)**2)**(-gamma_i/(gamma_i-1))
        Pe[i,j] = Pe_ij
        for k, Pb in enumerate(Patm):
            deltaP_ijk = np.abs(Pe_ij-Pb)
            deltaP[i,j,k] = deltaP_ijk


mean_deltaP = np.mean(deltaP, axis=2)

P0_labels = P0_array / 1e6

i_best, j_best = np.unravel_index(np.argmin(mean_deltaP), mean_deltaP.shape)
best_P0 = P0_array[i_best]
best_eps = epsilons[j_best]
best_deltaP = mean_deltaP[i_best, j_best]

print(f"Optimal Config: P₀ = {best_P0/1e6:.2f} MPa, ε = {best_eps:.2f}, Mean ΔP = {best_deltaP:.2f} Pa")

plt.figure(figsize=(10, 6))
c = plt.contourf(epsilons, P0_labels, mean_deltaP, levels=50, cmap='viridis')
plt.colorbar(c, label='Mean |Pe - Pa| [Pa]')
plt.xlabel('Expansion Ratio (ε)')
plt.ylabel('Chamber Pressure [MPa]')
plt.title('Mean Pressure Mismatch Across Altitudes')
plt.scatter(best_eps, best_P0 / 1e6, color='red', label='Best Match')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

plt.plot(np.linspace(0, 84, natm), deltaP[i_best, j_best, :])
plt.xlabel("Altitude (km)")
plt.ylabel("ΔP = |Pe - Pa| (Pa)")
plt.title("Pressure Mismatch vs Altitude at Best Configuration")
plt.grid(True)
plt.show()
