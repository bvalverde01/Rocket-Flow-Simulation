import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from atmosphere import grid_profile

gas = ct.Solution('gri30.yaml')
n = 200
nphis = 200
natms = 200
atm = grid_profile(n=natms, h_min=0, h_max=84000)
Patm = atm['P']
fuel = 'H2:2'
ox = 'O2:1'
P0 = 7e6
T0 = np.zeros(n)
phis = np.linspace(0.6, 5, nphis)
gamma = np.zeros(n)
Rspec = np.zeros(n)
OFact = np.zeros(n)
ve = np.zeros(n)    
Isp = np.zeros((nphis, natms))
g0 = 9.807
F = 8e6

for i, phi in enumerate(phis):
    gas.set_equivalence_ratio(phi = phi, fuel=fuel, oxidizer=ox)
    gas.TP = 298.15, P0
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
    for j, P_back in enumerate(Patm):
        exp_term = (P_back / P0)**((gamma_i - 1.0) / gamma_i)
        ve[i] = np.sqrt( (2 * gamma_i) / (gamma_i - 1.0) * R_i * T0_i * (1.0 - exp_term) )
        Isp[i, j] = ve[i] / g0

    

OFsto = gas.stoich_air_fuel_ratio('H2:2', 'O2:1')
OFact = OFsto/phis


import matplotlib.ticker as ticker

# Compute average Isp across pressures for each OF ratio
mean_Isp_per_OF = np.mean(Isp, axis=1)
best_i = np.argmax(mean_Isp_per_OF)
best_OF = OFact[best_i]
best_avg_isp = mean_Isp_per_OF[best_i]

# Create contour plot with OFact instead of phi
fig, ax = plt.subplots(figsize=(8, 6))
contour = ax.contourf(Patm, OFact, Isp, levels=50, cmap='viridis')    
cbar = plt.colorbar(contour, ax=ax)
cbar.set_label("Isp (s)")

# Add contour lines
lines = ax.contour(Patm, OFact, Isp, levels=np.arange(300, 550, 25), colors='k', linewidths=0.5)
ax.clabel(lines, inline=1, fontsize=8)

# Highlight best average-Isp OF
ax.axhline(best_OF, color='red', linestyle='--', linewidth=1.5,
           label=f"Best O/F ≈ {best_OF:.2f}")
ax.text(Patm[0], best_OF + 0.2,
        f"Best avg Isp ≈ {best_avg_isp:.1f} s",
        color='red', fontsize=9)

# Axis formatting
ax.set_xlabel("Atmospheric Pressure (Pa)")
ax.set_ylabel("Mixture Ratio (O/F actual)")
ax.set_title("Isp Contour Plot with Optimal O/F Highlighted")
ax.invert_xaxis()
ax.grid(True)
ax.legend()

# Format tick labels for scientific notation if needed
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x:.0f}'))

plt.tight_layout()
plt.show()



