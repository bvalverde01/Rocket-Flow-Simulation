import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ──────────────────────────────
# 1) Constants
g0 = 9.80665             # m/s² at sea level
RE = 6.371e6             # Earth radius (m)
M  = 0.0289644           # molar mass of air (kg/mol)
R0 = 8.314462618         # universal gas constant (J/mol·K)
P0 = 101325              # Pa at sea level

# ──────────────────────────────
# 2) Temperature profile (U.S. Standard Atmosphere 1976, up to 86 km)
def T(h):
    if   h < 11000:  return 288.15 - 0.0065 * h
    elif h < 20000:  return 216.65
    elif h < 32000:  return 216.65 + 0.0010 * (h - 20000)
    elif h < 47000:  return 228.65 + 0.0028 * (h - 32000)
    elif h < 51000:  return 270.65
    elif h < 71000:  return 270.65 - 0.0028 * (h - 51000)
    elif h <= 86000: return 214.65 - 0.0020 * (h - 71000)
    else:            return np.nan  # out of range

# ──────────────────────────────
# 3) Gravity vs. altitude
def g(h):
    return g0 * (RE / (RE + h))**2

# ──────────────────────────────
# 4) Pressure ODE: dP/dh = - P·(M·g(h)) / (R0·T(h))
def dPdh(h, P):
    return -P[0] * (M * g(h)) / (R0 * T(h))

# ──────────────────────────────
# 5) Build altitude arrays
h_temp     = np.linspace(0,  84000, 200)   # for Temp & g plots
h_pressure = np.linspace(0, 100000, 500)   # for pressure integration

# 6) Evaluate T and g for the first two plots
Temp_vals = np.array([T(h) for h in h_temp])
g_vals    = np.array([g(h) for h in h_temp])

# 7) Numerically integrate pressure from 0 to 100 km
sol = solve_ivp(
    dPdh,
    t_span=(0, 100000),
    y0=[P0],
    t_eval=h_pressure,
    rtol=1e-6,
    atol=1e-6
)
pressure = sol.y[0]  # Pa

# ──────────────────────────────
# 8) Plot Temperature vs. Altitude
fig1, ax1 = plt.subplots()
ax1.plot(h_temp/1000, Temp_vals, '-')
ax1.set(
    title='Atmospheric Temperature vs. Altitude',
    xlabel='Altitude (km)',
    ylabel='Temperature (K)'
)
ax1.grid()

# 9) Plot Gravity vs. Altitude (log-x)
fig2, ax2 = plt.subplots()
ax2.semilogx(h_temp, g_vals, '-')
ax2.set(
    title='Gravitational Acceleration vs. Altitude',
    xlabel='Altitude (m)',
    ylabel='g (m/s²)'
)
ax2.grid()

# 10) Plot Pressure vs. Altitude (log-y)
fig3, ax3 = plt.subplots()
ax3.semilogy(sol.t/1000, pressure, '-',)
ax3.set(
    title='Atmospheric Pressure vs. Altitude',
    xlabel='Altitude (km)',
    ylabel='Pressure (Pa)'
)
ax3.grid()

plt.tight_layout()
plt.show()
