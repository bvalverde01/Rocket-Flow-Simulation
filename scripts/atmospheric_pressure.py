import matplotlib.pyplot as plt
import numpy as np

n = 200
h_array = np.linspace(0, 84000, n)
g0 = 9.807
Re = 6378137
g = np.zeros(n)
Temp = np.zeros(n)

for i, h in enumerate(h_array):
    g[i] = g0*(Re/(Re+h))**2

def T(h):
    if h < 11000:
        return 288.15 - 0.0065 * h
    elif h < 20000:
        return 216.65
    elif h < 32000:
        return 216.65 + 0.0010 * (h - 20000)
    elif h < 47000:
        return 228.65 + 0.0028 * (h - 32000)
    elif h < 51000:
        return 270.65
    elif h < 71000:
        return 270.65 - 0.0028 * (h - 51000)
    elif h <= 86000:
        return 214.65 - 0.0020 * (h - 71000)
    else:
        raise ValueError("Altitude out of range for T(h): must be <= 86 km")

Temp = np.array([T(h) for h in h_array])


fig1, ax1 = plt.subplots()
ax1.grid()
ax1.plot(h_array, Temp)
ax1.set(title='Atmospheric Temperature vs Altitude', xlabel='Altitude [m]', ylabel='Atmospheric Temperature [K]')

fig2, ax2 = plt.subplots()
ax2.grid()
ax2.semilogx()
ax2.semilogy()
ax2.plot(h_array, g)
ax2.set(title='Gravitational Acceleration on Earth vs. Altitude', xlabel='Altitude [m]', ylabel='Gravitational Acceleration [m/s^2]')


plt.show()