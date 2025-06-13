# in your notebook cell
from frozen_combustion_sweep import sweep_loxlh2
import matplotlib.pyplot as plt

data = sweep_loxlh2(phi_min=0.6, phi_max=5.0, n=200)
phis  = data['phi']
OFs   = data['OF']
Isp  = data['Isp']
Ve    = data['Ve']
T0    = data['T0']
gamma = data['gamma']
Te = data['Te']
ae = data['ae']
Mach = data['Mach']

fig, ax = plt.subplots()
ax.grid()
h1 = ax.plot(OFs, Isp)
ax.set(xlabel='OF Ratio', ylabel='Isp [s]')
plt.show()
