# in your notebook cell
from combustion_sweep import sweep_loxlh2
import matplotlib.pyplot as plt

data = sweep_loxlh2(phi_min=0.6, phi_max=5.0, n=200)
phis  = data['phi']
OFs   = data['OF']
Isps  = data['Isp']
Ve    = data['Ve']
T0    = data['T0']
gamma = data['gamma']

F = 8e6
mdot = F/Ve

At = ((mdot* np.sqrt(T0))/P
# now plot or analyze as you like
'''plt.plot(OFs, Isps)
plt.xlabel("O/F ratio")
plt.ylabel("Isp (s)")
plt.show()'''
