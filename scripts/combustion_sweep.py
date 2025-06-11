# scripts/combustion_sweep.py
import cantera as ct
import numpy as np

def sweep_loxlh2(phi_min=0.6, phi_max=5.0, n=200,
                 P_chamber=7e6, P_exit=101325, g0=9.80665):
    gas = ct.Solution('gri30.yaml')
    phis   = np.linspace(phi_min, phi_max, n)
    T0     = np.zeros(n)
    gamma  = np.zeros(n)
    Rspec  = np.zeros(n)
    ve     = np.zeros(n)
    Isp    = np.zeros(n)

    # compute stoich O/F once:
    OFsto = gas.stoich_air_fuel_ratio('H2:2', 'O2:1')

    for i, phi in enumerate(phis):
        gas.set_equivalence_ratio(phi, 'H2', 'O2')
        gas.TP = 298.15, P_chamber
        gas.equilibrate('HP')

        cp = gas.cp_mass; cv = gas.cv_mass
        gamma[i]  = cp / cv
        T0[i]     = gas.T
        Rspec[i]  = cp - cv

        exp_term = (P_exit / P_chamber)**((gamma[i] - 1.0) / gamma[i])
        ve[i]     = np.sqrt((2 * gamma[i])/(gamma[i] - 1)*Rspec[i]*T0[i]*(1-exp_term))
        Isp[i]    = ve[i] / g0

    # compute actual O/F
    OFact = OFsto / phis

    return {
        'phi': phis,
        'OF': OFact,
        'T0': T0,
        'gamma': gamma,
        'Rspec': Rspec,
        'Ve': ve,
        'Isp': Isp
    }

if __name__ == '__main__':
    # simple test
    data = sweep_loxlh2()
    print("Peak Isp:", data['Isp'].max(), "at O/F =", data['OF'][data['Isp'].argmax()])
