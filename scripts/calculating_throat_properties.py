from frozen_combustion_sweep import sweep_loxlh2
import numpy as np
import matplotlib.pyplot as plt

data = sweep_loxlh2(phi_min=0.6, phi_max=5.0, n=400)
phis  = data['phi']
OFs   = data['OF']
Isp  = data['Isp']
Ve    = data['Ve']
T0    = data['T0']
gamma = data['gamma']
Te = data['Te']
ae = data['ae']
Mach = data['Mach']

max_isp_location = np.argmax(Isp)
corresponding_Isp = Isp[max_isp_location]
corresponding_ve = Ve[max_isp_location]
corresponding_phi = phis[max_isp_location]
corresponding_Of = OFs[max_isp_location]
corresponding_T0 = T0[max_isp_location]
corresponding_gamma = gamma[max_isp_location]
corresponding_Te = Te[max_isp_location]
corresponding_ae = ae[max_isp_location]
corresponding_Mach = Mach[max_isp_location]

print("Isp:", corresponding_Isp, "[s]\n" 
      "Ve:", corresponding_ve, '[m/s]\n' 
      'Phi:', corresponding_phi, '\n' 
      'OF ratio_actual:', corresponding_Of, '\n' 
      'Chamber Temperature:', corresponding_T0, '[K]\n' 
      'Gamma:', corresponding_gamma, '\n' 
      'Exit Temperature:', corresponding_Te, '[K]\n' 
      'Speed of Sound at the Exit:', corresponding_ae, '[m/s]\n' 
      'Mach Number at the exit:', corresponding_Mach)