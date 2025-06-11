import cantera as ct
import math
import numpy as np
import matplotlib.pyplot as plt

gas = ct.Solution('h2o2.yaml')
gas.TPX = 1200, 10*ct.one_atm, 'H2:1, N2:0.1'

s0 = gas.s
h0 = gas.h
p0 = gas.P
T0 = gas.T

mdot = 1
n_points = 200
data = np.zeros((n_points, 4))
