import cantera as ct

g = ct.Solution('gri30.yaml')
g.TPX = 300.0, ct.one_atm, 'CH4:0.95,O2:2, N2:7.52'
g.equilibrate('HP')

rf = g.forward_rates_of_progress
rr = g.reverse_rates_of_progress
for i in range(g.n_reactions):
    if g.reaction(i).reversible and rf[i] != 0.0:
        print(' %4i  %10.4g  ' % (i, (rf[i] - rr[i])/rf[i]))
