import cantera as ct

g = ct.Solution('gri30.yaml')
r = g.reaction(2)
print(r.rate.input_data)
