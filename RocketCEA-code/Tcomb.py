from rocketcea.cea_obj_w_units import CEA_Obj
from pylab import *

Pc = 0
MR = 0
MRArr = []
TArr = []

TObj = CEA_Obj(oxName='LOX', fuelName='LH2', temperature_units='K', pressure_units='MPa' )
while True:
    MR += 0.01
    Pc += 1/145
    MRArr.append(MR)
    T = TObj.get_Tcomb(Pc, MR)
    TArr.append(T)
    plot(MRArr, TArr)
    if MR >= 15:
        break

legend(loc='best')
grid(True)
title('Combustion Temperature vs. O/F Ratio')
xlabel('O/F Ratio')
ylabel('Temperature [K]')

show()