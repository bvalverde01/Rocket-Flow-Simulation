from rocketcea.cea_obj_w_units import CEA_Obj
import matplotlib.pyplot as plt

Pc = 0
MR = 0
PArr = []
MRArr = []
TArr = []
f, (ax1, ax2) = plt.subplots(1,2 , sharey=True)
ax1.set(xlabel='O/F Ratio', ylabel='Temperature [K]')
ax2.set(xlabel='Chamber Pressure')
ax1.grid()
ax2.grid()

TObj = CEA_Obj(oxName='LOX', fuelName='LH2', temperature_units='K', pressure_units='MPa' )
while True:
    MR += 0.01
    Pc += 1/145
    MRArr.append(MR)
    PArr.append(Pc)
    T = TObj.get_Tcomb(Pc, MR)
    TArr.append(T)
    ax1.plot(MRArr, TArr)
    ax2.plot(PArr, TArr)
    if MR >= 15:
        break

plt.show()