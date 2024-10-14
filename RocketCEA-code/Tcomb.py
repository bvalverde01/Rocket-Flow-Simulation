from rocketcea.cea_obj_w_units import CEA_Obj
import matplotlib.pyplot as plt
import numpy as np

figure, ax = plt.subplots()
ax.set(ylabel='Temperature [K]')
ax.set(xlabel='O/F Ratio')
ax.set_title('Combustion Temperature vs O/F Ratio at different Chamber Pressures')
ax.grid()

TObj = CEA_Obj(oxName='LOX', fuelName='LH2', temperature_units='K', pressure_units='MPa' )


for P in [0.101325, 1, 2, 3, 4, 5, 6, 7 , 8, 9, 10]:
    MR = 0 
    MRArr = []
    TArr = []
    
     # Reset MR for each new value of P
    while MR <= 15:
        MRArr.append(MR)
        T = TObj.get_Tcomb(P, MR)
        TArr.append(T)
        print(P, MR, T)
        MR += 0.5
    

    ax.plot(MRArr, TArr, label=f'{P} MPa')

ax.legend()

plt.show()