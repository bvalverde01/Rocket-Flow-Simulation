from rocketcea.cea_obj_w_units import CEA_Obj
import matplotlib.pyplot as plt
import numpy as np

IObj = CEA_Obj(oxName='LOX', fuelName="LH2", temperature_units='K', pressure_units='MPa', specific_heat_units='kJ/kg-K', isp_units='m/s')
figure, ax = plt.subplots()
ax.set(ylabel='Specific Impulse [s]')
ax.set(xlabel='O/F Ratio')
ax.set_title('Specific Impulse vs O/F Ratio at different Chamber Pressures')
ax.grid()

for P in [5, 6, 7 , 8, 9, 10]:
    MR = 1 
    MRArr = []
    IArr = []


    while MR <= 15: #I want to try the IvacCstrTc function from Rocketcea
        MRArr.append(MR)
        I = IObj.estimate_Ambient_Isp(Pc=P, MR=MR, Pamb=0.101325)
        A = I[0::2]
        IArr.append(A)
        MR += 1

    '''max_I = max(IArr)
    max_I_index = IArr.index(max_I)
    max_I = round(max_I, 2)
    corresponding_MR = MRArr[max_I_index]
    corresponding_MR = round(corresponding_MR, 2)
    print(P, max_I, corresponding_MR)'''
    
    ax.plot(MRArr, IArr, label=f'{P} MPa')

ax.legend()

plt.show()