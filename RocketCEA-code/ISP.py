from rocketcea.cea_obj_w_units import CEA_Obj
import matplotlib.pyplot as plt
import numpy as np

IObj = CEA_Obj(oxName='LOX', fuelName="LH2", pressure_units='MPa')
figure, ax = plt.subplots()
ax.set(ylabel='Specific Impulse [s]')
ax.set(xlabel='O/F Ratio')
ax.set_title('Specific Impulse vs O/F Ratio at different Chamber Pressures')
ax.grid()

for P in [0.101325, 1, 2, 3, 4, 5, 6, 7 , 8, 9, 10]:
    MR = 1 
    MRArr = []
    IArr = []
    
     # Reset MR for each new value of P
    while MR <= 15:
        MRArr.append(MR)
        I = IObj.get_Isp(Pc=P, MR=MR)
        IArr.append(I)
        MR += 0.1

    max_I = max(IArr)
    max_I_index = IArr.index(max_I)
    max_I = round(max_I, 2)
    corresponding_MR = MRArr[max_I_index]
    corresponding_MR = round(corresponding_MR, 2)
    print(P, max_I, corresponding_MR)
    
    ax.plot(MRArr, IArr, label=f'{P} MPa')

ax.legend()

plt.show()























'''PArr = []
MRArr = []
IPArr = []
IMRArr = []
figure, ax = plt.subplots()
ax.grid()

colorr = 'tab:red'
ax.set_xlabel('Chamber Pressure [MPa]', color=colorr)
ax.set_ylabel('Specific Impluse [m/s]')
ax.tick_params(axis='x', labelcolor=colorr)

ax2 = ax.twiny()

colorb = 'tab:blue'
ax2.set_xlabel('O/F Ratio', color=colorb)
ax2.tick_params(axis='x', labelcolor=colorb)

ispObj = CEA_Obj(oxName='LOX', fuelName="LH2", isp_units='m/s')

while True:
    MR += 0.01
    Pc += 1/145
    MRArr.append(MR)
    PArr.append(Pc)
    Ip = np.array(ispObj.get_Isp(Pc))
    IM = np.array(ispObj.get_Isp(MR))
    IPArr.append(Ip)
    IMRArr.append(IM)
    ax2.plot(MRArr, IMRArr, color=colorb)
    ax.plot(PArr, IPArr, color=colorr)
    if MR >= 15:
        break

plt.show()'''