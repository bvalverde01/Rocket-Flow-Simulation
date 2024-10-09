from rocketcea.cea_obj_w_units import CEA_Obj
import matplotlib.pyplot as plt
import numpy as np

Pc = 0
MR = 0
PArr = []
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

plt.show()