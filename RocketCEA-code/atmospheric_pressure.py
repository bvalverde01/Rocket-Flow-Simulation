import math
from pylab import *

P_sealevel = 1.225
P = 0
h = 0
g = 9.81
M = 0.0289652
R = 8.31446
T = 288.15
L = 0.0065
hArr = []
PArr = []

while True:
    h += 1000
    P = P_sealevel * (math.e**((((g*M)/(R*L))-1)*math.log(1-((L*h)/T))))
    hArr.append(h)
    PArr.append(P)
    plot(hArr, PArr)
    if h == 40000:
        print(h, P)
        break

legend(loc='best')
grid(True)
title('Atmospheric Density vs Height')
xlabel('Height [m]')
ylabel('Density [kg/m^3]')

show()


