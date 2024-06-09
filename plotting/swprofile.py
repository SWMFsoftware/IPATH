import math
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
plt.switch_backend('agg')
density_n = []
sw_spd = []
xdims = 1500
x1min = 0.05
x1max = 2
x1rat = 1.0
vo=52.483


x1a = [0.0]*xdims

if x1rat == 1.0:
    for i in xrange(0, xdims):
          x1a[i] = x1min+ i * (x1max-x1min)/(xdims-1)
else:
     x1a[0] = x1min
     deltax = (x1max - x1min) * (1 - x1rat) / (1 - x1rat**(xdims-2))
     x1a[1] = x1a[0] + deltax
     for i in xrange(2, xdims):
          x1a[i] = x1min + deltax * (1 - x1rat**(i-1)) / (1 - x1rat)
density = []
f1 = open('../CME/path_output/solar_wind_profile.dat', 'r')
for line in f1:
       line = line.strip()
       columns = line.split()
       density_n.append(float(columns[1]))
       sw_spd.append(float(columns[2]))
f1.close
for j in xrange(0,xdims):
       density_n[j] = density_n[j]*x1a[j]**2.
       sw_spd[j]=sw_spd[j]*52.483
f2 = open('../CME/path_output/shock_posn_comp.dat', 'r')

shock_loc     = []
for line in f2:
       line = line.strip()
       columns = line.split()
       
       shock_loc.append(float(columns[0]))

f2.close

x=np.linspace(0,2,xdims)

fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(16,12))
ax0.set_title('solar wind properties vs radial distance',fontsize=20)
ax0.plot(x,density_n)
ax0.set_ylabel('N ${R}^{2}$$({cm}^{-3}$ ${AU}^{2})$',fontsize=15)
ax1.plot(x,sw_spd)
ax1.set_ylabel('V_SW (km/s)',fontsize=15)
ax1.set_xlabel('R(AU)',fontsize=15)
plt.savefig('SW-profile.png')
