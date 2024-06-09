import math
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
plt.switch_backend('agg')

#min_e = 1e5
#max_e = 1e9
p_num = 20
phi_num = 25
#shell_num = 29

time_no  = 29
AU = 1.5e11
cme_center =88.
del_phi = 5.

pi = 3.14159265359


omega = 2.87e-6
r_bnd = 0.05*AU


t_o = 2858068.3
mp = 1.67e-27
co = 3e8
eo = 1.6e-19
vo = 52483.25

time_index = [5,9,13,17,21]


theta_bn_t    = []
shock_loc_t   = []
shk_spd_t     = []
comp_rat_t    = []
max_e_t       = []
inter_phi     = []

root_dir = '../CME/path_output/'

f1 = open(root_dir+'shock_posn_comp.dat', 'r')

shock_loc     = []
comp_rat      = []
shk_spd       = []
shk_begin     = []
shk_end       = []
theta_bn      = []
theta_vn      = []
max_e         = []
min_e         = []
real_spd      = []


line_no = 0
for line in f1:
       line = line.strip()
       columns = line.split()

       shock_loc.append(float(columns[0]))
       comp_rat.append(float(columns[1]))
       shk_spd.append(float(columns[2]))
       shk_begin.append(float(columns[3]))
       shk_end.append(float(columns[4]))
       theta_bn.append(float(columns[5]))
       theta_vn.append(float(columns[6]))
       
       real_spd.append(float(columns[2])*math.cos(float(columns[6])/180.*pi))
       line_no = line_no+1

f1.close
shell_num = line_no/phi_num



print shell_num

f2 = open(root_dir+'shock_momenta.dat', 'r')


time_all      =[]

for line in f2:
       line = line.strip()
       columns = line.split()
       time_all.append(float(columns[1]))
       min_e.append(float(columns[4]))
       max_e.append(float(columns[5])/1e6)       

f2.close

phi_0 = []
for i in xrange(0, phi_num):
       phi_0.append(int(cme_center +del_phi*(i-math.floor(phi_num/2.))))
print phi_0

foot_phi = [x for x in range(phi_num)]
# xtime = []



shell_time =[]
for i in xrange(0, shell_num):
       #current_time = (time_all[time_index[0]*phi_num] - 0.4)*t_o/3600
       shell_time.append('%(time).1F %(unit)s' %{"time":time_all[i*phi_num], "unit":"hrs"})


fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(2, 2, figsize=(16,12))

ax0.plot(phi_0, comp_rat[time_index[0]*phi_num:(time_index[0]+1)*phi_num],'b',  \
         phi_0, comp_rat[time_index[1]*phi_num:(time_index[1]+1)*phi_num],'r',  \
         phi_0, comp_rat[time_index[2]*phi_num:(time_index[2]+1)*phi_num],'g',  \
         phi_0, comp_rat[time_index[3]*phi_num:(time_index[3]+1)*phi_num],'y',  \
         phi_0, comp_rat[time_index[4]*phi_num:(time_index[4]+1)*phi_num],'k' )
ax0.set_ylabel('s', fontsize=18)
ax0.set_title('Compression ratio', fontsize=20)
ax0.set_xlim([33,143])
ax0.set_xlabel('longitude $\phi (^\circ)$', fontsize=15)
ax0.legend([shell_time[time_index[0]], shell_time[time_index[1]],\
            shell_time[time_index[2]],shell_time[time_index[3]],\
            shell_time[time_index[4]]], loc = 8)
ax0.tick_params(axis='both', which='major', labelsize=14)

ax1.plot(phi_0, shock_loc[time_index[0]*phi_num:(time_index[0]+1)*phi_num],'b',  \
         phi_0, shock_loc[time_index[1]*phi_num:(time_index[1]+1)*phi_num],'r',  \
         phi_0, shock_loc[time_index[2]*phi_num:(time_index[2]+1)*phi_num],'g',  \
         phi_0, shock_loc[time_index[3]*phi_num:(time_index[3]+1)*phi_num],'y',  \
         phi_0, shock_loc[time_index[4]*phi_num:(time_index[4]+1)*phi_num],'k')
ax1.set_ylabel('$R_{shk} (AU)$', fontsize=18)
ax1.set_title('Shock front location', fontsize=20)
ax1.set_xlim([33,143])
ax1.set_xlabel('longitude $\phi (^\circ)$', fontsize=15)
ax1.tick_params(axis='both', which='major', labelsize=14)



ax2.plot(phi_0, shk_spd[time_index[0]*phi_num:(time_index[0]+1)*phi_num],'b',  \
         phi_0, shk_spd[time_index[1]*phi_num:(time_index[1]+1)*phi_num],'r',  \
         phi_0, shk_spd[time_index[2]*phi_num:(time_index[2]+1)*phi_num],'g',  \
         phi_0, shk_spd[time_index[3]*phi_num:(time_index[3]+1)*phi_num],'y',  \
         phi_0, shk_spd[time_index[4]*phi_num:(time_index[4]+1)*phi_num],'k'  )
ax2.set_ylabel('$V_{shk} (km/s)$', fontsize=18)
ax2.set_title('shock speed', fontsize=20)
ax2.set_xlim([33,143])
ax2.set_xlabel('longitude $\phi (^\circ)$', fontsize=15)
ax2.tick_params(axis='both', which='major', labelsize=14)


ax3.plot(phi_0, theta_bn[time_index[0]*phi_num:(time_index[0]+1)*phi_num],'b',  \
         phi_0, theta_bn[time_index[1]*phi_num:(time_index[1]+1)*phi_num],'r',  \
         phi_0, theta_bn[time_index[2]*phi_num:(time_index[2]+1)*phi_num],'g',  \
         phi_0, theta_bn[time_index[3]*phi_num:(time_index[3]+1)*phi_num],'y',  \
         phi_0, theta_bn[time_index[4]*phi_num:(time_index[4]+1)*phi_num],'k')
ax3.set_ylabel('$\\theta_{BN} (^\circ)$', fontsize=18)
ax3.set_title('$\\theta_{BN}$', fontsize=20)
ax3.set_xlim([33,143])
ax3.set_xlabel('longitude $\phi (^\circ)$', fontsize=15)
ax3.tick_params(axis='both', which='major', labelsize=14)
ax3.legend([shell_time[time_index[0]], shell_time[time_index[1]],\
            shell_time[time_index[2]],shell_time[time_index[3]],\
            shell_time[time_index[4]]], loc = 2)
plt.savefig('Shock parameters.png')

'''
plt.figure(2,figsize=(16,12))
plt.plot(phi_0, real_spd[time_index[0]*phi_num:(time_index[0]+1)*phi_num],'b',  \
         phi_0, real_spd[time_index[1]*phi_num:(time_index[1]+1)*phi_num],'r',  \
         phi_0, real_spd[time_index[2]*phi_num:(time_index[2]+1)*phi_num],'g',  \
         phi_0, real_spd[time_index[3]*phi_num:(time_index[3]+1)*phi_num],'y', )
xlim(45,155)
plt.title('Real speed', fontsize=20)
plt.xlabel('longitude $\phi (^\circ)$', fontsize=15)
plt.ylabel('$V_{real} (km/s)$', fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.legend([shell_time[time_index[0]], shell_time[time_index[1]],\
            shell_time[time_index[2]],shell_time[time_index[3]],\
            shell_time[time_index[4]]])
plt.savefig('Real speed.png')
'''

fig = plt.figure(3,figsize=(16,12))
ax = fig.add_axes([0.15, 0.1, 0.75,0.8]) 
ax.plot(phi_0, max_e[time_index[0]*phi_num:(time_index[0]+1)*phi_num],'b',  \
         phi_0, max_e[time_index[1]*phi_num:(time_index[1]+1)*phi_num],'r',  \
         phi_0, max_e[time_index[2]*phi_num:(time_index[2]+1)*phi_num],'g',  \
         phi_0, max_e[time_index[3]*phi_num:(time_index[3]+1)*phi_num],'y', \
         phi_0, max_e[time_index[4]*phi_num:(time_index[4]+1)*phi_num],'k',)
ax.set_xlim([33,143])
#ax.set_ylim([0.01,1000])
ax.set_title('Maximum particle energy', fontsize=20)
ax.set_xlabel('longitude $\phi (^\circ)$', fontsize=15)
ax.set_ylabel('Max Energy (Mev)', fontsize=18)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.legend([shell_time[time_index[0]], shell_time[time_index[1]],\
            shell_time[time_index[2]],shell_time[time_index[3]],\
            shell_time[time_index[4]]])
ax.set_yscale('log')
plt.savefig('./Max energy.png')


#plt.show()















