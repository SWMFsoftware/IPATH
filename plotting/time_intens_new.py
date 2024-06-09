import math
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')

#plot for time intensity profiles for the paper

if_esp = 0
root_dir = './'
root_dir2 = '../transport/proton_80_1.09/'

plot_title = 'Observer'
anno = ' '
phi_e = 88.0        # Earth longitude
r_e = 1.0  #radius of the point of interest(AU)
f4 = open(root_dir2 + 'fp_total', 'r')


energy_index = [10, 14, 18, 22, 26, 30]
time_index = [3,5,8,11,14]

AU = 1.5e11
cme_center =88.
del_phi = 5.

pi = 3.14159265359
omega = 2.87e-6

time_start = 0.40
time_end = 0.50


#normalization factors
t_o = 2858068.3
mp = 1.67e-27
co = 3e8
eo = 1.6e-19
vo = 52483.25
no = 1e6

#=======================================================================
#n_time_trsp = 40 # transport output time number
t_num = 50      # extended time number with ESP    
p_num_trsp = 40

n_spectra = 9


p_num = 400
phi_num = 25




phi_0 = []
for i in xrange(0, phi_num):
       phi_0.append(int(cme_center +del_phi*(i-math.floor(phi_num/2.))))

phi_no = (int(phi_e) - phi_0[0])/5
print phi_no 


#energy_index = [3,5,7,9,11,13]
#energy_index = [2,5,8,11,14,17]



####################################################################################################  
# first read shock_momenta.dat to calculate shock arrival time
f1 = open(root_dir2 + 'shock_momenta.dat', 'r')

energy        =[]
p0            =[]
fp            =[]
f_hi          =[]
shell_t       =[]
shell_r       =[]
max_e         =[]

line_no = 0
for line in f1:
       line = line.strip()
       columns = line.split()

       shell_r.append(float(columns[0]))
       shell_t.append(float(columns[1]))
       max_e.append(float(columns[5])/1e6)       
       
       line_no = line_no + 1


t_zeus_no = line_no/phi_num

print 't_zeus_no', t_zeus_no

shell_time =[]

shock_loc = []
time_zeus = []
time_AU = 0


for i in range(0, t_zeus_no):
       shock_loc.append(shell_r[i*phi_num + phi_no])
       time_zeus.append(shell_t[i*phi_num + phi_no])
       shell_time.append('%(time).1F %(unit)s' %{"time":shell_t[i*phi_num], "unit":"hrs"})
       
       if shell_r[i*phi_num + phi_no] <= r_e:
              time_AU = i

print time_AU
if if_esp == 0:
       arriv_time = ((r_e - shock_loc[time_AU])* time_zeus[time_AU+1] + (shock_loc[time_AU+1]-r_e) \
              * time_zeus[time_AU]) / (shock_loc[time_AU+1]-shock_loc[time_AU])   # in hours
else:
       arriv_time = (time_end-0.4)*t_o/3600.

leave_time = (0.4589-0.4)*t_o/3600.
print 'time_AU, arrival time(h)', time_AU, arriv_time

print arriv_time

xtime = []

if if_esp == 0:
       for i in xrange(0,t_num):
              xtime.append((i+1)*arriv_time/t_num)

if if_esp == 1:
       for i in xrange(0,t_num):
              xtime.append((i *(time_end -time_start)/t_num)*t_o/3600.)


####################################################################################################       
# now read the fp
fp_t = []
p_num_a =[]
energy1 = []
energy1Mev = []
p_0 = []

for line in f4:
       line = line.strip()
       columns = line.split()

       p_num_a.append(columns[0])
       fp_t.append(float(columns[1])) #*(mp*vo*AU)**(-3.)) 
       fp.append(float(columns[1])*eo*100.*no*(mp*vo)**(-3.))     
#       fp.append(float(columns[1])*eo*100.*no*(mp*vo)**(-3.)/(4.* pi))   


fp_transport = []       
for i in xrange(0, p_num_trsp):
       emin = 1e5
       emax = 1e9	
       e_0 = emin * ((emax/emin)**(1./(p_num_trsp-1)))**i
       energy1.append(e_0)
       energy1Mev.append(e_0/1.e6)
       gamma0 = (e_0*eo + mp*co**2.0)/(mp*co**2.0)
       p_0.append(math.sqrt(((e_0*eo + mp*co**2.0)**2. -(mp*co**2.0)**2.0))/co)

#print energy1Mev

for time_no in xrange(0, t_num):  # convert f(p) into J(T) = f(p)*p^2
       fp_temp = []
       for i in xrange(0, p_num_trsp):
              fp_temp.append(fp_t[time_no*p_num_trsp+i]* p_0[i]*p_0[i] )
       fp_transport.append(fp_temp)
       
       
time_intensity = []


bg_ti_temp  = []
#===== background J(T), empirical

for i in xrange(0, p_num_trsp):
       coeff = 5.e-9  # base on observation
       vel = p_0[i]/mp*(1./math.sqrt(1.+(p_0[i]/(mp*co))**2.0))
       bg_ti_temp.append( coeff * energy1[i]**2. * (energy1[i]/energy1[0])**(-3.5) )
       print 'vel:', energy1[i], vel
'''
for i in xrange(0, p_num_trsp):
       vel = p_0[i]/mp*(1./math.sqrt(1.+(p_0[i]/(mp*co))**2.0))
       bg_ti_temp.append( 4.0*40000000./energy1[i]*1e6 /4./pi* (energy1[i]/energy1[0])**(-3.5) )
       print 'vel:', energy1[i], vel
'''

       
#=======================================================================
#calculate time intensity profiles  J(T,t)
e_legd =[]
for i in xrange(0, p_num_trsp):
       e_legd.append('%(energy).1F %(unit)s' %{"energy":energy1[i]/1e6, "unit":"MeV"})
       ti_temp = []
       for itime in xrange(0, t_num):
              ti_temp.append(fp_transport[itime][i]*eo*100.*no*(mp*vo)**(-3.)/(4.* pi)  + bg_ti_temp[i])#/1.5/1.5)
                                                 #{-------- normalization factors-----------------}
       time_intensity.append(ti_temp)
       
#=======================================================================
#calculate time integrated intensity (fluence, unit is #/(cm^2 MeV) )

time_int_spec = []
time_legend	=	[]
for i in xrange(0, n_spectra):
       time_fp_temp = []
       for j in xrange(0, p_num_trsp):
              aaa=(fp_t[4*i*p_num_trsp+j] + 2.*fp_t[(4*i+1)*p_num_trsp+j] + \
	              2.*fp_t[(4*i+2)*p_num_trsp+j] + fp_t[(4*i+3)*p_num_trsp+j])/2.*(xtime[1]-xtime[0])*3600. \
                     * p_0[j]*p_0[j] *eo*100.*no*(mp*vo)**(-3.)

              if aaa < 1e1:
                     aaa=np.nan  
              time_fp_temp.append(aaa)
       time_int_spec.append(time_fp_temp)
       print i, len(xtime)
       time_legend.append('%(time1).1F-%(time2).1F (hr)' %{"time1":xtime[4*i], "time2":xtime[4*(i+1)]})



#=======================================================================
total_fp=np.zeros(shape=(p_num_trsp,1))

for i in range(0, p_num_trsp):
       for j in range(0, time_no):
              total_fp[i] = total_fp[i]+ (float(fp[j*p_num_trsp+i])* p_0[i]*p_0[i])
       total_fp[i] = total_fp[i] * (arriv_time - time_start)*t_o/time_no
       if total_fp[i]< 100.:
              total_fp[i] = 0.0


	
#=======================================================================
# calculate energy-integrated intenisty   
 
flux10 = np.zeros(shape=(t_num,1))
flux50 = np.zeros(shape=(t_num,1))
flux100 = np.zeros(shape=(t_num,1))

for itime in xrange(0, t_num):

       for i in xrange(0, p_num_trsp):

              if energy1Mev[i] >= 10:
                  flux10[itime] = flux10[itime] + time_intensity[i][itime]  

              if energy1Mev[i] >= 50:
                  flux50[itime] = flux50[itime] + time_intensity[i][itime] 

              if energy1Mev[i] >= 100:
                  flux100[itime] = flux100[itime] + time_intensity[i][itime]    
  

#############################################################################################################################
#             PLOTTING
#############################################################################################################################
#------------------------------------------------
#plot Emax 
'''
fig = plt.figure(1,figsize=(16,12))
ax = fig.add_axes([0.15, 0.1, 0.75,0.8]) 

plt.plot(phi_0, max_e[time_index[0]*phi_num:(time_index[0]+1)*phi_num],'b',  \
         phi_0, max_e[time_index[1]*phi_num:(time_index[1]+1)*phi_num],'r',  \
         phi_0, max_e[time_index[2]*phi_num:(time_index[2]+1)*phi_num],'g',  \
         phi_0, max_e[time_index[3]*phi_num:(time_index[3]+1)*phi_num],'y', \
         phi_0, max_e[time_index[4]*phi_num:(time_index[4]+1)*phi_num],'k', linewidth = 1.3)
plt.ylabel('maximum energy(MeV)', fontsize=20)
plt.xlabel('longitude $\phi$', fontsize=20)
plt.xlim([31, 145])
#plt.ylim([0, 100])
plt.legend([shell_time[time_index[0]], shell_time[time_index[1]],\
            shell_time[time_index[2]],shell_time[time_index[3]],\
            shell_time[time_index[4]]],prop={'size':20}, loc=2)
ax.set_yscale('log')
'''
#------------------------------------------------
#plot time intensity


plt.figure(2, figsize=(16,12))
plt.plot(xtime, time_intensity[energy_index[0]], 'k-', xtime, time_intensity[energy_index[1]], 'g-', \
         xtime, time_intensity[energy_index[2]], 'r-', xtime, time_intensity[energy_index[3]], 'b-', \
         xtime, time_intensity[energy_index[4]], 'y-', xtime, time_intensity[energy_index[5]], 'c-', \
         linewidth = 2.5)
#print('***')
#print(time_intensity[2])
#print('???')
#print(time_intensity[:][2])
#print('???')

with open('p3', 'w') as f1:
    for itime in xrange(0, t_num):
        for i in xrange(0, p_num_trsp):
            f1.write('%15.6f\t %15.6f\t %15.6f\n' %(xtime[itime], energy1Mev[i], time_intensity[i][itime]))

plt.title(plot_title+" ($\phi$ ="+str(phi_e)+"$^\circ$)", fontsize=35)
plt.yscale('log')
plt.ylim([1e-4, 1e6])
plt.xlim([0, 40])
plt.xlabel('time (hours)', fontsize=30)
plt.axvline(arriv_time, color='black', linestyle='dashed', linewidth=2)
plt.annotate(anno, xy=(0.83,0.81), xycoords='figure fraction', color='red',fontsize = 70 )
#plt.axvline(leave_time, color='black', linestyle='dashed', linewidth=2)
plt.ylabel('$J_T(T)$ $(counts/(cm^2 s sr MeV))$', fontsize=30)
plt.legend([e_legd[energy_index[0]], e_legd[energy_index[1]], e_legd[energy_index[2]], \
            e_legd[energy_index[3]], e_legd[energy_index[4]], e_legd[energy_index[5]]], \
            loc=2, ncol = 3, borderaxespad=0., shadow = True, fontsize=25)
plt.tick_params(axis='both', which='major', labelsize=36)
plt.tick_params(axis='both', which='major', labelsize=24)
plt.savefig(root_dir + 'time-intensity.png')

#------------------------------------------------
#plot time integrated 

fig=plt.figure(3, figsize=(9,12))
ax = fig.add_axes([0.15, 0.1, 0.75,0.8]) 

ax.plot(energy1Mev, time_int_spec[1],'-o' ,linewidth=1.3)

#ax.plot(energy1Mev, total_fp,'-o' ,linewidth=3)

ax.set_xlim([3e-1,1e3])
ax.set_ylim([1e2,4e10])
ax.set_xlabel('Energy(MeV/nucleon)', fontsize=25)
ax.set_ylabel('differential fluence $(counts/(cm^2 MeV))$', fontsize=25)
ax.tick_params(axis='both', which='major', labelsize=22)

for i in [2,4,6,8]:
	ax.plot(energy1Mev, time_int_spec[i],'-o',linewidth=1.3)
	ax.set_yscale('log')
	ax.set_xscale('log')

ax.legend([time_legend[1],time_legend[2],time_legend[4],time_legend[6],time_legend[8]], prop={'size':20}, loc=3)
ax.set_title('time-integrated intensity at 1AU, '+str(phi_e)+'$^\circ$ \n', fontsize=26)
plt.savefig(root_dir + 'time-integrated-intensity.png')

#------------------------------------------------
#plot energy integrated 
plt.figure(4, figsize=(16,12))
plt.plot(xtime, flux10, 'r-',linewidth = 2.5, label=">= 10 MeV")
plt.plot(xtime, flux50, 'b-',linewidth = 2.5, label=">= 10 MeV")
plt.plot(xtime, flux100, 'g-',linewidth = 2.5, label=">= 10 MeV")
plt.title(plot_title+" ($\phi$ ="+str(phi_e)+"$^\circ$)", fontsize=35)
plt.yscale('log')
plt.ylim([1e-4, 1e6])
plt.xlim([0, 30])
plt.xlabel('time (hours)', fontsize=30)
#plt.axvline(arriv_time, color='black', linestyle='dashed', linewidth=2)
#plt.annotate(anno, xy=(0.83,0.81), xycoords='figure fraction', color='red',fontsize = 70 )
#plt.axvline(leave_time, color='black', linestyle='dashed', linewidth=2)
plt.ylabel('Protons $(counts/(cm^2 s sr))$', fontsize=30)

label = [">= 10 MeV", ">= 50 MeV", ">= 100 MeV"]  
plt.legend(label, loc = 0, ncol = 1, borderaxespad=0., shadow = True, fontsize=25)

plt.tick_params(axis='both', which='major', labelsize=36)
plt.tick_params(axis='both', which='major', labelsize=24)

plt.savefig(root_dir + 'proton_flux.png')



#plt.show()




f1.close

f4.close
