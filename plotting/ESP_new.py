import math
import numpy as np
import matplotlib.pyplot as plt

#plot for time intensity profiles for the paper

if_esp = 1
root_dir = '../transport/proton_80_1.09/'
plot_title = 'Observer'
file_name = 'fp_total'
phi_e =100.0         # Earth longitude
cme_center =100.
r_e = 1.0  #radius of the point of interest(AU)
time_start = 0.5
time_end = 0.57
t_num = 50      # extended time number with ESP    
p_num_trsp = 40
n_spectra = 9
p_num = 400
phi_num = 25
f4 = open(root_dir + file_name, 'r')


#energy_index = [8,11,14,17,20,23]
#energy_index = [8,10,12,14,16,18]
#energy_index = [4,6,8,10,12,14]

#energy_index = [8,9,10,11,12,13,14]
#energy_index = [5,7,9,11,13,15,18]
#energy_index = [5,6,7,8,9,10]   

#energy_index = [1,3,5, 7,9,11]
energy_index = [10, 14, 18, 21, 23, 27,31]



del_phi = 5.
#normalization factors
t_o = 2858068.3
mp = 1.67e-27
co = 3e8
eo = 1.6e-19
vo = 52483.25
no = 1e6
AU = 1.5e11
pi = 3.14159265359
omega = 2.87e-6

#=======================================================================
#n_time_trsp = 40 # transport output time number

phi_0 = []
for i in xrange(0, phi_num):
       phi_0.append(int(cme_center +del_phi*(i-math.floor(phi_num/2.))))

phi_no = (int(phi_e) - phi_0[0])/5
print phi_no 

fp_t = []
p_num_a =[]
energy1 = []
energy1Mev = []
p_0 = []
v_0 = []
time_rise = []


for line in f4:
       line = line.strip()
       columns = line.split()

       p_num_a.append(columns[0])
       fp_t.append(float(columns[1])) #*(mp*vo*AU)**(-3.))      
       

fp_transport = []       
for i in xrange(0, p_num_trsp):
       emin = 1e5
       emax = 1e9	
       e_0 = emin * ((emax/emin)**(1./(p_num_trsp-1)))**i
       energy1.append(e_0)
       energy1Mev.append(e_0/1.e6)
       gamma0 = (e_0*eo + mp*co**2.0)/(mp*co**2.0)
       p_0.append(math.sqrt(((e_0*eo + mp*co**2.0)**2. -(mp*co**2.0)**2.0))/co)
       v_0.append(p_0[i] / mp /gamma0 )
       time_rise.append(1.1*AU /v_0[i] /3600.)


####################################################################################################  
# first read shock_momenta.dat to calculate shock arrival time
f1 = open(root_dir + 'shock_momenta.dat', 'r')

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



f3 = open(root_dir+'all_shell_bndy.dat', 'r')
shell_bndy =[]
for line in f3:
       line = line.strip()
       columns = line.split()

       shell_bndy.append(float(columns[3]))

shl_front = []
shl_end   = []

shell_start = range(t_zeus_no)
shell_start[0] = 0
iii = 2
for i in xrange(1, t_zeus_no):
       shell_start[i] = shell_start[i-1] + iii*phi_num
       iii = iii+1

#print 'shell',shell_start[0], shell_start[1],shell_start[2],shell_start[3]

for i in xrange(0, t_zeus_no):
       #shell front location at each time step at this angle:
       shl_front.append(shell_bndy[shell_start[i]+ (i+2)*phi_no+ i+1])
       #shell end location at each time step at this angle:
       shl_end.append(shell_bndy[shell_start[i]+ (i+2)*phi_no])
'''
print 'phi_no', phi_no
print shl_front
print shl_end
'''
arrive_index = 0
leave_index  = 0
shock_loc = []
time_zeus = []


for i in xrange(0, t_zeus_no):
       shock_loc.append(shell_r[i*phi_num + phi_no])
       time_zeus.append(shell_t[i*phi_num + phi_no])

       if shl_front[i] <= r_e:
              arrive_index = i
       if shl_end[i] <= r_e:
              leave_index = i
arriv_time = (0.6591-0.6)*t_o/3600. # the arrival time of the shock complex

leav_time = (0.6872-0.6)*t_o/3600. # the time when the shock complex inner boundary leaves the observer

arrive_time = ((r_e - shl_front[arrive_index])* time_zeus[arrive_index+1] + (shl_front[arrive_index+1]-r_e) \
              * time_zeus[arrive_index]) / (shl_front[arrive_index+1]-shl_front[arrive_index])   # in hours
leave_time = ((r_e - shl_end[leave_index])* time_zeus[leave_index+1] + (shl_end[leave_index+1]-r_e) \
              * time_zeus[leave_index]) / (shl_end[leave_index+1]-shl_end[leave_index])   # in hours
print "arriv_time, leav_time", arriv_time, leav_time              

leave_index = 0
first_index = 0
first_time = 0.0

xtime = []

if if_esp == 0:
       for i in xrange(0,t_num):
              xtime.append((i+1)*arrive_time/t_num)

if if_esp == 1:
       for i in xrange(0,t_num):
              xtime.append((i *(time_end -time_start)/t_num)*t_o/3600.)


if if_esp == 1:
       ########################################################################
       # adding ESP part 
       ########################################################################
       f2 = open(root_dir+ 'dist_all_shl.dat','r')
       col1 = []
       col2 = []
       col3 = []
       col4 = []
       energy        =[]
       for line in f2:
              line = line.strip()
              columns = line.split()

              energy.append(float(columns[0]))
              col1.append(float(columns[1]))
              col2.append(columns[2])
              col3.append(float(columns[3]))
              col4.append(columns[4])

       shell_start = range(t_zeus_no)

       for i in xrange(0, t_zeus_no):
              for j in xrange(0, i):
                     shell_start[i] = shell_start[i]+j

       print 'shell', shell_start[1],shell_start[2],shell_start[3]

       #-------------------------------------------------
       # search the bin of p_0 in p_0_zeus
       energy_zeus = energy[0:p_num]
       ipleft = [0]*p_num_trsp
       ishlleft = [0]*p_num_trsp

       for i in xrange(0, p_num_trsp):
              for y in xrange(0, p_num):
                     if energy_zeus[y] <= energy1[i]:
                            ipleft[i] = y


       shell_bndy_start = [0]*t_zeus_no
       for i in xrange(0, t_zeus_no):
              if i == 0:
                     shell_bndy_start[i] = 0
              else:
                     for j in xrange(0, i):
                            shell_bndy_start[i] = shell_bndy_start[i] + phi_num*(j+2)


       fp_esp = []
       fp_10 = []
       time_esp = []
       for i in xrange(0,20):      # Seperate the ESP time period into 20 bins
              time_now = arrive_time + (leave_time-arrive_time)*i/19.
              time_esp.append(time_now)
              for j in xrange(0, t_zeus_no):
                     if time_zeus[j] <= time_now:
                            time_left = j
              time_right = time_left + 1

              ratio = (1.-i/19.) # when ratio =0, observer is in the innermost shell of the shock complex


              left_shl_no = int(np.floor(ratio*time_left)) # shell number of the observer
                                                      # total shell number here is time_left+1
              right_shl_no = int(np.floor(ratio*time_right))
                                                      # total shell number here is time_right+1

              fp_left = col1[shell_start[time_left]*p_num*phi_num+(phi_no-1)*(time_left+1)*p_num + \
                          left_shl_no*p_num: \
                             shell_start[time_left]*p_num*phi_num+(phi_no-1)*(time_left+1)*p_num + \
                          (left_shl_no+1)*p_num]
              
              fp_right = col1[shell_start[time_right]*p_num*phi_num+(phi_no-1)*(time_right+1)*p_num + \
                          right_shl_no*p_num: \
                             shell_start[time_right]*p_num*phi_num+(phi_no-1)*(time_right+1)*p_num + \
                          (right_shl_no+1)*p_num]

              fp_esp_temp = []

              for k in xrange(0, p_num_trsp):
                     fp_left_temp = ((fp_left[ipleft[k]]*((energy_zeus[ipleft[k]+1]) - (energy1[k]) ) + \
                                  fp_left[ipleft[k]+1] *((energy1[k]) - (energy_zeus[ipleft[k]]))) /  \
                                   ((energy_zeus[ipleft[k]+1]) - (energy_zeus[ipleft[k]])))* p_0[k]*p_0[k]
                     
                     fp_right_temp = ((fp_right[ipleft[k]]*((energy_zeus[ipleft[k]+1]) - (energy1[k]) ) + \
                                  fp_right[ipleft[k]+1] *((energy1[k]) - (energy_zeus[ipleft[k]]))) /  \
                                   ((energy_zeus[ipleft[k]+1]) - (energy_zeus[ipleft[k]])))* p_0[k]*p_0[k]
                     fp_esp_temp.append( abs(   ( fp_left_temp * (time_zeus[time_right] - time_now) + \
                                              fp_right_temp * (time_now - time_zeus[time_left]) )  \
                                            / (time_zeus[time_right] - time_zeus[time_left]) )   )

                     # test
                     # print time_zeus[time_left], '<', time_now,'<', time_zeus[time_right]

              fp_esp.append(fp_esp_temp)



####################################################################################################       
# now read the fp

for time_no in xrange(0, t_num):  # convert f(p) into J(T) = f(p)*p^2
       fp_temp = []
       for i in xrange(0, p_num_trsp):
              fp_temp.append((fp_t[time_no*p_num_trsp+i])* p_0[i]*p_0[i])

       fp_transport.append(fp_temp)
       
       
time_intensity = []
original_intens = []


bg_ti_temp  = []
#===== background J(T), empirical

for i in xrange(0, p_num_trsp):
       coeff = 5.e-9  # base on observation
       vel = p_0[i]/mp*(1./math.sqrt(1.+(p_0[i]/(mp*co))**2.0))
#       bg_ti_temp.append(coeff * energy1[i]**2. * (energy1[i]/energy1[0])**(-3.5) \
#                     * (r_e/1.0)**-2.0 )
       bg_ti_temp.append(0.0)
#       bg_ti_temp.append(50. * (energy1[i]/1.e6)**(-2.213) \
#                     * (r_e/1.0)**-2.0 )

       print 'vel:', energy1[i], vel


       
#=======================================================================
#calculate time intensity profiles  J(T,t)
e_legd =[]

p_time = []


for i in xrange(0, t_num):
       if xtime[i] <= arrive_time:
              t1 = i
       if xtime[i] <= leave_time:
              t2 = i

print "t1, t2:", t1, t2

if if_esp == 1:              
       for i in xrange(0, p_num_trsp):
              e_legd.append('%(energy).1F %(unit)s' %{"energy":energy1[i]/1e6, "unit":"MeV"})
              ti_temp = []
              oti_temp = []
              p_time_temp = []


              # before shock arrival
              '''
              p_time_temp.append(0.0)
              p_time_temp.append(time_rise[i])
              ti_temp.append(bg_ti_temp[i])
              ti_temp.append(bg_ti_temp[i])
              '''
              
              for itime in xrange(0, t1):
                     p_time_temp.append(xtime[itime])
                     ti_temp.append(fp_transport[itime][i]*eo*100.*no*(mp*vo)**(-3.)/(4.* pi) +bg_ti_temp[i])
                     
              # esp phase
              for itime in xrange(0,20):
                     p_time_temp.append(time_esp[itime]) 

                     t_zeus = 0
                     for ii in xrange(t1, t2):
                            if xtime[ii] <= time_esp[itime]:
                                   t_zeus = ii
                     
                     f_esp_temp = abs((fp_transport[t_zeus][i] * (xtime[t_zeus+1]-time_esp[itime]) + fp_transport[t_zeus+1][i]  \
                                   * (time_esp[itime] - xtime[t_zeus]) ) / (xtime[t_zeus+1] - xtime[t_zeus])) 
                     
                     ti_temp.append((fp_esp[itime][i] + f_esp_temp/2.718)*eo*100.*no*(mp*vo)**(-3.)/(4.* pi) +bg_ti_temp[i])
              
#                     ti_temp.append((fp_esp[itime][i])*eo*100.*no*(mp*vo)**(-3.)/(4.* pi))
              
              #after shock leaving
              for itime in xrange(t2+1, t_num):
                     p_time_temp.append(xtime[itime])
                     ti_temp.append(fp_transport[itime][i]*eo*100.*no*(mp*vo)**(-3.)/(4.* pi) )# +bg_ti_temp[i])#/1.5/1.5)
                                                        #{-------- normalization factors-----------------}
              for itime in xrange(0, t_num):
                     oti_temp.append(fp_transport[itime][i]*eo*100.*no*(mp*vo)**(-3.)/(4.* pi) )# +bg_ti_temp[i])

       #       print p_time_temp
              p_time.append(p_time_temp)
              time_intensity.append(ti_temp)
              original_intens.append(oti_temp)

if if_esp == 0:
       for i in xrange(0, p_num_trsp):
              e_legd.append('%(energy).1F %(unit)s' %{"energy":energy1[i]/1e6, "unit":"MeV"})

              oti_temp = []
              for itime in xrange(0, t_num):
                     oti_temp.append(fp_transport[itime][i]*eo*100.*no*(mp*vo)**(-3.)/(4.* pi) +bg_ti_temp[i])

       #       print p_time_temp
              time_intensity.append(oti_temp)
#print p_time[0]

       
#calculate event integrated spectra

#print len(p_time[0]),len(p_time[8])


event_spec = []
event_spec_no_esp = []
event_spec_small = []
b4_arrival_spec = []

ESP_spec = []
ESP_spec2 = []
ESP_spec3 = []
ESP_spec4 = []
for i in xrange(0, p_num_trsp):
       ESP_spec.append(time_intensity[i][t1+1])
       ESP_spec2.append(time_intensity[i][t1+6])
       ESP_spec3.append(time_intensity[i][t1+11])
       ESP_spec4.append(time_intensity[i][t1+16])

dtt = (time_esp[2] - time_esp[1])*3600.

print "dtt is", dtt

ESP_TOT_SPEC = []
for j in xrange(0, p_num_trsp):
       spec_temp = 0.0
       for itime in xrange(t1,t1+20):
              spec_temp = spec_temp + time_intensity[j][itime]* dtt
       ESP_TOT_SPEC.append(spec_temp)



print "time:", time_esp[1]-time_esp[0],time_esp[6]-time_esp[0],time_esp[11]-time_esp[0],time_esp[16]-time_esp[0]
for i in xrange(0, p_num_trsp):
       aaa = xtime[0]*  (fp_transport[0][i]*eo*100.*no*(mp*vo)**(-3.)/(4.* pi)  +bg_ti_temp[i]) *3600./2.
       bbb = xtime[0]*  time_intensity[i][0]*3600./2.      
       ccc = xtime[0]*  time_intensity[i][0]*3600./2.
       
       if if_esp == 1:
              for j in xrange(1, len(p_time[0])):
                     bbb = bbb + (p_time[i][j]-p_time[i][j-1]) * time_intensity[i][j]*3600.
              
              for j in xrange(1, t_num):
                     aaa = aaa + (xtime[j] - xtime[j-1])* (fp_transport[j][i]*eo*100.*no*(mp*vo)**(-3.)/(4.* pi)  +bg_ti_temp[i]) *3600.
                     
              for j in xrange(1, t1+1):
                     ccc = ccc + (xtime[j] - xtime[j-1])* time_intensity[i][j]*3600.
              
              event_spec.append( bbb )
              event_spec_no_esp.append( aaa )
              b4_arrival_spec.append( ccc )
       
       if if_esp == 0:
              for j in xrange(1, t_num):
                     aaa = aaa + (xtime[j] - xtime[j-1])* (fp_transport[j][i]*eo*100.*no*(mp*vo)**(-3.)/(4.* pi)  +bg_ti_temp[i]) *3600.
              event_spec.append( aaa )
              event_spec_small.append( aaa/10000. )



'''
f5=open(root_dir+file_name+'_event_spec_no_esp.dat', 'w')
f6=open(root_dir+file_name+'_event_spec.dat', 'w')
for i in xrange(0, p_num_trsp):
       f6.write(str(energy1Mev[i])+'      '+str(event_spec[i])+'\n')
       if if_esp == 1:
              f5.write(str(energy1Mev[i])+'      '+str(event_spec_no_esp[i])+'\n')

f5.close
f6.close       
'''         

'''
plt.figure(1, figsize=(16,12))
plt.plot(energy1,spec_Elli_Ram, 'k', \
         linewidth = 1.5)

plt.title('iPATH result for STEREO-A', fontsize=35)
plt.yscale('log')
plt.xscale('log')
plt.ylim([1e-15, 1e-5])
plt.xlim([1e5,1e9])
plt.xlabel('time after CME eruption (hours)', fontsize=30)
plt.axvline(arriv_time, color='black', linestyle='dashed', linewidth=2)
#plt.axvline(leave_time, color='black', linestyle='dashed', linewidth=2)
plt.axvline(xtime[4], color='red', linestyle='dashed', linewidth=2)
plt.ylabel('$J_T(T)$ $(counts/(cm^2 s sr MeV))$', fontsize=30)
plt.legend([e_legd[energy_index[0]], e_legd[energy_index[1]], e_legd[energy_index[2]], \
            e_legd[energy_index[3]], e_legd[energy_index[4]], e_legd[energy_index[5]]], \
            loc=2, ncol = 3, borderaxespad=0., shadow = True, fontsize=25)
plt.tick_params(axis='both', which='major', labelsize=36)
plt.tick_params(axis='both', which='major', labelsize=24)
'''

  
#############################################################################################################################
#             PLOTTING
#############################################################################################################################
#------------------------------------------------
#plot time intensity
plt.figure(2, figsize=(14,14))


if if_esp == 1:
       plt.plot(p_time[energy_index[0]], time_intensity[energy_index[0]], 'r', p_time[energy_index[1]], time_intensity[energy_index[1]], 'y', \
                p_time[energy_index[2]], time_intensity[energy_index[2]], 'g', p_time[energy_index[3]], time_intensity[energy_index[3]], 'b', \
                p_time[energy_index[4]], time_intensity[energy_index[4]], 'purple', p_time[energy_index[5]], time_intensity[energy_index[5]], 'k', \
                 linewidth = 1.5)
if if_esp == 0:
       plt.plot(xtime, time_intensity[energy_index[0]], 'r', xtime, time_intensity[energy_index[1]], 'y', \
                xtime, time_intensity[energy_index[2]], 'g', xtime, time_intensity[energy_index[3]], 'b', \
                xtime, time_intensity[energy_index[4]], 'purple', xtime, time_intensity[energy_index[5]], 'k', \
                linewidth = 1.5)

plt.title('Time-intensity Profiles', fontsize=35)
plt.yscale('log')
plt.ylim([1e-3, 1e6])
#plt.xlim([0,42])
plt.xlabel('time after CME eruption (hours)', fontsize=30)
plt.axvline(arrive_time, color='black', linestyle='dashed', linewidth=2)
plt.axvline(leave_time, color='black', linestyle='dashed', linewidth=2)
#plt.axvline(xtime[4], color='red', linestyle='dashed', linewidth=2)
plt.ylabel('$J_T(T)$ $(counts/(cm^2 s sr MeV))$', fontsize=30)
plt.legend([e_legd[energy_index[0]], e_legd[energy_index[1]], e_legd[energy_index[2]], \
            e_legd[energy_index[3]], e_legd[energy_index[4]], e_legd[energy_index[5]] ], \
            loc=2, ncol = 3, borderaxespad=0., shadow = True, fontsize=25)
plt.tick_params(axis='both', which='major', labelsize=36)
plt.tick_params(axis='both', which='major', labelsize=24)
plt.savefig('./time_intens'+file_name+'.png')



#------------------------------------------------
#plot  integrated 
'''
fig=plt.figure(4, figsize=(11,12))
ax = fig.add_axes([0.15, 0.1, 0.75,0.8]) 


if if_esp == 1:
       ax.plot(energy1Mev, event_spec,'-*', energy1Mev, b4_arrival_spec,'-*' , energy1Mev, event_spec_no_esp, '-*' ,linewidth=1.3)
       ax.legend(['event spec', 'pre-arrival spec', 'event spec no esp'])
if if_esp == 0:
       ax.plot(energy1Mev, event_spec, '-*' ,linewidth=1.3)
       ax.legend(['event spec'])

ax.plot(energy1Mev, ESP_spec,  energy1Mev, ESP_spec2, \
        energy1Mev, ESP_spec3, energy1Mev, ESP_spec4,linewidth=1.3)

labels=['During ESP', 'Pre Arrival', 'Total']

ax.plot(energy1Mev, ESP_TOT_SPEC, linewidth=2)
ax.plot(energy1Mev, b4_arrival_spec, linewidth=2)
ax.plot(energy1Mev, event_spec, linewidth=2)

ax.set_xscale('log')
ax.set_yscale('log')
plt.xlim([1, 1000])
plt.ylim([1e2, 1e9])
#ax.set_xlabel('Energy(MeV/nucleon)', fontsize=25)
#ax.set_ylabel('differential fluence $(counts/(cm^2 MeV))$', fontsize=25)
#ax.set_ylabel('$J_T(T)$ $(counts/(cm^2 s sr MeV))$', fontsize=25)

ax.set_ylabel('Proton fluence $(counts/(cm^2 sr MeV))$', fontsize=25)
ax.tick_params(axis='both', which='major', labelsize=22)

plt.legend(labels, ncol=5, loc='upper right', 
#           bbox_to_anchor=[0.5, 1.1], 
           columnspacing=1.0, labelspacing=0.0,
           handletextpad=0.0, handlelength=1.5,
           fancybox=True, shadow=True, fontsize=22)
#ax.set_title('event integrated intensity at 1AU, '+str(phi_e)+'$^\circ$ \n', fontsize=26)
ax.set_title('Time Integrated Energy spectra', fontsize=30)
plt.savefig('./spectrum_'+file_name+'.png')
'''
#--------------------------------------------------------------------------------------
# plot spectrum time evolution
'''
time_no=len(xtime)
print time_no 
All_spec = []
for i in xrange(0, time_no):
       temp =[]
       for j in xrange(0, p_num_trsp):
              temp.append(time_intensity[j][i])
       All_spec.append(temp)

print np.shape(All_spec), np.shape(time_intensity)

plt.figure(3, figsize=(10,12))
spec_index = [4, 11, 18, 25, 32, 39]
colormap = plt.cm.gist_ncar
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, len(spec_index))])


t0 = 0.
labels=[]
for i in xrange(0, len(spec_index)):
       plt.plot(energy1Mev, All_spec[spec_index[i]], linewidth=2)
       labels.append(str(int(xtime[spec_index[i]] - t0))+ ' hrs')

plt.plot(energy1Mev, event_spec_small, '-*' ,linewidth=2)
labels.append('Event Spectrum')

plt.legend(labels, ncol=5, loc='upper center', 
           bbox_to_anchor=[0.5, 1.1], 
           columnspacing=1.0, labelspacing=0.0,
           handletextpad=0.0, handlelength=1.5,
           fancybox=True, shadow=True, fontsize=20)

plt.legend(labels, fontsize=20)
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-2, 1e5])
plt.xlabel('Energy(MeV/nucleon)', fontsize=25)
plt.ylabel('Proton fluence $(counts/(cm^2 sr MeV))$', fontsize=25)
plt.tick_params(axis='both', which='major', labelsize=22)
plt.savefig(root_dir+'spectrum_'+file_name+'.png')
'''





f1.close

f4.close
if if_esp == 1:
       f2.close
       f3.close



