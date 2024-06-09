import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.projections import PolarAxes
from osgeo import gdal
import math
plt.switch_backend('agg')

#tr = PolarAxes.PolarTransform()

#root_dir = './Jul_23_event_2_closer_norm/'
root_dir = './'
n_start = 1
n_end   = 30


#fileName = './zeus_cycle_24/zhto098TW'
#gdal_dataset=gdal.Open(fileName)
#gdal_dataset.GetSubDatasets()

xdims = 1500
x1min = 0.05
x1max = 2
x1rat = 1.0



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


x3a = []

for i in xrange(0, 360):
       x3a.append(float(i)*361./360.)
azimuths = np.radians(x3a)

r, theta = np.meshgrid(x1a, azimuths)
#print r.shape, theta.shape


# for field lines:
num = 10  # number of field lines
theta_space = 360./num
dx1a = x1a[1] - x1a[0]
#dt = 0.00015/5
#dt = 0.0002
dt = 0.00015/15
maxt = 10000
pi = 3.14159265359

max_n = 0.0


################
# for each HDF file:
print 'Calculate maximnum for all HDF4 data...'
for ii in xrange(n_start, n_start+2):
       file_no = "../CME/zhto{:03d}TW".format(ii)
       # read the solar wind values from HDF files
       print 'HDF4_SDS:UNKNOWN:{:s}:3'.format(file_no)
       v1_data = gdal.Open('HDF4_SDS:UNKNOWN:{:s}:3'.format(file_no))  # (360,1,x)
       b1_data = gdal.Open('HDF4_SDS:UNKNOWN:{:s}:15'.format(file_no))
       b3_data = gdal.Open('HDF4_SDS:UNKNOWN:{:s}:23'.format(file_no))
       dens_data = gdal.Open('HDF4_SDS:UNKNOWN:{:s}:27'.format(file_no))
       inter_energy_data= gdal.Open('HDF4_SDS:UNKNOWN:{:s}:31'.format(file_no))
       
       v1 = v1_data.ReadAsArray()
       b1 = b1_data.ReadAsArray()
       b3 = b3_data.ReadAsArray()
       dens = dens_data.ReadAsArray()
       inter_energy = inter_energy_data.ReadAsArray()
       
       dens_norm = np.ndarray((360, xdims))

       for i in xrange(0, 360):
              for j in xrange(0,xdims):
                     dens_norm[i,j] = dens[i,0,j]*x1a[j]**2.
       if (np.max(dens_norm) > max_n):
              max_n = np.max(dens_norm)

max_n = int(max_n) * 0.8
kk = int(max_n/10)
max_n = kk * 10
################



# for each HDF file:
for ii in xrange(n_start, n_end+1):
       file_no = "../CME/zhto{:03d}TW".format(ii)
       # read the solar wind values from HDF files
       print 'HDF4_SDS:UNKNOWN:{:s}:3'.format(file_no)
       v1_data = gdal.Open('HDF4_SDS:UNKNOWN:{:s}:3'.format(file_no))  # (360,1,x)
       b1_data = gdal.Open('HDF4_SDS:UNKNOWN:{:s}:15'.format(file_no))
       b3_data = gdal.Open('HDF4_SDS:UNKNOWN:{:s}:23'.format(file_no))
       dens_data = gdal.Open('HDF4_SDS:UNKNOWN:{:s}:27'.format(file_no))
       inter_energy_data= gdal.Open('HDF4_SDS:UNKNOWN:{:s}:31'.format(file_no))
       
       v1 = v1_data.ReadAsArray()
       b1 = b1_data.ReadAsArray()
       b3 = b3_data.ReadAsArray()
       dens = dens_data.ReadAsArray()
       inter_energy = inter_energy_data.ReadAsArray()
       
       dens_norm = np.ndarray((360, xdims))

       for i in xrange(0, 360):
              for j in xrange(0,xdims):
                     dens_norm[i,j] = dens[i,0,j]*x1a[j]**2.
#                      dens_norm[i,j] = v1[i,0,j]
       # plot field lines

       fl0_r = []
       fl0_th = []
       for i in xrange(0,num):
              mfl_r = []
              mfl_th = []
              temp_r = x1min
              temp_th = (0.0 + i*theta_space)
       
              while temp_r <= x1max:
                     if x1rat == 1.0:
                           r_index = (temp_r-x1min)/dx1a
                     else:
                           r_index = math.log(1 - (temp_r-x1min) * (1 - x1rat) /dx1a) / math.log(x1rat)

                     dr  = b1[int(round(temp_th)-1), 0, int(r_index)] * dt
                     dth = b3[int(round(temp_th)-1), 0, int(r_index)] * dt / temp_r *180./pi
                     

                     temp_r = temp_r + dr
                     temp_th = temp_th + dth
                     
                     if temp_th < 0:
                            temp_th += 360.
                     if temp_th >= 360:
                            temp_th -= 360.

                     mfl_r.append(temp_r)
                     mfl_th.append(temp_th/180.*pi)
              fl0_r.append(mfl_r)
              fl0_th.append(mfl_th)


       all_fl_r =[]
       all_fl_th =[]

       mid_pos = [[1, 0], [1,120],[1,250]]
       target = []
       for i in xrange(0,3):
              mfl_r = []
              mfl_th = []

              temp_r = mid_pos[i][0]
              temp_th = mid_pos[i][1]
              target.append([temp_th/180.*pi, temp_r])

              while temp_r >= x1min:
                     mfl_r.append(temp_r)
                     mfl_th.append(temp_th/180.*pi)

                     if x1rat == 1.0:
                           r_index = (temp_r-x1min)/dx1a
                     else:
                           r_index = math.log(1 - (temp_r-x1min) * (1 - x1rat) /dx1a) / math.log(x1rat)

                     dr  = b1[int(round(temp_th)-1), 0, int(r_index)] * dt
                     dth = b3[int(round(temp_th)-1), 0, int(r_index)] * dt / temp_r *180./pi

                     temp_r = temp_r - dr
                     temp_th = temp_th - dth

                     if temp_th < 0:
                            temp_th += 360.
                     if temp_th >= 360:
                            temp_th -= 360.

              mfl_r.reverse()
              mfl_th.reverse()

              temp_r = mid_pos[i][0]
              temp_th = mid_pos[i][1]

              while temp_r <= x1max:

                     if x1rat == 1.0:
                           r_index = (temp_r-x1min)/dx1a
                     else:
                           r_index = math.log(1 - (temp_r-x1min) * (1 - x1rat) /dx1a) / math.log(x1rat)

                     dr  = b1[int(round(temp_th)-1), 0, int(r_index)] * dt
                     dth = b3[int(round(temp_th)-1), 0, int(r_index)] * dt / temp_r *180./pi

                     temp_r = temp_r + dr
                     temp_th = temp_th + dth
                     
                     if temp_th < 0:
                            temp_th += 360.
                     if temp_th >= 360:
                            temp_th -= 360.

                     mfl_r.append(temp_r)
                     mfl_th.append(temp_th/180.*pi)

              all_fl_r.append(mfl_r)
              all_fl_th.append(mfl_th)

       ticks = []
       

       for i in xrange(1,11):
              ticks.append(i*max_n/10.)

       print target

       print np.min(dens_norm), np.max(dens_norm)
       #fig = plt.subplots(subplot_kw=dict(projection='polar'), figsize=(12,12))
       
       fig = plt.figure(1, figsize=(12,12))
       ax = fig.add_axes([0.1, 0.03, 0.8, 0.8], projection='polar')


       ax.set_rgrids(np.arange(0.5,2.5,0.5))
       ax.tick_params(pad = 15)

       pcm = ax.contourf(theta, r, dens_norm, cmap='jet',  extend='max', \
              levels = np.linspace(0,max_n,150), vmin=0.0, vmax=max_n, yunits ='AU')



       ax.plot(np.linspace(0,2*pi,360), [1.0]*360, 'k', linewidth=2.0)

       for i in xrange(0,num):
              ax.plot(fl0_th[i], fl0_r[i], 'k')


#       for i in xrange(0,3):
#              ax.plot(all_fl_th[i], all_fl_r[i], 'w-', linewidth=1.5)
#              ax.plot(all_fl_th[i], all_fl_r[i], 'k--', linewidth=1.5)
#       for i in xrange(0,3):
#              ax.plot(target[i][0], target[i][1], 'ro',linewidth=2.0)

#       ax.plot(70./180.*pi, 1.5, 'ro',linewidth=2.0)

#       ax.annotate('Earth', xy=(mid_pos[0][1]/180.*pi, 1.1), color='w', fontsize=32)
#       ax.annotate('STA', xy=(132./180.*pi, 1.3), color='w', fontsize=32)
#       ax.annotate('STB', xy=(mid_pos[2][1]/180.*pi, 1.3), color='w', fontsize=32)
#       ax.annotate('D', xy=(70/180.*pi, 1.6), color='w', fontsize=32)


#       ax.annotate('E', xy=(mid_pos[0][1]/180.*pi, 1.1), color='r', fontsize=28)
#       ax.annotate('A', xy=(mid_pos[1][1]/180.*pi, 1.1), color='r', fontsize=28)
#       ax.annotate('B', xy=(mid_pos[2][1]/180.*pi, 1.1), color='r', fontsize=28)

            # ,  # theta, radius
            # xytext=(0.05, 0.05),    # fraction, fraction
            # textcoords='figure fraction',
            # arrowprops=dict(facecolor='black', shrink=0.05),
            # horizontalalignment='left',
            # verticalalignment='bottom',
            # )





       ax.grid(True,alpha= 0.4,linestyle='--',color='black')
       ax.set_rmin(0.0)
       ax.set_rmax(x1max)
       ax.tick_params(axis='both', labelsize=18)
	





       # magnetic field lines
       cbaxes = fig.add_axes([0.1, 0.9, 0.8, 0.03]) 
       cb = plt.colorbar(pcm, cax = cbaxes,orientation='horizontal',\
             ticks= ticks) 
       cb.ax.tick_params(direction='in',width=1,length=10, labelsize = 15)
       


       cbaxes.set_title('$R^2 N(AU^2cm^{-3})$', fontsize=18)

#       plt.text(0, -0.8, r'$mu=100, sigma=15$')
       plt.text(0.65, -2.3, 'Running: t = {:5.1f} hour'.format(ii * 1.98 * 2),fontsize=26)


       #cbar=fig.colorbar(pcm, orientation='horizontal', shrink=0.8, aspect=25,\
       #      ticks= ticks, pad=0.1, fraction=0.03 )
       #cbar.ax.set_title('$R^2 N(AU^2cm^{-3})$', fontsize=18)

       # plt.show()
       # break
       plt.savefig(root_dir+'CME-Density-{:03d}.png'.format(ii))


       plt.close(fig)
       
       



