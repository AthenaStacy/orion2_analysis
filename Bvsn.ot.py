import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as plt
import numpy as np
from yt.mods import *

import fields
import fields_bfield
import fft

pc_cm = 3.08567758e18
au_cm = 1.49597871e13
Msun = 1.98892e33
rho_nh = 1./(1.2195*1.6726e-24)

shrink_fac = 0.0625
my_level = 0

#datanum1 = '0010'  #--corresponding times and data steps for orion256
#datanum2 = '0500'
#datanum3 = '1100'
#datanum4 = '1700'

#datanum1 = '0020'  #--corresponding times and data steps for orion256
#datanum2 = '1100'
#datanum3 = '2300'
#datanum4 = '3300'


datanum1 = '0050'  #--corresponding times and data steps for orion256
datanum2 = '2266'
datanum3 = '4712'
datanum4 = '7097'

odir1 = '/work/00863/minerva/orion_Otest/lowBhighRes2/'
odir2 = '/work/00863/minerva/orion/ot_twodim/'
odir3 = '/work/00863/minerva/orion/ot_twodim_midres/'
odir4 = '/work/00863/minerva/orion/ot_twodim_highres/'

fname1 = odir1+"data."+ datanum1 +".3d.hdf5"

file_level = ''

file_dens = odir2+"gadget2dens_ot_" +file_level+ datanum2
file_vx = odir2+"gadget2vx_ot_"+file_level+datanum2
file_vy = odir2+"gadget2vy_ot_"+file_level+datanum2
file_vz = odir2+"gadget2vz_ot_"+file_level+datanum2
file_bx = odir2+"gadget2bx_ot_"+file_level+datanum2
file_by = odir2+"gadget2by_ot_"+file_level+datanum2
file_bz = odir2+"gadget2bz_ot_"+file_level+datanum2

file_dens3 = odir3+"gadget2dens_ot_"+file_level+datanum3
file_bx3 = odir3+"gadget2bx_ot_"+file_level+datanum3
file_by3 = odir3+"gadget2by_ot_"+file_level+datanum3
file_bz3 = odir3+"gadget2bz_ot_"+file_level+datanum3

file_dens4 = odir4+"gadget2dens_ot_"+file_level+datanum4
file_bx4 = odir4+"gadget2bx_ot_"+file_level+datanum4
file_by4 = odir4+"gadget2by_ot_"+file_level+datanum4
file_bz4 = odir4+"gadget2bz_ot_"+file_level+datanum4

pf = load(fname1)

bin_num = 100
bin_num_doub = 100.0

rmax = pf.domain_right_edge[0]/pc_cm

value, location = pf.h.find_max("Density")
data = pf.h.sphere(location, 1.5*rmax/pf['pc'])

def _logD(field, data):
    log_density = np.log10(data['Density'])
    return log_density
add_field("logD", function=_logD)

dims = list(pf.domain_dimensions)
ngrid_x = 1024
ngrid_y = 1024
ngrid_z = 1

proj = pf.h.proj(2, "Density", weight_field="Density", center=pf.domain_center)
#proj = pf.h.slice(2, 0.0, center=pf.domain_center)
w = (pf.h.domain_left_edge[0], pf.h.domain_right_edge[0], pf.h.domain_left_edge[1], pf.h.domain_right_edge[1])
frb1 = FixedResolutionBuffer(proj, w, (ngrid_x, ngrid_x), periodic=True)
rho = frb1["density"]
bfield = frb1['Bmag']
bx = frb1['X-magnfield']
by = frb1['Y-magnfield']
bz = frb1['Z-magnfield']

#############################################################################################
ngrid_x_gadget = 512
ngrid_y_gadget = 512
ngrid_z_gadget = 1
dim_arr = [ngrid_x_gadget,ngrid_y_gadget]

rho_gadget = np.fromfile(file_dens, dtype=float, count=-1)
rho_gadget = numpy.reshape(rho_gadget, dim_arr)

Bx_gadget = np.fromfile(file_bx, dtype=float, count=-1)
Bx_gadget = numpy.reshape(Bx_gadget, dim_arr)

By_gadget = np.fromfile(file_by, dtype=float, count=-1)
By_gadget = numpy.reshape(By_gadget, dim_arr)

Bz_gadget = np.fromfile(file_bz, dtype=float, count=-1)
Bz_gadget = numpy.reshape(Bz_gadget, dim_arr)

B_gadget = np.sqrt(Bx_gadget*Bx_gadget + By_gadget*By_gadget + Bz_gadget*Bz_gadget)

rho3 = np.fromfile(file_dens3, dtype=float, count=-1)
rho3 = numpy.reshape(rho3, dim_arr)

Bx3 = np.fromfile(file_bx3, dtype=float, count=-1)
Bx3 = numpy.reshape(Bx3, dim_arr)

By3 = np.fromfile(file_by3, dtype=float, count=-1)
By3 = numpy.reshape(By3, dim_arr)

Bz3 = np.fromfile(file_bz3, dtype=float, count=-1)
Bz3 = numpy.reshape(Bz3, dim_arr)

B3 = np.sqrt(Bx3*Bx3 + By3*By3 + Bz3*Bz3)

rho4 = np.fromfile(file_dens4, dtype=float, count=-1)
rho4 = numpy.reshape(rho4, dim_arr)

Bx4 = np.fromfile(file_bx4, dtype=float, count=-1)
Bx4 = numpy.reshape(Bx4, dim_arr)

By4 = np.fromfile(file_by4, dtype=float, count=-1)
By4 = numpy.reshape(By4, dim_arr)

Bz4 = np.fromfile(file_bz4, dtype=float, count=-1)
Bz4 = numpy.reshape(Bz4, dim_arr)

B4 = np.sqrt(Bx4*Bx4 + By4*By4 + Bz4*Bz4)


rho2_min = np.amin(rho_gadget)
rho2_max = np.amax(rho_gadget)

#############################################################################################
#create bin and calculate probability distribution function
print 'start the binning'

array1 = rho
array2 = rho_gadget
array3 = rho3
array4 = rho4

bfield = np.log10(bfield)
B_gadget = np.log10(B_gadget)
B3 = np.log10(B3)
B4 = np.log10(B4)

bx = np.log10(np.abs(bx))
Bx_gadget = np.log10(np.abs(Bx_gadget))
Bx3 = np.log10(np.abs(Bx3))
Bx4 = np.log10(np.abs(Bx4))

by = np.log10(np.abs(by))
By_gadget = np.log10(np.abs(By_gadget))
By3 = np.log10(np.abs(By3))
By4 = np.log10(np.abs(By4))

arrmin1 = np.amin(array1)
arrmin2 = np.amin(array2)
arrmin3 = np.amin(array3)
arrmin4 = np.amin(array4)
arrmin = [arrmin1, arrmin2, arrmin3, arrmin4]

arrmax1 = np.amax(array1)
arrmax2 = np.amax(array2)
arrmax3 = np.amax(array3)
arrmax4 = np.amax(array4)
arrmax = [arrmax1, arrmax2, arrmax3, arrmax4]

array_min = np.amin(arrmin)
array_max = np.amax(arrmax)

array_bin = []
for i in range(bin_num):
        i_doub = float(i)
        array_bin.append(array_min + (array_max-array_min)*i_doub/bin_num_doub)

#level=0

num =  np.zeros(bin_num)
num2 = np.zeros(bin_num)
num3 = np.zeros(bin_num)
num4 = np.zeros(bin_num)


################################################################################################

print 'array_min =', array_min, 'array_max =', array_max

plt.title('Black = Orion, Blue = SPH low res., Red = SPH med. res., Green = SPH high res', fontsize = 10)
plt.plot(rho, bfield, 'k')
plt.plot(rho_gadget, B_gadget)
plt.plot(rho3, B3, 'r')
plt.plot(rho4, B4, 'g')
ax = plt.gca()
ax.set_xlabel('density', fontsize=9)
ax.set_ylabel('log Bmag ', fontsize=9)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.axis((array_min, array_max, -22, -15))

plt.savefig("Bvsn.eps")





