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

#shrink_fac = 0.125
#my_level = 3

shrink_fac = 0.0625
my_level = 4

#datanum1 = '0000'
#datanum2 = '0000'
#datanum1 = '0020'
#datanum2 = '0746'
#datanum1 = '0023'
#datanum2 = '0827'

datanum1 = '0029'  #--corresponding times and data steps for orion256
datanum2 = '0679'

odir1 = '/work/00863/minerva/orion/bfield_comp_orion256/'
odir2 = '/work/00863/minerva/orion/'
odir3 = '/work/00863/minerva/orion//bfield_comp_vpot/'

fname1 = odir1+"data."+ datanum1 +".3d.hdf5"

file_level = 'lev4_'

file_dens = odir2+"gadget2dens_ideal_" +file_level+ datanum2
file_vx = odir2+"gadget2vx_ideal_"+file_level+datanum2
file_vy = odir2+"gadget2vy_ideal_"+file_level+datanum2
file_vz = odir2+"gadget2vz_ideal_"+file_level+datanum2
file_bx = odir2+"gadget2bx_ideal_"+file_level+datanum2
file_by = odir2+"gadget2by_ideal_"+file_level+datanum2
file_bz = odir2+"gadget2bz_ideal_"+file_level+datanum2 

file_bx_vp = odir3+"gadget2bx_ideal_"+file_level+datanum2
file_by_vp = odir3+"gadget2by_ideal_"+file_level+datanum2
file_bz_vp = odir3+"gadget2bz_ideal_"+file_level+datanum2


pf = load(fname1)

bin_num = 100
bin_num_doub = 100.0
gref = 256
oref = gref

rmax = pf.domain_right_edge[0]/pc_cm

value, location = pf.h.find_max("Density")
data = pf.h.sphere(location, 1.5*rmax/pf['pc'])

def _logD(field, data):
    log_density = np.log10(data['Density'])
    return log_density
add_field("logD", function=_logD)

dims = list(pf.domain_dimensions)
cube = pf.h.covering_grid(level=my_level,left_edge=pf.h.domain_left_edge*shrink_fac,dims=[oref,oref,oref])
dens = cube['logD']
bfield =  cube['Bmag']
vx =      cube['x-velocity'] / 1.e5
vy =      cube['y-velocity'] / 1.e5
vz =      cube['z-velocity'] / 1.e5

#############################################################################################
dens_gadget = np.fromfile(file_dens, dtype=float, count=-1)
dens_gadget = numpy.reshape(dens_gadget, [gref,gref,gref])
dens_gadget = np.log10(dens_gadget)

vx_gadget = np.fromfile(file_vx, dtype=float, count=-1)
vx_gadget = numpy.reshape(vx_gadget, [gref,gref,gref])
vx_gadget = vx_gadget / 1.e5

vy_gadget = np.fromfile(file_vy, dtype=float, count=-1)
vy_gadget = numpy.reshape(vy_gadget, [gref,gref,gref])
vy_gadget = vy_gadget / 1.e5

vz_gadget = np.fromfile(file_vz, dtype=float, count=-1)
vz_gadget = numpy.reshape(vz_gadget, [gref,gref,gref])
vz_gadget = vz_gadget / 1.e5

Bx_gadget = np.fromfile(file_bx, dtype=float, count=-1)
Bx_gadget = numpy.reshape(Bx_gadget, [gref,gref,gref])

By_gadget = np.fromfile(file_by, dtype=float, count=-1)
By_gadget = numpy.reshape(By_gadget, [gref,gref,gref])

Bz_gadget = np.fromfile(file_bz, dtype=float, count=-1)
Bz_gadget = numpy.reshape(Bz_gadget, [gref,gref,gref])

B_gadget = np.sqrt(Bx_gadget*Bx_gadget + By_gadget*By_gadget + Bz_gadget*Bz_gadget)

Bx_vp = np.fromfile(file_bx_vp, dtype=float, count=-1)
Bx_vp = numpy.reshape(Bx_vp, [gref,gref,gref])

By_vp = np.fromfile(file_by_vp, dtype=float, count=-1)
By_vp = numpy.reshape(By_vp, [gref,gref,gref])

Bz_vp = np.fromfile(file_bz, dtype=float, count=-1)
Bz_vp = numpy.reshape(Bz_vp, [gref,gref,gref])

B_vp = np.sqrt(Bx_vp*Bx_vp + By_vp*By_vp + Bz_vp*Bz_vp)

dens2_min = np.amin(dens_gadget)
dens2_max = np.amax(dens_gadget)


#############################################################################################
#create bin and calculate probability distribution function
print 'start the binning'

array1 = dens
array2 = dens_gadget
array3 = dens_gadget

do_bfield = 1
if do_bfield == 1:
	array1 = bfield
	array2 = B_gadget
	array3 = B_vp

print 'array1.shape =', array1.shape
print 'array2.shape =', array2.shape

array_min = np.amin(array1)
array_max = np.amax(array1)

array_bin = []
for i in range(bin_num):
        i_doub = float(i)
        array_bin.append(array_min + (array_max-array_min)*i_doub/bin_num_doub)

#level=0
#lmax = (level+1) * dims[0] - 1
lmax = oref 

num =  np.zeros(bin_num)
num2 = np.zeros(bin_num)
num3 = np.zeros(bin_num)

for i in range(lmax):
	for j in range(lmax):
		for k in range(lmax):
			for l in range(bin_num-1):
				if array1[i,j,k] > array_bin[l] and array1[i,j,k] < array_bin[l+1]:
					num[l] = num[l] + 1			
                                if array2[i,j,k] > array_bin[l] and array2[i,j,k] < array_bin[l+1]:
                                        num2[l] = num2[l] + 1
                                if array3[i,j,k] > array_bin[l] and array3[i,j,k] < array_bin[l+1]:
                                        num3[l] = num3[l] + 1
				if i%10 == 0 and j%10 == 0 and k == 0:
					print 'i =', i, 'j =', j, 'k =', k

################################################################################################

num_min = np.amin(num)
num_max = np.amax(num)
num2_min = np.amin(num2)
num2_max = np.amax(num2)

print 'num_min =', num_min, 'num_max =', num_max
print 'num2_min =', num_min, 'num2_max =', num2_max
print 'array_min =', array_min, 'array_max =', array_max


print 'total =', np.sum(num)
num  = num  / np.sum(num)
num2 = num2 / np.sum(num2)
num3 = num3 / np.sum(num3)

plt.title('Black=Orion, Blue=Roberts, Red=Vec. Pot.', fontsize = 10)
plt.plot(array_bin, num, 'k')
plt.plot(array_bin, num2)
plt.plot(array_bin, num3, 'r')
ax = plt.gca()
ax.set_yscale('log')
ax.set_xlabel('log density[g cm^-3]', fontsize=9)
ax.set_ylabel('N', fontsize=9)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.axis((array_min, array_max, 1.e-8, 1))

plt.savefig("pdf.eps")





