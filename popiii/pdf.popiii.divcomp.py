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

my_level = 0

datanum1 = '0029'  #--corresponding times and data steps for orion256
datanum2 = '0445'

odir1 = '/global/scratch/minerva/orion/bfield_comp_orion256/'
odir2 = '/global/scratch/minerva/orion/'
odir3 = '/global/scratch/minerva/orion/'
odir4 = '/global/scratch/minerva/orion/'
odir5 = '/global/scratch/minerva/orion/'
odir6 =  '/global/scratch/minerva/popiii_Bscope2/'

fname1 = odir1+"data."+ datanum1 +".3d.hdf5"

file_level = ''

file_dens = odir2+"gadget2dens_ideal_" +file_level+ datanum2
file_vx = odir2+"gadget2vx_ideal_"+file_level+datanum2
file_vy = odir2+"gadget2vy_ideal_"+file_level+datanum2
file_vz = odir2+"gadget2vz_ideal_"+file_level+datanum2
file_bx = odir2+"gadget2bx_ideal_"+file_level+datanum2
file_by = odir2+"gadget2by_ideal_"+file_level+datanum2
file_bz = odir2+"gadget2bz_ideal_"+file_level+datanum2 

file_bx_ns = odir3+"gadget2bx_nosmooth_"+file_level+datanum2
file_by_ns = odir3+"gadget2by_nosmooth_"+file_level+datanum2
file_bz_ns = odir3+"gadget2bz_nosmooth_"+file_level+datanum2

file_bx_ls = odir4+"gadget2bx_lowsmooth_"+file_level+datanum2
file_by_ls = odir4+"gadget2by_lowsmooth_"+file_level+datanum2
file_bz_ls = odir4+"gadget2bz_lowsmooth_"+file_level+datanum2

file_bx_hs = odir5+"gadget2bx_hismooth_"+file_level+datanum2
file_by_hs = odir5+"gadget2by_hismooth_"+file_level+datanum2
file_bz_hs = odir5+"gadget2bz_hismooth_"+file_level+datanum2

file_bx_cl = odir6+"gadget2bx_clean_"+file_level+datanum2
file_by_cl = odir6+"gadget2by_clean_"+file_level+datanum2
file_bz_cl = odir6+"gadget2bz_clean_"+file_level+datanum2

bin_num = 100
bin_num_doub = 100.0

extra=0
gref = 128+extra
oref = 128+extra

######################################################################################
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

Bx_ns = np.fromfile(file_bx_ns, dtype=float, count=-1)
Bx_ns = numpy.reshape(Bx_ns, [gref,gref,gref])

By_ns = np.fromfile(file_by_ns, dtype=float, count=-1)
By_ns = numpy.reshape(By_ns, [gref,gref,gref])

Bz_ns = np.fromfile(file_bz_ns, dtype=float, count=-1)
Bz_ns = numpy.reshape(Bz_ns, [gref,gref,gref])

B_ns = np.sqrt(Bx_ns*Bx_ns + By_ns*By_ns + Bz_ns*Bz_ns)

Bx_ls = np.fromfile(file_bx_ls, dtype=float, count=-1)
Bx_ls = numpy.reshape(Bx_ls, [gref,gref,gref])

By_ls = np.fromfile(file_by_ls, dtype=float, count=-1)
By_ls = numpy.reshape(By_ls, [gref,gref,gref])

Bz_ls = np.fromfile(file_bz_ls, dtype=float, count=-1)
Bz_ls = numpy.reshape(Bz_ls, [gref,gref,gref])

B_ls = np.sqrt(Bx_ls*Bx_ls + By_ls*By_ls + Bz_ls*Bz_ls)

Bx_hs = np.fromfile(file_bx_hs, dtype=float, count=-1)
Bx_hs = numpy.reshape(Bx_hs, [gref,gref,gref])

By_hs = np.fromfile(file_by_hs, dtype=float, count=-1)
By_hs = numpy.reshape(By_hs, [gref,gref,gref])

Bz_hs = np.fromfile(file_bz_hs, dtype=float, count=-1)
Bz_hs = numpy.reshape(Bz_hs, [gref,gref,gref])

B_hs = np.sqrt(Bx_hs*Bx_hs + By_hs*By_hs + Bz_hs*Bz_hs)

Bx_cl = np.fromfile(file_bx_cl, dtype=float, count=-1)
Bx_cl = numpy.reshape(Bx_cl, [gref,gref,gref])

By_cl = np.fromfile(file_by_cl, dtype=float, count=-1)
By_cl = numpy.reshape(By_cl, [gref,gref,gref])

Bz_cl = np.fromfile(file_bz_cl, dtype=float, count=-1)
Bz_cl = numpy.reshape(Bz_cl, [gref,gref,gref])

B_cl = np.sqrt(Bx_cl*Bx_cl + By_cl*By_cl + Bz_cl*Bz_cl)


dens2_min = np.amin(dens_gadget)
dens2_max = np.amax(dens_gadget)


#############################################################################################
#create bin and calculate probability distribution function
print 'start the binning'

array1 = dens_gadget
array2 = dens_gadget
array3 = dens_gadget

do_bfield = 1
if do_bfield == 1:
	array1 = np.log10(B_gadget)
	array2 = np.log10(B_gadget)
	array3 = np.log10(B_ns)
        array4 = np.log10(B_ls)
        array5 = np.log10(B_hs)
        array6 = np.log10(B_cl)

print 'array1.shape =', array1.shape
print 'array2.shape =', array2.shape

array_min = np.amin(array1)
array_max = np.amax(array1)

array_bin = []
for i in range(bin_num):
        i_doub = float(i)
        array_bin.append(array_min + (array_max-array_min)*i_doub/bin_num_doub)

print 'array_bin =', array_bin

#level=0
#lmax = (level+1) * dims[0] - 1
lmax = oref-1 

print 'lmax =', lmax

num =  np.zeros(bin_num)
num2 = np.zeros(bin_num)
num3 = np.zeros(bin_num)
num4 = np.zeros(bin_num)
num5 = np.zeros(bin_num)
num6 = np.zeros(bin_num)

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
                                if array4[i,j,k] > array_bin[l] and array4[i,j,k] < array_bin[l+1]:
                                        num4[l] = num4[l] + 1
                                if array5[i,j,k] > array_bin[l] and array5[i,j,k] < array_bin[l+1]:
                                        num5[l] = num5[l] + 1
                                if array6[i,j,k] > array_bin[l] and array6[i,j,k] < array_bin[l+1]:
                                        num6[l] = num6[l] + 1
				if i%10 == 0 and j%10 == 0 and k%10 == 0:
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
num4 = num4 / np.sum(num4)
num5 = num5 / np.sum(num5)
num6 = num6 / np.sum(num6)

plt.title('Red = no smoothing, Blue = infrequent smoothing, Green = frequent smoothing, Magenta = v. frequent smoothing', fontsize = 10)
plt.plot(array_bin, num, 'k')
plt.plot(array_bin, num2)
plt.plot(array_bin, num3, 'r')
plt.plot(array_bin, num4, 'g')
plt.plot(array_bin, num5, 'm')
plt.plot(array_bin, num6, 'b--')
ax = plt.gca()
ax.set_yscale('log')
ax.set_xlabel('log B [Gauss]', fontsize=9)
ax.set_ylabel('N', fontsize=9)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.axis((array_min, array_max, 1.e-8, 1))

plt.savefig("pdf.eps")





