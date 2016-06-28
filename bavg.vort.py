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

time = 16.0
t_str = '16'
datanum1  = '00'+t_str  #--corresponding times and data steps for orion256
datanum1a = '00'+t_str

datanum2 = ''+t_str+'00'
datanum3 = ''+t_str+'00'
datanum4 = ''+t_str+'00'

odir1  = '/work/00863/minerva/orion_Vtest/midres/'
odir1a = '/work/00863/minerva/orion_Vtest/hires/'
odir2 = '/work/00863/minerva/orion/vort/'
odir3 = '/work/00863/minerva/orion/vort_hires/'
odir4 = '/work/00863/minerva/orion/vort_anal/'

fname1  = odir1  + "data."+ datanum1  + ".3d.hdf5"
fname1a = odir1a + "data."+ datanum1a + ".3d.hdf5" 

file_level = ''

file_dens = odir2+"gadget2dens_ot_" +file_level+ datanum2
file_vx = odir2+"gadget2vx_ot_"+file_level+datanum2
file_vy = odir2+"gadget2vy_ot_"+file_level+datanum2
file_vz = odir2+"gadget2vz_ot_"+file_level+datanum2
file_bx = odir2+"gadget2bx_ot_"+file_level+datanum2
file_by = odir2+"gadget2by_ot_"+file_level+datanum2
file_bz = odir2+"gadget2bz_ot_"+file_level+datanum2

file_dens3 = odir3+"gadget2dens_ot_"+file_level+datanum3
file_vx3 = odir3+"gadget2vx_ot_"+file_level+datanum3
file_vy3 = odir3+"gadget2vy_ot_"+file_level+datanum3
file_vz3 = odir3+"gadget2vz_ot_"+file_level+datanum3
file_bx3 = odir3+"gadget2bx_ot_"+file_level+datanum3
file_by3 = odir3+"gadget2by_ot_"+file_level+datanum3
file_bz3 = odir3+"gadget2bz_ot_"+file_level+datanum3

file_dens4 = odir4+"gadget2dens_ot_"+file_level+datanum4
file_vx4 = odir4+"gadget2vx_ot_"+file_level+datanum4
file_vy4 = odir4+"gadget2vy_ot_"+file_level+datanum4
file_vz4 = odir4+"gadget2vz_ot_"+file_level+datanum4
file_bx4 = odir4+"gadget2bx_ot_"+file_level+datanum4
file_by4 = odir4+"gadget2by_ot_"+file_level+datanum4
file_bz4 = odir4+"gadget2bz_ot_"+file_level+datanum4

##################################################################################################
pf = load(fname1)

rmax = pf.domain_right_edge[0]/pc_cm

value, location = pf.h.find_max("Density")
data = pf.h.sphere(location, 1.5*rmax/pf['pc'])

def _logD(field, data):
    log_density = np.log10(data['Density'])
    return log_density
add_field("logD", function=_logD)

dims = list(pf.domain_dimensions)
ngrid_x = dims[0]
ngrid_y = dims[1]
ngrid_z = 1

#proj = pf.h.proj(2, "Density", weight_field="Density", center=pf.domain_center)
proj = pf.h.slice(2, 0.0, center=pf.domain_center)
w = (pf.h.domain_left_edge[0], pf.h.domain_right_edge[0], pf.h.domain_left_edge[1], pf.h.domain_right_edge[1])
frb1 = FixedResolutionBuffer(proj, w, (ngrid_x, ngrid_x), periodic=False)
rho = frb1["density"]
bfield = frb1['Bmag']
bx = frb1['X-magnfield']
by = frb1['Y-magnfield']
bz = frb1['Z-magnfield']
vel = frb1['z-velocity']
#####################################################################################################
pfA = load(fname1a)

rmax = pfA.domain_right_edge[0]/pc_cm

value, location = pfA.h.find_max("Density")
data = pfA.h.sphere(location, 1.5*rmax/pf['pc'])

def _logD(field, data):
    log_density = np.log10(data['Density'])
    return log_density
add_field("logD", function=_logD)

dims = list(pfA.domain_dimensions)
ngrid_xA = dims[0]
ngrid_yA = dims[1]
ngrid_zA = 1

#proj = pfA.h.proj(2, "Density", weight_field="Density", center=pf.domain_center)
proj = pfA.h.slice(2, 0.0, center=pfA.domain_center)
w = (pfA.h.domain_left_edge[0], pfA.h.domain_right_edge[0], pfA.h.domain_left_edge[1], pfA.h.domain_right_edge[1])
frb1 = FixedResolutionBuffer(proj, w, (ngrid_xA, ngrid_xA), periodic=False)
rhoA = frb1["density"]
bfieldA = frb1['Bmag']
bxA = frb1['X-magnfield']
byA = frb1['Y-magnfield']
bzA = frb1['Z-magnfield']
velA = frb1['z-velocity']
###########################################################################################################


#############################################################################################
bin_num = 100
bin_num_doub = 100.0

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
vx_gadget = np.fromfile(file_vx, dtype=float, count=-1)
vx_gadget = numpy.reshape(vx_gadget, dim_arr)
vy_gadget = np.fromfile(file_vy, dtype=float, count=-1)
vy_gadget = numpy.reshape(vy_gadget, dim_arr)
vz_gadget = np.fromfile(file_vz, dtype=float, count=-1)
vz_gadget = numpy.reshape(vz_gadget, dim_arr)
v_gadget = np.sqrt(vx_gadget*vx_gadget + vy_gadget*vy_gadget + vz_gadget*vz_gadget)


rho3 = np.fromfile(file_dens3, dtype=float, count=-1)
rho3 = numpy.reshape(rho3, dim_arr)
Bx3 = np.fromfile(file_bx3, dtype=float, count=-1)
Bx3 = numpy.reshape(Bx3, dim_arr)
By3 = np.fromfile(file_by3, dtype=float, count=-1)
By3 = numpy.reshape(By3, dim_arr)
Bz3 = np.fromfile(file_bz3, dtype=float, count=-1)
Bz3 = numpy.reshape(Bz3, dim_arr)
B3 = np.sqrt(Bx3*Bx3 + By3*By3 + Bz3*Bz3)
vx3 = np.fromfile(file_vx3, dtype=float, count=-1)
vx3 = numpy.reshape(vx3, dim_arr)
vy3 = np.fromfile(file_vy3, dtype=float, count=-1)
vy3 = numpy.reshape(vy3, dim_arr)
vz3 = np.fromfile(file_vz3, dtype=float, count=-1)
vz3 = numpy.reshape(vz3, dim_arr)
v3 = np.sqrt(vx3*vx3 + vy3*vy3 + vz3*vz3)

rho4 = np.fromfile(file_dens4, dtype=float, count=-1)
rho4 = numpy.reshape(rho4, dim_arr)
Bx4 = np.fromfile(file_bx4, dtype=float, count=-1)
Bx4 = numpy.reshape(Bx4, dim_arr)
By4 = np.fromfile(file_by4, dtype=float, count=-1)
By4 = numpy.reshape(By4, dim_arr)
Bz4 = np.fromfile(file_bz4, dtype=float, count=-1)
Bz4 = numpy.reshape(Bz4, dim_arr)
B4 = np.sqrt(Bx4*Bx4 + By4*By4 + Bz4*Bz4)
vx4 = np.fromfile(file_vx4, dtype=float, count=-1)
vx4 = numpy.reshape(vx4, dim_arr)
vy4 = np.fromfile(file_vy4, dtype=float, count=-1)
vy4 = numpy.reshape(vy4, dim_arr)
vz4 = np.fromfile(file_vz4, dtype=float, count=-1)
vz4 = numpy.reshape(vz4, dim_arr)
v4 = np.sqrt(vx4*vx4 + vy4*vy4 + vz4*vz4)


#############################################################################################
#create bin and calculate probability distribution function
print 'start the binning'

array1 = rho
array1A = rhoA
array2 = rho_gadget
array3 = rho3
array4 = rho4

do_bfield = 1
if do_bfield == 1:
#        array1 = np.log10(bfield)
#	array1A = np.log10(bfieldA)
#	array2 = np.log10(B_gadget)
#	array3 = np.log10(B3)
#        array4 = np.log10(B4)
       array1 = np.array(bfield)
       array1A = np.array(bfieldA)
       array2 = np.array(B_gadget)
       array3 = np.array(B3)
       array4 = np.array(B4)

do_bx = 0
if do_bx == 1:
        array1 = np.log10(np.abs(bx))
        array2 = np.log10(np.abs(Bx_gadget))
        array3 = np.log10(np.abs(Bx3))
        array4 = np.log10(np.abs(Bx4))

do_by = 0
if do_by == 1:
        array1 = np.log10(np.abs(by))
        array2 = np.log10(np.abs(By_gadget))
        array3 = np.log10(np.abs(By3))
        array4 = np.log10(np.abs(By4))

do_vel = 0
if do_vel == 1:
        array1 = vel*1.e5
        array2 = v_gadget
        array3 = v3
        array4 = v4

print 'array1.shape =', array1.shape
print 'array2.shape =', array2.shape

arrmin1 = np.amin(array1)
arrmin1A = np.amin(array1A)
arrmin2 = np.amin(array2)
arrmin3 = np.amin(array3)
arrmin4 = np.amin(array4)
arrmin = np.array([time, arrmin1, arrmin1A, arrmin2, arrmin3, arrmin4])
print 'arrmin', arrmin

arrmax1 = np.amax(array1)
arrmax1A = np.amax(array1A)
arrmax2 = np.amax(array2)
arrmax3 = np.amax(array3)
arrmax4 = np.amax(array4)
arrmax = np.array([time, arrmax1, arrmax1A, arrmax2, arrmax3, arrmax4])
print 'arrmax', arrmax

array_min = np.amin(arrmin)
array_max = np.amax(arrmax)

s = str(arrmin).strip('[]')
f = open('bmin_vort', 'a')
f.write(s)
f.write("\n")
f.close

s = str(arrmax).strip('[]')
f = open('bmax_vort', 'a')
f.write(s)
f.write("\n")
f.close

#############################################################################################
#calculate the rms value!

arr1_sq = array1*array1
arr1A_sq = array1A*array1A
arr2_sq = array2*array2
arr3_sq = array3*array3
arr4_sq = array4*array4

mean1 = np.mean(arr1_sq)
mean1A = np.mean(arr1A_sq)
mean2 = np.mean(arr2_sq)
mean3 = np.mean(arr3_sq)
mean4 = np.mean(arr4_sq)

rms1 = np.sqrt(mean1)
rms1A = np.sqrt(mean1A)
rms2 = np.sqrt(mean2)
rms3 = np.sqrt(mean3)
rms4 = np.sqrt(mean4)

rms = np.array([time, rms1, rms1A, rms2, rms3, rms4])
print 'rms =', rms

s = str(rms).strip('[]')
f = open('brms_vort', 'a')
f.write(s)
f.write("\n")
f.close

