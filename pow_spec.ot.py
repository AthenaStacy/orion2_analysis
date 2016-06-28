import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as plt
import numpy as np
from yt.mods import *

import fields
import fields_bfield
import fft

"""
Make a power spectrum.  
"""

pc_cm = 3.08567758e18
au_cm = 1.49597871e13
Msun = 1.98892e33
rho_nh = 1./(1.2195*1.6726e-24)

shrink_fac = 1.0
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


odir1 = '/work/00863/minerva/orion_Otest/lowBhighRes/'
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


pf = load(fname1)

# a FFT operates on uniformly gridded data.  We'll use the yt
# covering grid for this.
  
low = pf.h.domain_left_edge*shrink_fac
dims = list(pf.domain_dimensions) 

dims[0] = dims[0] * (my_level + 1)
dims[1] = dims[1] * (my_level + 1)
dims[2] = dims[2] * (my_level + 1)

ngrid_x = 1024
ngrid_y = 1024
ngrid_z = 1

#nx, ny, nz = dims
nx = 512
ny = 512
nz = 1
nx2 = nx
ny2 = ny
nz2 = nz
nindex_rho = 1./3.

nz = 0
nz2 = 0
print 'nx =', nx, 'ny = ', ny, 'nz = ', nz

Kk  = np.zeros( (nx/2+1, ny/2+1))
Kk2 = np.zeros( (nx2/2+1, ny2/2+1))
Kk3 = np.zeros( (nx2/2+1, ny2/2+1))
Kk4 = np.zeros( (nx2/2+1, ny2/2+1))

#cube = pf.h.covering_grid(level=my_level,left_edge=pf.h.domain_left_edge*shrink_fac,dims=[ngrid_x,ngrid_y,ngrid_z])
#rho = cube['density']
#x = cube['x']
#y = cube['y']
#z = cube['z']
#bfield =  cube['Bmag']

proj = pf.h.proj(2, "Density", weight_field="Density", center=pf.domain_center)
#proj = pf.h.slice(2, 0.0, center=pf.domain_center)
w = (pf.h.domain_left_edge[0], pf.h.domain_right_edge[0], pf.h.domain_left_edge[1], pf.h.domain_right_edge[1])
frb1 = FixedResolutionBuffer(proj, w, (ngrid_x, ngrid_x), periodic=True)
rho = frb1["density"]
bfield = frb1['Bmag']
vel = frb1['velocity']

cell_volume = ( (pf.h.domain_right_edge[0] - pf.h.domain_left_edge[0]) / ngrid_x)
cell_volume = pow(cell_volume, 2)
mass = rho * cell_volume 
KE = 0.5 * mass * pow(vel,2)

print "KE = ", np.sum(KE)
print "domain_center =", pf.domain_center
print "rho = ", rho
print "Bmag = ", bfield

#############################################################################################
dim_arr = [nx,ny]

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
do_bfield = 1
do_dens = 0
do_vel = 1

if do_bfield == 1:
	array1 = bfield
	array2 = B_gadget
	array3 = B3
	array4 = B4

if do_dens == 1:
	array1 = rho
	array2 = rho_gadget
	array3 = rho3
	array4 = rho4

if do_vel == 1:
        array1 = vel*1.e5
        array2 = v_gadget
        array3 = v3
        array4 = v4
#############################################################################################
#calculate the rms value!
level=0
lmax = (level+1) * dims[0] - 1

arr1_sq = array1*array1
arr2_sq = array2*array2
arr3_sq = array3*array3
arr4_sq = array4*array4

mean1 = np.mean(arr1_sq)
mean2 = np.mean(arr2_sq)
mean3 = np.mean(arr3_sq)
mean4 = np.mean(arr4_sq)

rms1 = np.sqrt(mean1)
rms2 = np.sqrt(mean2)
rms3 = np.sqrt(mean3)
rms4 = np.sqrt(mean4)

#for i in range(lmax):
#        for j in range(lmax):
#                for k in range(ngrid_z):
			#print 'i = ', i, 'j = ', j, 'k =', k, 'array2 = ', array2[i][j][k] 
print 'max1 = ', np.amax(array1)
print 'rms1 = ', rms1
print 'rms2 = ', rms2
print 'rms3 = ', rms3
print 'rms4 = ', rms4
##############################################################################################

#nx = ngrid_x
#ny = ngrid_y
nz = 1
print 'nx = ', nx, 'ny =', ny, 'nz =', nz
ru = np.fft.fftn(array1)[0:nx/2+1,0:ny/2+1]
ru = 8.0*ru/(nx*ny*nz)
Kk += np.abs(ru)**2

#nx2 = ngrid_x
#ny2 = ngrid_y
nz2 = 1
print 'nx2 = ', nx2, 'ny2 =', ny2, 'nz2 =', nz2
ru = np.fft.fftn(array2)[0:nx2/2+1,0:ny2/2+1]
ru = 8.0*ru/(nx2*ny2*nz2)
Kk2 += np.abs(ru)**2

ru = np.fft.fftn(array3)[0:nx2/2+1,0:ny2/2+1]
ru = 8.0*ru/(nx2*ny2*nz2)
Kk3 += np.abs(ru)**2

ru = np.fft.fftn(array4)[0:nx2/2+1,0:ny2/2+1]
ru = 8.0*ru/(nx2*ny2*nz2)
Kk4 += np.abs(ru)**2

# wavenumbers
L = (pf.h.domain_right_edge - pf.h.domain_left_edge) / pc_cm
print 'L[0] = ', L[0], 'L[1] =', L[1], 'L[2] =', L[2] 

kx = fft.rfftfreq(nx)*nx/L[0]
ky = fft.rfftfreq(ny)*ny/L[1]
kz = fft.rfftfreq(nz)*nz/L[2]

kx2 = fft.rfftfreq(nx2)*nx2/L[0]
ky2 = fft.rfftfreq(ny2)*ny2/L[1]
kz2 = fft.rfftfreq(nz2)*nz2/L[2]
    
# physical limits to the wavenumbers
dims_d = float(ngrid_x)
kmin = np.min(1.0/L[0])
kmax = np.max(0.5*dims_d/L[0])
print 'kmin =', kmin, 'kmax = ', kmax
 
kbins = np.arange(kmin, kmax, kmin)
N = len(kbins)

print 'N =', N 
print 'kbins =', kbins

# bin the Fourier KE into radial kbins
kx3d, ky3d = np.meshgrid(kx, ky, indexing="ij")
k = np.sqrt(kx3d**2 + ky3d**2)
print 'k =', k

whichbin = np.digitize(k.flat, kbins)
ncount = np.bincount(whichbin)    

E_spectrum  = np.zeros(len(ncount)-1)
E2_spectrum = np.zeros(len(ncount)-1)
E3_spectrum = np.zeros(len(ncount)-1)
E4_spectrum = np.zeros(len(ncount)-1)

for n in range(1,len(ncount)):
	E_spectrum[n-1]  = np.sum(Kk.flat[whichbin==n])
	E2_spectrum[n-1] = np.sum(Kk2.flat[whichbin==n])
	E3_spectrum[n-1] = np.sum(Kk3.flat[whichbin==n])
	E4_spectrum[n-1] = np.sum(Kk4.flat[whichbin==n])

k = 0.5*(kbins[0:N-1] + kbins[1:N])


E_spectrum  = E_spectrum[1:N]
E2_spectrum = E2_spectrum[1:N]
E3_spectrum = E3_spectrum[1:N]
E4_spectrum = E4_spectrum[1:N]

index = np.argmax(E_spectrum)
kmax = k[index]
Emax = E_spectrum[index]

print 'len(k) =', len(k), 'k =', k
print 'len(E3) =', len(E3_spectrum), 'E3_spectrum =', E3_spectrum

plt.loglog(k, E_spectrum, 'k')
plt.loglog(k, E2_spectrum)
plt.loglog(k, E3_spectrum, 'r')
plt.loglog(k, E4_spectrum, 'g')
#plt.loglog(k, Emax*(k/kmax)**(-5./3.), ls=":", color="0.5")
#plt.axis((1.e0, 1.e2, 1.e-38, 1.e-35))

plt.xlabel(r"$k$")
plt.ylabel(r"$E(k)dk$")

plt.savefig("spectrum.eps")

