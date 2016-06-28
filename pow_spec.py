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

#shrink_fac = 1
#my_level = 0
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

# a FFT operates on uniformly gridded data.  We'll use the yt
# covering grid for this.
  
low = pf.h.domain_left_edge*shrink_fac
dims = list(pf.domain_dimensions) 

dims[0] = dims[0] * (my_level + 1)
dims[1] = dims[1] * (my_level + 1)
dims[2] = dims[2] * (my_level + 1)

gref = 256
oref = gref

nx, ny, nz = dims
nx2 = gref
ny2 = gref
nz2 = gref
nindex_rho = 1./3.

Kk  = np.zeros( (nx/2+1, ny/2+1, nz/2+1))
Kk2 = np.zeros( (nx2/2+1, ny2/2+1, nz2/2+1))

cube = pf.h.covering_grid(level=my_level,left_edge=pf.h.domain_left_edge*shrink_fac,dims=[oref,oref,oref])
rho = cube['density']
x = cube['x']
y = cube['y']
z = cube['z']
bfield =  cube['Bmag']

#############################################################################################
rho_gadget = np.fromfile(file_dens, dtype=float, count=-1)
rho_gadget = numpy.reshape(rho_gadget, [gref,gref,gref])

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
#############################################################################################

"""
level=0
lmax = (level+1) * dims[0] - 1

for i in range(lmax):
        for j in range(lmax):
                for k in range(lmax):
 			rho[i,j,k] = i*j*k
			rho[i,j,k] = np.sqrt(x[i,j,k]*x[i,j,k] + y[i,j,k]*y[i,j,k] + z[i,j,k]*z[i,j,k])
			rho[i,j,k] = np.random.random()
#rho[10,10,10] = 1.e10
"""

array1 = rho
array2 = rho_gadget
array3 = B_vp

nx, ny, nz = array1.shape
print 'nx = ', nx, 'ny =', ny, 'nz =', nz
ru = np.fft.fftn(array1)[0:nx/2+1,0:ny/2+1,0:nz/2+1]
ru = 8.0*ru/(nx*ny*nz)
Kk += np.abs(ru)**2

nx2, ny2, nz2 = array2.shape
print 'nx2 = ', nx2, 'ny2 =', ny2, 'nz2 =', nz2
ru = np.fft.fftn(array2)[0:nx2/2+1,0:ny2/2+1,0:nz2/2+1]
ru = 8.0*ru/(nx2*ny2*nz2)
Kk2 += np.abs(ru)**2

# wavenumbers
L = (pf.h.domain_right_edge - pf.h.domain_left_edge) / pc_cm
print 'L[0] = ', L[0], 'L[1] =', L[1], 'L[2] =', L[2], 'nx = ', nx, 'ny = ', ny, 'nz =', nz 

kx = fft.rfftfreq(nx)*nx/L[0]
ky = fft.rfftfreq(ny)*ny/L[1]
kz = fft.rfftfreq(nz)*nz/L[2]

kx2 = fft.rfftfreq(nx2)*nx2/L[0]
ky2 = fft.rfftfreq(ny2)*ny2/L[1]
kz2 = fft.rfftfreq(nz2)*nz2/L[2]
    
# physical limits to the wavenumbers
dims_d = float(dims[0])
kmin = np.min(1.0/L)
kmax = np.max(0.5*dims_d/L)
    
kbins = np.arange(kmin, kmax, kmin)
N = len(kbins)

# bin the Fourier KE into radial kbins
kx3d, ky3d, kz3d = np.meshgrid(kx, ky, kz, indexing="ij")
k = np.sqrt(kx3d**2 + ky3d**2 + kz3d**2)

whichbin = np.digitize(k.flat, kbins)
ncount = np.bincount(whichbin)    
E_spectrum  = np.zeros(len(ncount)-1)

kx3d, ky3d, kz3d = np.meshgrid(kx2, ky2, kz2, indexing="ij")
k = np.sqrt(kx3d**2 + ky3d**2 + kz3d**2)

whichbin = np.digitize(k.flat, kbins)
ncount = np.bincount(whichbin)
E2_spectrum = np.zeros(len(ncount)-1)

for n in range(1,len(ncount)):
	E_spectrum[n-1]  = np.sum(Kk.flat[whichbin==n])
	E2_spectrum[n-1] = np.sum(Kk2.flat[whichbin==n])

k = 0.5*(kbins[0:N-1] + kbins[1:N])


E_spectrum  = E_spectrum[1:N]
E2_spectrum = E2_spectrum[1:N]

index = np.argmax(E_spectrum)
kmax = k[index]
Emax = E_spectrum[index]

plt.loglog(k, E_spectrum, 'k')
plt.loglog(k, E2_spectrum)
#plt.loglog(k, Emax*(k/kmax)**(-5./3.), ls=":", color="0.5")
#plt.axis((1.e0, 1.e2, 1.e-38, 1.e-35))

plt.xlabel(r"$k$")
plt.ylabel(r"$E(k)dk$")

plt.savefig("spectrum.eps")

