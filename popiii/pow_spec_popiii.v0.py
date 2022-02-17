import yt as yt
import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as plt
import numpy as np
from yt.mods import *
import fields_simp2
#import fields_bfield
#import tracer_def
import fft
import h5py

"""
Make a power spectrum.  
"""

pc_to_cm = 3.08567758e18
au_to_cm = 1.49597871e13
Msun = 1.98892e33
rho_nh = 1./(1.2195*1.6726e-24)

odir1 = '/work/00863/minerva/stampede2/popiii_bfieldA/'
datanum1 = '0050'  #--corresponding times and data steps for orion256
fname1 = odir1+"data."+ datanum1 +".3d.hdf5"
#fname1 = '/nobackupp12/astacy/to_share/data.0004.3d.hdf5'

pf = yt.load(fname1)

for i in sorted(pf.h.field_list):
  print(i)

# a FFT operates on uniformly gridded data.  We'll use the yt
# covering grid for this.

max_level = 0
#max_level = pf.h.index.max_level

#dims = pf.h.domain_dimensions
#dims = [2048,2048,2048]
#dims = [1024,1024,1024]
#dims = [512,512,512]
#dims = [256,256,256]
dims = [128,128,128]
#dims = dims * pow(2,max_level)

fac = 1.0
left_edge = fac*pf.h.domain_left_edge
right_edge = pf.h.domain_right_edge

left_corner =  fac*left_edge
right_corner = fac*right_edge
region = pf.h.box(left_corner, right_corner)

cube = pf.smoothed_covering_grid(level=max_level, left_edge=left_edge, dims=dims, fields = ['density', 'X-momentum', 'Y-momentum', 'Z-momentum'], field_parameters = {'center':[0,0,0]})
#cube = pf.h.arbitrary_grid(left_edge = left_corner, right_edge = right_corner, dims = dims)

#cube = region

print('max level = ', pf.h.index.max_level)

'''
vturb = cube['vturb']
temp = cube['Temp']
cs = np.power(1.38e-16 * temp / 1.22 / 1.67e-24, 0.5) /1.e5
mturb = vturb / cs
print('avg mach = ', np.mean(mturb))
'''

nx, ny, nz = np.asarray(dims)

# wavenumbers
L = (pf.h.domain_right_edge - pf.h.domain_left_edge) / (pf.h.domain_right_edge[0] - pf.h.domain_left_edge[0]) 
print ('L[0] = ', L[0], 'L[1] =', L[1], 'L[2] =', L[2], 'nx = ', nx, 'ny = ', ny, 'nz =', nz) 

fac_here = 4
nx2, ny2, nz2 = int(nx/fac_here), int(ny/fac_here), int(nz/fac_here)

kx = np.fft.fftfreq(nx2)*nx2/L[0]
ky = np.fft.fftfreq(ny2)*ny2/L[1]
kz = np.fft.fftfreq(nz2)*nz2/L[2]
print ('kx shape = ', kx.shape)

print (' cube shape = ', cube['x'].shape)

'''
print('x0 = ', cube['x'].reshape((nx,nx,nx))[0,0,0], 'y0 = ', cube['y'].reshape((nx,nx,nx))[0,0,0], 'z0 = ', cube['z'].reshape((nx,nx,nx))[0,0,0])
print('x1/2 = ', cube['x'].reshape((nx,nx,nx))[nx/2,ny/2,nz/2], 'y0 = ', cube['y'].reshape((nx,nx,nx))[nx/2,nx/2,nz/2], 'z0 = ', cube['z'].reshape((nx,nx,nx))[nx/2,ny/2,nz/2])
print('x1 = ', cube['x'].reshape((nx,nx,nx))[nx-1,ny-1,nz-1], 'y0 = ', cube['y'].reshape((nx,nx,nx))[nx-1,nx-1,nz-1], 'z0 = ', cube['z'].reshape((nx,nx,nx))[nx-1,ny-1,nz-1])
'''

print('max xmom = ', np.max(cube['X-momentum'])) 
print('min xmom = ', np.min(cube['X-momentum'])) 
print('avg xmom = ', np.mean(cube['X-momentum']))
print('max ymom = ', np.max(cube['Y-momentum'])) 
print('min ymom = ', np.min(cube['Y-momentum'])) 
print('avg ymom = ', np.mean(cube['Y-momentum']))
print('max zmom = ', np.max(cube['Z-momentum'])) 
print('min zmom = ', np.min(cube['Z-momentum'])) 
print('avg zmom = ', np.mean(cube['Z-momentum']))

print('x0 = ', cube['x'][0,0,0], 'y0 = ', cube['y'][0,0,0], 'z0 = ', cube['z'][0,0,0])
print('x1/2 = ', cube['x'][int(nx/2),int(ny/2),int(nz/2)], 'y1/2 = ', cube['y'][int(nx/2),int(nx/2),int(nz/2)], 'z1/2 = ', cube['z'][int(nx/2),int(ny/2),int(nz/2)])
print('x1 = ', cube['x'][nx-1,ny-1,nz-1], 'y1 = ', cube['y'][nx-1,nx-1,nz-1], 'z1 = ', cube['z'][nx-1,ny-1,nz-1])


# physical limits to the wavenumbers
'''
if fac_here == 1:
	factor = 0.5
	l2 = 0
	l3 = nx
elif fac_here == 2:
	factor = 0.25
	l2 = int(nx/4)
	l3 = int(3*nx/4)
elif fac_here == 4:
'''
factor = 1.0 / fac_here / 2.0
l2 = int(nx/2 - nx/fac_here/2)
l3 = int(nx/2 + nx/fac_here/2)
	
kmin = np.min(1.0/L) 
kmax = np.max(factor*dims[0]/L) + 1


kbins = np.arange(kmin, kmax, kmin)
N = len(kbins)

print('N = ', N)

# bin the Fourier KE into radial kbins
kx3d, ky3d, kz3d = np.meshgrid(kx, ky, kz, indexing="ij")
k = np.sqrt(kx3d**2 + ky3d**2 + kz3d**2) 

#fields = ['X-momentum', 'Y-momentum', 'Z-momentum']
fields = ['x-velocity', 'y-velocity', 'z-velocity']
#fields = ['x_vel_norm', 'y_vel_norm', 'z_vel_norm']
#fields = ['vturb']
#fields = ['vrad']
#fields = ['Temp']

Kk  = np.zeros((l3-l2, l3-l2, l3-l2))
for i in range(len(fields)):
	array = cube[fields[i]].d #/ cube['density'].d
	#array = array.reshape((nx,nx,nx), order='A')
	array = array[l2:l3, l2:l3, l2:l3]
	array2 = array
	#array[0,0,0] = array2[64,64,64]
	#array[64,64,64] = array2[0,0,0]
	ru = np.fft.fftn(array)
	#ru = ru[0:kx,0:ky,0:kz]
	ru = np.abs(ru)
	ru = np.power((ru),2)
	ru = ru/(nx*ny*nz)/8
	Kk += ru

#Kk = np.fft.fftshift(Kk)

print('k shape = ', k.shape, 'kmin = ', kmin, 'kmax = ', kmax)
print('KK shape = ', Kk.shape)
print('kbins = ', kbins)

print('Kk_tot = ', np.sum(Kk))
print('Kk_min = ', np.min(Kk))
print('Kk_max = ', np.max(Kk))
print('Kk_mean = ', np.mean(Kk))

whichbin = np.digitize(k.flat, kbins)
ncount = np.bincount(whichbin)    
E_spectrum  = np.zeros(len(ncount)-1)
print('len(ncount)-1 ', len(ncount)-1)

#shift_num = (l3-l2)/2 -1
shift_num=0
rad_grid = np.ones((l3-l2, l3-l2, l3-l2)) * 1.e3

'''
for xx in range(l3-l2):
	for yy in range(l3-l2):
		for zz in range(l3-l2):
			rad_grid[xx][yy][zz] = np.sqrt(np.power((xx - shift_num),2) + np.power((yy - shift_num),2) + np.power((zz - shift_num),2))
'''
'''	
			k_here = np.sqrt((xx - shift_num)**2 + (yy - shift_num)**2 + (zz - shift_num)**2)
			for nn in range(0,len(ncount)-1):
				if(k_here > kbins[nn] - 0.5 and k_here < kbins[nn] + 0.5):
					E_spectrum[nn] = E_spectrum[nn] + Kk[xx][yy][zz] 
'''

'''
print('rad_grid = ', rad_grid.shape)
print('min rad_grid = ', np.amin(rad_grid), 'max rad_grid = ', np.amax(rad_grid))
'''

'''
for nn in range(0,len(ncount)-1):
	ind = np.where((rad_grid > kbins[nn] - 0.5) & (rad_grid < kbins[nn] + 0.5) ) 
	E_spectrum[nn] = np.sum(Kk[ind])
print('test A')
'''

for nn in range(1,len(ncount)):
	E_spectrum[nn-1]  = np.sum(Kk.flat[whichbin==nn])
print('test B')

E_spectrum  = E_spectrum[0:N-2]
#k = 0.5*(kbins[0:N-1] + kbins[1:N])
k = kbins[0:N-2]

index = np.argmax(E_spectrum) 
kmax = k[index]
Emax = E_spectrum[index]

print('E spectrum shape = ', E_spectrum.shape)
print('k shape = ', k.shape)

print('\n Here is k!')
for nn in range(len(k)):
	print(k[nn])

print('\n Here is E!')
for nn in range(len(E_spectrum)):
        print(E_spectrum[nn])

print('\n Here is E*k5/3')
for nn in range(len(E_spectrum)):
        print(E_spectrum[nn]*pow(k[nn],5./3.))

ax = plt.gca()
plt.loglog(k, E_spectrum*pow(k,5./3.) , 'k')
#plt.loglog(k, E_spectrum , 'k')
plt.xlabel(r"k")
#plt.ylabel(r"k$^{5/3}$P(k)")
plt.ylabel(r"P(k)")
ax.set_xlim(1.e0, 3.e2)
#ax.set_ylim(1.e3, 1.e7)

plt.savefig("spectrum.eps", bbox_inches='tight')

