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
datanum1 = '0000'  #--corresponding times and data steps for orion256
fname1 = odir1+"data."+ datanum1 +".3d.hdf5"
#fname1 = '/nobackupp12/astacy/to_share/data.0004.3d.hdf5'

pf = yt.load(fname1)

for i in sorted(pf.field_list):
  print(i)

# a FFT operates on uniformly gridded data.  We'll use the yt
# covering grid for this.

max_level = 0
#max_level = pf.index.max_level

#dims = [512,512,512]
#dims = [256,256,256]
dims = [128,128,128]
#dims = dims * pow(2,max_level)

fac = 1.0
left_edge = fac*pf.domain_left_edge
right_edge = pf.domain_right_edge

left_corner =  fac*left_edge
right_corner = fac*right_edge
region = pf.box(left_corner, right_corner)

cube = pf.smoothed_covering_grid(level=max_level, left_edge=left_edge, dims=dims, fields = ['density', 'X-momentum', 'Y-momentum', 'Z-momentum'], field_parameters = {'center':[0,0,0]})

#cube = region

print('max level = ', pf.index.max_level)

'''
vturb = cube['vturb']
temp = cube['Temp']
cs = np.power(1.38e-16 * temp / 1.22 / 1.67e-24, 0.5) /1.e5
mturb = vturb / cs
print('avg mach = ', np.mean(mturb))
'''

nx, ny, nz = np.asarray(dims)

# wavenumbers
L = (pf.domain_right_edge - pf.domain_left_edge) / (pf.domain_right_edge[0] - pf.domain_left_edge[0]) 
print ('L[0] = ', L[0], 'L[1] =', L[1], 'L[2] =', L[2], 'nx = ', nx, 'ny = ', ny, 'nz =', nz) 

fac_here = 1
nx2, ny2, nz2 = int(nx/fac_here), int(ny/fac_here), int(nz/fac_here)

kx = np.fft.fftfreq(nx2)*nx2/L[0]
ky = np.fft.fftfreq(ny2)*ny2/L[1]
kz = np.fft.fftfreq(nz2)*nz2/L[2]
print ('kx shape = ', kx.shape)

print (' cube shape = ', cube['x'].shape)


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

Kk  = np.zeros((l3-l2, l3-l2, l3-l2))
for i in range(len(fields)):
	array = cube[fields[i]].d #/ cube['density'].d
	array = array[l2:l3, l2:l3, l2:l3]
	ru = np.fft.fftn(array)
	#ru = np.abs(ru)
	if i == 0:
		vk_x = ru
	if i == 1:
		vk_y = ru
	if i == 2:
		vk_z = ru
	ru = np.abs(ru)
	ru = np.power((ru),2)
	ru = ru/(nx*ny*nz)/8
	Kk += ru

fac = np.sqrt(3.0)

vc = vk_x + vk_y + vk_z
vc = vc / fac
vc_x = vc / fac
vc_y = vc / fac
vc_z = vc / fac

vs_x1 = (vk_z - vk_y) / fac
vs_y1 = (vk_x - vk_z) / fac
vs_z1 = (vk_y - vk_x) / fac

vs_x = (vs_y1 - vs_z1) / fac
vs_y = (vs_z1 - vs_x1) / fac
vs_z = (vs_x1 - vs_y1) / fac

'''
Pvc = np.sqrt(3 * np.fft.fftn(vc)*np.fft.fftn(vc))

Pvs = np.sqrt(np.fft.fftn(vs_x)*np.fft.fftn(vs_x) + np.fft.fftn(vs_y)*np.fft.fftn(vs_y) + np.fft.fftn(vs_z)*np.fft.fftn(vs_z))

Pv = np.sqrt(np.fft.fftn(vk_x)*np.fft.fftn(vk_x) + np.fft.fftn(vk_y)*np.fft.fftn(vk_y) + np.fft.fftn(vk_z)*np.fft.fftn(vk_z))
'''

Pvc = abs(vc_x*vc_x + vc_y*vc_y + vc_z*vc_z)

Pvs = abs(vs_x*vs_x + vs_y*vs_y + vs_z*vs_z)

Pv  = abs(vk_x*vk_x + vk_y*vk_y + vk_z*vk_z)

print('vc shape = ', vc.shape)
print('vk_x shape = ', vk_x.shape, 'vk_y shape = ', vk_y.shape, 'vk_z shape = ', vk_z.shape,)
print('vs_x shape = ', vs_x.shape, 'vs_y shape = ', vs_y.shape, 'vs_z shape = ', vs_z.shape,)
print('k shape = ', k.shape, 'kmin = ', kmin, 'kmax = ', kmax)
print('KK shape = ', Kk.shape)
print('kbins = ', kbins)

whichbin = np.digitize(k.flat, kbins)
ncount = np.bincount(whichbin)    
E_spectrum  = np.zeros(len(ncount)-1)
Pvc_spec = np.zeros(len(ncount)-1)
Pvs_spec = np.zeros(len(ncount)-1)
Pv_spec = np.zeros(len(ncount)-1)
print('len(ncount)-1 ', len(ncount)-1)

#shift_num = (l3-l2)/2 -1
shift_num=0
rad_grid = np.ones((l3-l2, l3-l2, l3-l2)) * 1.e3

for nn in range(1,len(ncount)):
	E_spectrum[nn-1]  = np.sum(Kk.flat[whichbin==nn])
	Pvs_spec[nn-1]  = np.sum(Pvs.flat[whichbin==nn])
	Pvc_spec[nn-1]  = np.sum(Pvc.flat[whichbin==nn])
	Pv_spec[nn-1]  = np.sum(Pv.flat[whichbin==nn])

E_spectrum  = E_spectrum[0:N-2]
Pvs_spec = Pvs_spec[0:N-2]
Pvc_spec = Pvc_spec[0:N-2]
Pv_spec  = Pv_spec[0:N-2]
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
ax.set_xlim(1.e0, 1.e2)
plt.savefig("spectrum.eps", bbox_inches='tight')

plt.clf()
ax = plt.gca()
plt.loglog(k, Pv_spec , 'k', label = 'overall')
#plt.loglog(k, Pvs_spec + Pvc_spec, 'g:', label = 'sum')
plt.loglog(k, Pvs_spec , 'b:', label = 'solenoidal')
plt.loglog(k, Pvc_spec , 'r--', label= 'compressive')
plt.xlabel(r"k")
plt.ylabel(r"P$_v$(k)")
ax.set_xlim(1.e0, 1.e2)
plt.legend(loc='upper right', fontsize=10)
plt.savefig("spectrum_vel.eps", bbox_inches='tight')

plt.clf()
ax = plt.gca()
plt.plot(k, Pvs_spec/Pvc_spec , 'k', label = 'soleniodal / compressive')
ax.set_xscale('log')
plt.xlabel(r"k")
plt.ylabel(r"solenoidal to compressive ratio")
ax.set_xlim(1.e0, 1.e2)
plt.legend(loc='upper right', fontsize=10)
plt.savefig("spectrum_ratio.eps", bbox_inches='tight')

print('vk_x=', vk_x[50][50][50])
print('sum = ',vc_x[50][50][50]+vs_x[50][50][50])
print('vc_x= ', vc_x[50][50][50])
print('vs_x = ', vs_x[50][50][50])

print('vk_y=', vk_y[50][50][50])
print('sum = ',vc_y[50][50][50]+vs_y[50][50][50])
print('vc_y= ', vc_y[50][50][50])
print('vs_y = ', vs_y[50][50][50])

print('vk_z=', vk_z[50][50][50])
print('sum = ',vc_z[50][50][50]+vs_z[50][50][50])
print('vc_z= ', vc_z[50][50][50])
print('vs_z = ', vs_z[50][50][50])

print('Pv  =', Pv[50][50][50])
print('sum = ',Pvc[50][50][50]+Pvs[50][50][50])
print('Pvc = ', Pvc[50][50][50])
print('Pvy = ', Pvs[50][50][50])

'''
plt.clf()
ax = plt.gca()
plt.loglog(k, vk_x , 'k', label = 'overall')
plt.loglog(k, vc_x+vs_x , 'g:', label = 'sum')
plt.loglog(k, vc_x , 'b:', label = 'solenoidal')
plt.loglog(k, vs_x , 'r--', label= 'compressive')
plt.xlabel(r"k")
plt.ylabel(r"P$_v$(k)")
ax.set_xlim(1.e0, 1.e2)
plt.legend(loc='upper right', fontsize=10)
plt.savefig("spectrum_check.eps", bbox_inches='tight')
'''

