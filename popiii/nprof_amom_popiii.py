import yt
from yt.units import dimensions
import matplotlib
matplotlib.use('ps')
from yt.mods import *
import numpy as na
#import pylab as pl
import matplotlib.pyplot as pl
#from string import rstrip
import fields_bfield
import tracer_def
from yt import YTArray

MRH = 0.76
h2conv = (2.0/1.3158)
heconv = (4.0/1.3158)
pc_to_cm = 3.08567758e18

def load_file(dir, datanum):

	pf = load(dir + 'data.' + datanum + ".3d.hdf5")
	value, location = pf.find_max("density")
	data = pf.sphere(location, 3.0*pc_to_cm)
	print ('value = ', value, 'location =', location)
	data.set_field_parameter("center", location)

	center = data.get_field_parameter('center')
	print("center main =", center)

	bulk_vel0 = data.get_field_parameter("bulk_velocity")	
	print('bul_vel0 = ', bulk_vel0)

	'''
	ind = np.where(data['Radius'] < 1.e17)
	vx = np.sum(data['x-velocity'][ind] *data['cell_mass'][ind])/ np.sum(data['cell_mass'][ind])
	vy = np.sum(data['y-velocity'][ind] *data['cell_mass'][ind])/ np.sum(data['cell_mass'][ind])
	vz = np.sum(data['z-velocity'][ind] *data['cell_mass'][ind])/ np.sum(data['cell_mass'][ind])

	bulk_vel = YTArray([vx, vy, vz], 'cm/s')
	data.set_field_parameter("bulk_velocity", bulk_vel)
	bulk_vel0 = data.get_field_parameter("bulk_velocity").in_units("cm/s")
	print('bul_vel0_cm_s = ', bulk_vel0)
	'''

	bin_num2 = 100
	nmin = 1.e-20
	nmax = 1.e-9

	plot = ProfilePlot(data, "n_h", ['vrot', 'amom', 'vrad2', 'vrad', 'radial_velocity', 'Temp'], n_bins = bin_num2, weight_field="cell_mass")

	profile = plot.profiles[0]
	den = profile.x
	vrot = np.sqrt(profile['vrot'])
	amom = profile['amom']
	vrad = np.sqrt(profile['vrad2'])      
	vrad_yt = profile['radial_velocity']  
	temp = profile['Temp'] 

	plot = ProfilePlot(data, "radius", ['cell_mass', 'amom'], n_bins = bin_num2, weight_field=None, accumulation=True)
	
	profile = plot.profiles[0]
	rad = profile.x
	mass_enc = profile['cell_mass']
	amom_enc = profile['amom']

	return rad, den, mass_enc, vrot, vrad, vrad_yt/1.e5, amom, amom_enc

def read_gadget(dir_name,fname):

	fpath = dir_name + fname
	gdat = [] 
	with open(fpath, "r") as f:
		gdat= [line.split() for line in f]

	nh = []
	mass = []
	vrot = []
	rad = []
	
	for i in range(len(gdat)):
		line_arr = gdat[i]
		nh.append(float(line_arr[6])*1.22*.76)
		rad.append(float(line_arr[11]))
		mass.append(float(line_arr[15]))
		vrot.append(float(line_arr[13]))

	nh = np.asarray(nh)
	rad = np.asarray(rad)
	mass = np.asarray(mass)
	vrot = np.asarray(vrot)

	mass = mass*1.e10/.7*1.989e33

	nbins = np.logspace(4, 8, 200)
	print('nbins = ', nbins)

	rbins = np.logspace(1,5,200)
	print('rbins = ', rbins)

	print('rad = ', rad[0:500])

	amom_enc = []
	mass_enc = []

	ind = np.where(nh < nbins[0])
	amom_enc.append(np.mean(vrot[ind] * rad[ind]))
	mass_enc.append(np.sum(mass[ind]))

	for i in range(1,len(nbins)):
		#ind = np.where(np.logical_and(rad > 0, rad < rbins[i]))	
		ind = np.where(np.logical_and(nh > nbins[i-1], nh < nbins[i]) == True)
		amom_enc.append(np.mean(vrot[ind] * rad[ind]))
		mass_enc.append(np.sum(mass[ind]))

	amom_enc = np.asarray(amom_enc)
	mass_enc = np.asarray(mass_enc)

	print('amom_enc = ', amom_enc)
	print('mass_enc = ', mass_enc)
	
	return nbins, rbins, mass_enc, amom_enc

dir_name = '/work/00863/minerva/stampede2/popiii_bfieldA/'
datanum = '0000'
rad, den, mass_enc, vrot, vrad, vrad_yt, amom, amom_enc = load_file(dir_name, datanum)

dir_name = '/work/00863/minerva/stampede2/'
fname = 'snapbin2_zoom10_ref4_0537'
n_gadget, rad_gadget, mass_gadget, amom_gadget = read_gadget(dir_name, fname) 

nmin = 1.e5
nmax = 1.e8

fsize = 12

print('let us plot!')

'''
print('mass_enc = ', mass_enc)
print('amom_enc = ', amom_enc)
print('amom = ', amom)

amin = 1.e38
amax = 1.e42
mmin = 1.e0
mmax = 1.e3
ax = pl.gca()
ax.plot(mass_enc/1.989e33, amom, linewidth=2.0, label=r'Orion')
ax.plot(mass_gadget/1.989e33, amom_gadget,'k:', linewidth=2.0, label=r'Gadget')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'M$_{\rm enc}$ [M$_{\odot}$]', fontsize=fsize)
ax.set_ylabel('ang. mom. [g cm s$^{-1}$]', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((mmin, mmax, amin, amax))
pl.legend(loc=1, fontsize=fsize)

pl.savefig('nprof_amom.eps')
'''

print('den = ', den)
print('amom = ', amom)

amin = 1.e21
amax = 1.e24
ax = pl.gca()
ax.plot(den, amom, linewidth=2.0, label=r'Orion')
ax.plot(n_gadget, amom_gadget*1.e5*1.5e13, 'k:', linewidth=2.0, label=r'Gadget')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'n$_{\rm H}$ [cm$^{-3}$]', fontsize=fsize)
ax.set_ylabel('ang. mom. [cm$^{2}$ s$^{-1}$]', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((nmin, nmax, amin, amax))
pl.legend(loc=1, fontsize=fsize)

pl.savefig('nprof_amom.eps')

