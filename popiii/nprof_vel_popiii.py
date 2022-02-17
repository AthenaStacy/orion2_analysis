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
	value, location = pf.h.find_max("density")
	#data = pf.h.sphere(location, 3.0*pc_to_cm)
	data = pf.h.sphere(location * 0, 3.0*pc_to_cm)
	print ('value = ', value, 'location =', location)
	data.set_field_parameter("center", location)

	center = data.get_field_parameter('center')
	print("center main =", center)

	bulk_vel0 = data.get_field_parameter("bulk_velocity")	
	print('bul_vel0 = ', bulk_vel0)

	bulk_vel0 = data.get_field_parameter("bulk_velocity").in_units("cm/s")
	print('bul_vel0_cm_s = ', bulk_vel0)
	#calc_angmom_0(data, bulk_vel0, 1.e1, 1.e17)
	calc_angmom_0(data, bulk_vel0, 0, 1.e15)

	normal = data.get_field_parameter("normal")
	print('normal = ', normal)

	bin_num2 = 100
	nmin = 1.e-20
	nmax = 1.e-9

	plot = ProfilePlot(data, "n_h", ['velocity', 'vturb2', 'vrad2', 'vrad', 'radial_velocity', 'Temp'], n_bins = bin_num2, weight_field="cell_mass")

	profile = plot.profiles[0]

	#print "profile data = ", profile

	den = profile.x
	vel = profile['velocity']
	vturb = np.sqrt(profile['vturb2'])
	vrad = np.sqrt(profile['vrad2'])      
	vrad_yt = profile['radial_velocity']  
	temp = profile['Temp'] 

	return den, vel, vturb, vrad, temp, vrad_yt/1.e5

def calc_angmom_0(data, bulk_vel0, rad_in, rad_out):

	print('rad_out = ', rad_out)

	angmom_x_avg = 0	
	angmom_y_avg = 0
	angmom_z_avg = 1

	vel_x_avg = 0
	vel_y_avg = 0
	vel_z_avg = 0

	#inds = np.where(na.logical_and(data['density'] > den_in, data['density'] < den_out))
	inds2 = np.where(data['Radius'] < rad_out)
	inds = np.where(na.logical_and(data['Radius'] > rad_in, data['Radius'] < rad_out))

	mass_tot = na.sum(data['cell_mass'][inds])
	marr = data['cell_mass'][inds]

	mass_tot2 = na.sum(data['cell_mass'][inds2])
	marr2 = data['cell_mass'][inds2]

	vel_x_avg = na.sum(data['x-velocity'][inds2].value * marr2.value) / mass_tot2.value
	vel_y_avg = na.sum(data['y-velocity'][inds2].value * marr2.value) / mass_tot2.value
	vel_z_avg = na.sum(data['z-velocity'][inds2].value * marr2.value) / mass_tot2.value

	print('vel_x = ', vel_x_avg, 'vel_y = ', vel_y_avg, 'vel_z = ', vel_z_avg)
	print('marr2 = ', marr2)
	print('mass_tot2  = ', mass_tot2)

	bulk_vel = YTArray([vel_x_avg, vel_y_avg, vel_z_avg], 'cm/s')

	angmom_x_avg = na.sum(data['angular_momentum_x'][inds2].value * marr2) / mass_tot2 / 1.e50
	angmom_y_avg = na.sum(data['angular_momentum_y'][inds2].value * marr2) / mass_tot2 / 1.e50
	angmom_z_avg = na.sum(data['angular_momentum_z'][inds2].value * marr2) / mass_tot2 / 1.e50

	data.set_field_parameter("bulk_velocity", bulk_vel)
	#data.set_field_parameter('normal', [angmom_x_avg, angmom_y_avg, angmom_z_avg])

	bulk_vel = data.get_field_parameter("bulk_velocity")
	print('bul_vel = ', bulk_vel)

	normal = data.get_field_parameter("normal")
	print('normal = ', normal)

dir_name = '/nobackupp12/astacy/popiii_bfieldA/'
datanum = '0245'
den, vel, vturb, vrad, temp, vrad_yt = load_file(dir_name, datanum)
cs = np.power(1.38e-16 * temp / 1.22 / 1.67e-24, 0.5) /1.e5
mturb = vturb / cs

'''
print('den = ', den)
print('vrad = ', vrad)
print('vrad_yt = ', vrad_yt)
print('cs = ', cs)
print('temp = ', temp)
print('vturb = ', vturb)
print('mturb = ', mturb)
'''

'''
dir_name = '/nobackupp12/astacy/popiii_bfieldA0/'
datanum = '0240'
den0, vel0, vturb0, vrad0, temp0 = load_file(dir_name, datanum)
cs0 = np.power(1.38e-16 * temp0 / 1.e-24, 0.5) /1.e5
mturb0 = vturb0 / cs0
'''

'''
den = [1.e4, 1.e6, 1.e9, 1.e10]
vrad = [-2, -1, -2, -1]
vturb = [2, 3, 3, 4]
mturb = [1, 2, 1, 2]
'''

nmin = 1.e4
nmax = 1.e14

fsize = 12

vmin = 0
vmax = 5
#pl.subplot(121)
ax = pl.gca()
ax.plot(den, vrad,'k', linewidth=2.0, label=r'|v$_{\rm rad}|$')
#ax.plot(den, abs(vrad_yt),'r', linewidth=2.0, label=r'|v$_{\rm rad yt}|$')
ax.plot(den, vturb,'k:', linewidth=2.0, label=r'v$_{\rm turb}$')
#ax.plot(den, vrad_yt,'r:', linewidth=2.0, label=r'v$_{\rm rad-yt}$')
pl.plot([1e5], [0],'b--', linewidth=2.0, label=r'M$_{\rm turb}$')
ax.set_xscale('log')
ax.set_xlabel(r'n$_{\rm H}$ [cm$^{-3}$]', fontsize=fsize)
ax.set_ylabel('velocity [km s$^{-1}$]', fontsize=fsize)
#ax.xaxis.set_ticklabels([])
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((nmin, nmax, vmin, vmax))
pl.legend(loc=1, fontsize=fsize)
#pl.title('Orion')
pl.text(1.e8, 4.5, 'Orion', fontsize=16)
ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
ax2.set_ylabel('M$_{turb}$', color='b')  # we already handled the x-label with ax1
ax2.plot(den, mturb, 'b--', linewidth=2.0, label=r'M$_{\rm turb}$')
ax2.tick_params(axis='y', labelcolor='b')
ax2.set_ylim(0, 2.0)
ax2.set_xlim(nmin, nmax)

'''
pl.subplot(122)
pl.plot(den, vturb,'k', linewidth=2.0, label='v$_{turb}$')
pl.plot(den, vrad,'g--', linewidth=2.0, label='v$_{rad}$')
pl.plot(den, mturb,'b:', linewidth=2.0, label='M$_{turb}$')
ax = pl.gca()
ax.set_xscale('log')
ax.set_xlabel('n$_{H}$ [cm$^{-3}$]', fontsize=10)
ax.set_ylabel('[km s$^{-1}$]', fontsize=10)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((nmin, nmax, vmin, vmax))
ax.set_title("MHD")
'''

pl.savefig('nprof_vel.eps')

'''
pl.clf()
pl.gca()
tmin = 0
tmax = 1200
ax = pl.gca()
ax.plot(den, temp,'k', linewidth=2.0, label=r'Temp [K]')
ax.set_xscale('log')
ax.set_xlabel(r'n$_{\rm H}$ [cm$^{-3}$]', fontsize=fsize)
ax.set_ylabel('Temp [K]', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((nmin, nmax, tmin, tmax))
#pl.legend(loc=1, fontsize=fsize)

pl.savefig('nprof_temp.eps')
'''
