import yt
from yt.units import dimensions
import matplotlib
matplotlib.use('ps')
from yt.mods import *
import numpy as na
import pylab as pl
#from string import rstrip
import fields_bfield
import tracer_def

MRH = 0.76
h2conv = (2.0/1.3158)
heconv = (4.0/1.3158)
pc_to_cm = 3.08567758e18

def load_file(dir, datanum):

	pf = load(dir + 'data.' + datanum + ".3d.hdf5")
	value, location = pf.find_max("density")
	data = pf.sphere(location, 0.05*pc_to_cm)
	print ('location =', location)

	bin_num2 = 100
	nmin = 1.e-20
	nmax = 1.e-9

	bulk_vel0 = data.get_field_parameter("bulk_velocity").in_units("cm/s")
	calc_angmom_0(data, bulk_vel0, 0, 1.e15)

	plot = ProfilePlot(data, "n_h", ['density', 'Bmag', 'MagEnergy', 'ThermalEnergy', 'KinEnergy', 'vturb2', 'cs2'], n_bins = bin_num2, weight_field="cell_mass")

	profile = plot.profiles[0]

	n_h = profile.x
	den = profile['density']
	bmag = profile['Bmag']
	U_B  = profile['MagEnergy'] / den
	U_K_alt = profile['KinEnergy']
	U_T_alt = profile['ThermalEnergy']
	v_turb = np.sqrt(profile['vturb2'])
	cs = np.sqrt(profile['cs2'])

	U_K = 0.5 * v_turb * v_turb
	U_T = 1.5 * cs * cs 

	print('nh = ', n_h)
	print('bmag = ', bmag)
	print('U_B = ', U_B)
	print('v_turb = ', v_turb)
	print('cs = ', cs)

	return n_h, U_B, U_K, U_T

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

        bulk_vel = YTArray([vel_x_avg, vel_y_avg, vel_z_avg], 'cm/s')

        angmom_x_avg = na.sum(data['angular_momentum_x'][inds2].value * marr2) / mass_tot2 / 1.e50
        angmom_y_avg = na.sum(data['angular_momentum_y'][inds2].value * marr2) / mass_tot2 / 1.e50
        angmom_z_avg = na.sum(data['angular_momentum_z'][inds2].value * marr2) / mass_tot2 / 1.e50

        data.set_field_parameter("bulk_velocity", bulk_vel)


dir_name = '/work2/00863/minerva/stampede2/popiii_bfieldA/'
datanum = '0245'
n_h, U_B, U_K, U_T = load_file(dir_name, datanum)

nmin = 1.e4
nmax = 1.e14
bmin = 0.9*min( [min(U_B.value), min(U_T.value)] )
bmax = 2.0*max( [max(U_B.value), max(U_T.value)] )

pl.plot(n_h, U_B,'k',   linewidth=2.0, label='magnetic specific energy')
pl.plot(n_h, U_K,'c--', linewidth=2.0, label='turbulent specific energy')
pl.plot(n_h, U_T,'g:',  linewidth=2.0, label='thermal specific energy')
fsize = 14
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'n$_{\rm H}$ [cm$^{-3}$]', fontsize=fsize)
ax.set_ylabel('specific energy [erg g$^{-1}$]', fontsize=fsize)
pl.title('t = 0 years')
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((nmin, nmax, bmin, bmax))
pl.legend(loc="lower right", fontsize=10)
pl.savefig('nprof_erat.eps', bbox_inches='tight')



