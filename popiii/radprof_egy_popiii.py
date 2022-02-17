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
G = YTQuantity(6.67e-8, 'cm**3/g/s**2')

#my_type = '0yr'
my_type = '1000yr'

def load_file(dir, datanum):

	pf = load(dir + 'data.' + datanum + ".3d.hdf5")

	if my_type == '0yr':
		value, location = pf.find_max("density")

	with open(dir + 'data.' + datanum + '.3d.sink', "r") as f:
		sinkdat = [line.split() for line in f]

	mass_sink = 0
	if my_type != '0yr':
		for i in range(1,2):
			line_arr = sinkdat[i]
			mass_sink = float(line_arr[0])
			xpos_sink = float(line_arr[1])
			ypos_sink = float(line_arr[2])
			zpos_sink = float(line_arr[3])

		location = [xpos_sink, ypos_sink, zpos_sink]

	data = pf.sphere(location, 0.1*pc_to_cm)
	print ('location =', location)

	bulk_vel0 = data.get_field_parameter("bulk_velocity").in_units("cm/s")
	print('bul_vel0_cm_s = ', bulk_vel0)
	calc_angmom_0(data, bulk_vel0, 0, 1.e15)

	bin_num2 = 100

	plot = ProfilePlot(data, "Radius", ['cs','vrot'], n_bins = bin_num2, weight_field="cell_mass")
	profile = plot.profiles[0]
	rad = profile.x
	cs = profile['cs']	
	vrot = profile['vrot']

	plot = ProfilePlot(data, "Radius", ['density','density2', 'Bmag', 'Bmag2'], n_bins = bin_num2, weight_field="cell_volume")
	profile = plot.profiles[0]
	rho = profile['density']
	bmag = profile['Bmag']
	bmag_rms = np.sqrt(profile['Bmag2']) 
	rho_rms = np.sqrt(profile['density2'])

	plot2 = ProfilePlot(data, "Radius", ['cell_mass', 'cs_mass'], n_bins = bin_num2, weight_field=None, accumulation=True)
	profile2 = plot2.profiles[0]
	menc = profile2['cell_mass']
	cs_avg = profile2['cs_mass']

	egy_grav = G * (menc + mass_sink* YTQuantity(1.0, 'g')) / rad
	egy_therm = 1.5 * cs * cs
	egy_rot = 0.5 * vrot * vrot * 1.e10
	egy_mag = bmag * bmag / 8.0 / 3.14159 / rho 

	print('egy_grav = ', egy_grav)
	print('egy_therm = ', egy_therm)
	print('egy_rot = ', egy_rot)
	print('egy_mag = ', egy_mag)

	return rad/1.5e13, egy_grav, egy_therm, egy_rot, egy_mag

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

dir_name = '/work2/00863/minerva/stampede2/popiii_bfieldA/'
if my_type == '0yr':
	datanum = '0245'
else:
	datanum = '0400'
rad, egy_grav, egy_therm, egy_rot, egy_mag = load_file(dir_name, datanum)

dir_name = '/work2/00863/minerva/stampede2/popiii_bfieldA0/'
if my_type == '0yr':
	datanum = '0240'
else:
	datanum = '0400'
rad0, egy_grav0, egy_therm0, egy_rot0, egy_mag0 = load_file(dir_name, datanum)

fsize = 14

rmin = 1.e1
rmax = 1.e4

mmin = 1.e8
mmax = 1.e14
pl.clf()
pl.cla()
pl.plot(rad, egy_grav,'k', linewidth=3.5, label='gravitational specific energy (MHD)')
pl.plot(rad, egy_therm,'b:', linewidth=3.5, label='thermal specific energy (MHD)')
pl.plot(rad, egy_rot,'r--', linewidth=3.5, label='rotational specific energy (MHD)')
pl.plot(rad, egy_mag,'g-.', linewidth=3.5, label='magnetic specific energy (MHD)')

ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel('mass ratio', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((rmin, rmax, mmin, mmax))
pl.legend(loc='upper right', fontsize=10)
if my_type == '0yr':
        pl.title('t = 0 years (MHD)')
else:
        pl.title('t = 1000 years (MHD)')
if my_type == '0yr':
        pl.savefig('radprof_egy_t240.eps')
else:
        pl.savefig('radprof_egy_t400.eps')

pl.clf()
pl.cla()

pl.plot(rad0, egy_grav0,'k', linewidth=3.5, label='gravitational specific energy (Hydro)')
pl.plot(rad0, egy_therm0,'b:', linewidth=3.5, label='thermal specific energy (Hydro)')
pl.plot(rad0, egy_rot0,'r--', linewidth=3.5, label='rotational specific energy (Hydro)')

ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel('mass ratio', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((rmin, rmax, mmin, mmax))
pl.legend(loc='upper right', fontsize=10)
if my_type == '0yr':
        pl.title('t = 0 years (Hydro)')
else:
        pl.title('t = 1000 years (Hydro)')
if my_type == '0yr':
        pl.savefig('radprof_egyH_t240.eps')
else:
        pl.savefig('radprof_egyH_t400.eps')

