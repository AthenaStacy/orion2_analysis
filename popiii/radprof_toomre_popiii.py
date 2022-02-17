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

G = YTQuantity(6.67e-8, 'cm**3/g/s**2')
MRH = 0.76
h2conv = (2.0/1.3158)
heconv = (4.0/1.3158)
pc_to_cm = 3.08567758e18

my_type = '0yr'
#my_type = '1000yr'

def load_file(dir, datanum):

	pf = load(dir + 'data.' + datanum + ".3d.hdf5")

	if my_type == '0yr':
		value, location = pf.find_max("density")

	with open(dir + 'data.' + datanum + '.3d.sink', "r") as f:
		sinkdat = [line.split() for line in f]

	mass_sink = 1 *  YTQuantity(1.0, 'g') 
	if my_type != '0yr':
		for i in range(1,2):
			line_arr = sinkdat[i]
			mass_sink = float(line_arr[0]) *  YTQuantity(1.0, 'g')
			xpos_sink = float(line_arr[1])
			ypos_sink = float(line_arr[2])
			zpos_sink = float(line_arr[3])

		location = [xpos_sink, ypos_sink, zpos_sink]

	data = pf.sphere(location, 0.5*pc_to_cm)
	print ('location =', location)

	#bulk_vel0 = data.get_field_parameter("bulk_velocity").in_units("cm/s")
	#print('bul_vel0_cm_s = ', bulk_vel0)
	#calc_angmom_0(data, bulk_vel0, 0, 1.e18)

	bin_num1 = 5000
	bin_num2 = 100
	nmin = 1.e-20
	nmax = 1.e-9

	plot = ProfilePlot(data, "Radius", ['mdot', 'density', 'Temp', 'vrad', 'Omega', 'Toomre'], n_bins = bin_num2, weight_field="cell_mass")

	profile = plot.profiles[0]

	print ("profile data = ", profile)

	rad = profile.x
	mdot = profile['mdot']
	den = profile['density']
	toomre = profile['Toomre']        
	temp = profile['Temp'] 
	vrad = profile['vrad']
	omega = profile['Omega']

	plot = ProfilePlot(data, "Radius", ['n_h'], n_bins = bin_num2, weight_field="cell_volume")
	profile = plot.profiles[0]
	volume = (4./3.) * 3.14159 * np.power(rad,3)
	den_2 = mass_sink / volume
	nh_rad = profile['n_h']

	plot = ProfilePlot(data, "n_h", ['Temp'], n_bins = bin_num2, weight_field="cell_volume")
	profile = plot.profiles[0]
	nh = profile.x
	temp_nh = profile['Temp']

	toomre_alt = omega * omega / 3.14159 / G / den
	toomre_alt_2 = omega * omega / 3.14159 / G / (den + den_2)

	mdot = 4.0 * 3.14159 * rad * rad * den * vrad

	print('rad = ', rad)
	print('vrad = ',vrad)
	print('den = ', den)
	print('mdot = ', mdot)
	print('omega= ', omega)
	print('toomre = ', toomre)
	print('toomre_alt = ',toomre)
	print('toomre_alt_2 = ', toomre_alt_2)

	mfac = 3.15e7 / 2.e33

	return rad/1.5e13, mdot*mfac, den, temp, nh_rad, nh, temp_nh, toomre, toomre_alt, toomre_alt_2, vrad, omega

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
rad, mdot, den, temp, nh_rad, nh, temp_nh, toomre, toomre_alt, toomre_alt_2, vrad, omega = load_file(dir_name, datanum)

rad_alt = [10.3029, 12.8750, 16.7430, 22.4757, 28.5357, 35.0944, 44.5567, 56.1202, 71.8194, 90.4616, 111.253, 137.913, 168.266, 202.053, 240.704, 289.023, 358.296, 444.181, 563.968, 682.684, 794.174, 969.047, 1163.84, 1408.96, 1678.0, 1937.56, 2200.95, 2520.10, 2885.55, 3330.31, 3813.19, 4507.38, 5079.34, 5815.66, 6819.60, 7933.47, 8869.07, 9994.23]
rad_alt = np.asarray(rad_alt)
rad_alt = rad_alt*1.5e13

temp_alt = [931.692, 925.816, 910.049, 910.134, 916.143, 922.143, 928.152, 940.099, 950.071, 954.097, 962.077, 974.020, 980.018, 991.951, 1001.90, 1023.74, 1027.76, 1027.82, 1023.93, 1016.06, 1002.24, 990.421, 962.750, 935.082, 903.449, 881.708, 854.022, 828.318, 800.634, 778.893, 755.169, 727.495, 707.727, 689.944, 672.168, 654.389, 640.560, 626.733]

dir_name = '/work2/00863/minerva/stampede2/popiii_bfieldA0/'
if my_type == '0yr':
	datanum = '0240'
else:
	datanum = '0400'
rad0, mdot0, den0, temp0, nh_rad0, nh0, temp_nh0, toomre0, toomre_alt0, toomre_alt_20, vrad0, omega0 = load_file(dir_name, datanum)


rmin = 1.e1
rmax = 1.e5
#rmax = 1.e6

print ('rad0', rad0)
print ('den0', den0)

fsize = 14

dmin = 1.e5
dmax = 1.e13
pl.clf()
pl.gca()
pl.plot(rad0, nh_rad0,'k', linewidth=2.0, label='hydro')
pl.plot(rad, nh_rad,'b--', linewidth=2.5, label='MHD')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('radius [cm]', fontsize=fsize)
ax.set_ylabel(r'n$_{\rm H}$ [cm$^{-3}$]', fontsize=fsize)
pl.legend(loc='upper right', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((rmin, rmax, dmin, dmax))
pl.savefig('radprof_nh.eps', bbox_inches='tight')

tmin = 300
tmax = 1300
pl.clf()
pl.gca()
pl.plot(rad0, temp0,'k', linewidth=2.0, label='hydro')
pl.plot(rad, temp,'b--', linewidth=2.0, label='MHD')
#pl.plot(rad_alt, temp_alt,'b--', linewidth=2.0, label='MHD')
ax = pl.gca()
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('radius [cm]', fontsize=fsize)
ax.set_ylabel('Temp [K]', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((rmin, rmax, tmin, tmax))
pl.legend(loc=1, fontsize=fsize)
pl.savefig('radprof_temp.eps')

tmin = 300
tmax = 1300
pl.clf()
pl.gca()
pl.plot(nh0, temp_nh0,'k', linewidth=2.0, label='hydro')
pl.plot(nh, temp_nh,'b--', linewidth=2.0, label='MHD')
ax = pl.gca()
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel(r'n$_{\rm H}$ [cm$^{-3}$]', fontsize=fsize)
ax.set_ylabel(r'Temp [K]', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((dmin, dmax, tmin, tmax))
pl.legend(loc=0, fontsize=fsize)
pl.savefig('nprof_temp.eps', bbox_inches='tight')

mmin = 1.e-3
mmax = 1.e0
pl.clf()
pl.gca()
pl.plot(rad0, mdot0,'k', linewidth=2.0, label='hydro')
pl.plot(rad, mdot,'b:', linewidth=2.0, label='MHD')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel('dM/dt [M$_{\odot}$ yr$^{-1}$]', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.legend(loc='upper right', fontsize=10)
pl.axis((rmin, rmax, mmin, mmax))
if my_type == '0yr':
	pl.title('t = 0 years')
	pl.savefig('radprof_mdot_t240.eps')
else:
	pl.title('t = 750 years')
	pl.savefig('radprof_mdot_t360.eps')


qmin = 0.5
#qmax = 3
qmax = 10
pl.clf()
pl.gca()
pl.title('t = 1000 years') 
pl.plot(rad0, toomre_alt0,'k', linewidth=2.0, label='hydro')
pl.plot(rad, toomre_alt,'b--', linewidth=2.0, label='MHD')
ax = pl.gca()
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel(r'Q$_{\rm local}$', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.legend(loc=1, fontsize=fsize)
pl.axis((rmin, rmax, qmin, qmax))
pl.savefig('radprof_toomre_alt.eps')

qmin = 0.1
qmax = 5.0
pl.gca()
pl.title('t = 1000 years')
pl.plot(rad0, toomre_alt_20,'k', linewidth=2.0, label='hydro')
pl.plot(rad, toomre_alt_2,'b--', linewidth=2.0, label='MHD')
ax = pl.gca()
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel(r'Q$_{\rm local}$', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.legend(loc=1, fontsize=fsize)
pl.axis((rmin, rmax, qmin, qmax))
pl.savefig('radprof_toomre_alt_2.eps')

omin=0
omax=2
pl.clf()
pl.gca()
pl.plot(rad0, omega/omega0,'k', linewidth=2.0)
#pl.plot(rad, omega,'b--', linewidth=2.0, label='MHD')
ax = pl.gca()
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel(r'Omega_MHD / Omega_Hydro', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.legend(loc=1, fontsize=fsize)
pl.axis((rmin, rmax, omin, omax))
pl.savefig('radprof_omega_t=0.eps')

omin=0
omax=4.0
pl.clf()
pl.gca()
pl.plot(rad0, den/den0,'k', linewidth=2.0)
ax = pl.gca()
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel(r'rho_MHD / rho_Hydro', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.legend(loc=1, fontsize=fsize)
pl.axis((rmin, rmax, omin, omax))
pl.savefig('radprof_rho_t=0.eps')

vmin=-5
vmax=10
pl.clf()
pl.gca()
pl.plot(rad0, vrad0,'k', linewidth=2.0, label = 'hydro')
pl.plot(rad, vrad,'b:', linewidth=2.0, label='MHD')
ax = pl.gca()
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel(r'v$_{rad}$', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.legend(loc=1, fontsize=fsize)
pl.axis((rmin, rmax, vmin, vmax))
if my_type == '0yr':
	pl.savefig('radprof_vrad_t240.eps')
else:
	pl.savefig('radprof_vrad_t360.eps')
