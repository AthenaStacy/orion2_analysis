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

def load_file(dir, datanum, my_type):

	pf = load(dir + 'data.' + datanum + ".3d.hdf5")

	if my_type == '0yr':
		value, location = pf.find_max("density")

	with open(dir + 'data.' + datanum + '.3d.sink', "r") as f:
		sinkdat = [line.split() for line in f]

	if my_type != '0yr':
        	for i in range(1,2):
                	line_arr = sinkdat[i]
                	xpos_sink = float(line_arr[1])
                	ypos_sink = float(line_arr[2])
                	zpos_sink = float(line_arr[3])

       		location = [xpos_sink, ypos_sink, zpos_sink]

	data = pf.sphere(location, 0.1*pc_to_cm)
	print ('location =', location)

	bin_num2 = 100

	plot = ProfilePlot(data, "Radius", ['mdot', 'cs', 'gamma', 'mu'], n_bins = bin_num2, weight_field="cell_mass")
	profile = plot.profiles[0]
	rad = profile.x
	mdot = profile['mdot']
	cs = profile['cs']
	gamma = profile['gamma']
	mu = profile['mu']

	plot = ProfilePlot(data, "Radius", ['density','Temp'], n_bins = bin_num2, weight_field="cell_volume")
	profile = plot.profiles[0]
	rad = profile.x
	rho = profile['density']
	temp = profile['Temp']

	mbe = 1.182 * cs**3 * np.sqrt(1.0 / G**3 / rho) / 2.e33

	print('rad = ', rad)
	print('cs = ', cs)
	print('rho = ', rho)
	print('temp =', temp)

	return rad/1.5e13, rho, cs/1.e5, gamma, mu

dir_name = '/work2/00863/minerva/stampede2/popiii_bfieldA/'
datanum = '0245'
my_type = '0yr'
rad0, rho0, cs0, gamma0, mu0 = load_file(dir_name, datanum, my_type)

datanum = '0400'
my_type = '1000yr'
rad1, rho1, cs1, gamma1, mu1 = load_file(dir_name, datanum, my_type)

dir_name = '/work2/00863/minerva/stampede2/popiii_bfieldA0/'
datanum = '0240'
my_type = '0yr'
rad0a, rho0a, cs0a, gamma0a, mu0a = load_file(dir_name, datanum, my_type)

datanum = '0400'
my_type = '1000yr'
rad1a, rho1a, cs1a, gamma1a, mu1a = load_file(dir_name, datanum, my_type)

fsize = 14

pl.clf()
pl.cla()
rmin = 1.e1
rmax = 1.e4
dmin = 1.e-17
dmax = 1.e-11
pl.plot(rad0, rho0,'k', linewidth=1.5, label='MHD t=0yr')
pl.plot(rad1, rho1,'b', linewidth=1.5, label='MHD t=1000yr')
pl.plot(rad0a, rho0a,'k:', linewidth=3.5, label='Hydro t=0yr')
pl.plot(rad1a, rho1a,'b:', linewidth=3.5, label='Hydro t=1000yr')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel('rho [cm$^{-3}$]', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((rmin, rmax, dmin, dmax))
pl.legend(loc=1, fontsize=10)
pl.savefig('radprof_rho_t0_t1000.eps', bbox_inches='tight')

pl.clf()
pl.cla()
rmin = 1.e1
rmax = 1.e4
dmin = 1
dmax = 6
pl.plot(rad0, cs0,'k', linewidth=1.5, label='MHD t=0yr')
pl.plot(rad1, cs1,'b', linewidth=1.5, label='MHD t=1000yr')
pl.plot(rad0a, cs0a,'k:', linewidth=3.5, label='Hydro t=0yr')
pl.plot(rad1a, cs1a,'b:', linewidth=3.5, label='Hydro t=1000yr')
ax = pl.gca()
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel('cs [km/s]', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((rmin, rmax, dmin, dmax))
pl.legend(loc=1, fontsize=10)
pl.savefig('radprof_cs_t0_t1000.eps', bbox_inches='tight')

pl.clf()
pl.cla()
rmin = 1.e1
rmax = 1.e4
dmin = 1
dmax = 2
pl.plot(rad0, gamma0,'k', linewidth=1.5, label='MHD t=0yr')
pl.plot(rad1, gamma1,'b', linewidth=1.5, label='MHD t=1000yr')
pl.plot(rad0a, gamma0a,'k:', linewidth=3.5, label='Hydro t=0yr')
pl.plot(rad1a, gamma1a,'b:', linewidth=3.5, label='Hydro t=1000yr')
ax = pl.gca()
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel('gamma', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((rmin, rmax, dmin, dmax))
pl.legend(loc=1, fontsize=10)
pl.savefig('radprof_gamma_t0_t1000.eps', bbox_inches='tight')

pl.clf()
pl.cla()
rmin = 1.e1
rmax = 1.e4
dmin = 1
dmax = 3
pl.plot(rad0, mu0,'k', linewidth=1.5, label='MHD t=0yr')
pl.plot(rad1, mu1,'b', linewidth=1.5, label='MHD t=1000yr')
pl.plot(rad0a, mu0a,'k:', linewidth=3.5, label='Hydro t=0yr')
pl.plot(rad1a, mu1a,'b:', linewidth=3.5, label='Hydro t=1000yr')
ax = pl.gca()
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel('mu', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((rmin, rmax, dmin, dmax))
pl.legend(loc=1, fontsize=10)
pl.savefig('radprof_mu_t0_t1000.eps', bbox_inches='tight')
