import yt
from yt.units import dimensions
import matplotlib
matplotlib.use('ps')
from yt.mods import *
import numpy as np
import pylab as pl
#from string import rstrip
import fields_bfield
import tracer_def

MRH = 0.76
h2conv = (2.0/1.3158)
heconv = (4.0/1.3158)
pc_to_cm = 3.08567758e18
G = YTQuantity(6.67e-8, 'cm**3/g/s**2')
cm = YTQuantity(1, 'cm')

my_type = '0yr'
#my_type = '1000yr'

def load_file(dir, datanum):

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

	#plot = ProfilePlot(data, "Radius", ['mdot', 'cs', 'density','Bmag'], n_bins = bin_num2, weight_field="cell_mass")
	plot = ProfilePlot(data, "Radius", ['cs'], n_bins = bin_num2, weight_field="cell_mass")
	profile = plot.profiles[0]
	rad = profile.x
	cs = profile['cs']	

	plot = ProfilePlot(data, "Radius", ['density','density2', 'Bmag', 'Bmag2'], n_bins = bin_num2, weight_field="cell_volume")
	profile = plot.profiles[0]
	rho = profile['density']
	bmag = profile['Bmag']
	bmag_rms = np.sqrt(profile['Bmag2']) 
	rho_rms = np.sqrt(profile['density2'])

	mbe = 1.182 * cs**3 * np.sqrt(1.0 / G**3 / rho) / 2.e33
	mcrit = bmag * rad * rad / 2. / np.sqrt(G) / 2.e33

	print('rad = ', rad)
	print('cs = ', cs)
	print('rho = ', rho)
	print('rho_rms ', rho_rms)
	print('bmag = ', bmag)
	print('bmag_rms', bmag_rms)

	plot2 = ProfilePlot(data, "Radius", ['cell_mass', 'cs_mass'], n_bins = bin_num2, weight_field=None, accumulation=True)
	profile2 = plot2.profiles[0]
	menc = profile2['cell_mass']
	cs_avg = profile2['cs_mass']
	cs_avg = cs_avg / menc

	beta = 8 * 3.14159 * rho * cs*cs / (bmag_rms*bmag_rms)
	beta_fac_avg = 1.0 + 0.74 / beta

	delta_x = 1.543e+18 / 128. / (np.power(2.,8.)) *  YTQuantity(1.0, 'cm') 
	jmax = (1.0/4.0)

	cs2_alt = 1.e-14 * YTQuantity(1.0, 'dyne / cm**2') / data['density']
	cs2 = data['cs'] * data['cs']
	cs2_max = np.maximum(cs2, cs2_alt)

	print('delta_x = ', delta_x)
	value, location = pf.find_min("dx")
	print('min dx = ', value)

	#rho_max = 3.14159 * jmax * jmax * cs * cs / G / delta_x / delta_x * (1.0 + 0.74 / beta) 
	rho_max = np.zeros(np.shape(beta))
	ratio_max = np.zeros(np.shape(beta))	
	beta_maxi = np.zeros(np.shape(beta))
	beta_fac = np.zeros(np.shape(beta))
	beta_fac_maxi = np.zeros(np.shape(beta))

	ind = np.where(data['Radius'] < rad[0])

	rho_max_here = 3.14159 * jmax * jmax * cs2_max[ind] / G / delta_x / delta_x * (1.0 + 0.74 / data['beta'][ind])
	beta_here = data['beta'][ind]
	beta_fac_here = 1.0 + 0.74 / data['beta'][ind]	

	rho_here = data['density'][ind]
	rho_max[0] = np.amax(rho_here)

	rho_here = rho_here.flatten()
	ratio = rho_here / rho_max_here
	ratio = ratio.flatten() 
	ratio_max[0] =  np.amax(ratio)
	beta_fac[0] = np.amax(beta_fac_here)

	beta_flat= beta_here.flatten()
	beta_fac_flat = beta_fac_here.flatten()
	imax = np.argmax(rho_here)
	beta_maxi[0] = beta_flat[imax]
	beta_fac_maxi[0] = beta_fac_flat[imax]	

	for i in range(1,len(rad)):
		ind =  np.where(np.logical_and(data['Radius'] > rad[i-1], data['Radius'] < rad[i]) == True)
		if (len(data['beta'][ind]) > 0):
			rho_max_here = 3.14159 * jmax * jmax * cs2_max[ind] / G / delta_x / delta_x * (1.0 + 0.74 / data['beta'][ind])
			beta_here = data['beta'][ind]
			beta_fac_here = 1.0 + 0.74 / data['beta'][ind]

			rho_here = data['density'][ind]
			rho_max[i] = np.amax(rho_here)

			rho_here = rho_here.flatten()
			ratio = rho_here / rho_max_here		
			ratio = ratio.flatten()
			ratio_max[i] =  np.amax(ratio)
			beta_fac[i] = np.amax(beta_fac_here)

			beta_flat = beta_here.flatten()
			beta_fac_flat = beta_fac_here.flatten()
			imax = np.argmax(rho_here)
			beta_maxi[i] = beta_flat[imax]
			beta_fac_maxi[i] = beta_fac_flat[imax]
			print('ratio shape', ratio.shape)
			print('beta_flat shape', beta_flat.shape)
			print('imax = ', imax, 'beta_maxi =',beta_maxi[i], 'beta_fac_maxi = ', beta_fac_maxi[i])

	print('rad = ', rad)
	print('ratio_max = ', ratio_max)
	print('rho = ', rho)
	print('cs_avg = ', cs_avg)
	print('mbe =', mbe)
	print('beta = ', beta)

	return rad/1.5e13, mbe, menc/2.e33, mcrit, bmag_rms, rho_rms, rho, rho_max, ratio_max, cs/1.e5, cs_avg/1.e5, beta, beta_maxi, beta_fac,beta_fac_avg, beta_fac_maxi


dir_name = '/work2/00863/minerva/stampede2/popiii_bfieldA/'
if my_type == '0yr':
	datanum = '0245'
else:
	datanum = '0420'
rad, mbe, menc, mcrit, bmag_rms, rho_rms, rho, rho_max, ratio_max, cs, cs_avg, beta, beta_maxi, beta_fac, beta_fac_avg, beta_fac_maxi = load_file(dir_name, datanum)

dir_name = '/work2/00863/minerva/stampede2/popiii_bfieldA0/'
if my_type == '0yr':
        datanum = '0240'
else:
        datanum = '0420'
rad0, mbe0,  menc0, mcrit0, bmag_rms0, rho_rms0, rho0, rho_max0, ratio_max0, cs0, cs_avg0, beta0, beta_maxi0, beta_fac0, beta_fac_avg0, beta_fac_maxi0 = load_file(dir_name, datanum)

fsize = 14

rmin = 1.e1
rmax = 1.e4

mmin = 1.e-1
mmax = 1.e2
pl.plot(rad, beta,'k', linewidth=1.5, label='beta')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel('beta', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((rmin, rmax, mmin, mmax))
pl.legend(loc=2, fontsize=10)
if my_type == '0yr':
	pl.title('t = 0 years')
else:
	pl.title('t = 1000 years')
if my_type == '0yr':
	pl.savefig('radprof_beta_t240.eps')
else:
	pl.savefig('radprof_beta_t420.eps')

pl.gca()
pl.clf()
print('rho/rho_max = ', rho_max)
mmin = 1.e-5
mmax = 1.e1
pl.plot(rad, ratio_max,'k', linewidth=2.5, label='MHD')
pl.plot(rad0, ratio_max0,'b:', linewidth=2.5, label='Hydro')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel(r'max $\mathrm{\rho}$ / $\mathrm{\rho_{max}}$', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((rmin, rmax, mmin, mmax))
pl.legend(loc=2, fontsize=10)
if my_type == '0yr':
        pl.title('t = 0 years')
else:
        pl.title('t = 1000 years')
if my_type == '0yr':
        pl.savefig('radprof_ratio_t240.eps')
else:
        pl.savefig('radprof_ratio_t420.eps')

pl.gca()
pl.clf()
bmin = 1.e-1
bmax = 1.e3
pl.plot(rad, beta_fac_maxi,'k:', linewidth=2.5, label=r'$\mathrm{\beta}$ factor at max $\mathrm{\rho}$')
pl.plot(rad, beta_fac_avg,'k', linewidth=2.5, label=r'avg $\mathrm{\beta}$ factor')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel(r'$\mathrm{\beta}$ factor', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((rmin, rmax, bmin, bmax))
pl.legend(loc=2, fontsize=10)
if my_type == '0yr':
        pl.title('t = 0 years')
else:
        pl.title('t = 1000 years')
if my_type == '0yr':
        pl.savefig('radprof_betafac_t240.eps')
else:
        pl.savefig('radprof_betafac_t420.eps')

pl.gca()
pl.clf()
bmin = 1.e-3
bmax = 1.e4
pl.plot(rad, beta_maxi,'k', linewidth=2.5, label=r'$\mathrm{\beta}$ at max $\mathrm{\rho}$')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel(r'$\mathrm{\beta}$', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((rmin, rmax, bmin, bmax))
pl.legend(loc=2, fontsize=10)
if my_type == '0yr':
        pl.title('t = 0 years')
else:
        pl.title('t = 1000 years')
if my_type == '0yr':
        pl.savefig('radprof_betamax_t240.eps')
else:
        pl.savefig('radprof_betamax_t420.eps')

pl.gca()
pl.clf()
dmin = 1.e-17
dmax = 1.e-10
pl.plot(rad, rho_max,'k', linewidth=2.5, label=r'MHD')
pl.plot(rad0, rho_max0,'b:', linewidth=2.5, label=r'Hydro')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel(r'max $\mathrm{\rho}$', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((rmin, rmax, dmin, dmax))
pl.legend(loc='upper right', fontsize=10)
if my_type == '0yr':
        pl.title('t = 0 years')
else:
        pl.title('t = 1000 years')
if my_type == '0yr':
        pl.savefig('radprof_rhomax_t240.eps')
else:
        pl.savefig('radprof_rhomax_t420.eps')
