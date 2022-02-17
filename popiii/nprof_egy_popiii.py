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

	if my_type != '0yr':
        	for i in range(1,2):
                	line_arr = sinkdat[i]
                	xpos_sink = float(line_arr[1])
                	ypos_sink = float(line_arr[2])
                	zpos_sink = float(line_arr[3])

       		location = [xpos_sink, ypos_sink, zpos_sink]

	data = pf.sphere(location, 0.1*pc_to_cm)
	print ('location =', location)

	bin_num2 = 25

	plot = ProfilePlot(data, "n_h", ['mdot', 'cs'], n_bins = bin_num2, weight_field="cell_mass")
	profile = plot.profiles[0]
	nh = profile.x
	mdot = profile['mdot']
	cs = profile['cs']	

	plot = ProfilePlot(data, "n_h", ['Radius', 'density','density2', 'Bmag', 'Bmag2'], n_bins = bin_num2, weight_field="cell_volume")
	profile = plot.profiles[0]
	rad = profile['Radius']
	rho = profile['density']
	bmag = profile['Bmag']
	bmag_rms = np.sqrt(profile['Bmag2']) 
	rho_rms = np.sqrt(profile['density2'])

	mbe = 1.182 * cs**3 * np.sqrt(1.0 / G**3 / rho) / 2.e33
	mcrit = bmag_rms * rad * rad / 2. / np.sqrt(G) / 2.e33

	print('rad = ', rad)
	print('cs = ', cs)
	print('rho = ', rho)
	print('rho_rms ', rho_rms)
	print('bmag = ', bmag)
	print('bmag_rms', bmag_rms)


	menc = []
	for nh_here in nh:
		inds = np.where(data['n_h'] > nh_here)
		mass_tot = na.sum(data['cell_mass'][inds])
		menc.append(mass_tot)
	menc = np.asarray(menc)	

	volume = (4./3.) * 3.14159 * np.power(rad,3)
	rho_avg = menc / volume

	beta = 8i * 3.14159 * rho * cs*cs2 / (bmag*bmag).

	print('volume = ', volume)
	print('rho_avg = ', rho_avg)
	print('mbe =', mbe)
	print('menc = ', menc)
	print('beta = ', beta)

	return nh, rad/1.5e13, mdot, mbe, menc/2.e33, mcrit, bmag_rms, rho_rms, rho, rho_avg, cs/1.e5


#dir_name = '/nobackupp7/astacy/popiii_bfac3/'
dir_name = '/work2/00863/minerva/stampede2/popiii_bfieldA/'
if my_type == '0yr':
	datanum = '0245'
else:
	datanum = '0400'
nh, rad, mdot, mbe, menc, mcrit, bmag_rms, rho_rms, rho, rho_avg, cs = load_file(dir_name, datanum)

dir_name = '/work2/00863/minerva/stampede2/popiii_bfieldA0/'
if my_type == '0yr':
	datanum = '0240'
else:
	datanum = '0400'
nh0, rad0, mdot0, mbe0, menc0, mcrit0, bmag_rms0, rho_rms0, rho0, rho_avg0, cs0 = load_file(dir_name, datanum)

tgrowth = []
tgrowth0 = []
tfrag = []
tfrag0 = []
for i in range(len(mbe0)):
	tgrowth0.append(2.e33*menc0[i] /mdot0[i])
	tfrag0.append(2.e33*mbe0[i]/mdot0[i])
for i in range(len(mbe)):
        tgrowth.append(2.e33*menc[i] /mdot[i])
        tfrag.append(2.e33*mbe[i]/mdot[i])

'''
print('menc = ', menc)
print('mbe = ', mbe)
print('mdot = ', mdot)
print('tgrowth = ', tgrowth)
print('tfrag = ',tfrag)
'''

fsize = 14

nmin = 1.e4
nmax = 1.e14

mmin = 1.e-2
mmax = 1.e3
pl.plot(nh0, menc0,'k', linewidth=1.5, label='M (hydro)')
pl.plot(nh, menc,'k', linewidth=3.5, label='M (MHD)')
pl.plot(nh0, mbe0,'b:', linewidth=1.5, label='M$_\mathrm{BE}$ (hydro)')
pl.plot(nh, mbe,'b:', linewidth=3.5, label='M$_\mathrm{BE}$ (MHD)')
pl.plot(nh, mcrit,'r--', linewidth=2.0, label='M$_\mathrm{\Phi}$ (MHD)')
#pl.plot(nh, mcrit+mbe,'g--', linewidth=2.0, label='M$_\mathrm{\Phi}$+M$_\mathrm{BE}$ (MHD)')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('nh [cm$^{-3}$]', fontsize=fsize)
ax.set_ylabel('mass [M$_{\odot}$]', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((nmin, nmax, mmin, mmax))
pl.legend(loc='upper right', fontsize=10)
if my_type == '0yr':
	pl.title('t = 0 years')
else:
	pl.title('t = 1000 years')
if my_type == '0yr':
	pl.savefig('nprof_menc_t240.eps')
else:
	pl.savefig('nprof_menc_t400.eps')


mmin = 1.e-2
mmax = 1.e1
pl.clf()
pl.cla()
pl.plot(nh0, mbe0/menc0,'k', linewidth=1.5, label='M$_\mathrm{BE}$/M (hydro)')
pl.plot(nh, mbe/menc,'b:', linewidth=3.5, label='M$_\mathrm{BE}$/M (MHD)')
pl.plot(nh, mcrit/menc,'r--', linewidth=2.0, label='M$_\mathrm{\Phi}$/M (MHD)')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('nh [cm$^{-3}$]', fontsize=fsize)
ax.set_ylabel('mass ratio', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((nmin, nmax, mmin, mmax))
pl.legend(loc='upper right', fontsize=10)
if my_type == '0yr':
        pl.title('t = 0 years')
else:
        pl.title('t = 1000 years')
if my_type == '0yr':
        pl.savefig('nprof_mrat_t240.eps')
else:
        pl.savefig('nprof_mrat_t400.eps')

tmin = 1.e0
tmax = 1.e4
pl.clf()
pl.cla()
pl.plot(nh0, tgrowth0,'k', linewidth=1.5, label='t$_\mathrm{growth}$ hydro')
pl.plot(nh, tgrowth,'k', linewidth=3.5, label='t$_\mathrm{growth}$ MHD')
pl.plot(nh0, tfrag0,'b:', linewidth=1.5, label='t$_\mathrm{frag}$ hydro')
pl.plot(nh, tfrag,'b:', linewidth=3.5, label='t$_\mathrm{frag}$ MHD')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('nh [cm$^{-3}$]', fontsize=fsize)
ax.set_ylabel('time [yr]', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((nmin, nmax, tmin, tmax))
pl.legend(loc=2, fontsize=10)
if my_type == '0yr':
        pl.title('t = 0 years')
        pl.savefig('nprof_time_t240.eps')
else:
        pl.title('t = 1000 years')
        pl.savefig('nprof_time_t400.eps')


