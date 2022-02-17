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
	plot = ProfilePlot(data, "Radius", ['mdot', 'cs'], n_bins = bin_num2, weight_field="cell_mass")
	profile = plot.profiles[0]
	rad = profile.x
	mdot = profile['mdot']
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

	volume = (4./3.) * 3.14159 * np.power(rad,3)
	rho_avg = menc / volume 	
	mbe_alt = 1.86 * cs_avg**3 * np.sqrt(1.0 / G**3 / rho_avg)/2.e33

	print('volume = ', volume)
	print('rho_avg = ', rho_avg)
	print('cs_avg = ', cs_avg)
	print('mbe =', mbe)
	print('mbe_alt = ', mbe_alt)
	print('mbe_alt_alt = ', mbe_alt / 1.86 * 1.182)

	rad_out = []
	mdot_out = []
	mbe_out = []
	mbe_alt_out = []
	menc_out = []
	mcrit_out = []
	mtot = []
	rho_out = []

	for i in range(len(mbe)):
		if (mbe[i] == mbe[i]):
			rad_out.append(rad[i])
			mdot_out.append(mdot[i])
			mbe_alt_out.append(mbe_alt[i])
			mbe_out.append(mbe[i])
			menc_out.append(menc[i])
			mcrit_out.append(mcrit[i])
			mtot.append(mcrit[i] + mbe[i])
			print('rad = ', rad[i]/1.5e13, 'menc = ', menc[i]/2.e33, 'mbe = ', mbe[i], 'mphi = ', mcrit[i],)

	return rad/1.5e13, mdot, mbe, mbe_alt, menc/2.e33, mcrit, bmag_rms, rho_rms, mtot, rho, rho_avg, cs/1.e5, cs_avg/1.e5


#dir_name = '/nobackupp7/astacy/popiii_bfac3/'
dir_name = '/work2/00863/minerva/stampede2/popiii_bfieldA/'
if my_type == '0yr':
	datanum = '0245'
else:
	datanum = '0400'
rad, mdot, mbe, mbe_alt, menc, mcrit, bmag_rms, rho_rms, mtot, rho, rho_avg, cs, cs_avg = load_file(dir_name, datanum)

dir_name = '/work2/00863/minerva/stampede2/popiii_bfieldA0/'
if my_type == '0yr':
	datanum = '0240'
else:
	datanum = '0400'
rad0, mdot0, mbe0, mbe_alt0, menc0, mcrit0, bmag_rms0, rho_rms0, mtot0, rho0, rho_avg0, cs0, cs_avg0 = load_file(dir_name, datanum)

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

rmin = 1.e1
rmax = 1.e4

mmin = 1.e-2
mmax = 1.e3
#pl.subplot(222)
pl.plot(rad0, menc0,'k', linewidth=1.5, label='M (hydro)')
pl.plot(rad, menc,'k', linewidth=3.5, label='M (MHD)')
pl.plot(rad0, mbe0,'b:', linewidth=1.5, label='M$_\mathrm{BE}$ (hydro)')
pl.plot(rad, mbe,'b:', linewidth=3.5, label='M$_\mathrm{BE}$ (MHD)')
pl.plot(rad, mcrit,'r--', linewidth=2.0, label='M$_\mathrm{\Phi}$ (MHD)')
#pl.plot(rad, mtot,'g-.', linewidth=2.0, label='M$_\mathrm{\Phi}$+M$_\mathrm{BE}$ (MHD)')
pl.plot(rad0, mbe_alt0,'g-.', linewidth=1.5, label='M$_\mathrm{BE}$ (hydro, rhobar)')
pl.plot(rad, mbe_alt,'g-.', linewidth=3.5, label='M$_\mathrm{BE}$ (MHD, rhobar)')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel('mass [M$_{\odot}$]', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((rmin, rmax, mmin, mmax))
pl.legend(loc=2, fontsize=10)
if my_type == '0yr':
	pl.title('t = 0 years')
else:
	pl.title('t = 1000 years')
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
if my_type == '0yr':
	pl.savefig('radprof_menc_t240.eps')
else:
	pl.savefig('radprof_menc_t400.eps')

rad_sink1 = [50,50,50]
rad_sink1 = np.asarray(rad_sink1)
rad_sink2 = [275,275,275]
rad_sink2 = np.asarray(rad_sink2)
m_sink = [1.e-2,1.e0,1.e1]
m_sink= np.asarray(m_sink)

mmin = 1.e-2
mmax = 1.e1
pl.clf()
pl.cla()
pl.plot(rad0, mbe0/menc0,'k', linewidth=1.5, label='M$_\mathrm{BE}$/M (hydro)')
pl.plot(rad, mbe/menc,'b:', linewidth=3.5, label='M$_\mathrm{BE}$/M (MHD)')
pl.plot(rad, mcrit/menc,'r--', linewidth=2.0, label='M$_\mathrm{\Phi}$/M (MHD)')
pl.plot(rad_sink1, m_sink,'g:', linewidth=2.0, label='closest sink (Hydro)')
pl.plot(rad_sink2, m_sink,'g--', linewidth=2.0, label='farthest sink (Hydro)')
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
        pl.title('t = 0 years')
else:
        pl.title('t = 1000 years')
if my_type == '0yr':
        pl.savefig('radprof_mrat_t240.eps')
else:
        pl.savefig('radprof_mrat_t400.eps')


tmin = 1.e0
tmax = 1.e4
pl.clf()
pl.cla()
#pl.subplot(221)
pl.plot(rad0, tgrowth0,'k', linewidth=1.5, label='t$_\mathrm{growth}$ hydro')
pl.plot(rad, tgrowth,'k', linewidth=3.5, label='t$_\mathrm{growth}$ MHD')
pl.plot(rad0, tfrag0,'b:', linewidth=1.5, label='t$_\mathrm{frag}$ hydro')
pl.plot(rad, tfrag,'b:', linewidth=3.5, label='t$_\mathrm{frag}$ MHD')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel('time [yr]', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((rmin, rmax, tmin, tmax))
pl.legend(loc=2, fontsize=10)
if my_type == '0yr':
        pl.title('t = 0 years')
        pl.savefig('radprof_time_t240.eps')
else:
        pl.title('t = 1000 years')
        pl.savefig('radprof_time_t400.eps')

cmin = 1.0
cmax = 6.0
tmin = 1.e0
tmax = 1.e4
pl.clf()
pl.cla()
pl.plot(rad0, cs0,'k', linewidth=1.5, label='cs Hydro')
pl.plot(rad, cs,'k', linewidth=3.5, label='cs MHD')
pl.plot(rad0, cs_avg0,'b:', linewidth=1.5, label='cs_bar Hydro')
pl.plot(rad, cs_avg,'b:', linewidth=3.5, label='cs_bar MHD')
ax = pl.gca()
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel('cs [km/s]', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((rmin, rmax, cmin, cmax))
pl.legend(loc='upper right', fontsize=10)
if my_type == '0yr':
        pl.title('t = 0 years')
        pl.savefig('radprof_cs_t240.eps')
else:
        pl.title('t = 1000 years')
        pl.savefig('radprof_cs_t400.eps')

dmin = 1.e-17
dmax = 1.e-11
tmin = 1.e0
tmax = 1.e4
pl.clf()
pl.cla()
pl.plot(rad0, rho0,'k', linewidth=1.5, label='rho Hydro')
pl.plot(rad, rho,'k', linewidth=3.5, label='rho MHD')
pl.plot(rad0, rho_avg0,'b:', linewidth=1.5, label='rho_bar Hydro')
pl.plot(rad, rho_avg,'b:', linewidth=3.5, label='rho_bar MHD')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel('density [g cm$^{-3}$', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((rmin, rmax, dmin, dmax))
pl.legend(loc='upper right', fontsize=10)
if my_type == '0yr':
        pl.title('t = 0 years')
        pl.savefig('radprof_rho_t240.eps')
else:
        pl.title('t = 1000 years')
        pl.savefig('radprof_rho_t400.eps')
