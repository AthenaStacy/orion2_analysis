import yt
from yt.units import dimensions
import matplotlib
matplotlib.use('ps')
from yt.mods import *
import numpy as na
import pylab as pl
from string import rstrip
import fields_bfield
import tracer_def

MRH = 0.76
h2conv = (2.0/1.3158)
heconv = (4.0/1.3158)
pc_to_cm = 3.08567758e18

def load_file(dir, datanum):

	pf = load(dir + 'data.' + datanum + ".3d.hdf5")
	#value, location = pf.h.find_max("density")
	#value, location = pf.h.find_max("cell_mass")

	#print data['Bmag'].units, data['Bmag'].units.dimensions
	#print data['MagEnergy'].units, data['MagEnergy'].units.dimensions

        with open(dir + 'data.' + datanum + '.3d.sink', "r") as f:
                sinkdat = [map(float, line.split()) for line in f]

        for i in range(1,2):
                line_arr = sinkdat[i]
                xpos_sink = (line_arr[1])
                ypos_sink = (line_arr[2])
                zpos_sink = (line_arr[3])

        location = [xpos_sink, ypos_sink, zpos_sink]
        data = pf.h.sphere(location, 0.1*pc_to_cm)
        print ('location =', location)

	bin_num1 = 5000
	bin_num2 = 500
	nmin = 1.e-20
	nmax = 1.e-9

        plot = ProfilePlot(data, "Radius", ['mdot', 'mbe', 'mcrit'], n_bins = bin_num2, weight_field="cell_mass")

        profile = plot.profiles[0]

        print "profile data = ", profile

	rad = profile.x
	mdot = profile['mdot']
	mbe = profile['mbe']
	mcrit = profile['mcrit']
	
	plot2 = ProfilePlot(data, "Radius", ['cell_mass'], n_bins = bin_num2, weight_field=None, accumulation=True)
	profile2 = plot2.profiles[0]
	menc = profile2['cell_mass']

        rad_out = []
        mdot_out = []
        mbe_out = []
        menc_out = []
        mcrit_out = []

        for i in range(len(mbe)):
        	if (mbe[i] == mbe[i]):
                	rad_out.append(rad[i]/1.5e13)
                	mdot_out.append(mdot[i])
                	mbe_out.append(mbe[i])
                	menc_out.append(menc[i]/2.e33)
                	mcrit_out.append(mcrit[i])

	return rad_out, mdot_out, mbe_out, menc_out, mcrit_out


#dir_name = '/nobackupp7/astacy/popiii_bfac3/'
dir_name = '/nobackupp7/astacy/popiii_bfieldA/'
#datanum = '0245'
datanum = '0400'
rad, mdot, mbe, menc, mcrit = load_file(dir_name, datanum)

dir_name = '/nobackupp7/astacy/popiii_bfieldA0/'
#datanum = '0240'
datanum = '0400'
rad0, mdot0, mbe0, menc0, mcrit0 = load_file(dir_name, datanum)

tgrowth = []
tgrowth0 = []
tfrag = []
tfrag0 = []
for i in range(len(mbe0)):
	tgrowth0.append(menc0[i] /mdot0[i])
	tfrag0.append(mbe0[i]/mdot0[i])
for i in range(len(mbe)):
        tgrowth.append(menc[i] /mdot[i])
        tfrag.append(mbe[i]/mdot[i])


fsize = 14

rmin = 1.e1
rmax = 1.e4

mmin = 1.e-2
mmax = 1.e3
#pl.subplot(222)
pl.plot(rad0, menc0,'k', linewidth=1.5, label='M$_\mathrm{enc}$ hydro')
pl.plot(rad, menc,'k', linewidth=3.5, label='M$_\mathrm{enc}$ MHD')
pl.plot(rad0, mbe0,'b:', linewidth=1.5, label='M$_\mathrm{BE}$ hydro')
pl.plot(rad, mbe,'b:', linewidth=3.5, label='M$_\mathrm{BE}$ MHD')
pl.plot(rad, mcrit,'r--', linewidth=2.0, label='M$_\mathrm{crit}$ MHD')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel('mass [M$_{\odot}$]', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((rmin, rmax, mmin, mmax))
pl.legend(loc=2, fontsize=10)
#pl.title('t$_\mathrm{acc}$ = 0 years')
pl.title('t$_\mathrm{acc}$ = 300 years')
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
pl.savefig('radprof_menc2.eps')

