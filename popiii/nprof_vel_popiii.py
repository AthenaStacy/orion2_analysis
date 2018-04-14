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
	value, location = pf.h.find_max("density")
	data = pf.h.sphere(location, 3.0*pc_to_cm)
	print ('location =', location)

	bin_num1 = 5000
	bin_num2 = 200
	nmin = 1.e-20
	nmax = 1.e-9

        plot = ProfilePlot(data, "density", ['velocity', 'vrot', 'vrad', 'Temp'], n_bins = bin_num2, weight_field="cell_mass")

        profile = plot.profiles[0]

        print "profile data = ", profile

	den = profile.x
	vel = profile['velocity']
	vrot = profile['vrot']
	vrad = profile['vrad']        
	temp = profile['Temp'] 

	return den, vel, vrot, vrad, temp


#dir_name = '/nobackupp7/astacy/popiii_bfac3/'
dir_name = '/nobackupp7/astacy/popiii_bfieldA/'
datanum = '0245'
den, vel, vrot, vrad, temp = load_file(dir_name, datanum)

dir_name = '/nobackupp7/astacy/popiii_bfieldA0/'
datanum = '0240'
den0, vel0, vrot0, vrad0, temp0 = load_file(dir_name, datanum)

nmin = 1.e-16
nmax = 1.e-10


vmin = -10
vmax = 10
pl.subplot(211)
pl.plot(den0, vel0,'k', linewidth=2.0, label='total velocity')
pl.plot(den0, vrad0,'g--', linewidth=2.0, label='v$_{rad}$')
pl.plot(den0, vrot0,'b:', linewidth=2.0, label='v$_{rot}$')
ax = pl.gca()
ax.set_xscale('log')
#ax.set_xlabel('density [g / cm$^3$]', fontsize=10)
ax.set_ylabel('[km s$^{-1}$]', fontsize=10)
ax.xaxis.set_ticklabels([])
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((nmin, nmax, vmin, vmax))
pl.legend(loc=3, fontsize=10)
ax.set_title("Hydro")

pl.subplot(212)
pl.plot(den, vel,'k', linewidth=2.0, label='total velocity')
pl.plot(den, vrad,'g--', linewidth=2.0, label='v$_{rad}$')
pl.plot(den, vrot,'b:', linewidth=2.0, label='v$_{rot}$')
ax = pl.gca()
ax.set_xscale('log')
ax.set_xlabel('density [g / cm$^3$]', fontsize=10)
ax.set_ylabel('[km s$^{-1}$]', fontsize=10)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((nmin, nmax, vmin, vmax))
ax.set_title("MHD")

pl.savefig('nprof_vel.eps')
