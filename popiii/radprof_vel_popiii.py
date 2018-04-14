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
	data = pf.h.sphere(location, 0.5*pc_to_cm)
	print ('location =', location)

	bin_num1 = 5000
	bin_num2 = 100
	nmin = 1.e-20
	nmax = 1.e-9

        plot = ProfilePlot(data, "Radius", ['velocity', 'vrot', 'vrad', 'Temp'], n_bins = bin_num2, weight_field="cell_mass")

        profile = plot.profiles[0]

        print "profile data = ", profile

	rad = profile.x
	vel = profile['velocity']
	vrot = profile['vrot']
	vrad = profile['vrad']        
	temp = profile['Temp'] 

	return rad/1.5e13, vel, vrot, vrad, temp


#dir_name = '/nobackupp7/astacy/popiii_bfac3/'
dir_name = '/nobackupp7/astacy/popiii_bfieldA/'
datanum = '0245'
rad, vel, vrot, vrad, temp = load_file(dir_name, datanum)

dir_name = '/nobackupp7/astacy/popiii_bfieldA0/'
datanum = '0240'
rad0, vel0, vrot0, vrad0, temp0 = load_file(dir_name, datanum)

fsize = 14

rmin = 1.e1
rmax = 1.e4

vmin = -5
vmax = 7
#pl.subplot(221)
pl.plot(rad0, vel0,'k', linewidth=2.5, label='v$_\mathrm{{tot}}$, hydro')
pl.plot(rad, vel,'k', linewidth=1.5, label='v$_\mathrm{tot}$, MHD')
pl.plot(rad0, vrad0,'g--', linewidth=2.5, label='v$_\mathrm{rad}$, hydro')
pl.plot(rad, vrad,'g--', linewidth=1.5, label='v$_\mathrm{rad}$, MHD')
pl.plot(rad0, vrot0,'b:', linewidth=2.5, label='v$_\mathrm{rot}$, hydro')
pl.plot(rad, vrot,'b:', linewidth=1.5, label='v$_\mathrm{rot}$, MHD')
ax = pl.gca()
ax.set_xscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel('[km s$^{-1}$]', fontsize=fsize)
#ax.xaxis.set_ticklabels([])
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((rmin, rmax, vmin, vmax))
#pl.legend(loc=3, fontsize=10)
#ax.set_title("Hydro")
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

'''
pl.subplot(222)
pl.plot(rad, vel,'k', linewidth=2.0, label='total velocity')
pl.plot(rad, vrad,'g--', linewidth=2.0, label='v$_{rad}$')
pl.plot(rad, vrot,'b:', linewidth=2.0, label='v$_{rot}$')
ax = pl.gca()
ax.set_xscale('log')
ax.set_xlabel('radius [AU]', fontsize=10)
ax.set_ylabel('[km s$^{-1}$]', fontsize=10)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((rmin, rmax, vmin, vmax))
ax.set_title("MHD")
'''
pl.savefig('radprof_vel.eps')
