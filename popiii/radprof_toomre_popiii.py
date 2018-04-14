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

        plot = ProfilePlot(data, "Radius", ['mdot', 'density', 'Toomre', 'Temp'], n_bins = bin_num2, weight_field="cell_mass")

        profile = plot.profiles[0]

        print "profile data = ", profile

	rad = profile.x
	mdot = profile['mdot']
	den = profile['density']
	toomre = profile['Toomre']        
	temp = profile['Temp'] 

	return rad/1.5e13, mdot, den, toomre, temp


#dir_name = '/nobackupp7/astacy/popiii_bfac3/'
dir_name = '/nobackupp7/astacy/popiii_bfieldA/'
datanum = '0245'
rad, mdot, den, toomre, temp = load_file(dir_name, datanum)

dir_name = '/nobackupp7/astacy/popiii_bfieldA0/'
datanum = '0240'
rad0, mdot0, den0, toomre0, temp0 = load_file(dir_name, datanum)


rmin = 1.e1
rmax = 1.e4

dmin = 1.e-18
dmax = 1.e-10
pl.subplot(221)
pl.plot(rad0, den0,'k', linewidth=2.0, label='hydro')
pl.plot(rad, den,'b--', linewidth=2.0, label='MHD')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=10)
ax.set_ylabel('density [g cm$^{-3}$]', fontsize=10)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((rmin, rmax, dmin, dmax))
#pl.legend(loc=2, fontsize=10)

tmin = 5.e2
tmax = 1.8e3
pl.subplot(222)
pl.plot(rad0, temp0,'k', linewidth=2.0, label='hydro')
pl.plot(rad, temp,'b--', linewidth=2.0, label='MHD')
ax = pl.gca()
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=10)
ax.set_ylabel('temperature [K]', fontsize=10)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((rmin, rmax, tmin, tmax))
pl.legend(loc=1, fontsize=10)

mmin = 1.e-3
mmax = 1.e0
pl.subplot(223)
pl.plot(rad0, mdot0,'k', linewidth=2.0, label='hydro')
pl.plot(rad, mdot,'b--', linewidth=2.0, label='MHD')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=10)
ax.set_ylabel('dM/dt [M$_{\odot}$ yr$^{-1}$]', fontsize=10)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((rmin, rmax, mmin, mmax))

qmin = 0.5
qmax = 3
pl.subplot(224)
pl.plot(rad0, toomre0,'k', linewidth=2.0, label='hydro')
pl.plot(rad, toomre,'b--', linewidth=2.0, label='MHD')
ax = pl.gca()
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=10)
ax.set_ylabel(r'Q$_{\rm local}$', fontsize=10)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((rmin, rmax, qmin, qmax))

pl.tight_layout()
pl.savefig('radprof_toomre.eps')


