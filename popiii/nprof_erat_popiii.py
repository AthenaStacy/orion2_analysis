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

	#p = pc.add_profile_sphere(2., "pc", ['density', 'MagEnergy', 'ThermalEnergy', 'KinEnergy'], weight='Density', x_bins = bin_num2, center = location, x_bounds = [nmin, nmax])
        plot = ProfilePlot(data, "density", ['MagEnergy', 'ThermalEnergy', 'KinEnergy'], n_bins = bin_num2, weight_field="cell_mass")

        profile = plot.profiles[0]

        print "profile data = ", profile

	den = profile.x
	U_B  = profile['MagEnergy']
        U_K = profile['KinEnergy']
        U_T = profile['ThermalEnergy']

	return den, U_B, U_K, U_T


#dir_name = '/nobackupp7/astacy/popiii_bfac3/'
dir_name = '/nobackupp7/astacy/popiii_bfieldA/'
datanum = '0245'
den, U_B, U_K, U_T = load_file(dir_name, datanum)

nmin = 1.e-20
nmax = 1.e-10
bmin = 0.9*min( [min(U_B), min(U_K), min(U_T)] )
bmax = 1.1*max( [max(U_B), max(U_K), max(U_T)] )

pl.plot(den, U_B,'c', linewidth=2.0, label='Magnetic Energy')
pl.plot(den, U_K,'k--', linewidth=2.0, label='Kinetic Energy')
pl.plot(den, U_T,'g:', linewidth=2.0, label='Thermal Energy')

fsize = 12

ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('density [g / cm$^3$]', fontsize=fsize)
ax.set_ylabel('energy density [erg cm$^{-3}$]', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((nmin, nmax, bmin, bmax))
pl.legend(loc=2, fontsize=10)
pl.savefig('nprof_erat.eps')

