import yt
from yt.units import dimensions
import matplotlib
matplotlib.use('ps')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

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
	#value, location = pf.h.find_max("cell_mass")
	data = pf.h.sphere(location, 0.1*pc_to_cm)
	#print ('location =', location)
	
	all_gas = pf.all_data()

        dens_gas1 = all_gas.cut_region(["obj['density'] < 2e-16"])
        dens_gas2 = dens_gas1.cut_region(["obj['density'] >  1e-16"])

	bmag = dens_gas2['Bmag']
	mass = dens_gas2['cell_mass']
	bmag_weighted = bmag*mass
	mass_tot = na.sum(mass)
	bmag_weighted = bmag_weighted/mass_tot
	bmag_avg = na.sum(bmag_weighted)

	time_here = pf.current_time / 3.14e7
	
	return time_here.value, bmag_avg.value


dir_name = '/nobackupp7/astacy/popiii_bfieldA/'

file_start = 150
file_end = 250
#file_start = 230
#file_end = 570

time_arr = []
b_arr = []
rad_arr = []
menc_arr = []

tacc0 = 0

for i in range(file_start, file_end, 10):
	datanum = '0' + str(i)
	time, bmag = load_file(dir_name, datanum)

        if (i == 240):
                tacc0 = time

	time_arr.append(time)
	b_arr.append(bmag)

for i in range(len(time_arr)):
        time_arr[i] = time_arr[i] - tacc0

print 'time_arr= ', time_arr
print 'b_arr = ', b_arr


time_arr_shift = time_arr + min(time_arr) + 100 #make sure time values are positive

lnx = np.log(time_arr_shift)
lny = np.log(b_arr)

n = float(len(lnx))
B = (n * np.sum(lnx*lny) - np.sum(lnx)*np.sum(lny))/(n*np.sum(lnx*lnx) - (np.sum(lnx))**2)

A = (np.sum(lny) - B*np.sum(lnx))/n
A = 2.71818**A

print 'A = ', A, 'B = ', B

#bfit_exp = 

bfit = A * (time_arr_shift**B)
print 'bfit = ', bfit 

fsize = 14


bmin = (min(b_arr) - (min(b_arr) % 10))/10
bmax = (max(b_arr) - (max(b_arr) % 10)) * 10
tmin = min(time_arr) - (min(time_arr) % 100)
tmax = max(time_arr) + (100 - max(time_arr) % 100)

#pl.subplot(222)
pl.plot(time_arr, b_arr,'k', linewidth=2.5, label=r'B$_{\rm mag}$')
pl.plot(time_arr, bfit, '--', linewidth=2.5, label='powerlaw fit')
ax = pl.gca()
#ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('time [yr]', fontsize=fsize)
ax.set_ylabel(r'B$_{\rm mag}$', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((tmin, tmax, bmin, bmax))
pl.legend(loc=2, fontsize=10)
#pl.title('t$_\mathrm{acc}$ = 0 years')
#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
pl.savefig('brho_vs_t.eps')

