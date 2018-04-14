import matplotlib
matplotlib.use('ps')

from yt.mods import *
import numpy as na
import pylab as pl
from string import rstrip
#import fields
import fields_bfield
from tracer_def import *

MRH = 0.76
pf_pc = 3e-18


###########################################################################

def load_file(dir_name, datanum):

	pf = load(dir_name + 'data.' + datanum + ".3d.hdf5")

	all_gas = pf.all_data()
	dens_gas1 = all_gas.cut_region(["obj['density'] > 1e-14"])
	dens_gas2 = dens_gas1.cut_region(["obj['cell_mass'].in_units('g') < 1e32"])
	dens_gas = dens_gas2.cut_region(["obj['Temp'] < 1500"])

	den = dens_gas['density']
	mass = dens_gas['cell_mass']

	mass_tot = na.sum(mass) / 2e33

	time_here = pf.current_time / 3.14e7

	return mass_tot.value, time_here.value
##############################################################################

dir_name = '/nobackupp7/astacy/popiii_bfieldA/'
dmass_arr = []
time_arr = []

file_start = 230
file_end = 570 
tacc0 = 0
time =0

for i in range (file_start, file_end, 10):
	filenum = '0' + str(i)

	dmass, time = load_file(dir_name, filenum)

	if (i == 240):
		tacc0 = time

	print 'i = ', i, 'dmass = ', dmass

	dmass_arr.append(dmass)
	time_arr.append(time)

for i in range(len(time_arr)):
	time_arr[i] = time_arr[i] - tacc0

print 'dmass_arr =', dmass_arr
print 'time_arr =', time_arr

###############################################

dir_name = '/nobackupp7/astacy/popiii_bfieldA0/'
dmass0_arr = []
time0_arr = []

file_start = 230
file_end = 620
tacc0 = 0
time0 =0

for i in range (file_start, file_end, 10):
        filenum = '0' + str(i)

        dmass, time = load_file(dir_name, filenum)

        if (i == 240):
                tacc0 = time

        print 'i = ', i, 'dmass = ', dmass

        dmass0_arr.append(dmass)
        time0_arr.append(time)

for i in range(len(time0_arr)):
        time0_arr[i] = time0_arr[i] - tacc0

print 'dmass0_arr =', dmass0_arr
print 'time0_arr =', time0_arr

dmin = int(0.9 * min([min(dmass0_arr), min(dmass_arr)]))
dmax = int(1.1 * max([max(dmass0_arr), max(dmass_arr)]))
tmin = min(time_arr) - (min(time_arr) % 100)
tmax = max(time_arr) + (100 - max(time_arr) % 100)

fsize = 14

pl.plot(time0_arr, dmass0_arr, label='Hydro')
pl.plot(time_arr, dmass_arr,'k', label='MHD')
ax = pl.gca()
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('time [yr]', fontsize=fsize)
ax.set_ylabel('mass [M$_{\odot}]$', fontsize=fsize)
pl.xticks(fontsize=fsize)
#ax.xaxis.set_ticklabels([])
pl.yticks(fontsize=fsize)
pl.axis((tmin, tmax, dmin, dmax))
pl.legend(fontsize=10, loc=2)
pl.show()
pl.savefig('dmass_vs_t.eps')


#dmass0_arr = [17.7999874554 g, 18.055258389 g, 20.4843278741 g, 19.2687218869 g, 19.6947671768 g]
#time0_arr = [-218.740810177 code_time, -111.329463302 code_time, 0.0 code_time, 124.014765974 code_time, 241.135596593 code_time]

#dmass_arr = [19.7346528358 g, 22.205534941 g, 23.6702473247 g, 21.4141403638 g, 26.5017569079 g]
#time_arr = [-224.562274413 code_time, -109.402133688 code_time, 0.0 code_time, 115.160140725 code_time, 230.320281449 code_time]
