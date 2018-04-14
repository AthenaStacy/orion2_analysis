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

	#print data['Bmag'].units, data['Bmag'].units.dimensions
	#print data['MagEnergy'].units, data['MagEnergy'].units.dimensions

	bin_num1 = 5000
	bin_num2 = 500
	nmin = 1.e-20
	nmax = 1.e-9

        plot = ProfilePlot(data, "Radius", ['mdot', 'mbe', 'mcrit'], n_bins = bin_num2, weight_field="cell_mass")

        profile = plot.profiles[0]

	rad = profile.x
	mdot = profile['mdot']
	mbe = profile['mbe']
	mcrit = profile['mcrit']
	
	all_gas = pf.all_data()
	dens_gas = all_gas.cut_region(["obj['cell_mass'].in_units('g') < 1e32"])
	
	plot2 = ProfilePlot(dens_gas, "Radius", ['cell_mass'], n_bins = bin_num2, weight_field=None, accumulation=True)

	profile2 = plot2.profiles[0]

	rad = profile2.x
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

	for i in range(len(menc_out)):
		if menc_out[i] > 50.0:
			imass = i
			break
	rad_here = rad[imass] 

	rad_here = np.ones_like(data['x']) * rad_here.value
	rad_here = na.mean(rad_here)
	#rad_here = rad_here.value
	#print 'rad_here = ', rad_here
	
	#dens_gas2 = pf.h.sphere(location, 0.1*pc_to_cm)	
	dens_gas2 = pf.h.sphere(location, rad_here)
	#dens_gas2 = dens_gas.cut_region(["obj['Radius'].in_units('cm') < 1.5e15"])
	bmag = dens_gas2['Bmag']
	mass = dens_gas2['cell_mass']
	rho = dens_gas2['density']
	bmag_weighted = bmag*mass
	rho_weighted = rho*mass	

	mass_tot = na.sum(mass)
	
	rho_weighted = rho_weighted/mass_tot
	rho_avg = na.sum(rho_weighted)
        bmag_weighted = bmag_weighted/mass_tot
        bmag_avg = na.sum(bmag_weighted)

	time_here = pf.current_time / 3.14e7
	
	return time_here.value, bmag_avg.value, rho_avg.value,  rad_here.value, mdot_out, mbe_out, menc_out[imass], mcrit_out


dir_name = '/nobackupp7/astacy/popiii_bfieldA/'

#file_start = 560
#file_end = 570
#file_start = 230
file_start = 150
file_end = 570

time_arr = []
b_arr = []
rho_arr = []
rad_arr = []
menc_arr = []

tacc0 = 0

for i in range(file_start, file_end, 10):
	datanum = '0' + str(i)
	time, bmag, rho, rad, mdot, mbe, menc, mcrit = load_file(dir_name, datanum)

        if (i == 240):
                tacc0 = time

	time_arr.append(time)
	b_arr.append(bmag)
	rho_arr.append(rho)
	rad_arr.append(rad)
	menc_arr.append(menc)

for i in range(len(time_arr)):
        time_arr[i] = time_arr[i] - tacc0

print 'rad_arr = ', rad_arr
print 'rho_arr = ', rho_arr
print 'menc_arr = ', menc_arr
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

bmin = 1.e-2
bmax = 1.e0
tmin = min(time_arr) - (min(time_arr) % 100)
tmax = max(time_arr) + (100 - max(time_arr) % 100)

#pl.subplot(222)
pl.plot(time_arr, b_arr,'k', linewidth=2.5, label=r'B$_{\rm mag}$')
pl.plot(time_arr, bfit, '--', linewidth=2.5, label='powerlaw fit')
ax = pl.gca()
#ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('radius [AU]', fontsize=fsize)
ax.set_ylabel(r'B$_{\rm mag}$ [M$_{\odot}$]', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((tmin, tmax, bmin, bmax))
pl.legend(loc=2, fontsize=10)
#pl.title('t$_\mathrm{acc}$ = 0 years')
#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
pl.savefig('bmenc_vs_t.eps')

