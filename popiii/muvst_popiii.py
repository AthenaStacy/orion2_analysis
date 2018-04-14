import matplotlib
matplotlib.use('ps')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

from yt.mods import *
import numpy as na
import pylab as pl
from string import rstrip
#import fields
import fields_bfield
from tracer_def import *

my_rho = YTQuantity(1, 'g / cm**3')
MRH = 0.76
pf_pc = 3e-18


###########################################################################

def load_file(dir_name, datanum):

	pf = load(dir_name + 'data.' + datanum + ".3d.hdf5")

	all_gas = pf.all_data()

	#dens_gas1 = all_gas.cut_region(["obj['cell_mass'].in_units('g') < 1e32"])
        dens_gas1 = all_gas.cut_region(["obj['Bmag'].in_units('gauss') < obj['bline']"])
        dens_gas2 = all_gas.cut_region(["obj['Bmag'].in_units('gauss') < obj['bline']"])

	dens_gas18  = dens_gas2.cut_region(["obj['density'] > 1.15e-18"])
	dens_gas18a = dens_gas2.cut_region(["obj['density'] < 1.1e-18"])
        dens_gas18b = dens_gas18a.cut_region(["obj['density'] >  1e-18"])

	dens_gas16  = dens_gas2.cut_region(["obj['density'] > 1.5e-16"])
        dens_gas16a = dens_gas2.cut_region(["obj['density'] < 2e-16"])
        dens_gas16b = dens_gas16a.cut_region(["obj['density'] >  1e-16"])

	dens_gas14  = dens_gas2.cut_region(["obj['density'] > 1.5e-14"])
        dens_gas14a = dens_gas2.cut_region(["obj['density'] < 2e-14"])
        dens_gas14b = dens_gas14a.cut_region(["obj['density'] >  1e-14"])

	den = dens_gas1['density']
	bmag = dens_gas1['Bmag']
	mass = dens_gas1['cell_mass']

	bmag_weighted = bmag*mass
	mass_tot = na.sum(mass)

	bmag_weighted = bmag_weighted/mass_tot

	bmag_avg = na.sum(bmag_weighted)

        mass18 = dens_gas18['cell_mass']
        mass18_tot = na.sum(mass18)
        mass16 = dens_gas16['cell_mass']
        mass16_tot = na.sum(mass16)
        mass14 = dens_gas14['cell_mass']
        mass14_tot = na.sum(mass14)

        mbe18 = dens_gas18b['mbe']
        mcrit18 = dens_gas18b['mcrit']
        mass18 = dens_gas18b['cell_mass']
        mbe16 = dens_gas16b['mbe']
        mcrit16 = dens_gas16b['mcrit']
        mass16 = dens_gas16b['cell_mass']
        mbe14 = dens_gas14b['mbe']
        mcrit14 = dens_gas14b['mcrit']
        mass14 = dens_gas14b['cell_mass']

        mbe18_weighted = mbe18*mass18
        mcrit18_weighted = mcrit18*mass18
        mass18_bin = na.sum(mass18)
        mbe16_weighted = mbe16*mass16
        mcrit16_weighted = mcrit16*mass16
        mass16_bin = na.sum(mass16)
        mbe14_weighted = mbe14*mass14
        mcrit14_weighted = mcrit14*mass14
        mass14_bin = na.sum(mass14)

        mbe18_weighted = mbe18_weighted/mass18_bin
        mcrit18_weighted = mcrit18_weighted/mass18_bin
        mbe16_weighted = mbe16_weighted/mass16_bin
        mcrit16_weighted = mcrit16_weighted/mass16_bin
        mbe14_weighted = mbe14_weighted/mass14_bin
        mcrit14_weighted = mcrit14_weighted/mass14_bin

        mbe18_avg = np.sum(mbe18_weighted)
        mcrit18_avg = np.sum(mcrit18_weighted)
        mbe16_avg = np.sum(mbe16_weighted)
        mcrit16_avg = np.sum(mcrit16_weighted)
        mbe14_avg = np.sum(mbe14_weighted)
        mcrit14_avg = np.sum(mcrit14_weighted)

	time_here = pf.current_time / 3.14e7

	return bmag_avg.value, mbe18_avg.value, mcrit18_avg.value, mass18_tot.value/2.e33, mbe16_avg.value, mcrit16_avg.value, mass16_tot.value/2.e33, mbe14_avg.value, mcrit14_avg.value, mass14_tot.value/2.e33, time_here.value
##############################################################################

dir_name = '/nobackupp7/astacy/popiii_bfieldA/'
bmag_arr = []

mbe18_arr = []
mcrit18_arr = []
mu_phi18_arr = []
menc18_arr = []

mbe16_arr = []
mcrit16_arr = []
mu_phi16_arr = []
menc16_arr = []

mbe14_arr = []
mcrit14_arr = []
mu_phi14_arr = []
menc14_arr = []

time_arr = []

file_start = 230
file_end = 570 
tacc0 = 0
time =0

for i in range (file_start, file_end, 10):
	filenum = '0' + str(i)

	bmag, mbe18, mcrit18, menc18, mbe16, mcrit16, menc16, mbe14, mcrit14, menc14, time = load_file(dir_name, filenum)

	if (i == 240):
		tacc0 = time

	print 'i = ', i, 'bmag = ', bmag, 'mbe18 = ', mbe18, 'menc18 = ', menc16, 'mcrit18 = ', mcrit18

	print 'i = ', i, 'bmag = ', bmag, 'mbe14 = ', mbe14, 'menc14 = ', menc14, 'mcrit14 = ', mcrit14

	bmag_arr.append(bmag)

        mbe18_arr.append(mbe18)
        mcrit18_arr.append(mcrit18)
        mu_phi18_arr.append(menc18/mcrit18)
        menc18_arr.append(menc18)

        mbe16_arr.append(mbe16)
        mcrit16_arr.append(mcrit16)
        mu_phi16_arr.append(menc16/mcrit16)
        menc16_arr.append(menc16)

        mbe14_arr.append(mbe14)
        mcrit14_arr.append(mcrit14)
        mu_phi14_arr.append(menc14/mcrit14)
        menc14_arr.append(menc14)

	time_arr.append(time)

for i in range(len(time_arr)):
	time_arr[i] = time_arr[i] - tacc0


print 'bmag_arr =', bmag_arr
print 'mu_phi16 = ', menc18/mcrit18, 'mu_phi16 = ', menc16/mcrit16, 'mu_phi14 = ', menc14/mcrit14
print 'menc = ', menc18, 'menc16 = ', menc16, 'menc14 = ', menc14
print 'mcrit = ', mcrit18, 'mcrit16 = ', mcrit16, 'mcrit14 = ', mcrit14
print 'mbe = ', mbe18
print 'time_arr =', time_arr

tmin = min(time_arr) - (min(time_arr) % 100) 
tmax = max(time_arr) + (100 - max(time_arr) % 100)

mmin = max(min(mu_phi14_arr) - 0.3, 0.02)
mmax = max(max(mu_phi18_arr),max(mu_phi16_arr)) + 4.0

fsize = 14

pl.clf()
pl.plot(time_arr, mu_phi18_arr, 'c', linewidth=2.0, label = r'10$^{-18}$ g cm$^{-3}$')
pl.plot(time_arr, mu_phi16_arr, 'k--', linewidth=2.0, label = r'10$^{-16}$ g cm$^{-3}$')
pl.plot(time_arr, mu_phi14_arr, 'g:', linewidth=2.0, label = r'10$^{-14}$ g cm$^{-3}$')
ax = pl.gca()
ax.set_yscale('log')
ax.set_xlabel('time [yr]', fontsize=fsize)
ax.set_ylabel(r'$\mu_{\phi}$', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.legend(loc=4, fontsize=10)
pl.yticks(fontsize=fsize)
pl.axis((tmin, tmax, mmin, mmax))

pl.savefig('muphi_vs_t.eps')

