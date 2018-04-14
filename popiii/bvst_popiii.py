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

	dens_gas1 = all_gas.cut_region(["obj['density'] > 1e-16"])
	dens_gas2 = dens_gas1.cut_region(["obj['cell_mass'].in_units('g') < 1e32"])
	#dens_gas = dens_gas2.cut_region(["obj['Bmag'].in_units('gauss') < 1e2"])
        dens_gas = dens_gas2.cut_region(["obj['Bmag'].in_units('gauss') < obj['bline']"])

        dens_gas_alt = dens_gas.cut_region(["obj['density'] < 1.1e-16"])

	den = dens_gas['density']
	bmag = dens_gas['Bmag']
        U_B  = dens_gas['MagEnergy']
        U_K = dens_gas['KinEnergy']
        U_T = dens_gas['ThermalEnergy']
	mass = dens_gas['cell_mass']
        mbe = dens_gas_alt['mbe']
        mcrit = dens_gas_alt['mcrit']
        mass_alt = dens_gas_alt['cell_mass']

	bmag_weighted = bmag*mass
	U_B_weighted = U_B*mass
	U_K_weighted = U_K*mass
	U_T_weighted = U_T*mass
        mbe_weighted = mbe*mass_alt
        mcrit_weighted = mcrit*mass_alt
	mass_tot = na.sum(mass)
        mass_tot_alt = na.sum(mass_alt)

	bmag_weighted = bmag_weighted/mass_tot
	U_B_weighted = U_B_weighted/mass_tot
        U_K_weighted = U_K_weighted/mass_tot
        U_T_weighted = U_T_weighted/mass_tot
        mbe_weighted = mbe_weighted/mass_tot_alt
        mcrit_weighted = mcrit_weighted/mass_tot_alt

	#bmag_avg = na.mean(bmag)
	bmag_avg = na.sum(bmag_weighted)
	ub_avg = np.sum(U_B_weighted)
        uk_avg = np.sum(U_K_weighted)
        ut_avg = np.sum(U_T_weighted)
        mbe_avg = np.sum(mbe_weighted)
        mcrit_avg = np.sum(mcrit_weighted)
        #mbe_avg = np.max(mbe)
        #mcrit_avg = np.max(mcrit)
	time_here = pf.current_time / 3.14e7

	return bmag_avg.value, ub_avg.value, uk_avg.value, ut_avg.value, mbe_avg.value, mcrit_avg.value, mass_tot.value/2.e33, time_here.value
##############################################################################

dir_name = '/nobackupp7/astacy/popiii_bfieldA/'
bmag_arr = []
ub_arr = []
uk_arr = []
ut_arr = []
mbe_arr = []
mcrit_arr = []
mu_phi_arr = []
menc_arr = []
time_arr = []

file_start = 230
file_end = 570 
tacc0 = 0
time =0

for i in range (file_start, file_end, 10):
	filenum = '0' + str(i)

	bmag, ub, uk, ut, mbe, mcrit, menc, time = load_file(dir_name, filenum)

	if (i == 240):
		tacc0 = time

	print 'i = ', i, 'bmag = ', bmag, 'mbe = ', mbe, 'mcrit = ', mcrit

	bmag_arr.append(bmag)
	ub_arr.append(ub)
        uk_arr.append(uk)
        ut_arr.append(ut)
        mbe_arr.append(mbe)
        mcrit_arr.append(mcrit)
        mu_phi_arr.append(menc/mcrit)
        menc_arr.append(menc)
	time_arr.append(time)

for i in range(len(time_arr)):
	time_arr[i] = time_arr[i] - tacc0


print 'bmag_arr =', bmag_arr
print 'ub_arr =', ub_arr
print 'uk_arr =', uk_arr
print 'ut_arr =', ut_arr
print 'mu_phi = ', menc/mcrit
print 'menc = ', menc
print 'mcrit = ', mcrit
print 'mbe = ', mbe
print 'time_arr =', time_arr

bmin = 1.e-4
bmax = 2.e0
tmin = min(time_arr) - (min(time_arr) % 100) 
tmax = max(time_arr) + (100 - max(time_arr) % 100)

mmin = min(mu_phi_arr) - 0.3
mmax = max(mu_phi_arr) + 0.3

fsize = 14
#pl.subplot(121)
#pl.tight_layout()
#pl.plot(time_arr, bmag_arr,'k')
pl.plot(time_arr, ub_arr,'c', linewidth=2.0, label='Magnetic Energy')
pl.plot(time_arr, uk_arr,'k--', linewidth=2.0, label='Kinetic Energy')
pl.plot(time_arr, ut_arr,'g:', linewidth=2.0, label='Thermal Energy')
ax = pl.gca()
#ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('time [yr]', fontsize=fsize)
ax.set_ylabel('energy density [erg cm$^{-3}$]', fontsize=fsize)
pl.xticks(fontsize=fsize)
#ax.xaxis.set_ticklabels([])
pl.legend(loc=2, fontsize=10)
pl.yticks(fontsize=fsize)
pl.axis((tmin, tmax, bmin, bmax))

pl.savefig('bfield_vs_t.eps')

pl.clf()
#pl.subplot(122)
pl.plot(time_arr, mu_phi_arr, linewidth=2.0)
ax = pl.gca()
#ax.set_yscale('log')
ax.set_xlabel('time [yr]', fontsize=fsize)
ax.set_ylabel(r'$\mu_{\phi}$', fontsize=fsize)
pl.xticks(fontsize=fsize)
#pl.legend(loc=2, fontsize=10)
pl.yticks(fontsize=fsize)
pl.axis((tmin, tmax, mmin, mmax))

pl.savefig('muphi_vs_t.eps')

#bmag_arr = [0.520534593083 gauss, 0.721367132429 gauss, 0.894877066294 gauss, 0.824716082478 gauss, 0.968238755129 gauss]
#ub_arr = [0.380665076417 erg/cm**3, 0.808329743378 erg/cm**3, 1.15098349644 erg/cm**3, 0.867740079737 erg/cm**3, 1.24570058491 erg/cm**3]
#uk_arr = [0.926342298059 erg/cm**3, 1.53808668001 erg/cm**3, 1.59200436889 erg/cm**3, 1.22684228904 erg/cm**3, 1.3440641372 erg/cm**3]
#ut_arr = [0.129430245508 erg/cm**3, 0.145671992685 erg/cm**3, 0.145619471621 erg/cm**3, 0.13988701511 erg/cm**3, 0.124565690316 erg/cm**3]
#time_ar: = [-224.562274413 code_time, -109.402133688 code_time, 0.0 code_time, 115.160140725 code_time, 230.320281449 code_time]
