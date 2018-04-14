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
	value, location = pf.h.find_max("Bmag")
	data = pf.h.sphere(location, 3.0*pc_to_cm)
	print ('location =', location)

	bin_num1 = 5000
	bin_num2 = 200
	nmin = 1.e-20
	nmax = 1.e-9

	#p = pc.add_profile_sphere(2., "pc", ['Density', 'Bmag'], weight='CellMassMsun', x_bins = bin_num2, center = location, x_bounds = [nmin, nmax])
	plot = ProfilePlot(data, "density", ['Bmag'], n_bins = bin_num2, weight_field="cell_mass")

	profile = plot.profiles[0]
	den = profile.x
	bmag = profile['Bmag']

	return den, bmag

fac = 1.e9

x_init = 1.e-20
x = []
y = []
y2 = []
for i in range (0,12):
        x_here = x_init * 10**i
        x.append(x_here)
        y.append(fac*1.e-16 * (x_here/x_init) ** 0.6666)
        y2.append(fac*1.e-14 * (x_here/x_init) ** 0.9)

print ('xline =', x)
print ('yline =', y)
print ('y2line =', y2)

dir_name =  '/nobackupp7/astacy/popiii_Bscope12/'
datanum = '0020'
den0, bmag0 = load_file(dir_name, datanum)

#dir_name = '/nobackupp7/astacy/popiii_bfac3/'
dir_name = '/nobackupp7/astacy/popiii_bfieldA/'
datanum  = '0100'
den, bmag = load_file(dir_name, datanum)

datanum = '0200'
den2, bmag2 = load_file(dir_name, datanum)

datanum = '0240'  #first sink forms
den3, bmag3 = load_file(dir_name, datanum)

datanum = '0250'
den4, bmag4 = load_file(dir_name, datanum)

datanum = '0350'
den5, bmag5 = load_file(dir_name, datanum)

nmin = 1.e-20
nmax = 1.e-10
bmin = 1.e-16*fac
bmax = 1.e-5*fac

fsize = 12

pl.plot(den0, bmag0*fac,'c', linewidth=2.0, label='initial B-field, no div. clean')
line1, = pl.plot(den, bmag,'k', linewidth=2.0, label="initial B-field")
pl.plot(den2, bmag2,'g--', linewidth=2.0, label='1000 yr before first sink forms')
line2, = pl.plot(den3, bmag3,'r:', linewidth=2.0, label='0 yr, first sink forms')
pl.plot(den4, bmag4,'-.', linewidth=2.0, label='100 yr')
line3, = pl.plot(den5, bmag5,'m--', linewidth=2.0, label='1000 yr')
pl.plot(x, y, 'b:')
pl.plot(x, y2, 'b:')

ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('density [g / cm$^3$]', fontsize=fsize)
ax.set_ylabel('B-field magnitude [G]', fontsize=fsize)
ax.annotate(r'$\rho^{2/3}$', xy=(1e-10,1e-10*fac),xytext=(3e-12, 1e-11*fac), fontsize=12)
ax.annotate(r'$\rho^{0.9}$', xy=(1e-10,1e-10*fac),xytext=(3e-14, 3e-9*fac), fontsize=12)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((nmin, nmax, bmin, bmax))

pl.legend(loc=2, fontsize=10)

pl.show()
pl.savefig('nprof_bfield.eps')

