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

######################################################################

def _nh_alt(field,data):
   return  (data["density"] * MRH / 1.67e-24 )
add_field("nh_alt",function=_nh_alt,units=r"cm^{-3}",take_log=True)

###########################################################################

def load_file(dir_name, datanum):

	pf = load(dir_name + 'data.' + datanum + ".3d.hdf5")
	value, location = pf.h.find_max("density")
	data = pf.h.sphere(location, 3.0/pf_pc)
	print 'location =', location
	pc = PlotCollection(pf, center = location)

	x=location[0]
	y=location[1]
	z=location[2]

	data['x'] = data['x'] - x
	data['y'] = data['y'] - y
	data['z'] = data['z'] - z	

	#data['Radius'] = data['x']

	rad_min = min(data['Radius'])
	rad_max = max(data['Radius'])
	rad_min = 1.e-5*rad_max
	rad_arr = []
	print 'rad_arr length = ', len(rad_arr), 'min_r =', rad_min, 'max_r = ', rad_max

	bin_num = 200


	p3 = pc.add_profile_sphere(2., "pc", ['Radius', 'nh_alt', 'Density', 'Bmag'], weight='Radius', x_bins = bin_num, center = location, x_bounds = [rad_min, rad_max])

	rad = p3.data['Radius']
	den = p3.data['Density']
	nh = p3.data['nh_alt']
	bmag = p3.data['Bmag']

	dmin = min(data['Density'])
	dmax = max(data['Density'])

	p4 = pc.add_profile_sphere(1., "pc", ['Density', 'CellMassMsun'], weight=None, accumulation=True, x_bins = bin_num, center = location, x_bounds = [dmax, dmin])

	n_menc = p4.data['Density']
	menc = p4.data['CellMassMsun']

	return rad, den, nh, bmag, n_menc, menc
##############################################################################

dir_name = '/nobackupp7/astacy/popiii_Bscope17/'

datanum = '0300'
rad, den, nh, bmag, d_menc, menc = load_file(dir_name, datanum)

datanum2 = '0500'
rad2, den2, nh2, bmag2, d_menc2, menc2 = load_file(dir_name, datanum2)

print 'rad =', rad
print 'den = ', den
print 'nh = ', nh
print 'bmag = ', bmag
print 'd_menc = ', d_menc
print 'menc = ', menc

nmin = 1.e6
nmax = 1.e14
bmin = 1.e-15
bmax = 1.e-5

dmin = 0.95*min([min(d_menc), min(d_menc2)])
dmax = 1.05*max([max(d_menc), max(d_menc2)])

#mmin = 0.9*min([min(menc), min(menc2)])
mmin = 1.e-2
mmax = 1.1*max([max(menc), max(menc2)])

pl.plot(nh, bmag,'k')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('number density [cm^{-3}]', fontsize=9)
ax.set_ylabel('B [G]', fontsize=9)
pl.xticks(fontsize=6)
#ax.xaxis.set_ticklabels([])
pl.yticks(fontsize=10)
pl.axis((nmin, nmax, bmin, bmax))

pl.savefig('radprof_bfield.eps')


pl.gca()
pl.plot(d_menc, menc,'k', label='0 yr, first sink forms')
pl.plot(d_menc2, menc2, label='6000 yr')
ax.set_xlabel(r'density [g cm^{-3}]', fontsize=9)
ax.set_ylabel(r'enclosed mass [M_{\odot}]', fontsize=9)
pl.xticks(fontsize=10)
#ax.xaxis.set_ticklabels([])
pl.yticks(fontsize=10)
pl.legend(loc=3, fontsize=10)
pl.axis((dmin, dmax, mmin, mmax))

pl.savefig('nprof_menc.eps')



