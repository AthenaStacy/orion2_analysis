import yt
from yt.units import dimensions
import matplotlib
matplotlib.use('ps')
from yt.mods import *
import numpy as na
import pylab as pl
#from string import rstrip
import fields_bfield
import tracer_def

MRH = 0.76
h2conv = (2.0/1.3158)
heconv = (4.0/1.3158)
pc_to_cm = 3.08567758e18
G = YTQuantity(6.67e-8, 'cm**3/g/s**2')

my_type = '0yr'
#my_type = '1000yr'

def load_file(dir, datanum):

	pf = load(dir + 'data.' + datanum + ".3d.hdf5")

	if my_type == '0yr':
		value, location = pf.find_max("density")

	with open(dir + 'data.' + datanum + '.3d.sink', "r") as f:
		sinkdat = [line.split() for line in f]

	if my_type != '0yr':
        	for i in range(1,2):
                	line_arr = sinkdat[i]
                	xpos_sink = float(line_arr[1])
                	ypos_sink = float(line_arr[2])
                	zpos_sink = float(line_arr[3])

       		location = [xpos_sink, ypos_sink, zpos_sink]

	data = pf.sphere(location, 0.1*pc_to_cm)
	print ('location =', location)

	bin_num2 = 100

	#plot = ProfilePlot(data, "Radius", ['mdot', 'cs', 'density','Bmag'], n_bins = bin_num2, weight_field="cell_mass")
	plot = ProfilePlot(data, "n_h", ['cs', 'density','Temp'], n_bins = bin_num2, weight_field="cell_volume")

	profile = plot.profiles[0]

	nh = profile.x
	cs = profile['cs']	
	rho = profile['density']
	temp = profile['Temp']

	mbe = 1.182 * cs**3 * np.sqrt(1.0 / G**3 / rho) / 2.e33

	print('nh = ', nh)
	print('cs = ', cs)
	print('rho = ', rho)
	print('temp =', temp)

	return nh, cs, rho, temp

dir_name = '/work2/00863/minerva/stampede2/popiii_bfieldA/'
if my_type == '0yr':
	datanum = '0245'
else:
	datanum = '0400'
nh, cs, rho, temp = load_file(dir_name, datanum)

dir_name = '/work2/00863/minerva/stampede2/popiii_bfieldA0/'
if my_type == '0yr':
	datanum = '0240'
else:
	datanum = '0400'
nh0, cs0, rho0, temp0  = load_file(dir_name, datanum)

tmin = 300
tmax = 1300
nmin = 1.e5
nmax = 1.e13

fsize = 14.0

pl.clf()
pl.gca()
pl.plot(nh0, temp0,'k', linewidth=2.0, label='hydro')
pl.plot(nh, temp,'b--', linewidth=2.0, label='MHD')
ax = pl.gca()
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('n$_{\rm H}$ [cm$^{-3}$]', fontsize=fsize)
ax.set_ylabel('Temp [K]', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((nmin, nmax, tmin, tmax))
pl.legend(loc=1, fontsize=fsize)
pl.savefig('nprof_temp.eps' , bbox_inches='tight')
