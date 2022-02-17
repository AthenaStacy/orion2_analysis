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
cm = YTQuantity(1, 'cm')

#my_type = '0yr'
my_type = '1000yr'

def load_file(dir, datanum):

	pf = load(dir + 'data.' + datanum + ".3d.hdf5")

	time = pf.current_time

	with open(dir + 'data.' + datanum + '.3d.sink', "r") as f:
		sinkdat = [line.split() for line in f]

	if int(datanum) > 245:
        	for i in range(1,2):
                	line_arr = sinkdat[i]
                	xpos_sink = float(line_arr[1])
                	ypos_sink = float(line_arr[2])
                	zpos_sink = float(line_arr[3])

       		location = YTArray([xpos_sink, ypos_sink, zpos_sink], "cm")
	else:
		value, location = pf.find_max("density")

	data_0 = pf.sphere(location, 0.1*pc_to_cm)
	print ('location =', location)

	#dx on finest grid is 47083703369140.625 cm
	#data = data_0#.cut_region(["obj['dx'] > 9.5e13"]) 

	x_left = location[0] - 1.496e+16*cm
	y_left = location[1] - 1.496e+16*cm
	z_left = location[2] - 1.496e+16*cm

	data  = pf.covering_grid(level=6, left_edge=[x_left,y_left,z_left], dims=pf.domain_dimensions * 2)
	#data = data_0

	data.set_field_parameter("center", location)

	jmax = (1.0/4.0)

	ind = np.where(data['Temp'] > 0)
	mrat_here = data['cell_mass'][ind] / data['mbe'][ind]
	mrat_here = mrat_here.flatten()
	mcrit_here = data['cell_mass'][ind] /(data['mphi'][ind] + data['mbe'][ind])
	mcrit_here = mcrit_here.flatten()
	mphi_here = data['cell_mass'][ind] /(data['mphi'][ind] )
	mphi_here = mphi_here.flatten()
	rho_crit = 3.14159 * jmax * jmax * data['cs'][ind] * data['cs'][ind]  / G / data['dx'][ind]  / data['dx'][ind]  * (1.0 + 0.74 / data['beta'][ind] )
	rho_crit_here = rho_crit.flatten()
	rho_here = data['density'][ind] 
	rho_here = rho_here.flatten()
	temp_here = data['Temp'][ind] 
	temp_here = temp_here.flatten()
	beta_here = data['beta'][ind] 
	beta_here = beta_here.flatten()

	mcrit_max = np.amax(mcrit_here)
	imax = np.argmax(mcrit_here)
	mrat_max = mrat_here[imax]
	mphi_max = mphi_here[imax]
	rho_crit_max = rho_crit_here[imax]
	rho_max = rho_here[imax]
	temp_max = temp_here[imax]
	beta_max = beta_here[imax]

	print('mrat_max = ', mrat_max)
	print('mcrit_max = ', mcrit_max)

	return time, mcrit_max.value, mrat_max.value, mphi_max.value, rho_crit_max.value, rho_max.value, temp_max.value, beta_max.value

fsize = 14

rmin = 1.e1
rmax = 1.e4

xline = [1e1, 1e2, 1e3, 1e4]
xline = np.asarray(xline)
yline = [1, 1, 1, 1]
yline = np.asarray(yline)

mmin = 1.e-1
mmax = 2.0

file_start = 570
file_end = 620


for i in range (file_start, file_end, 10):
	datanum = '0' + str(i)

	dir_name = '/work2/00863/minerva/stampede2/popiii_bfieldA/'
	time, mcrit_max, mrat_max, mphi_max, rho_crit_max, rho_max, temp_max, beta_max = load_file(dir_name, datanum)

	dir_name = '/work2/00863/minerva/stampede2/popiii_bfieldA0/'
	time0, mcrit_max0, mrat_max0, mphi_max0, rho_crit_max0, rho_max0, temp_max0, beta_max0 = load_file(dir_name, datanum)

	f = open("mrat.txt", "a")
	line = str(i) + ' ' + str(time.value) + ' ' + str(mrat_max) + ' ' + str(mcrit_max) + ' ' + str(mphi_max) + ' ' + str(rho_crit_max) + ' ' + str(rho_max) + ' ' + str(temp_max) + ' ' + str(beta_max)
	f.write(line)
	f.write('\n')
	f.close()

	f = open("mrat0.txt", "a")
	line = str(i) + ' ' + str(time0.value) + ' ' + str(mrat_max0) + ' ' + str(mcrit_max0) + ' ' + str(mphi_max0) + ' ' + str(rho_crit_max0) + ' ' + str(rho_max0) + ' ' + str(temp_max0) + ' ' + str(beta_max0)
	f.write(line)
	f.write('\n')
	f.close()

