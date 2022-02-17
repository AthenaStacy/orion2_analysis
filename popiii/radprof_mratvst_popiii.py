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
	data.set_field_parameter("center", location)

	bin_num2 = 100

	#plot = ProfilePlot(data, "Radius", ['mdot', 'cs', 'density','Bmag'], n_bins = bin_num2, weight_field="cell_mass")
	plot = ProfilePlot(data_0, "Radius", ['mdot', 'cs'], n_bins = bin_num2, weight_field="cell_mass")
	profile = plot.profiles[0]
	rad = profile.x
	mdot = profile['mdot']
	cs = profile['cs']	

	plot = ProfilePlot(data_0, "Radius", ['density','density2', 'Bmag', 'Bmag2'], n_bins = bin_num2, weight_field="cell_volume")
	profile = plot.profiles[0]
	rho = profile['density']
	bmag = profile['Bmag']
	bmag_rms = np.sqrt(profile['Bmag2']) 
	rho_rms = np.sqrt(profile['density2'])

	mbe = 1.182 * cs**3 * np.sqrt(1.0 / G**3 / rho) / 2.e33
	mcrit = bmag * rad * rad / 2. / np.sqrt(G) / 2.e33

	plot2 = ProfilePlot(data_0, "Radius", ['cell_mass', 'cs_mass'], n_bins = bin_num2, weight_field=None, accumulation=True)
	profile2 = plot2.profiles[0]
	menc = profile2['cell_mass']
	cs_avg = profile2['cs_mass']
	cs_avg = cs_avg / menc

	volume = (4./3.) * 3.14159 * np.power(rad,3)
	rho_avg = menc / volume 	

	mrat_max = np.zeros(np.shape(rad))
	mcrit_max = np.zeros(np.shape(rad))
	mphi_max = np.zeros(np.shape(rad))
	rho_max = np.zeros(np.shape(rad))
	temp_max = np.zeros(np.shape(rad))
	beta_max = np.zeros(np.shape(rad))

	ind = np.where(np.logical_and(data['Radius'] < rad[0], data['Temp'] > 0))
	if (len(data['Radius'][ind] > 0)):
		mrat_here = data['cell_mass'][ind] / data['mbe'][ind]
		mrat_here = mrat_here.flatten()
		mcrit_here = data['cell_mass'][ind] /(data['mphi'][ind] + data['mbe'][ind])
		mcrit_here = mcrit_here.flatten()
		mphi_here = data['cell_mass'][ind] /(data['mphi'][ind])
		mphi_here = mphi_here.flatten()
		rho_here = data['density'][ind]
		rho_here = rho_here.flatten()
		temp_here = data['Temp'][ind]
		temp_here = temp_here.flatten()
		beta_here = data['beta'][ind]
		beta_here = beta_here.flatten()

		mcrit_max[0] = np.amax(mcrit_here)
		imax = np.argmax(mcrit_here)
		mrat_max[0] = mrat_here[imax]
		mphi_max[0] = mphi_here[imax]
		rho_max[0] = rho_here[imax]
		temp_max[0] = temp_here[imax]
		beta_max[0] = beta_here[imax]

	for i in range(1,len(rad)):
		ind =  np.where(np.logical_and(np.logical_and(data['Radius'] > rad[i-1], data['Temp'] > 0), data['Radius'] < rad[i]) == True)
		if (len(data['Radius'][ind] > 0)):
			mrat_here = data['cell_mass'][ind] / data['mbe'][ind]
			mrat_here = mrat_here.flatten()
			mcrit_here = data['cell_mass'][ind] / (data['mphi'][ind] + data['mbe'][ind])
			mcrit_here = mcrit_here.flatten()
			mphi_here = data['cell_mass'][ind] /(data['mphi'][ind])
			mphi_here = mphi_here.flatten()
			rho_here = data['density'][ind]
			rho_here = rho_here.flatten()
			temp_here = data['Temp'][ind]
			temp_here = temp_here.flatten()
			beta_here = data['beta'][ind]
			beta_here = beta_here.flatten()

			mcrit_max[i] = np.amax(mcrit_here)
			imax = np.argmax(mrat_here)
			mrat_max[i] = mrat_here[imax]
			mphi_max[i] = mphi_here[imax]
			rho_max[i] = rho_here[imax]
			temp_max[i] = temp_here[imax]
			beta_max[i] = beta_here[imax]

	print('rad = ',rad)
	print('mrat_max = ', mrat_max)
	print('mcrit_max = ', mcrit_max)
	print('menc = ',  menc/2.e33)

	return rad/1.5e13, time, mbe, menc/2.e33, mcrit, mcrit_max, mrat_max, mphi_max, rho, rho_avg, rho_max, temp_max, beta_max

fsize = 14

rmin = 1.e1
rmax = 1.e4

xline = [1e1, 1e2, 1e3, 1e4]
xline = np.asarray(xline)
yline = [1, 1, 1, 1]
yline = np.asarray(yline)

mmin = 1.e-1
mmax = 2.0

file_start = 510
file_end = 530


for i in range (file_start, file_end, 10):
	datanum = '0' + str(i)

	'''
	dir_name = '/work2/00863/minerva/stampede2/popiii_bfieldA/'
	rad, time, mbe, menc, mcrit, mcrit_max, mrat_max, mphi_max, rho, rho_avg, rho_max, temp_max, beta_max0 = load_file(dir_name, datanum)
	'''

	dir_name = '/work2/00863/minerva/stampede2/popiii_bfieldA0/'
	rad0, time0, mbe0, menc0, mcrit0, mcrit_max0, mrat_max0, mphi_max0, rho0, rho_avg0, rho_max0, temp_max0, beta_max0 = load_file(dir_name, datanum)

	imax = np.argmax(mcrit_max)
	mrat_max_max = mrat_max[imax]
	rho_max_max = rho_max[imax]
	temp_max_max = temp_max[imax]
	mcrit_max_max = mcrit_max[imax]
	mphi_max_max = mphi_max[imax]
	beta_max_max = beta_max[imax]

	imax = np.argmax(mcrit_max0)
	mrat_max_max0 = mrat_max0[imax]
	rho_max_max0 = rho_max0[imax]
	temp_max_max0 = temp_max0[imax]
	mcrit_max_max0 = mcrit_max0[imax]
	mphi_max_max0 = mphi_max0[imax]
	beta_max_max0 = beta_max0[imax]

	f = open("mrat.txt", "a")
	line = str(i) + ' ' + str(time.value) + ' ' + str(mrat_max_max) + ' ' + str(mcrit_max_max) + ' ' + str(mphi_max_max)+ ' ' + str(rho_max_max) + ' ' + str(temp_max_max) + ' ' + str(beta_max_max)
	f.write(line)
	f.write('\n')
	f.close()

	f = open("mrat0.txt", "a")
	line = str(i) + ' ' + str(time0.value) + ' ' + str(mrat_max_max0) + ' ' + str(mcrit_max_max0) + ' ' + str(mphi_max_max0) + ' ' + str(rho_max_max0) + ' ' + str(temp_max_max0) + ' ' + str(beta_max_max0)
	f.write(line)
	f.write('\n')
	f.close()

	'''
	time_acc = time.value/3.14e7 - 8830.0
	time_str = str(int(time_acc))
	title = 't = {0} years'.format(time_str)
	figname = 'radprof_mrat_{0}_l6.eps'.format(str(i))

	pl.clf()
	pl.cla()
	pl.plot(rad0, mrat_max0,'k', linewidth=1.5, label='M/M$_\mathrm{crit}$ (hydro)')
	pl.plot(rad, mrat_max,'b:', linewidth=3.5, label='M/M$_\mathrm{BE}$ (MHD)')
	pl.plot(rad, mcrit_max,'r--', linewidth=2.0, label='M/M$_\mathrm{crit}$ (MHD)')
	pl.plot(xline, yline,'g:', linewidth=2.0, label='M/M$_\mathrm{crit}$ = 1')
	ax = pl.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel('radius [AU]', fontsize=fsize)
	ax.set_ylabel('max M/M$_\mathrm{crit}$', fontsize=fsize)
	pl.xticks(fontsize=fsize)
	pl.yticks(fontsize=fsize)
	pl.axis((rmin, rmax, mmin, mmax))
	pl.legend(loc='upper right', fontsize=10)
	pl.title(title)
	pl.savefig(figname)
	'''
