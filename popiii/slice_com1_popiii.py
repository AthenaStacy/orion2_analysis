import yt
from yt.units import dimensions
import matplotlib
matplotlib.use('pdf')

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

from yt.mods import *
import numpy as na
import fields_bfield
from tracer_def import *

rdir1 = '/work/00863/minerva/stampede2/popiii_bfieldA/'
#rdir = '/work/00863/minerva/stampede2/popiii_bfieldA0/'

#jlist = [245]
#jlist = [246]
jlist = [250,400,560]
#jlist = [10, 20, 30]

'''
field_name = 'density'
fmin = 1.e-16
fmax = 1.e-14
'''

'''
field_name= 'Temp'
fmin = 500.
fmax = 3000.
'''

field_name = 'log(B)'
fmin = -5
fmax = -1

'''
field_name = "Toomre"
fmin = 0.1
fmax = 10
'''
#######################################################################
##Read in sink information
res = 128.0
nlev = 6.0
pc_to_cm = 3.08567758e18
######################################################################################

#length = 0.25 #length of box we want in pc
length = 0.025

fig = plt.figure()

grid = AxesGrid(fig, (0.075,0.075,0.8,0.4),
		nrows_ncols = (1, 3),
		axes_pad = 0.05,
		label_mode = "L",
		share_all = True,
		cbar_location="right",
		cbar_mode="single",
		cbar_size="3%",
		cbar_pad="0%")

for en,j in enumerate(jlist): 

	jstr = str(j)

	if j < 100:
		datanum = '00' + jstr
	if j < 10:
		datanum = '000' + jstr
	if j >= 100:
		datanum = '0' + jstr

	if (en < 3):
		rdir = rdir1
	else:
		rdir = rdir2

	pf = load(rdir + 'data.' + datanum + '.3d.hdf5')

	rad  = 0.5*(pf.domain_right_edge[0]-pf.domain_left_edge[0])/pc_to_cm
	sink_size =  pf.domain_right_edge[0] / res / (2**(nlev-1))
	sink_size = sink_size * 2.0* 1.0
	print ('sink_size =', sink_size)

	xpos_sink = []
	ypos_sink = []
	zpos_sink = []
	mass_sink = []

	if j<250:
		read_sink = 0
	else:
		read_sink = 1

	if read_sink == 1:
		with open(rdir + 'data.' + datanum + '.3d.sink', "r") as f:
			sinkdat = [line.split() for line in f]
		line_arr = sinkdat[0]
		num_sink = int(line_arr[0])
		nmerge = int(line_arr[1])
		for i in range(0,num_sink):
			line_arr = sinkdat[i+1]
			xpos_sink.append(float(line_arr[1]))
			ypos_sink.append(float(line_arr[2]))
			zpos_sink.append(float(line_arr[3]))
			mass_sink.append(float(line_arr[0]))
		xpos_sink0 = xpos_sink[0]
		ypos_sink0 = ypos_sink[0]
		zpos_sink0 = zpos_sink[0]
		for i in range(0,num_sink):
			dis = (xpos_sink[i]-xpos_sink0)**2 + (ypos_sink[i]-ypos_sink0)**2 +  (zpos_sink[i]-zpos_sink0)**2
			dis = dis**0.5
			print ('mass =', mass_sink[i], 'x =', xpos_sink[i], ' y =', ypos_sink[i], ' z= ', zpos_sink[i], ' dis = ', dis)

	#value, location = pf.h.find_max("Bmag")
  	#value, location = pf.h.find_max("density")

	value, location = pf.find_max("density")

	p=SlicePlot(pf,'x', field_name, center=location, fontsize=14)
	p.hide_axes()
	p.set_width(length,'pc')
	p.set_zlim(field_name, fmin, fmax)
	#p.annotate_grids()
	if read_sink == 1:
		for i in range(0,num_sink):
			center = [xpos_sink[i], ypos_sink[i], zpos_sink[i]]
			p.annotate_sphere(center, sink_size, circle_args={'color':'red'})

	p.set_axes_unit('AU')

	if en == 0:
		p.annotate_title("t = 0 yr, central 5000 AU")
	if en == 1:
		p.annotate_title("t = 1000 yr, central 5000 AU")
	if en == 2:
		p.annotate_title("t = 2000 yr, central 5000 AU")

	plot = p.plots[field_name]
	plot.figure = fig
	plot.axes = grid[en].axes
	plot.cax = grid.cbar_axes[en]

	# Finally, this actually redraws the plot.
	p._setup_plots()


plt.savefig('slice_com.png')
