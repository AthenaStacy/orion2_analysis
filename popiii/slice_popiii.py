import yt
from yt.units import dimensions
import matplotlib
matplotlib.use('pdf')

from yt.mods import *
import numpy as na
from string import rstrip
#import fields_iso
#import fields
import fields_bfield
import fields_wpot
from tracer_def import *

#rdir = '/nobackupp7/astacy/popiii_bfieldA/'
rdir = '/nobackupp7/astacy/popiii_bfieldA0/'
#rdir =  '/nobackupp7/astacy/popiii_chemtest/'

jstart = 335
jlist = [250,247]

#dmin = 1.e-18
#dmax = 1.e-13
#dmin = 1.e-16
#dmax = 1.e-15
dmin = 1.e-16
dmax = 1.e-14

gmin = 1.4
gmax = 1.67

tmin = 200.
tmax = 5000.

bmin = -6
bmax = -1

qmin = 0.1
qmax = 10

#######################################################################
##Read in sink information
read_sink = 1
res = 128.0
nlev = 4.0
pc_to_cm = 3.08567758e18
######################################################################################

#length = 0.25 #length of box we want in pc
length = 0.05

for x in xrange(2,3):

        j = jstart
 	#j = jlist[x-1]
	jstr = str(j)

	if j < 100:
		datanum = '00' + jstr
	if j < 10:
		datanum = '000' + jstr
        if j >= 100:
                datanum = '0' + jstr

	pf = load(rdir + 'data.' + datanum + '.3d.hdf5')
	
	rad  = 0.5*(pf.domain_right_edge[0]-pf.domain_left_edge[0])/pc_to_cm
	sink_size =  pf.domain_right_edge[0] / res / (2**(nlev-1))
	sink_size = sink_size * 2.0* 1.0
	print ('sink_size =', sink_size)

	xpos_sink = []
	ypos_sink = []
	zpos_sink = []
	mass_sink = []

	if read_sink == 1:
        	with open(rdir + 'data.' + datanum + '.3d.sink', "r") as f:
                	sinkdat = [map(float, line.split()) for line in f]
        	line_arr = sinkdat[0]
        	num_sink = int(line_arr[0])
        	nmerge = int(line_arr[1])
        	for i in range(0,num_sink):
                	line_arr = sinkdat[i+1]
                	xpos_sink.append(line_arr[1])
                	ypos_sink.append(line_arr[2])
                	zpos_sink.append(line_arr[3])
			mass_sink.append(line_arr[0])
		xpos_sink0 = xpos_sink[0]
		ypos_sink0 = ypos_sink[0]
		zpos_sink0 = zpos_sink[0]
		for i in range(0,num_sink):
			dis = (xpos_sink[i]-xpos_sink0)**2 + (ypos_sink[i]-ypos_sink0)**2 +  (zpos_sink[i]-zpos_sink0)**2
			dis = dis**0.5
			print ('mass =', mass_sink[i], 'x =', xpos_sink[i], ' y =', ypos_sink[i], ' z= ', zpos_sink[i], ' dis = ', dis)
		

	value, location = pf.h.find_max("density")
        #value, location = pf.h.find_max("Temp")
	#value, location = pf.h.find_max("Bmag")
	#value, location = pf.h.find_max("cell_mass")
	print ('max density location = ', location)

	data = pf.h.sphere(location, length*pc_to_cm)
        edgeL = [ -length*pc_to_cm, -length*pc_to_cm, -length*pc_to_cm]
        edgeR = [  length*pc_to_cm,  length*pc_to_cm,  length*pc_to_cm]
        data = pf.h.region(location, edgeL, edgeR)

        dims = list(pf.domain_dimensions)
        ngrid_x = dims[0]
        ngrid_y = dims[1]
        proj = pf.h.slice(2, location[2], center=location)
        w = (pf.h.domain_left_edge[0], pf.h.domain_right_edge[0], pf.h.domain_left_edge[1], pf.h.domain_right_edge[1])
        frb1 = FixedResolutionBuffer(proj, w, (ngrid_x, ngrid_x), periodic=False)

	p=SlicePlot(pf,'x', 'Temp', center=location)
	p.set_width(length,'pc')
        p.set_buff_size(256)
	p.set_zlim('Temp', tmin, tmax)
	p.set_cmap('Temp', 'Rainbow18')
        if read_sink == 1:
		for i in range(0,num_sink):
			center = [xpos_sink[i], ypos_sink[i], zpos_sink[i]]
			print('sink_size = ', sink_size, 'sink loc = ', center)
			p.annotate_sphere(center, sink_size, {'fill':True})
	p.save('slice_t'+datanum)

        p=SlicePlot(pf,'x', 'density', center=location)
        p.set_width(length,'pc')
        p.set_zlim('density', dmin, dmax)
	#p.annotate_grids()
        if read_sink == 1:
                for i in range(0,num_sink):
                        center = [xpos_sink[i], ypos_sink[i], zpos_sink[i]]
                        p.annotate_sphere(center, sink_size, {'fill':True})
	p.save('slice_t'+datanum)

        p=SlicePlot(pf,'x', 'Bmag-log', center=location)
        p.set_width(length,'pc')
        p.set_zlim('Bmag-log', bmin, bmax)
	p.set_cmap('Bmag-log', 'Rainbow18')
        #p.annotate_grids()
        if read_sink == 1:
                for i in range(0,num_sink):
                        center = [xpos_sink[i], ypos_sink[i], zpos_sink[i]]
                        p.annotate_sphere(center, sink_size, {'fill':True})
        p.save('slice_t'+datanum)

        p=SlicePlot(pf,'x', 'Toomre', center=location)
        p.set_width(length,'pc')
        p.set_zlim('Toomre', qmin, qmax)
        #p.set_cmap('Toomre', 'Rainbow18')
        #p.annotate_grids()
        if read_sink == 1:
                for i in range(0,num_sink):
                        center = [xpos_sink[i], ypos_sink[i], zpos_sink[i]]
                        p.annotate_sphere(center, sink_size, {'fill':True})
        p.save('slice_t'+datanum)

	'''
	print ('len density =', len(data['density']))
	print ('min density=', min(data['density']), 'max density=', max(data['density']))
	print ('min Temp=', min(data['Temp']), 'max Temp=', max(data['Temp']))
	print ('min Q=', min(data['Toomre']), 'max Q=', max(data['Toomre']))
	'''
