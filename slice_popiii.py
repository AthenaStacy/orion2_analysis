import matplotlib
matplotlib.use('pdf')

from yt.mods import *
import numpy as na
from string import rstrip
#import fields_iso
#import fields
#import fields_bfield
#import fields_wpot
#from tracer_def import *

#rdir = "/work/00863/minerva/popiii_Bscope6/"
rdir = "/work/00863/minerva/popiii_Bscope10/"

#jstart = 209
jstart = 187

dmin = 1.e-18
dmax = 1.e-13

gmin = 1.4
gmax = 1.67

tmin = 200.
tmax = 5000.

bmin = -14
bmax = -10
########################################################get B-field logs###############
#######################################################################
##Read in sink information
read_sink = 1
res = 16
nlev = 10
pc_to_cm = 3.08567758e18
######################################################################################

######################################################################################


length = 0.01 #length of box we want in pc

for x in xrange(1,2):

        j = jstart
#	j = 10*x + 90
	jstr = str(j)

	if j < 100:
		datanum = '00' + jstr
        if j >= 100:
                datanum = '0' + jstr

	pf = load(rdir + 'data.' + datanum + '.3d.hdf5')

        value, location = pf.h.find_max("Density")
        #value, location = pf.h.find_max("Temp.")

        print 'LINE 57'
	
	rad  = 0.5*(pf.domain_right_edge[0]-pf.domain_left_edge[0])*pf['pc']
	sink_size =  pf.domain_right_edge[0] / res / (2**(nlev-1))
	sink_size = sink_size * 4.5 / pc_to_cm
	print 'sink_size =', sink_size

	xpos_sink = []
	ypos_sink = []
	zpos_sink = []

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

		center1 = [xpos_sink[0], ypos_sink[0], zpos_sink[0]]
                location = center1

        print 'LINE 82'

	print 'max density location = ', location

	#data = pf.h.sphere(location, length/pf['pc'])
        #edgeL = [ -length/pf['pc'], -length/pf['pc'], -length/pf['pc']]
        #edgeR = [  length/pf['pc'],  length/pf['pc'],  length/pf['pc']]
        #data = pf.h.region(location, edgeL, edgeR)

        #ad = pf.h.all_data()
        #cr = ad.cut_region(["grid['density'] > 1e-15"])
        #data.clear_cache()
        #cube = pf.h.covering_grid(0, left_edge=edgeL, dims=[128, 128, 128], fields = ['Temp.', 'density'])

	#pc=PlotCollection(pf, center=[0,0,0])
	#pc=PlotCollection(pf, center=location)

        #dims = list(pf.domain_dimensions)
        #ngrid_x = dims[0]
        #ngrid_y = dims[1]
        #proj = pf.h.slice(2, location[2], center=location)
        #w = (pf.h.domain_left_edge[0], pf.h.domain_right_edge[0], pf.h.domain_left_edge[1], pf.h.domain_right_edge[1])
        #frb1 = FixedResolutionBuffer(proj, w, (ngrid_x, ngrid_x), periodic=False)

#	p=SlicePlot(pf,'x','Temp.', center=location)
#	p.set_width(length,'pc')
        #p.set_buff_size(256)
#	p.set_zlim('Temp.', tmin, tmax)
#	p.set_cmap('Temp.', 'Rainbow18')
#        if read_sink == 1:
#		p.annotate_sphere(center1, sink_size/pf['pc'], {'fill':True})
#	p.save('data_'+datanum)

        p=SlicePlot(pf,'z', 'density', center=location)
        #p.set_width(length,'pc')
        p.zoom(10)
        p.set_zlim('density', dmin, dmax)
	#p.annotate_grids()
        if read_sink == 1:
                p.annotate_sphere(center1, sink_size/pf['pc'], {'fill':True})
	p.save('data_'+datanum)

#	p=pc.add_slice('Temp_alt', 'z')
#	p.set_width(length,'pc')
#	p.set_zlim(tmin, tmax)
#	p.set_cmap("Rainbow18")

#	p=pc.add_slice(iH2, 'x')
#	p.set_width(length,'pc')
#	p=pc.add_slice(iH2, 'y')
#	p.set_width(length,'pc')
#	p=pc.add_slice(iH2, 'z')
#	p.set_width(length,'pc')

#	p=pc.add_slice('density', 'z')
#	p.set_width(length,'pc')
#	pc.plots[-1].modify["grids"]()
#	p.set_zlim(dmin, dmax)

	#print 'min density=', min(data['density']), 'max density=', max(data['density'])
	#print 'min tracer1=', min(data['tracer1']), 'max tracer1=', max(data['tracer1'])
	#print 'min tracer2=', min(data['tracer2']), 'max tracer2=', max(data['tracer2'])
	#print 'min tracer3=', min(data['tracer3']), 'max tracer3=', max(data['tracer3'])
	#print 'min tracer5=', min(data['tracer5']), 'max tracer5=', max(data['tracer5'])
	#print 'min tracer6=', min(data['tracer6']), 'max tracer6=', max(data['tracer6'])
	#print 'min tracer15=', min(data['tracer15']), 'max tracer15=', max(data['tracer15'])
	#print 'min tracer17=', min(data['tracer17']), 'max tracer17=', max(data['tracer17'])
	#print 'min gravpot=', min(data['gravitational-potential']), 'max gravpot=', max(data['gravitational-potential'])
	#print 'min Temp=', min(data['Temp.']), 'max Temp=', max(data['Temp.'])

	#print 'B average = ', np.mean(data['Bmag'])
	#print 'dens average =', np.mean(data['density']) 

	#pc.save('data_'+ datanum)


