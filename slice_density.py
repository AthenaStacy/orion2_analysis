import matplotlib
matplotlib.use('pdf')

from yt.mods import *
import numpy as na
from string import rstrip
import fields_iso
import fields_wpot

#datanum = '0021'
datanum = '0000'

#pf = load("/work/00863/minerva/orion/bfield_ref2_maptest1_2lev_nosink/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion/bfield_comp1/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion/bfield_maptest2/" + 'data.' + datanum + '.3d.hdf5')
pf = load("/work/00863/minerva/orion/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion/bfield_ref3_maptest1/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion/gravpot_2pc_256Nest_HR10/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion_Gtest/" + 'data.' + datanum + '.3d.hdf5')
print pf.h.field_list

value, location = pf.h.find_max("Density")
data = pf.h.sphere(location, 1.e0/pf['pc'])
print 'max density location = ', location

#SlicePlot(pf, 'x', "z-velocity", width = (1, 'pc')).save()

#print 'min tracer1=', min(data['tracer1']), 'max tracer1=', max(data['tracer1'])

pc=PlotCollection(pf, center=[0,0,0])
#pc=PlotCollection(pf, center=location)

p1=SlicePlot(pf, 0, "Density")
p2=SlicePlot(pf, 1, "Density")
p3=SlicePlot(pf, 2, "Density")

#p.set_zlim(1.e-23, 1.e-20)
p1.set_cmap('all','Paired')
p2.set_cmap('all','Paired')
p3.set_cmap('all','Paired')

print_extremes = 0

if print_extremes == 1:
	print 'min density=', min(data['density']), 'max density=', max(data['density'])
	print 'min gravpot=', min(data['gravitational-potential']), 'max gravpot=', max(data['gravitational-potential'])
	print 'vx_min=', min(data['x-velocity']), 'vx_max=', max(data['x-velocity'])
	print 'vy_min=', min(data['y-velocity']), 'vy_max=', max(data['y-velocity'])
	print 'vz_min=', min(data['z-velocity']), 'vz_max=', max(data['z-velocity'])
	print 'xmom_min=', min(data['X-momentum']), 'zmom_max=', max(data['X-momentum'])
	print 'ymom_min=', min(data['Y-momentum']), 'ymom_max=', max(data['Y-momentum'])
	print 'zmom_min=', min(data['Z-momentum']), 'xmom_max=', max(data['Z-momentum'])


p1.save('density_x'+ datanum)
p2.save('density_y'+ datanum)
p3.save('density_z'+ datanum)
