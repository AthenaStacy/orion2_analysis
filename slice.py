import matplotlib
matplotlib.use('pdf')

from yt.mods import *
import numpy as na
from string import rstrip
#import fields_iso
import fields
import fields_wpot

#datanum = '0012'
#datanum = '0010'
datanum = '0023'

#pf = load("/work/00863/minerva/orion/bfield_ref2_maptest1_2lev_nosink/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion/bfield_comp1/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion/bfield_maptest2/" + 'data.' + datanum + '.3d.hdf5')
pf = load("/work/00863/minerva/orion/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion_Gtest/" + 'data.' + datanum + '.3d.hdf5')
print pf.h.field_list

value, location = pf.h.find_max("Density")
data = pf.h.sphere(location, 1.e0/pf['pc'])
print 'max density location = ', location

#SlicePlot(pf, 'x', "z-velocity", width = (1, 'pc')).save()

#print 'min tracer1=', min(data['tracer1']), 'max tracer1=', max(data['tracer1'])

#pc=PlotCollection(pf, center=[0,0,0])
pc=PlotCollection(pf, center=location)

p=pc.add_slice('density', 'x')
p.set_zlim(1.e-21, 1.e-14)
p=pc.add_slice('density', 'y')
p.set_zlim(1.e-21, 1.e-14)
p=pc.add_slice('density', 'z')
p.set_zlim(1.e-21, 1.e-14)


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

vmin = -2.e5
vmax = 2.e5

print_vrad = 0
print_vel = 0
print_cs = 0
print_mach = 0

if print_vrad == 1:
	print 'min vrad=', min(data['Vrad']), 'max vrad=', max(data['Vrad'])
	p=pc.add_slice('Vrad', 'x')
	p.set_zlim(vmin, vmax)
	p=pc.add_slice('Vrad', 'y')
	p.set_zlim(vmin, vmax)
	p=pc.add_slice('Vrad', 'z')
	p.set_zlim(vmin, vmax)

vmin = 0
vmax = 2.e5
if print_vel == 1:
	print 'min velocity=', min(data['velocity']), 'max velocity=', max(data['velocity'])
	p=pc.add_slice('velocity', 'x')
	p.set_zlim(vmin, vmax)
	p=pc.add_slice('velocity', 'y')
	p.set_zlim(vmin, vmax)
	p=pc.add_slice('velocity', 'z')
	p.set_zlim(vmin, vmax)

vmin = 0
vmax = 2.e5
if print_cs == 1:
	print 'min cs', min(data['cs']), 'max cs=', max(data['cs'])
	p=pc.add_slice('cs', 'x')
	p.set_zlim(vmin, vmax)
	p=pc.add_slice('cs', 'y')
	p.set_zlim(vmin, vmax)
	p=pc.add_slice('cs', 'z')
	p.set_zlim(vmin, vmax)

vmin = 0
vmax = 2
if print_mach == 1:
	print 'min mach', min(data['mach']), 'max mach=', max(data['mach'])
	p=pc.add_slice('mach', 'x')
	p.set_zlim(vmin, vmax)
	p=pc.add_slice('mach', 'y')
	p.set_zlim(vmin, vmax)
	p=pc.add_slice('mach', 'z')
	p.set_zlim(vmin, vmax)

tmin = 1.e1
tmax = 1.e4
print 'min temp', min(data['Temperature']), 'max temp=', max(data['Temperature'])
print 'min pres', min(data['Pressure']), 'max pres=', max(data['Pressure'])
p=pc.add_slice('Temperature', 'x')
p.set_zlim(tmin, tmax)
p=pc.add_slice('Temperature', 'y')
p.set_zlim(tmin, tmax)
p=pc.add_slice('Temperature', 'z')
p.set_zlim(tmin, tmax)
#ProjectionPlot(pf, 'x', 'Temperature', width=(1.3,'pc')).save()

p=pc.add_slice('energy-density', 'x')
p=pc.add_slice('x-velocity', 'x')

pc.save('data_'+ datanum)
