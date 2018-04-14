import matplotlib
matplotlib.use('pdf')

from yt.mods import *
import numpy as na
from string import rstrip
#import fields_iso
import fields
import fields_bfield
import fields_wpot
from tracer_def import *

datanum = '0409'

pf = load("/nobackupp7/astacy/popiii_Bscope17/" + 'data.' + datanum + '.3d.hdf5')

print  pf.h.field_list
print 'dimensions =', pf.dimensionality

#value, location = pf.h.find_max("Density")
value, location = pf.h.find_max("Bmag")
data = pf.h.sphere(location, 0.1/pf['pc'])
print 'max density location = ', location


#pc=PlotCollection(pf, center=[0,0,0])
pc=PlotCollection(pf, center=location)

dims = list(pf.domain_dimensions)
ngrid_x = dims[0]
ngrid_y = dims[1]
#proj = pf.h.slice(2, 0.0, center=pf.domain_center)
proj = pf.h.slice(2, location[2], center=location)
w = (pf.h.domain_left_edge[0], pf.h.domain_right_edge[0], pf.h.domain_left_edge[1], pf.h.domain_right_edge[1])
frb1 = FixedResolutionBuffer(proj, w, (ngrid_x, ngrid_x), periodic=False)
bx  = frb1['X-magnfield']
by  = frb1['Y-magnfield']
bz  = frb1['Z-magnfield']
bx_face  = frb1['Bx_face']
by_face  = frb1['By_face']
bz_face  = frb1['Bz_face']
rho = frb1['Density']
velx = frb1['x-velocity']
x1  = frb1['x']
y1  = frb1['y']
z1  = frb1['z']

########################################################get B-field logs###############
def _DivB_log(field,data):
    return np.log10((data["absDivB"]))
add_field("DivB-log",function=_DivB_log, take_log=False,
          units=r'Gauss / cm')


dims = list(pf.domain_dimensions)
DeltaX = pf.domain_right_edge[0] / dims[0]
print 'DeltaX =', DeltaX
def _errDivB(field,data):
    return np.log10(data["absDivB"] * DeltaX / data["Bmag"])
add_field("errDivB",function=_errDivB, take_log=False,
          units=r' ')

######################################################################################

print 'min density=', min(data['density']), 'max density=', max(data['density'])
print 'min Bfield-X=', min(data['X-magnfield']), 'max Bfield-X=', max(data['X-magnfield'])
print 'min Bfield-Y=', min(data['Y-magnfield']), 'max Bfield-Y=', max(data['Y-magnfield'])
print 'min Bfield-Z=', min(data['Z-magnfield']), 'max Bfield-Z=', max(data['Z-magnfield'])
print 'min log-Bfield-X=', min(data['Bx-log']), 'max log-Bfield-X=', max(data['Bx-log'])
print 'min log-Bfield-Y=', min(data['By-log']), 'max log-Bfield-Y=', max(data['By-log'])
print 'min log-Bfield-Z=', min(data['Bz-log']), 'max log-Bfield-Z=', max(data['Bz-log'])

bmin = -1.e-14
bmax =  1.e-14

p=pc.add_slice('X-magnfield', 'x')
p.set_zlim(bmin, bmax)
p=pc.add_slice('Y-magnfield', 'x')
p.set_zlim(bmin, bmax)
p=pc.add_slice('Z-magnfield', 'x')
p.set_zlim(bmin, bmax)

p=pc.add_slice('X-magnfield', 'y')
p.set_zlim(bmin, bmax)
p=pc.add_slice('Y-magnfield', 'y')
p.set_zlim(bmin, bmax)
p=pc.add_slice('Z-magnfield', 'y')
p.set_zlim(bmin, bmax)

p=pc.add_slice('X-magnfield', 'z')
p.set_zlim(bmin, bmax)
p=pc.add_slice('Y-magnfield', 'z')
p.set_zlim(bmin, bmax)
p=pc.add_slice('Z-magnfield', 'z')
p.set_zlim(bmin, bmax)


#bmin = -15
#bmax = -11
bmin = -18
bmax = -8
p=pc.add_slice('Bx-log', 'z')
p.set_zlim(bmin, bmax)
p.set_cmap("Rainbow18")

p=pc.add_slice('By-log', 'z')
p.set_zlim(bmin, bmax)
p.set_cmap("Rainbow18")

p=pc.add_slice('Bz-log', 'z')
p.set_zlim(bmin, bmax)
p.set_cmap("Rainbow18")


p=pc.add_slice('Bmag-log', 'x')
p.set_zlim(bmin, bmax)
p.set_cmap("Rainbow18")

p=pc.add_slice('Bmag-log', 'y')
p.set_zlim(bmin, bmax)
p.set_cmap("Rainbow18")

p=pc.add_slice('Bmag-log', 'z')
p.set_zlim(bmin, bmax)
p.set_cmap("Rainbow18")


dmin = 1.e-25 
dmax = 1.e-10
p=pc.add_slice('Density', 'x')
p.set_zlim(dmin, dmax)
p.set_cmap("Rainbow18")
p=pc.add_slice('Density', 'y')
p.set_zlim(dmin, dmax)
p.set_cmap("Rainbow18")
p=pc.add_slice('Density', 'z')
p.set_zlim(dmin, dmax)
p.set_cmap("Rainbow18")

tmin = 100.
tmax = 1.e4
p=pc.add_slice('Temperature', 'x')
p.set_zlim(tmin, tmax)
p.set_cmap("Rainbow18")
p=pc.add_slice('Temperature', 'y')
p.set_zlim(tmin, tmax)
p.set_cmap("Rainbow18")
p=pc.add_slice('Temperature', 'z')
p.set_zlim(tmin, tmax)
p.set_cmap("Rainbow18")

dmin = 1.e-45
dmax = 1.e-37
p=pc.add_slice('DivB', 'x')
p.set_zlim(dmin, dmax)
p.set_cmap("Rainbow18")
p=pc.add_slice('DivB', 'y')
p.set_zlim(dmin, dmax)
p.set_cmap("Rainbow18")
p=pc.add_slice('DivB', 'z')
p.set_zlim(dmin, dmax)
p.set_cmap("Rainbow18")

gmin = 1e-25
gmax = 1e-24
p=pc.add_slice('tracer1', 'x')
p.set_zlim(gmin, gmax)
p=pc.add_slice('tracer1', 'y')
p.set_zlim(gmin, gmax)
p=pc.add_slice('tracer1', 'z')
p.set_zlim(gmin, gmax)

p=pc.add_slice('tracer6', 'x')
p.set_zlim(gmin, gmax)
p=pc.add_slice('tracer6', 'y')
p.set_zlim(gmin, gmax)
p=pc.add_slice('tracer6', 'z')
p.set_zlim(gmin, gmax)


print 'total divB = ', data['DivB'].sum()
print 'max divB = ', max(data['DivB'])
print 'min density=', min(data['density']), 'max density=', max(data['density'])
print 'min Temp=', min(data['Temperature']), 'max Temp=', max(data['Temperature'])
print 'min tracer1=', min(data['tracer1']), 'max tracer1=', max(data['tracer1'])
print 'min tracer6=', min(data['tracer6']), 'max tracer6=', max(data['tracer6'])
print 'min Bx=', min(data['X-magnfield']), 'max Bx=', max(data['X-magnfield'])
print 'min By=', min(data['Y-magnfield']), 'max By=', max(data['Y-magnfield'])
print 'min Bz=', min(data['Z-magnfield']), 'max Bz=', max(data['Z-magnfield'])

print 'B average = ', np.mean(data['Bmag'])
print 'dens average =', np.mean(data['density']) 

pc.save('data_'+ datanum)


