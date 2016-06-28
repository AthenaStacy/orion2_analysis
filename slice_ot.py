import matplotlib
matplotlib.use('pdf')

from yt.mods import *
import numpy as na
from string import rstrip
import fields
import fields_iso
import fields_bfield
import fields_wpot

datanum = '0020'

#pf = load("/work/00863/minerva/orion_Otest/lowBlowRes/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion_Otest/lowB/" + 'data.' + datanum + '.3d.hdf5')
pf = load("/work/00863/minerva/orion_Otest/lowBhighRes/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion_Otest/" + 'data.' + datanum + '.3d.hdf5')

print pf.h.field_list

value, location = pf.h.find_max("Density")
data = pf.h.sphere(location, 1.e3/pf['pc'])
print 'max density location = ', location

#SlicePlot(pf, 'x', "z-velocity", width = (1, 'pc')).save()

#print 'min tracer1=', min(data['tracer1']), 'max tracer1=', max(data['tracer1'])

pc=PlotCollection(pf, center=[0.0,0.0,0.0])
#pc=PlotCollection(pf, center=location)


########################################################get B-field logs###############
def _Bx_log(field,data):
    return np.log10(np.abs(data["X-magnfield"]))
add_field("Bx-log",function=_Bx_log, take_log=False,
          units=r'Gauss')

def _By_log(field,data):
    return np.log10(np.abs(data["Y-magnfield"]))
add_field("By-log",function=_By_log, take_log=False,
          units=r'Gauss')

def _Bz_log(field,data):
    return np.log10(np.abs(data["Z-magnfield"]))
add_field("Bz-log",function=_Bz_log, take_log=False,
          units=r'Gauss')

def _Bmag_log(field,data):
    return np.log10(np.abs(data["Bmag"]))
add_field("Bmag-log",function=_Bmag_log, take_log=False,
          units=r'Gauss')

def _DivB_log(field,data):
    return np.log10((data["absDivB"]))
add_field("DivB-log",function=_DivB_log, take_log=False,
          units=r'Gauss / cm')
######################################################################################

print 'min Bfield-X=', min(data['X-magnfield']), 'max Bfield-X=', max(data['X-magnfield'])
print 'min Bfield-Y=', min(data['Y-magnfield']), 'max Bfield-Y=', max(data['Y-magnfield'])
print 'min Bfield-Z=', min(data['Z-magnfield']), 'max Bfield-Z=', max(data['Z-magnfield'])
print 'min log-Bfield-X=', min(data['Bx-log']), 'max log-Bfield-X=', max(data['Bx-log'])
print 'min log-Bfield-Y=', min(data['By-log']), 'max log-Bfield-Y=', max(data['By-log'])
print 'min log-Bfield-Z=', min(data['Bz-log']), 'max log-Bfield-Z=', max(data['Bz-log'])
print 'min vlog=', min(data['vlog']), 'max vlog=', max(data['vlog'])

#bmin =  min(data['Z-magnfield'])
#bmax = max(data['Z-magnfield'])
#bmin =  -1.e-18
#bmax =  1.e-18
bmin = -1  
bmax = 1

#p=pc.add_slice('X-magnfield', 'z')
#p.set_zlim(bmin, bmax)
#p=pc.add_slice('Y-magnfield', 'z')
#p.set_zlim(bmin, bmax)
#p=pc.add_slice('Z-magnfield', 'z')
#p.set_zlim(bmin, bmax)

bmin = -21
bmax = -18
p=pc.add_slice('Bx-log', 'z')
p.set_zlim(bmin, bmax)
p.set_cmap( "Rainbow18")

bmin = -40
p=pc.add_slice('By-log', 'z')
p.set_zlim(bmin, bmax)
p.set_cmap( "Rainbow18")

p=pc.add_slice('Bz-log', 'x')
p.set_zlim(bmin, bmax)
p.set_cmap( "Rainbow18")

bmin = -21
bmax = -18
#bmin = -11 
#bmax = -8
p=pc.add_slice('Bmag-log', 'z')
p.set_zlim(bmin, bmax)
p.set_cmap( "Rainbow18")

#emin = 1e-18
#emax = 1e-16
emin = 2.1e-7 
emax = 0.3

dmin = 1.e-1
dmax = 5.e-1
#dmin = 1.e4
#dmax = 1.e5
p=pc.add_slice('density', 'z')
p.set_zlim(dmin, dmax)

vmin = 0.1
vmax = 1.1
p=pc.add_slice('x-velocity', 'z')
#p.set_zlim(vmin, vmax)

p=pc.add_slice('y-velocity', 'z')
#p.set_zlim(vmin, vmax)

p=pc.add_slice('z-velocity', 'x')
p.set_zlim(vmin, vmax)

p=pc.add_slice('vlog', 'z')
#p.set_zlim(vmin, vmax)

p=pc.add_slice('vlog', 'x')

print 'total divB = ', data['DivB'].sum()

ntot = np.mean(data['density']*data['density'])
nrms = np.sqrt(ntot)
print 'nmin=', min(data['density']), 'nmax=', max(data['density']), 'nrms = ', nrms


Btot = np.mean(data['Bmag']*data['Bmag'])
Brms = np.sqrt(Btot)
print 'Bmin =', min(data['Bmag']), 'Bmax = ', max(data['Bmag']), 'Brms = ', Brms


print 'B average = ', np.mean(data['Bmag'])
print 'dens average =', np.mean(data['density']) 

pc.save('data_'+ datanum)


