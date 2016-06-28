import matplotlib
matplotlib.use('pdf')

from yt.mods import *
import numpy as na
from string import rstrip
import fields_iso
import fields_bfield
import fields_wpot

datanum = '0000'

#pf = load("/work/00863/minerva/orion/bfield_comp1/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion/bfield_maptest2/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion/" + 'data.' + datanum + '.3d.hdf5')
pf = load("/work/00863/minerva/orion_Btest/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion_Btest/bfield_homolog/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion_Btest/bfield_new/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion_Btest/bfield_analyt/" + 'data.' + datanum + '.3d.hdf5')

print pf.h.field_list

value, location = pf.h.find_max("Density")
data = pf.h.sphere(location, 1.e3/pf['pc'])
print 'max density location = ', location

#SlicePlot(pf, 'x', "z-velocity", width = (1, 'pc')).save()

#print 'min tracer1=', min(data['tracer1']), 'max tracer1=', max(data['tracer1'])

pc=PlotCollection(pf, center=[0,0,0])
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

def _DivB_log(field,data):
    return np.log10((data["absDivB"]))
add_field("DivB-log",function=_DivB_log, take_log=False,
          units=r'Gauss / cm')

def _absBx(field,data):
    return (np.abs(data["X-magnfield"]))
add_field("absBx",function=_Bx_log, take_log=False,
          units=r'Gauss')

def _absBy(field,data):
    return (np.abs(data["Y-magnfield"]))
add_field("absBy",function=_By_log, take_log=False,
          units=r'Gauss')

def _absBz(field,data):
    return (np.abs(data["Z-magnfield"]))
add_field("absBz",function=_Bz_log, take_log=False,
          units=r'Gauss')
######################################################################################

#bmin = -2.e-9
#bmax =  2.e-9
bmin =  1.e-8
bmax =  1.e-5

print 'min Bfield-X=', min(data['X-magnfield']), 'max Bfield-X=', max(data['X-magnfield'])
print 'min Bfield-Y=', min(data['Y-magnfield']), 'max Bfield-Y=', max(data['Y-magnfield'])
print 'min Bfield-Z=', min(data['Z-magnfield']), 'max Bfield-Z=', max(data['Z-magnfield'])
print 'min log-Bfield-X=', min(data['Bx-log']), 'max log-Bfield-X=', max(data['Bx-log'])
print 'min log-Bfield-Y=', min(data['By-log']), 'max log-Bfield-Y=', max(data['By-log'])
print 'min log-Bfield-Z=', min(data['Bz-log']), 'max log-Bfield-Z=', max(data['Bz-log'])

p=pc.add_slice('absBx', 'z')
#p.set_zlim(bmin, bmax)
p=pc.add_slice('absBy', 'z')
#p.set_zlim(bmin, bmax)
p=pc.add_slice('absBz', 'z')
#p.set_zlim(bmin, bmax)

p=pc.add_slice('absBx', 'x')
#p.set_zlim(bmin, bmax)
p=pc.add_slice('absBy', 'x')
#p.set_zlim(bmin, bmax)
p=pc.add_slice('absBz', 'x')
#p.set_zlim(bmin, bmax)

p=pc.add_slice('absBx', 'y')
#p.set_zlim(bmin, bmax)
p=pc.add_slice('absBy', 'y')
#p.set_zlim(bmin, bmax)
p=pc.add_slice('absBz', 'y')
#p.set_zlim(bmin, bmax)


emin = 1e-18
emax = 1e-16
p=pc.add_slice('MagneticEnergy', 'x')
p.set_zlim(emin, emax)

p=pc.add_slice('VorticityX', 'x')
#p=pc.add_slice('VorticityMagnitude', 'x')

p=pc.add_slice('DivB', 'x')
p=pc.add_slice('DivB', 'y')
p=pc.add_slice('DivB', 'z')
#p=pc.add_slice('DivB-log', 'x')

dmin = 1.e-27
dmax = 1.e-26
p=pc.add_slice('density', 'x')
#p.set_zlim(dmin, dmax)
p=pc.add_slice('density', 'y')
#p.set_zlim(dmin, dmax)
p=pc.add_slice('density', 'z')
#p.set_zlim(dmin, dmax)

print 'total divB = ', data['DivB'].sum()
print 'min density=', min(data['density']), 'max density=', max(data['density'])
print 'min gravpot=', min(data['gravitational-potential']), 'max gravpot=', max(data['gravitational-potential'])

print 'B average = ', np.mean(data['Bmag'])
print 'dens average =', np.mean(data['density']) 

pc.save('data_'+ datanum)


