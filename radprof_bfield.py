import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as pl
#import pylab as pl

from yt.mods import *
import numpy as na
from string import rstrip
import fields
import fields_bfield

au_cm = 1.49597871e13
Msun = 1.98892e33

#datanum = '0005'
datanum = '0000'

#pf = load("/work/00863/minerva/orion/bfield_ref2_maptest1_2lev_nosink/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion/bfield_comp1/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion/bfield_maptest2/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion/bfield_ref3_maptest1/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion/bfield_ref3_maptest1_2lev/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion_Btest/bfield_orig/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion_Btest/bfield_homolog/" + 'data.' + datanum + '.3d.hdf5')
pf = load("/work/00863/minerva/orion_Btest/" + 'data.' + datanum + '.3d.hdf5')

print pf.h.field_list

value, location = pf.h.find_max("Density")
data = pf.h.sphere(location, 1.e0/pf['pc'])
print 'max density location = ', location

pc = PlotCollection(pf, center = location)

print pf.h.field_list

bin_num_square = 32
bin_doub = 32.0
DeltaX = pf.domain_right_edge[0] / bin_doub

x=location[0]
y=location[1]
z=location[2]

def _rVelocity(field, data):
    '''
    The infall velocity. In this problem the center is at density peak.
    '''
    vr = data['x-velocity']*data['x'] + data['y-velocity']*data['y'] + data['z-velocity']*data['z']
    vr = (vr) / na.sqrt(data['x']**2 + data['y']**2 + data['z']**2)
    vr = vr/1.e5
    return vr
add_field("radial-velocity", function=_rVelocity, take_log=False,
          units=r'\rm{km}/\rm{s}')

def _Radius(field, data):
    '''
    Distance from central point. In this problem the center is at density peak.
    '''
    rad =  na.sqrt((data['x'] )**2 + (data['y'])**2 + (data['z'])**2)
    rad = rad/1.5e13
    return rad
add_field("Radius", function=_Radius, take_log=False,
          units=r'\rm{AU}')

print 'temp_min=', min(data['Temperature']), 'temp_max=', max(data['Temperature'])
print 'vrad_min=', min(data['radial-velocity']), 'vrad_max=', max(data['radial-velocity'])
print 'dens_min=', min(data['Density']), 'dens_max=', max(data['Density'])

p1 = pc.add_profile_sphere(2., "pc", ['Radius', 'number-density'], weight='CellVolume')
p2 = pc.add_profile_sphere(2., "pc", ['Radius', 'radial-velocity'], weight='CellVolume')
p3 = pc.add_profile_sphere(2., "pc", ['density', 'Temperature', 'Bmag', 'ufield', 'MagneticEnergy', 'DivB', 'absDivB'], weight='CellVolume')

rad = p1.data['Radius']
nh = p1.data['number-density'] 
rad2 = p2.data['Radius']
vrad = p2.data['radial-velocity'] 

dens_arr = p3.data['density']
temp = p3.data['Temperature']
bfield = p3.data['Bmag']
ufield = p3.data['MagneticEnergy']
ufield_norm = ufield / (dens_arr**1.3333333)
DivB = p3.data['DivB']
absDivB = p3.data['absDivB']

errDivB = absDivB * DeltaX / bfield

print 'rad =', rad
print 'nh =', nh
print 'rad2 =', rad2
print 'vrad = ', vrad
print 'bfield =', bfield
print 'ufield_norm =', ufield_norm
print 'error divB =', errDivB
print 'total divB = ', data['DivB'].sum()

dens_min = 1.e-24
dens_max = 1.e-14

bmin = 1.e-15
bmax = 1.e-5

umin = 1.e5
umax = 1.e15

DivBmin = -1.e-22
DivBmax = 1.e-22


pl.subplot(221)
pl.plot(dens_arr, bfield, 'k')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel('B-field [G]', fontsize=9)
ax.set_xlabel('density [g cm^{-3}]', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((dens_min, dens_max, bmin, bmax))

pl.subplot(222)
pl.plot(dens_arr, ufield_norm,'k')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('density [g cm^{-3}]', fontsize=9)
ax.set_ylabel('U / rho^(4/3) [cgs]', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((dens_min, dens_max, umin, umax))

pl.subplot(223)
pl.plot(dens_arr, errDivB,'k')
ax = pl.gca()
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('density [g cm^{-3}]', fontsize=9)
ax.set_ylabel('error DivB', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((dens_min, dens_max, -1.e0, 1.e0))

pl.subplot(224)
pl.plot(dens_arr, absDivB,'k')
ax = pl.gca()
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('density [g cm^{-3}]', fontsize=9)
ax.set_ylabel('|DivB| [Gauss / cm]', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((dens_min, dens_max, -DivBmax, DivBmax))
#pl.axis((dens_min, dens_max, 0, DivBmax))

pl.savefig('bfield.eps')

M_i = data["CellMassMsun"].sum()
print "total mass = ", M_i
