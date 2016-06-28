import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as pl
#import pylab as pl

from yt.mods import *
import numpy as na
from string import rstrip
import fields_iso
import fields_wpot

au_cm = 1.49597871e13

wdir = '/work/00863/minerva/orion/gravpot_0.8pc/'
datanum = '0000'
setup = 3

pf = load(wdir + "data." + datanum + ".3d.hdf5")
value, location = pf.h.find_max("Density")
data = pf.h.sphere(location, 2.0/pf['pc'])
pc = PlotCollection(pf, center = location)

print pf.h.field_list

x=location[0]
y=location[1]
z=location[2]

massenc_com = []
vx = []
vy = []
vz = []
sub_sphere =  pf.h.sphere(location, 6.e16/pf['cm'])
massenc_com.append(sub_sphere["CellMassMsun"].sum())
#vx.append((sub_sphere["x-velocity"]*sub_sphere["CellMassMsun"]).sum()/(sub_sphere["CellMassMsun"].sum()))
#vy.append((sub_sphere["y-velocity"]*sub_sphere["CellMassMsun"]).sum()/(sub_sphere["CellMassMsun"].sum()))
#vz.append((sub_sphere["z-velocity"]*sub_sphere["CellMassMsun"]).sum()/(sub_sphere["CellMassMsun"].sum()))

vx.append(0)
vy.append(0)
vz.append(0)

print 'location =', location, 'x=', x, 'y=', y, 'z=', z, 'vx =', vx, 'vy =', vy, 'vz =', vz

def _rVelocity(field, data):
    '''
    The infall velocity. In this problem the center is at density peak.
    '''
    vr = (data['x-velocity'] - vx[0])*(data['x'] - x) + (data['y-velocity'] - vy[0])*(data['y'] - y) + (data['z-velocity'] - vz[0])*(data['z'] - z)
    vr = (vr) / na.sqrt((data['x'] - x)**2 + (data['y'] - y)**2 + (data['z'] - z)**2)
    vr = vr/1.e5
    return vr
add_field("radial-velocity", function=_rVelocity, take_log=False,
          units=r'\rm{km}/\rm{s}')

def _Radius(field, data):
    '''
    Distance from central point. In this problem the center is at density peak.
    '''
    rad =  na.sqrt((data['x'] - x)**2 + (data['y'] - y)**2 + (data['z'] - z)**2)
    rad = rad/1.5e13
    return rad
add_field("Radius", function=_Radius, take_log=False,
          units=r'\rm{AU}')

def _rotVelocity(field, data):
    '''
    Rotational velocity. In this problem the center is at density peak.
    '''
    vrotx = (data['z-velocity'] - vz[0])*(data['y'] - y) - (data['y-velocity'] - vy[0])*(data['z'] - z)
    vroty = (data['x-velocity'] - vx[0])*(data['z'] - z) - (data['z-velocity'] - vz[0])*(data['x'] - x)
    vrotz = (data['y-velocity'] - vy[0])*(data['x'] - x) - (data['x-velocity'] - vx[0])*(data['y'] - y)
    vrot = na.sqrt(vrotx**2 + vroty**2 + vrotz**2)
    vrot = (vrot) / na.sqrt((data['x'] - x)**2 + (data['y'] - y)**2 + (data['z'] - z)**2)
    vrot = vrot/1.e5
    return vrot
add_field("rot-velocity", function=_rotVelocity, take_log=False,
          units=r'\rm{km}/\rm{s}')

##############################################################################

rmin = 100.
rmax = 3.e5
nmin = 1e1
nmax = 1e16
bin_num1 = 1000

p1 = pc.add_profile_sphere(1., "pc", ['Radius', 'density', 'number-density', 'radial-velocity', 'rot-velocity'], weight='CellMassMsun', x_bins=bin_num1, x_log=False, x_bounds = [rmin,rmax], center = location)

p2 = pc.add_profile_sphere(1., "pc", ['Radius', 'gravitational-potential', 'tracer4', 'dGmag', 'dGmag_G2', 'dGdx', 'dGdy', 'dGdz', 'dGdx_G2', 'dGdy_G2', 'dGdz_G2'], weight='CellMassMsun', x_bins=bin_num1, x_log=False, x_bounds = [rmin,rmax], center = location) 

p3 = pc.add_profile_sphere(1., "pc", ['Radius', 'CellMassMsun'], weight=None, accumulation=True, x_bins=bin_num1, x_log=False, x_bounds = [rmin,rmax], center = location)

pc.save('radprof')

rad1 = p1.data['Radius']
dens = p1.data['density'] 
nh = p1.data['number-density']

rad = p2.data['Radius']
gpot_orion = -p2.data['gravitational-potential']
gpot_gadget = p2.data['tracer4']
gacc_orion = p2.data['dGmag']
gacc_gadget = p2.data['dGmag_G2']
dGdx = p2.data['dGdx']
dGdy = p2.data['dGdy']
dGdz = p2.data['dGdz']
dGdx_G2 = -p2.data['dGdx_G2']
dGdy_G2 = -p2.data['dGdy_G2']
dGdz_G2 = -p2.data['dGdz_G2']

rad_arr = p3.data['Radius']
massenc_arr = p3.data['CellMassMsun']

print 'min dGmag=', min(data['dGmag']), 'max dGmag=', max(data['dGmag'])

####################################################################################
#calculate estimated gravitational potential and gravitational acceleration
gpot = []
gacc = []
rad_pot = []
for i in range(bin_num1):
        rad_pot.append(rad_arr[i])
        rad_in = 1.5e13*rad_arr[i]
        gpot.append(6.67e-8*2.e33*(massenc_arr[i])/rad_in)
        gacc.append(6.67e-8*2.e33*(massenc_arr[i])/(rad_in*rad_in))

#####################################################################################

#####################################################################################

rmin = 1.e2
rmax = 3.e5
nmax=1.e12

pl.subplot(221)
pl.plot(rad, gpot_gadget)
pl.plot(rad, gpot_orion, 'r--')
pl.plot(rad_pot, gpot, 'k')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Radius [AU]', fontsize=9)
ax.set_ylabel('Gravitational Potential [cgs]', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((rmin, rmax, 1e11, 1e13))

pl.subplot(222)
pl.plot(rad, gacc_gadget)
pl.plot(rad, gacc_orion, 'r--')
pl.plot(rad_pot, gacc, 'k')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Radius [AU]', fontsize=9)
ax.set_ylabel('mag(dG/dr)', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((rmin, rmax, 1e-7, 1e-3))

pl.subplot(223)
pl.plot(rad, dGdx_G2)
pl.plot(rad, dGdx, 'r--')
pl.plot(rad_pot, gacc, 'k')
ax = pl.gca()
ax.set_xscale('log')
ax.set_xlabel('Radius [AU]', fontsize=9)
ax.set_ylabel('dGdx [cgs]', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((rmin, rmax, -1e-5, 1e-5))

pl.subplot(224)
pl.plot(rad_arr, massenc_arr,'k')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Radius [AU]', fontsize=9)
ax.set_ylabel('Enclosed Mass [M_{\odot}', fontsize=9)
pl.xticks(fontsize=9)
pl.yticks(fontsize=9)
pl.axis((rmin, rmax, 1.e1, 1e4))

pl.savefig('plot.ps')

M_i = data["CellMassMsun"].sum()
print "total mass = ", M_i
