import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as pl
#import pylab as pl

from yt.mods import *
import numpy as na
from string import rstrip
import fields
import fields_bfield

pc_cm = 3.08567758e18
au_cm = 1.49597871e13
Msun = 1.98892e33
rho_nh = 1./(1.2195*1.6726e-24)

datanum = '0000'

#pf = load("/work/00863/minerva/orion/bfield_ref2_maptest1_2lev_nosink/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion/bfield_comp1/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion/bfield_maptest2/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion/bfield_ref3_maptest1/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion/bfield_ref3_maptest1_2lev/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion_Btest/bfield_orig/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion_Btest/bfield_new/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion_Btest/bfield_homolog/" + 'data.' + datanum + '.3d.hdf5')
pf = load("/work/00863/minerva/orion_Btest/" + 'data.' + datanum + '.3d.hdf5')

print pf.h.field_list

bin_num_square = 32
bin_doub = 32.0
DeltaX = pf.domain_right_edge[0] / bin_doub

rmax = pf.domain_right_edge[0]/pc_cm

value, location = pf.h.find_max("Density")
data = pf.h.sphere([0,0,0], rmax/pf['pc'])
print 'max density location = ', location

pc = PlotCollection(pf, center = location)

print pf.h.field_list


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
    rad = rad/pc_cm
    return rad
add_field("Radius", function=_Radius, take_log=False,
          units=r'\rm{AU}')

print 'temp_min=', min(data['Temperature']), 'temp_max=', max(data['Temperature'])
print 'vrad_min=', min(data['radial-velocity']), 'vrad_max=', max(data['radial-velocity'])
print 'dens_min=', min(data['Density']), 'dens_max=', max(data['Density'])

p1 = pc.add_profile_sphere(rmax, "pc", ['Radius', 'number-density'], weight='CellVolume')
p2 = pc.add_profile_sphere(rmax, "pc", ['Radius', 'radial-velocity'], weight='CellVolume')
p3 = pc.add_profile_sphere(rmax, "pc", ['Radius', 'density', 'Temperature', 'Bmag', 'ufield', 'MagneticEnergy', 'DivB', 'absDivB'], weight='CellVolume')

rad_arr = p3.data['Radius']

#x_gadget = [ -133.51959, -124.90541,  -116.29128,  -107.67710,  -99.062920,  -90.448737,  -81.834614,      -73.220432,      -64.606249,      -55.992067,
#      -47.377944,      -38.763762,      -30.149579,      -21.535397,      -12.921274,      -4.3070912,       4.3070912,       12.921274,       21.535397,       30.149579,
#       38.763762,       47.377944,       55.992067,       64.606249,       73.220432,       81.834614,       90.448737,       99.062920,       107.67710,       116.29128,
#       124.90541,       133.51959] 
#B_gadget = [1.3208950e-09,   1.3207000e-09,   1.3207000e-09,   1.3207000e-09,   1.3207000e-09,   1.3207000e-09,   1.3207000e-09,   1.3207000e-09,   1.3207000e-09,   1.3207000e-09,
#   1.3207000e-09,   1.3207000e-09,   1.3207000e-09,   1.3206850e-09,   1.3186650e-09,   1.2565650e-09,   1.1498600e-09,   8.3262900e-10,   3.0788000e-10,   1.6989800e-10,
#   1.6529500e-10,   1.6508600e-10,   1.6508600e-10,   1.6508600e-10,   1.6508300e-10,   1.6508300e-10,   1.6508300e-10,  1.6508300e-10,   1.6508300e-10,    1.6508300e-10,
#   1.6508300e-10,   1.6527400e-10]
#n_gadget = [0.00052736300,   0.00052736300,   0.00052736300,   0.00052736300,   0.00052736300,   0.00052736300,   0.00052736300,   0.00052736300,   0.00052736300,   0.00052736300,
#   0.00052736300,   0.00052736300,   0.00052736300,   0.00052737150,   0.00052670050,   0.00051934950,   0.00042450600,   0.00024767300,   8.3125150e-05,   6.5921200e-05,
#   6.5921200e-05,   6.5918200e-05,   6.5921200e-05,   6.5921200e-05,   6.5921200e-05,   6.5921200e-05,   6.5921200e-05,   6.5921200e-05,   6.5921200e-05,   6.5921200e-05,
#   6.5921200e-05,   6.5921200e-05]

x_gadget = [-267.03916,      -249.81081,      -232.58251,      -215.35414,      -198.12584,      -180.89747,      -163.66917,      -146.44080,      -129.21250,      -111.98413,      -94.755828, -77.527463,      -60.299158,      -43.070793,      -25.842488,      -8.6141527,       8.6141825,       25.842547,       43.070912,       60.299158,       77.527523,       94.755888,       111.98425,   129.21250,       146.44086,       163.66923,       180.89759,       198.12584,       215.35420,       232.58257,       249.81087,       267.03918]

B_gadget = [1.3207000e-09,   1.3207000e-09,   1.3207000e-09,   1.3207000e-09,   1.3207000e-09,   1.3207000e-09,   1.3207000e-09,   1.3207000e-09,   1.3207000e-09,   1.3207000e-09,   1.3207000e-09,  1.3207000e-09,   1.3207000e-09,   1.3207000e-09,   1.3206925e-09,   1.2876150e-09,   9.9124450e-10,   2.6188600e-10,   1.6519050e-10,   1.6509400e-10,   1.6508300e-10,   1.6508300e-10,   1.6508300e-10,   1.6508300e-10,   1.6508300e-10,   1.6508300e-10,   1.6508300e-10,   1.6508300e-10,  1.6508300e-10,   1.6508300e-10,   1.6508300e-10,   1.6508300e-10]

n_gadget = [0.00052736300,   0.00052736300,   0.00052736300,   0.00052736300,   0.00052736300,   0.00052736300,   0.00052736300,   0.00052736300,   0.00052736300,   0.00052736300,   0.00052736300, 0.00052736300,   0.00052736300,   0.00052736300,   0.00052736835,   0.00051949084,   0.00036580593,   7.2425333e-05,   6.5920200e-05,   6.5921200e-05,   6.5921200e-05,   6.5921200e-05,   6.5921200e-05, 6.5921200e-05,   6.5921200e-05,   6.5921200e-05,   6.5921200e-05,   6.5921200e-05,   6.5921200e-05,   6.5921200e-05,   6.5921200e-05,   6.5921200e-05]

#n_gadget = [0.00052021465,   0.00052021465,   0.00052021465,   0.00052021465,   0.00052021465,   0.00052021465,   0.00052021465,   0.00052021465,   0.00052021465,   0.00052021465,   0.00052021465, 0.00052021465,   0.00052021465,   0.00052021465,   0.00052021465,   0.00052021549,   0.00016229686,   9.6730332e-05,   6.4487112e-05,   6.4487112e-05,   6.4486888e-05,   6.4486888e-05,  6.4487112e-05,   6.4487112e-05,   6.4486888e-05,   6.4486888e-05,   6.4487112e-05,   6.4487112e-05,   6.4486888e-05,   6.4477161e-05,   6.4477161e-05,   6.4477161e-05]

#n_gadget = [0.00039825796,   0.00039825796,   0.00039825796,   0.00039825796,   0.00039825796,   0.00039825796,   0.00039825796,   0.00039825796,   0.00039825796,   0.00039825796,   0.00039825796,   0.00039825796,   0.00039825796,   0.00039825796,   0.00039826200,   0.00039426737,   0.00025594247,   5.8388772e-05,   4.9774585e-05,   4.9778752e-05,  4.9778997e-05,   4.9778997e-05,   4.9778997e-05,   4.9778997e-05,   4.9778997e-05,   4.9778997e-05,   4.9778997e-05,   4.9778997e-05,   4.9778997e-05,   4.9778997e-05,   4.9778997e-05,   4.9778997e-05]


# Chose center of data (selecting highest density point)
c = pf.h.find_max('Density')[1]
rayx = pf.h.ortho_ray(0,(c[1],c[2]))
#rayx = pf.h.ortho_ray(0,(0,0))
# Save the distance axes since they are used often
unit = 'pc'
axX = rayx['x']*pf[unit]

density = rayx['density']
nh = density * rho_nh
bfield = rayx['Bmag']
ufield = rayx['ufield'] 
ufield_norm = ufield / (density**1.33333)
absDivB = rayx['absDivB']

errDivB = absDivB * DeltaX / bfield

print 'bfield =', bfield
print 'ufield_norm =', ufield_norm
print 'error divB =', errDivB
print 'total divB = ', data['DivB'].sum()

dens_min = 1.e-28
dens_max = 1.e-26

nmin = 1.e-5
nmax = 1.e-3

bmin = 1.e-10
bmax = 1.e-8

rmin = -rmax
rmax = rmax

umin = 1.e16
umax = 1.e18

DivBmin = -1.e-22
DivBmax = 1.e-22

pl.subplot(221)
#pl.plot(axX, density,'k')
pl.plot(axX, nh,'k')
pl.plot(x_gadget, n_gadget)
ax = pl.gca()
#ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('X [pc]', fontsize=9)
ax.set_ylabel('nh [cm^-3]', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
#pl.axis((rmin, rmax, dens_min, dens_max))
pl.axis((rmin, rmax, nmin, nmax))

pl.subplot(222)
pl.plot(axX, bfield, 'k')
pl.plot(x_gadget, B_gadget)
ax = pl.gca()
#ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel('B-field [G]', fontsize=9)
ax.set_xlabel('X [pc]', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((rmin, rmax, bmin, bmax))

#pl.subplot(222)
#pl.plot(axX, ufield_norm,'k')
#ax = pl.gca()
##ax.set_xscale('log')
#ax.set_yscale('log')
#ax.set_xlabel('axX [pc]', fontsize=9)
#ax.set_ylabel('U / rho^(4/3) [cgs]', fontsize=9)
#pl.xticks(fontsize=10)
#pl.yticks(fontsize=10)
#pl.axis((rmin, rmax, umin, umax))

pl.subplot(223)
pl.plot(nh, bfield, 'k')
pl.plot(n_gadget, B_gadget)
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('nh [cm^-3]', fontsize=9)
ax.set_ylabel('B-field [G]', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((nmin, nmax, bmin, bmax))

pl.subplot(224)
pl.plot(axX, errDivB,'k')
ax = pl.gca()
#ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('X [pc]', fontsize=9)
ax.set_ylabel('error DivB (divB * cell_size / B)', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((rmin, rmax, 1.e-10, 1.e-2))

pl.savefig('bfield.eps')

M_i = data["CellMassMsun"].sum()
print "total mass = ", M_i
