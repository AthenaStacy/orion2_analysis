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

rdir = '/home1/00863/minerva/gadget_runs/'

datanum = '0000'
gnum = '0000'

pf = load("/work/00863/minerva/orion/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion_Btest/bfield_homolog/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion_Btest/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion_Gtest/" + 'data.' + datanum + '.3d.hdf5')

print pf.h.field_list

bin_num = 128
bin_doub = 128.0
DeltaX = pf.domain_right_edge[0] / bin_doub

rmax = pf.domain_right_edge[0]/pc_cm

value, location = pf.h.find_max("Density")
data = pf.h.sphere(location, rmax/pf['pc'])
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

p1 = pc.add_profile_sphere(rmax, "pc", ['Radius', 'number-density'], weight='CellVolume', x_bins=bin_num)
p2 = pc.add_profile_sphere(rmax, "pc", ['Radius', 'radial-velocity', 'x-velocity', 'y-velocity', 'z-velocity', 'Temperature'], weight='CellVolume', x_log=False, x_bins=bin_num)
p3 = pc.add_profile_sphere(rmax, "pc", ['Radius', 'density', 'Temperature', 'Bmag', 'ufield', 'MagneticEnergy', 'DivB', 'absDivB'], weight='CellVolume', x_log=False, x_bins=bin_num)

rad_arr = p3.data['Radius']

# Chose center of data (selecting highest density point)
c = pf.h.find_max('Density')[1]
rayx = pf.h.ortho_ray(0,(c[1],c[2]))
#rayx = pf.h.ortho_ray(0,(0,0))
# Save the distance axes since they are used often
unit = 'pc'

take_ray = 1
take_profile = 0

if take_ray == 1:
	axX = rayx['x']*pf[unit]
	density = rayx['density']
	bfield =  rayx['Bmag']
	ufield =  rayx['ufield'] 
	absDivB = rayx['absDivB']
        vx =      rayx['x-velocity'] / 1.e5
        vy =      rayx['y-velocity'] / 1.e5
        vz =      rayx['z-velocity'] / 1.e5
        temp =    rayx['Temperature']
	rmin = -rmax

if take_profile == 1:
	axX =     p3.data['Radius']
	density = p3.data['density']
        bfield =  p3.data['Bmag']
        ufield =  p3.data['ufield']
        absDivB = p3.data['absDivB']	
	vx =      p2.data['x-velocity'] / 1.e5
        vy =      p2.data['y-velocity'] / 1.e5
        vz =      p2.data['z-velocity'] / 1.e5
	temp =    p2.data['Temperature'] 
	rmin = 0

nh = density * rho_nh
ufield_norm = ufield / (density**1.33333)
errDivB = absDivB * DeltaX / bfield


print 'bfield =', bfield
print 'ufield_norm =', ufield_norm
print 'error divB =', errDivB
print 'total divB = ', data['DivB'].sum()
print 'velx = ', vx
print 'velz = ', vz

###################################################################
with open(rdir+"bfield_comp_"+gnum, "r") as f:
        Bdat = [map(float, line.split()) for line in f]

narr1 = bin_num
if take_profile == 1: 
	narr1 = narr1*2

line_arr = []
x_gadget = []
n_gadget = []
B_gadget = []
vx_gadget = []
vy_gadget = []
vz_gadget = []
temp_gadget = []

for i in range(narr1):
        line_arr = Bdat[i]
        x_gadget.append(line_arr[0])
        n_gadget.append(line_arr[1])
        B_gadget.append(line_arr[5])
	vx_gadget.append(line_arr[6])
	vy_gadget.append(line_arr[7])
	vz_gadget.append(line_arr[8])
	temp_gadget.append(line_arr[9])
#####################################################################

dens_min = 1.e-28
dens_max = 1.e-26

nmin = 1.e3
nmax = 1.e13

bmin = 1.e-10
bmax = 1.e-5

umin = 1.e8
umax = 1.e15

DivBmin = -1.e-22
DivBmax = 1.e-22

tmin = 0
tmax = 2000

vmin = -10
vmax = 10

plot_B = 0
plot_vel = 1

size = 8.5

if plot_vel == 1:

        pl.subplot(321)
        #pl.plot(axX, density,'k')
        pl.plot(axX, nh,'k')
        pl.plot(x_gadget, n_gadget)
        ax = pl.gca()
        ax.set_yscale('log')
        ax.set_xlabel('Radius [pc]', fontsize=size)
        ax.set_ylabel('nh [cm^-3]', fontsize=size)
        pl.xticks(fontsize=size)
        pl.yticks(fontsize=size)
        #pl.axis((rmin, rmax, dens_min, dens_max))
        pl.axis((rmin, rmax, nmin, nmax))

        pl.subplot(322)
        pl.plot(axX, temp, 'k')
        pl.plot(x_gadget, temp_gadget)
        ax = pl.gca()
        ax.set_ylabel('Temp [K]', fontsize=size)
        ax.set_xlabel('Radius [pc]', fontsize=size)
        pl.xticks(fontsize=size)
        pl.yticks(fontsize=size)
        pl.axis((rmin, rmax, tmin, tmax))

        pl.subplot(323)
        pl.plot(axX, vx, 'k')
        pl.plot(x_gadget, vx_gadget)
        ax = pl.gca()
        ax.set_ylabel('x-velocity [km / s]', fontsize=size)
        ax.set_xlabel('Radius [pc]', fontsize=size)
        pl.xticks(fontsize=size)
        pl.yticks(fontsize=size)
        pl.axis((rmin, rmax, vmin, vmax))

        pl.subplot(324)
        pl.plot(axX, vy, 'k')
        pl.plot(x_gadget, vy_gadget)
        ax = pl.gca()
        ax.set_ylabel('y-velocity [km / s]', fontsize=size)
        ax.set_xlabel('Radius [pc]', fontsize=size)
        pl.xticks(fontsize=size)
        pl.yticks(fontsize=size)
        pl.axis((rmin, rmax, vmin, vmax))

        pl.subplot(325)
        pl.plot(axX, vz, 'k')
        pl.plot(x_gadget, vz_gadget)
        ax = pl.gca()
        ax.set_ylabel('z-velocity [km / s]', fontsize=size)
        ax.set_xlabel('Radius [pc]', fontsize=size)
        pl.xticks(fontsize=size)
        pl.yticks(fontsize=size)
        pl.axis((rmin, rmax, vmin, vmax))



if plot_B == 1:

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
	pl.axis((rmin, rmax, 1.e-10, 1.e0))

pl.savefig('bfield.eps')

M_i = data["CellMassMsun"].sum()
print "total mass = ", M_i
