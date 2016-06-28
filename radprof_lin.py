import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as pl
#import pylab as pl

from yt.mods import *
import numpy as na
from string import rstrip
import fields
#import fields_iso

au_cm = 1.49597871e13

#wdir = '/work/00863/minerva/orion//maptest_iso_ref_newvel/'
wdir = '/work/00863/minerva/orion/maptest_ideal/'
rdir = '/home1/00863/minerva/research_programs/'
snapnum = '000'
setup = 3

#pf = load("/nobackupp7/astacy/orion/chemtest_iso/data.0021.3d.hdf5")
pf = load(wdir+"data.0000.3d.hdf5")
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

######x vel's
def _xVelocity(field, data):
    velx = (data['x-velocity'] - vx[0])
    return velx/1.e5
add_field("velx", function=_xVelocity, take_log=False,
          units=r'\rm{km}/\rm{s}')

def _xvrad(field, data):
    vradx = data['radial-velocity']*(data['x'] - x)/(au_cm*data['Radius'])
    return vradx
add_field("vradx", function=_xvrad, take_log=False,
          units=r'\rm{km}/\rm{s}')

def _xvrot(field, data):
    vrotx = (data['z-velocity'] - vz[0])*(data['y'] - y) - (data['y-velocity'] - vy[0])*(data['z'] - z)
    return vrotx/1.e5/(au_cm*data['Radius'])
add_field("vrotx", function=_xvrot, take_log=False,
          units=r'\rm{km}/\rm{s}')

#######y vel's
def _yVelocity(field, data):
    vely = (data['y-velocity'] - vy[0])
    return vely/1.e5
add_field("vely", function=_yVelocity, take_log=False,
          units=r'\rm{km}/\rm{s}')

def _yvrad(field, data):
    vrady = data['radial-velocity']*(data['y'] - y)/(au_cm*data['Radius'])
    return vrady
add_field("vrady", function=_yvrad, take_log=False,
          units=r'\rm{km}/\rm{s}')

def _yvrot(field, data):
    vroty = (data['x-velocity'] - vx[0])*(data['z'] - z) - (data['z-velocity'] - vz[0])*(data['x'] - x)
    return vroty/1.e5/(au_cm*data['Radius'])
add_field("vroty", function=_yvrot, take_log=False,
          units=r'\rm{km}/\rm{s}')

######z vel's
def _zVelocity(field, data):
    velz = (data['z-velocity'] - vz[0])
    return velz/1.e5
add_field("velz", function=_zVelocity, take_log=False,
          units=r'\rm{km}/\rm{s}')

def _zvrad(field, data):
    vradz = data['radial-velocity']*(data['z'] - z)/(au_cm*data['Radius'])
    return vradz
add_field("vradz", function=_zvrad, take_log=False,
          units=r'\rm{km}/\rm{s}')

def _zvrot(field, data):
    vrotz = (data['y-velocity'] - vy[0])*(data['x'] - x) - (data['x-velocity'] - vx[0])*(data['y'] - y)
    return vrotz/1.e5/(au_cm*data['Radius'])
add_field("vrotz", function=_zvrot, take_log=False,
          units=r'\rm{km}/\rm{s}')

def _rotVelocity(field, data):
    '''
    Rotational velocity. In this problem the center is at density peak.
    '''
    vel_rot = na.power(data['vrotx']*data['vrotx'] + data['vroty']*data['vroty'] + data['vrotz']*data['vrotz'], 0.5)
    return vel_rot
add_field("rot-velocity", function=_rotVelocity, take_log=False,
          units=r'\rm{km}/\rm{s}')


################read in other data files#######################################

with open(rdir+"bin_HR10_map_gas_000.dat.lin", "r") as f:
	gasdatA = [map(float, line.split()) for line in f]

with open(rdir+'bin_HR10_map_vel_000.dat.lin', "r") as f:
        veldatA = [map(float, line.split()) for line in f]

with open(rdir+'bin_HR10_map_gas_'+snapnum+'.dat.lin', "r") as f:
        gasdatB = [map(float, line.split()) for line in f]

with open(rdir+'bin_HR10_map_vel_'+snapnum+'.dat.lin', "r") as f:
        veldatB = [map(float, line.split()) for line in f]

with open(rdir+'bin_HR10_map_gas_nh_000.dat', "r") as f:
        gasdat_nhA = [map(float, line.split()) for line in f]

with open(rdir+'bin_HR10_map_gas_nh_'+snapnum+'.dat.lin', "r") as f:
        gasdat_nhB = [map(float, line.split()) for line in f]

radA = []
nhA = []
densA = []
tempA = []
rmencA = []
csA = []
vradA = []
vrotA = []
turbA = []
mencA = []
presA = []
nh_nhA = []
temp_nhA = []

radB = []
nhB = []
densB = []
tempB = []
rmencB = []
csB = []
vradB = []
vrotB = []
turbB = []
mencB = []
presB = []
nh_nhB = []
temp_nhB = []

narr1 = 200
line_arr = []
fac = 1.2195*1.67e-24
for i in range(narr1):
        line_arr = gasdatA[i]
	radA.append(line_arr[0])
	nhA.append(line_arr[1])
        densA.append(line_arr[1]*fac)
	tempA.append(line_arr[2])
	mencA.append(line_arr[6])
	presA.append(1.38e-16*nhA[i]*tempA[i])

for i in range(narr1):
        line_arr = veldatA[i]
	rmencA.append(line_arr[0])
        csA.append(line_arr[1])
        turbA.append(line_arr[2])
        vradA.append(line_arr[3])
        vrotA.append(line_arr[4])
        turbA[i] = turbA[i]/csA[i]


narr = 200
for i in range(narr):
	line_arr = gasdatB[i]
	radB.append(line_arr[0])
	nhB.append(line_arr[1])
        densB.append(line_arr[1]*fac)
	tempB.append(line_arr[2])
	mencB.append(line_arr[6])
	presB.append(1.38e-16*nhB[i]*tempB[i])

for i in range(narr):
	line_arr = veldatB[i]
        rmencB.append(line_arr[0])
        csB.append(line_arr[1])
        turbB.append(line_arr[2])
        vradB.append(line_arr[3])
        vrotB.append(line_arr[4])
        turbB[i] = turbB[i]/csB[i]

for i in range(narr1):
        line_arr = gasdat_nhA[i]
        nh_nhA.append(line_arr[0])
        temp_nhA.append(line_arr[1])

for i in range(narr):
        line_arr = gasdat_nhB[i]
        nh_nhB.append(line_arr[0])
        temp_nhB.append(line_arr[1])

print 'radA =', radA
print 'vradA = ', vradA
##############################################################################

rmin = 100.
rmax = 3.e5
nmin = 1e2
nmax = 1e16
#bin_num1 = 7500
bin_num1 = 200
bin_num2 = 200

p1 = pc.add_profile_sphere(1., "pc", ['Radius', 'density', 'number-density', 'radial-velocity', 'rot-velocity', 'Temperature', 'Pressure', 'velx', 'vely', 'velz', 'vradx', 'vrady', 'vradz', 'vrotx', 'vroty', 'vrotz', 'cs'], weight='CellMassMsun', x_bins=bin_num1, x_log=False, x_bounds = [rmin,rmax], center = location)

p2 = pc.add_profile_sphere(1., "pc", ['number-density', 'Temperature'], weight='CellMassMsun')

p3 = pc.add_profile_sphere(1., "pc", ['Radius', 'CellMassMsun'], weight=None, accumulation=True, x_bins=bin_num1, x_log=False, x_bounds = [rmin,rmax], center = location)

pc.save('radprof')

rad = p1.data['Radius']
dens = p1.data['density'] 
nh = p1.data['number-density']
temp = p1.data['Temperature']

vrad = p1.data['radial-velocity'] 
vrot = p1.data['rot-velocity']
velx = p1.data['velx']
vely = p1.data['vely']
velz = p1.data['velz']
vradx = p1.data['vradx']
vrady = p1.data['vrady']
vradz = p1.data['vradz']
vrotx = p1.data['vrotx']
vroty = p1.data['vroty']
vrotz = p1.data['vrotz']
cs = p1.data['cs']
cs = cs/1.e5

pres = p1.data['Pressure']

nh_alt = p2.data['number-density']
temp_alt = p2.data['Temperature']

rad_arr = p3.data['Radius']
massenc_arr = p3.data['CellMassMsun']

vrot = na.power(vrotx*vrotx + vroty*vroty + vrotz*vrotz, 0.5)

turb = []
rad_turb = []
amom = []
amom.append(massenc_arr[i]*vrot[i]*rad[i]*1.5e13*1.e5*2.e33)
for i in range(bin_num1-1):
        turb_arr = (velx[i] - vradx[i] - vrotx[i])**2 + (vely[i] - vrady[i] - vroty[i])**2 + (velz[i] - vradz[i] - vrotz[i])**2;
        turb.append(na.sqrt(turb_arr)/cs[i])
        rad_turb.append(rad[i])
        if i > 0:
                amom.append((massenc_arr[i] - massenc_arr[i-1])*vrot[i]*rad[i]*1.5e13*1.e5*2.e33)
####################################################################################
#alternate calculation for density
dens2 = []
nh2 = []
rad2 = []
for i in range(bin_num1 - 1):
        rad_out = 1.5e13*rad[i+1]
        rad_in = 1.5e13*rad[i]
        rad2.append(rad_in/1.5e13)
        volume = (4./3.)*3.14159*(rad_out**3 - rad_in**3)
        dens2.append(2.e33*(massenc_arr[i+1] - massenc_arr[i])/volume)
        nh2.append(dens2[i]/fac)

#####################################################################################


#print 'radA =', radA
#print 'vradA = ', vradA
#print 'nhA =', nhA
#print 'tempA =', tempA
#print 'rad =', rad
#print 'nh =', nh

rmin = 1.e2
rmax = 3.e5
nmax=1.e12

pl.subplot(221)
#pl.plot(rad_nh, nh, 'k')
pl.plot(rad2, nh2, 'k')
#pl.plot(radA, nhA)
pl.plot(radB, nhB, 'r--')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Radius [AU]', fontsize=9)
ax.set_ylabel('number density [cm^{-3}]', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((rmin, rmax, 1e3, nmax))

if setup== 3:
        pl.subplot(222)
        pl.plot(rad, temp,'k')
        #pl.plot(radA, tempA)
        pl.plot(radB, tempB, 'r--')
        ax = pl.gca()
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Radius [AU]', fontsize=9)
        ax.set_ylabel('Temperature [K]', fontsize=9)
        pl.xticks(fontsize=10)
        pl.yticks(fontsize=10)
        pl.axis((rmin, rmax, 1e2, 1e4))

if setup== 4:
	pl.subplot(222)
	pl.plot(nh2, temp2,'k')
	#pl.plot(nh_nhA, temp_nhA)
	pl.plot(nh_nhB, temp_nhB, 'r--')
	ax = pl.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel('number density [cm^{-3}]', fontsize=9)
	ax.set_ylabel('Temperature [K]', fontsize=9)
	pl.xticks(fontsize=10)
	pl.yticks(fontsize=10)
	pl.axis((1e2, 1e12, 1e2, 1e4))

if setup== 3:
        pl.subplot(223)
        pl.plot(rad, vrad, 'k')
        #pl.plot(rmencA, vradA)
        pl.plot(rmencB, vradB, 'r--')
        ax = pl.gca()
        ax.set_xscale('log')
        ax.set_xlabel('Radius [AU]', fontsize=9)
        ax.set_ylabel('radial-velocity [km s^{-1}]', fontsize=9)
        pl.xticks(fontsize=10)
        pl.yticks(fontsize=10)
        pl.axis((rmin, rmax, -6, 6))


if setup== 4:
	pl.subplot(223)
	pl.plot(rad, vrot, 'k')
	#pl.plot(rmencA, vrotA)
	pl.plot(rmencB, vrotB, 'r--')
	ax = pl.gca()
	ax.set_xscale('log')
	ax.set_xlabel('Radius [AU]', fontsize=9)
	ax.set_ylabel('rotat-velocity [km s^{-1}]', fontsize=9)
	pl.xticks(fontsize=10)
	pl.yticks(fontsize=10)
	pl.axis((rmin, rmax, 0, 10))


if setup == 1:
	pl.subplot(224)
	pl.plot(radc, pres,'k')
	#pl.plot(radA, presA)
	pl.plot(radB, presB, 'r--')
	ax = pl.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel('Radius [AU]', fontsize=9)
	ax.set_ylabel('Pressure [dyn cm^{-2}]', fontsize=9)
	pl.xticks(fontsize=10)
	pl.yticks(fontsize=10)
	pl.axis((rmin, rmax, 1e-10, 1e-3))

if setup == 2 or setup == 3:
	pl.subplot(224)
	pl.plot(rad_arr, massenc_arr,'k')
	#pl.plot(rmencA, mencA)
	pl.plot(rmencB, mencB, 'r--')
	ax = pl.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel('Radius [AU]', fontsize=9)
	ax.set_ylabel('Enclosed Mass [M_{\odot}', fontsize=9)
	pl.xticks(fontsize=9)
	pl.yticks(fontsize=9)
	pl.axis((rmin, rmax, 1.1, 3e3))

pl.savefig('plot.ps')

M_i = data["CellMassMsun"].sum()
print "total mass = ", M_i
