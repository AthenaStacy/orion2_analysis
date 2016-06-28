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

#wdir = '/work/00863/minerva/orion/maptest_iso_ref_newvel/'
wdir = '/work/00863/minerva/orion/maptest_ideal32/'
#wdir = '/work/00863/minerva/orion/'
rdir = '/home1/00863/minerva/research_programs/'
snapnum = '000'
setup = 2

#pf = load("/nobackupp7/astacy/orion/chemtest_iso/data.0021.3d.hdf5")
pf = load(wdir+"data.0000.3d.hdf5")
value, location = pf.h.find_max("Density")
data = pf.h.sphere(location, 3.0/pf['pc'])
pc = PlotCollection(pf, center = location)

print pf.h.field_list

#location[0] = location[1] = location[2] = 0

x=location[0]
y=location[1]
z=location[2]

massenc_com = []
vx = []
vy = []
vz = []
sub_sphere =  pf.h.sphere(location, 3.e18/pf['cm'])
massenc_com.append(sub_sphere["CellMassMsun"].sum())

#vx.append((sub_sphere["x-velocity"]*sub_sphere["CellMassMsun"]).sum()/(sub_sphere["CellMassMsun"].sum()))
#vy.append((sub_sphere["y-velocity"]*sub_sphere["CellMassMsun"]).sum()/(sub_sphere["CellMassMsun"].sum()))
#vz.append((sub_sphere["z-velocity"]*sub_sphere["CellMassMsun"]).sum()/(sub_sphere["CellMassMsun"].sum()))

vx.append(0)
vy.append(0)
vz.append(0)

#x=0
#y=0
#z=0

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
    rad = rad/au_cm
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
    vrotx = vrotx / na.sqrt((data['x'] - x)**2 + (data['y'] - y)**2 + (data['z'] - z)**2)
    vroty = vroty / na.sqrt((data['x'] - x)**2 + (data['y'] - y)**2 + (data['z'] - z)**2)
    vrotz = vrotz / na.sqrt((data['x'] - x)**2 + (data['y'] - y)**2 + (data['z'] - z)**2)
    vrot = na.sqrt(vrotx**2 + vroty**2 + vrotz**2)
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


#print 'pres_min=', min(data['Pressure']), 'pres_max=', max(data['Pressure'])
print 'temp_min=', min(data['Temperature']), 'temp_max=', max(data['Temperature'])
print 'vrad_min=', min(data['radial-velocity']), 'vrad_max=', max(data['radial-velocity'])
print 'dens_min=', min(data['density']), 'dens_max=', max(data['density'])
print 'rad_min', min(data['Radius']), 'rad_max', max(data['Radius'])
print 'velz_min', min(data['velz']), 'velz_max', max(data['velz'])

################read in other data files#######################################

with open(rdir+"bin_HR10_map_vel_000.dat.lin", "r") as f:
	veldatA = [map(float, line.split()) for line in f]

with open(rdir+'bin_HR10_map_vel_'+snapnum+'.dat'+'.lin', "r") as f:
        veldatB = [map(float, line.split()) for line in f]

rmencA = []
radA = []
csA = []
turbA = []
velxA = []
velyA = []
velzA = []
vradA = []
vrotA = []
vradxA = []
vradyA = []
vradzA = []
vrotxA = []
vrotyA = []
vrotzA = []
amomA = []

rmencB = []
radB = []
csB = []
turbB = []
velxB = []
velyB = []
velzB = []
vradB = []
vrotB = []
vradxB = []
vradyB = []
vradzB = []
vrotxB = []
vrotyB = []
vrotzB = []
amomB = []

narr = 200
narr2 = 200
line_arr = []
fac = 1.2195*1.67e-24

for i in range(narr):
        line_arr = veldatA[i]
	radA.append(line_arr[0])
        csA.append(line_arr[1])
        turbA.append(line_arr[2])
        vradA.append(line_arr[3])
        vrotA.append(line_arr[4])
        velxA.append(line_arr[5])
        velyA.append(line_arr[6])
        velzA.append(line_arr[7])
        vradxA.append(line_arr[8])
        vradyA.append(line_arr[9])
        vradzA.append(line_arr[10])
        vrotxA.append(line_arr[11])
        vrotyA.append(line_arr[12])
        vrotzA.append(line_arr[13])
	amomB.append(line_arr[14])
	turbA[i] = turbA[i]/csA[i]

for i in range(narr2):
	line_arr = veldatB[i]
	radB.append(line_arr[0])
        csB.append(line_arr[1])
        turbB.append(line_arr[2])
        vradB.append(line_arr[3])
        vrotB.append(line_arr[4])
        velxB.append(line_arr[5])
        velyB.append(line_arr[6])
        velzB.append(line_arr[7])
        vradxB.append(line_arr[8])
        vradyB.append(line_arr[9])
        vradzB.append(line_arr[10])
        vrotxB.append(line_arr[11])
        vrotyB.append(line_arr[12])
        vrotzB.append(line_arr[13])
	amomB.append(line_arr[14])
	turbB[i] = turbB[i]/csB[i]

#print 'radA =', radA
#print 'vradA = ', vradA
##############################################################################

rmin = 100
rmax = 3e5
nmin = 1e2
nmax = 1e16
bin_num1 = 200
bin_num2 = 200

p1 = pc.add_profile_sphere(1., "pc", ['Radius', 'number-density', 'radial-velocity', 'rot-velocity', 'velx', 'vely', 'velz', 'vradx', 'vrady', 'vradz', 'vrotx', 'vroty', 'vrotz', 'cs'], weight='CellMassMsun', x_bins=bin_num1, x_log=False, x_bounds = [rmin,rmax], center = location)
print 'location =', location

p3 = pc.add_profile_sphere(1., "pc", ['Radius', 'CellMassMsun'], weight=None, accumulation=True, x_bins=bin_num1, x_log=False, x_bounds = [rmin,rmax], center = location)

pc.save('radprof')

rad = p1.data['Radius']
nh = p1.data['number-density']
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

rad_arr = p3.data['Radius']
massenc_arr = p3.data['CellMassMsun']

vrot = na.power(vrotx*vrotx + vroty*vroty + vrotz*vrotz, 0.5)

turb = []
rad_turb = []
amom = []
amom.append(massenc_arr[i]*vrot[i]*rad[i]*au_cm*1.e5*2.e33)
for i in range(bin_num1-1):
        turb_arr = (velx[i] - vradx[i] - vrotx[i])**2 + (vely[i] - vrady[i] - vroty[i])**2 + (velz[i] - vradz[i] - vrotz[i])**2;
        turb.append(na.sqrt(turb_arr)/cs[i])
        rad_turb.append(rad[i])
        if i > 0:
                amom.append((massenc_arr[i] - massenc_arr[i-1])*vrot[i]*rad[i]*au_cm*1.e5*2.e33)

##############################################################################
#calculate errors
rad_err = []
vrot_err = []
velx_err = []
vely_err = []
velz_err = []
turb_err = []
cs_err = []

for i in range(narr-1):
        print rad[i], radB[i], vrot[i], vrotB[i], velx[i], velxB[i], turb[i], turbB[i]
        if nh[i] > 0:
                rad_err.append(rad[i])
                vrot_err.append((vrot[i] - vrotB[i])/vrotB[i])
                velx_err.append((velx[i] - velxB[i])/velxB[i])
                vely_err.append((vely[i] - velyB[i])/velyB[i])
                velz_err.append((velz[i] - velzB[i])/velzB[i])
                turb_err.append((turb[i] - turbB[i])/turbB[i])
                cs_err.append((cs[i] - csB[i])/csB[i])
##############################################################################

rmin = 1.e2
rmax = 3.e5

vmax = 5.0
vmin = -5.0

errmin = -0.1
errmax = 0.1

if setup == 0:
	pl.subplot(221)
	pl.plot(rad, vrad, 'k')
	pl.plot(rmencA, vradA)
	pl.plot(rmencB, vradB, 'r--')
	ax = pl.gca()
	ax.set_xscale('log')
	ax.set_xlabel('Radius [AU]', fontsize=9)
	ax.set_ylabel('radial-velocity [km s^{-1}]', fontsize=9)
	pl.xticks(fontsize=10)
	pl.yticks(fontsize=10)
	pl.axis((rmin, rmax, vmin, vmax))

if setup == 1:
        pl.subplot(221)
        pl.plot(rad_err, vrot_err)
        ax = pl.gca()
        ax.set_xscale('log')
        ax.set_xlabel('Radius [AU]', fontsize=9)
        ax.set_ylabel('rot-velocity error', fontsize=9)
        pl.xticks(fontsize=10)
        pl.yticks(fontsize=10)
        pl.axis((rmin, rmax, errmin,errmax))

if setup == 0 or setup == 1:
	pl.subplot(222)
	pl.plot(rad_err, velx_err)
	ax = pl.gca()
	ax.set_xscale('log')
	ax.set_xlabel('Radius [AU]', fontsize=9)
	ax.set_ylabel('x-velocity error', fontsize=9)
	pl.xticks(fontsize=10)
	pl.yticks(fontsize=10)
	pl.axis((rmin, rmax, errmin, errmax))

	pl.subplot(223)
	pl.plot(rad_err, vely_err)
	ax = pl.gca()
	ax.set_xscale('log')
	ax.set_xlabel('Radius [AU]', fontsize=9)
	ax.set_ylabel('y-velocity error', fontsize=9)
	pl.xticks(fontsize=10)
	pl.yticks(fontsize=10)
	pl.axis((rmin, rmax, errmin, errmax))

	pl.subplot(224)
	pl.plot(rad_err, velz_err)
	ax = pl.gca()
	ax.set_xscale('log')
	ax.set_xlabel('Radius [AU]', fontsize=9)
	ax.set_ylabel('z-velocity error', fontsize=9)
	pl.xticks(fontsize=10)
	pl.yticks(fontsize=10)
	pl.axis((rmin, rmax, errmin, errmax))


if setup == 2:
        pl.subplot(221)
        pl.plot(rad_err, vrot_err)
        ax = pl.gca()
        ax.set_xscale('log')
        ax.set_xlabel('Radius [AU]', fontsize=9)
        ax.set_ylabel('rot-velocity error', fontsize=9)
        pl.xticks(fontsize=10)
        pl.yticks(fontsize=10)
        pl.axis((rmin, rmax, errmin,errmax))

        pl.subplot(222)
        pl.plot(rad_err, turb_err)
        ax = pl.gca()
        ax.set_xscale('log')
        ax.set_xlabel('Radius [AU]', fontsize=9)
        ax.set_ylabel('turbulence error', fontsize=9)
        pl.xticks(fontsize=10)
        pl.yticks(fontsize=10)
        pl.axis((rmin, rmax, errmin,errmax))

        pl.subplot(223)
        pl.plot(rad_err, cs_err)
        ax = pl.gca()
        ax.set_xscale('log')
        ax.set_xlabel('Radius [AU]', fontsize=9)
        ax.set_ylabel('soundspeed error', fontsize=9)
        pl.xticks(fontsize=10)
        pl.yticks(fontsize=10)
        pl.axis((rmin, rmax, errmin,errmax))

pl.savefig('plot.ps')

M_i = data["CellMassMsun"].sum()
print "total mass = ", M_i
