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
Msun = 1.98892e33

#wdir = '/work/00863/minerva/orion//maptest_iso_ref_newvel/'
#wdir = '/work/00863/minerva/orion/gravpot_0.8pc_128/'
wdir = '/work/00863/minerva/orion/'
rdir = '/home1/00863/minerva/research_programs/'
file_prefix = 'bin_HR10'
snapnum = '0007'
setup = 3

#pf = load("/nobackupp7/astacy/orion/chemtest_iso/data.0021.3d.hdf5")
pf = load(wdir+"data.0000.3d.hdf5")
value, location = pf.h.find_max("Density")
data = pf.h.sphere(location, 2.0/pf['pc'])
pc = PlotCollection(pf, center = location)

print pf.h.field_list

bin_num_square = 64
bin_doub = 64.0
DeltaX = pf.domain_right_edge[0] / bin_doub

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

def _xmomentum(field, data):
    return (data["X-momentum"]*data["CellMassMsun"]/data["density"])
add_field("xmomentum", function=_xmomentum, take_log=False, units=r'\rm{km}/\rm{s}')

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

def _ymomentum(field, data):
    return (data["Y-momentum"]*data["CellMassMsun"]/data["density"])
add_field("ymomentum", function=_ymomentum, take_log=False, units=r'\rm{km}/\rm{s}')

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

def _zmomentum(field, data):
    return (data["Z-momentum"]*data["CellMassMsun"]/data["density"])
add_field("zmomentum", function=_zmomentum, take_log=False, units=r'\rm{km}/\rm{s}')

################read in other data files#######################################

with open(rdir+ 'bin_HR10_gas_0000.dat', "r") as f:
	gasdatA = [map(float, line.split()) for line in f]

with open(rdir+'bin_HR10_vel_0000.dat.lin', "r") as f:
        veldatA = [map(float, line.split()) for line in f]

with open(rdir + file_prefix + '_gas_' + snapnum +'.dat', "r") as f:
        gasdatB = [map(float, line.split()) for line in f]

with open(rdir + file_prefix + '_vel_' + snapnum + '.dat.lin', "r") as f:
        veldatB = [map(float, line.split()) for line in f]

with open(rdir+'bin_HR10_gas_nh_0000.dat', "r") as f:
        gasdat_nhA = [map(float, line.split()) for line in f]

with open(rdir + file_prefix + '_gas_nh_' + snapnum +'.dat', "r") as f:
        gasdat_nhB = [map(float, line.split()) for line in f]

with open(rdir+'bin_HR10_energy_0000.dat', "r") as f:
        edatA = [map(float, line.split()) for line in f]

with open(rdir + file_prefix + '_energy_' + snapnum +'.dat', "r") as f:
        edatB = [map(float, line.split()) for line in f]


radA = []
nhA = []
densA = []
tempA = []
h2A = []
elecA = []
hdA = []
rmencA = []
csA = []
vradA = []
vrotA = []
turbA = []
mencA = []
presA = []
nh_nhA = []
temp_nhA = []
xmomA = []
ymomA = []
zmomA = []
eintA = []

radB = []
nhB = []
densB = []
tempB = []
h2B = []
elecB = []
hdB = []
rmencB = []
csB = []
vradB = []
vrotB = []
turbB = []
mencB = []
presB = []
nh_nhB = []
temp_nhB = []
xmomB = []
ymomB = []
zmomB = []
eintB = []
gpotB = []
gaccB = []

narr1 = 200
line_arr = []
fac = 1.2195*1.67e-24
for i in range(narr1):
        line_arr = gasdatA[i]
	radA.append(line_arr[0])
	nhA.append(line_arr[1])
        densA.append(line_arr[1]*fac)
	tempA.append(line_arr[2])
        h2A.append(line_arr[3])
        elecA.append(line_arr[4])
        hdA.append(line_arr[5])
	mencA.append(line_arr[6])
	presA.append(1.38e-16*nhA[i]*tempA[i])

narr_vel = 200
for i in range(narr_vel):
        line_arr = veldatA[i]
	rmencA.append(line_arr[0])
        csA.append(line_arr[1])
        turbA.append(line_arr[2])
        vradA.append(line_arr[3])
        vrotA.append(line_arr[4])
        turbA[i] = turbA[i]/csA[i]


narr = 64
for i in range(narr):
	line_arr = gasdatB[i]
	radB.append(line_arr[0])
	nhB.append(line_arr[1])
        densB.append(line_arr[1]*fac)
	tempB.append(line_arr[2])
        h2B.append(line_arr[3])
        elecB.append(line_arr[4])
        hdB.append(line_arr[5])
	mencB.append(line_arr[6])
	presB.append(1.38e-16*nhB[i]*tempB[i])

for i in range(narr_vel):
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

for i in range(narr1):
        line_arr = edatA[i]
        eintA.append(line_arr[0])
        xmomA.append(line_arr[1])
	ymomA.append(line_arr[2])
	zmomA.append(line_arr[3])

for i in range(narr):
        line_arr = edatB[i]
        eintB.append(line_arr[0])
        xmomB.append(line_arr[1])
        ymomB.append(line_arr[2])
        zmomB.append(line_arr[3])
	gpotB.append(line_arr[4])
        if i < 1:
                gaccB.append(0)
        if i >= 1:
                gaccB.append( (gpotB[i] - gpotB[i-1])/DeltaX )


print 'xmomB =', xmomB
print 'ymomB = ', ymomB
##############################################################################

rmin = 100.
rmax = 3.e5
rmax_arr = pf.domain_right_edge[0]
nmin = 1e2
nmax = 1e16
#bin_num1 = 7500
bin_num1 = 200
bin_num2 = 200

i_doub = 1.0
bin_num_square = 64
bin_doub = 64.0

rad_arr = []
xmom_arr = []
ymom_arr = []
zmom_arr = []
eint_arr = []
gpot_arr = []
gacc_arr = []
massenc_arr = []

dfac = 0.99
location = [0.0,0.0,0.0]

for i in range(bin_num_square):
        Lval = -rmax_arr  * (i_doub/bin_doub)
        Rval =  rmax_arr  * (i_doub/bin_doub)
        if i > 0:
                data_in = data

        data = pf.h.region(location, [dfac*Lval, dfac*Lval, dfac*Lval], [dfac*Rval, dfac*Rval, dfac*Rval])
        xmom_arr.append(data["xmomentum"].sum()*Msun)
	ymom_arr.append(data["ymomentum"].sum()*Msun)
	zmom_arr.append(data["zmomentum"].sum()*Msun) 
	eint_arr.append((data["InternalEnergy"]*data["CellMassMsun"]).sum()*Msun)
        massenc_arr.append(data["CellMassGrams"].sum() / 1.98892e33 )
        if i == 0:
                gpot_arr.append( (data['gravitational-potential'] * data['CellMassMsun']).sum() / (massenc_arr[i]) )
                #gpot_arr[i] = np.abs(gpot_arr[i])
        if i > 0:
                data_ann = pf.h.boolean([data, "NOT", data_in])
                gpot_arr.append( (data_ann['gravitational-potential'] * data_ann['CellMassMsun']).sum() / data_ann['CellMassMsun'].sum() )
                #gpot_arr[i] = np.abs(gpot_arr[i])
        if i == 0:
                gacc_arr.append(0)
        if i > 0:
                gacc_arr.append( np.abs(gpot_arr[i] - gpot_arr[i-1]) / DeltaX )
        rad_arr.append(Rval / au_cm )
        i_doub = i_doub + 1.0

#####################################################################################

gpotB_min = 1.e30
gpot_arr_min = 1.e30

for i in range(bin_num_square):
        if gpot_arr[i] < gpot_arr_min:
                gpot_arr_min = gpot_arr[i]
        if gpotB[i] < gpotB_min:
                gpotB_min = gpotB[i]

for i in range(bin_num_square):
        gpot_arr[i] = gpot_arr[i] - gpot_arr_min
        gpotB[i] = gpotB[i] - gpotB_min

print 'gpot_arr =', gpot_arr
print 'gacc_arr =', gacc_arr
print 'eint_arr =', eint_arr
print 'eintB =', eintB
#####################################################################################
#calculate errors
xmom_err = []
ymom_err = []
zmom_err = []
eint_err = []
gpot_err = []
gacc_err = []

for i in range(bin_num_square):
	xmom_err.append(np.abs((xmom_arr[i] - xmomB[i])/xmomB[i]))
	ymom_err.append(np.abs((ymom_arr[i] - ymomB[i])/ymomB[i]))	
	zmom_err.append(np.abs((zmom_arr[i] - zmomB[i])/zmomB[i]))
	eint_err.append(np.abs((eint_arr[i] - eintB[i])/eintB[i]))
        if i == 0:
                gpot_err.append(1.e-10)
                gacc_err.append(1.e-10)
	if i > 0:
        	gpot_err.append(np.abs(gpot_arr[i] - gpotB[i])/gpotB[i])
        	gacc_err.append(np.abs(gacc_arr[i] - gaccB[i])/gaccB[i])

print 'eint_err', eint_err
print 'gacc_err', gacc_err
print 'gpot_err', gpot_err
#####################################################################################

rmin = 1.e3
rmax = 1.e5

nmax=1.e12

mmin = -2.e41
#mmin = 0
mmax = 1.e41

err_min = 1.e-9
err_max = 1.e-1

pl.subplot(321)
pl.plot(rad_arr, xmom_err, 'k')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Radius [AU]', fontsize=9)
ax.set_ylabel('total x-momentum error', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((rmin, rmax, err_min, err_max))

pl.subplot(322)
pl.plot(rad_arr, ymom_err,'k')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Radius [AU]', fontsize=9)
ax.set_ylabel('total y-momentum error', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((rmin, rmax, err_min, err_max))

pl.subplot(323)
pl.plot(rad_arr, zmom_err,'k')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Radius [AU]', fontsize=9)
ax.set_ylabel('total z-momentum error', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((rmin, rmax, err_min, err_max))

pl.subplot(324)
pl.plot(rad_arr, eint_err,'k')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Radius [AU]', fontsize=9)
ax.set_ylabel('internal energy error', fontsize=9)
pl.xticks(fontsize=9)
pl.yticks(fontsize=9)
pl.axis((rmin, rmax, err_min, err_max))

pl.subplot(325)
pl.plot(rad_arr, gpot_err,'k')
ax = pl.gca()
#ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Radius [AU]', fontsize=9)
ax.set_ylabel('grav. pot. error', fontsize=9)
pl.xticks(fontsize=9)
pl.yticks(fontsize=9)
pl.axis((rmin, rmax, 1.e-3, 1.e0))

pl.subplot(326)
pl.plot(rad_arr, gacc_err,'k')
ax = pl.gca()
#ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Radius [AU]', fontsize=9)
ax.set_ylabel('grav. acc. error', fontsize=9)
pl.xticks(fontsize=9)
pl.yticks(fontsize=9)
pl.axis((rmin, rmax, 1.e-3, 1.e0))

pl.savefig('plot.ps')

M_i = data["CellMassMsun"].sum()
print "total mass = ", M_i
