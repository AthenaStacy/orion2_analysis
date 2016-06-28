import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as pl
#import pylab as pl

from yt.mods import *
import numpy as na
from string import rstrip
#import fields
import fields_iso

wdir = '/work/00863/minerva/orion/maptest_isosphere/'
#wdir = '/work/00863/minerva/orion/'
rdir = '/home1/00863/minerva/research_programs/'
snapnum = '001'
onum = '0000'
setup = 3

#pf = load("/nobackupp7/astacy/orion/chemtest_iso/data.0021.3d.hdf5")
pf = load(wdir+"data."+onum+".3d.hdf5")
#value, location = pf.h.find_max("Density")
value, location = pf.h.find_max("CellMassMsun")
data = pf.h.sphere(location, 1.e4/pf['pc'])
pc = PlotCollection(pf, center = location)

print pf.h.field_list

x=location[0]
y=location[1]
z=location[2]

massenc_com = []
vx = []
vy = []
vz = []
sub_sphere =  pf.h.sphere(location, 6.e19/pf['cm'])
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

#print 'pres_min=', min(data['Pressure']), 'pres_max=', max(data['Pressure'])
print 'temp_min=', min(data['Temperature']), 'temp_max=', max(data['Temperature'])
print 'vrad_min=', min(data['radial-velocity']), 'vrad_max=', max(data['radial-velocity'])
print 'dens_min=', min(data['density']), 'dens_max=', max(data['density'])
print 'rad_min', min(data['Radius']), 'rad_max', max(data['Radius'])

################read in other data files#######################################

with open(rdir+"iso_map_gas_002.dat", "r") as f:
	gasdatA = [map(float, line.split()) for line in f]

with open(rdir+'iso_map_turb_002.dat', "r") as f:
        turbdatA = [map(float, line.split()) for line in f]

with open(rdir+'iso_map_gas_'+snapnum+'.dat', "r") as f:
        gasdatB = [map(float, line.split()) for line in f]

with open(rdir+'iso_map_turb_'+snapnum+'.dat', "r") as f:
        turbdatB = [map(float, line.split()) for line in f]

with open(rdir+'iso_map_gas_nh_002.dat', "r") as f:
        gasdat_nhA = [map(float, line.split()) for line in f]

with open(rdir+'iso_map_gas_nh_'+snapnum+'.dat', "r") as f:
        gasdat_nhB = [map(float, line.split()) for line in f]

radA = []
nhA = []
densA = []
tempA = []
rmencA = []
vradA = []
vrotA = []
mencA = []
presA = []
nh_nhA = []
temp_nhA = []

radB = []
nhB = []
densB = []
tempB = []
rmencB = []
vradB = []
vrotB = []
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

for i in range(narr1):
        line_arr = turbdatA[i]
	rmencA.append(line_arr[0])
#        if line_arr[2] < 0:
#        	line_arr[2] = line_arr[2]*(-1)
	vradA.append(line_arr[2])
	vrotA.append(line_arr[3])
	mencA.append(line_arr[5])

narr = 200
for i in range(narr):
	line_arr = gasdatB[i]
	radB.append(line_arr[0])
	nhB.append(line_arr[1])
        densB.append(line_arr[1]*fac)
	tempB.append(line_arr[2])

for i in range(narr):
	line_arr = turbdatB[i]
	rmencB.append(line_arr[0])
#        if line_arr[2] < 0:
#                line_arr[2] = line_arr[2]*(-1)
	vradB.append(line_arr[2])
	vrotB.append(line_arr[3])
	mencB.append(line_arr[5])

for i in range(narr1):
        line_arr = gasdat_nhA[i]
        nh_nhA.append(line_arr[0])
        temp_nhA.append(line_arr[1])

for i in range(narr1):
        line_arr = gasdat_nhB[i]
        nh_nhB.append(line_arr[0])
        temp_nhB.append(line_arr[1])

for i in range(narr1):
	presA.append(1.38e-16*nhA[i]*tempA[i])
for i in range(narr):
	presB.append(1.38e-16*nhB[i]*tempB[i])      

print 'radA =', radA
print 'vradA = ', vradA
##############################################################################
##############################################################################
#read in sink mass

mass_sink_arr = []

f=open(wdir+"data."+onum+".3d.sink")
f.readline()
mass_sink_arr.append(float(f.readline().split(' ')[0])/2.e33)

mass_sink = mass_sink_arr[0]
###################################################################################

rmin = 1e2
rmax = 1e6
nmin = 1e-2
nmax = 1e10
bin_num1 = 500
bin_num2 = 100

p1 = pc.add_profile_sphere(1000., "pc", ['Radius', 'density', 'number-density', 'radial-velocity', 'rot-velocity', 'Temperature', 'Pressure'], weight='CellMassMsun', x_bins=bin_num1, x_log=True, x_bounds = [rmin,rmax], center = location)

p1b = pc.add_profile_sphere(1000., "pc", ['Radius','radial-velocity','rot-velocity'], weight='CellMassMsun', x_bins=bin_num2, x_log=True, x_bounds = [rmin,rmax], center = location)

p1c = pc.add_profile_sphere(1000., "pc", ['Radius', 'Pressure'], weight='CellMassMsun', x_bins=bin_num1, x_log=False, x_bounds = [rmin,rmax], center = location)

p2 = pc.add_profile_sphere(1000., "pc", ['number-density', 'Temperature'], weight='CellMassMsun')

p3 = pc.add_profile_sphere(1000., "pc", ['Radius', 'CellMassMsun'], weight=None, accumulation=True, x_bins=bin_num1, x_log=True, x_bounds = [rmin,rmax], center = location)

pc.save('radprof')

rad = p1.data['Radius']
dens = p1.data['density'] 
nh = p1.data['number-density']
temp = p1.data['Temperature']
vrada = p1.data['radial-velocity']
vrota = p1.data['rot-velocity']

radb = p1b.data['Radius']
vrad = p1b.data['radial-velocity'] 
vrot = p1b.data['rot-velocity']

radc = p1c.data['Radius']
pres = p1c.data['Pressure']

nh2 = p2.data['number-density']
temp2 = p2.data['Temperature']

rad_arr = p3.data['Radius']
massenc_arr = p3.data['CellMassMsun']

print 'rad =', rad
print 'nh=', nh
print 'mass_enc=', massenc_arr

################################################################################################
#print/write sim data to file!
masstot_arr = []

f = open("maptest_isosphere_"+onum+".dat", "w")
#with open('maptest_isosphere_'+onum+'.dat', 'w') as f:
#f.write("Oh hi\n");
#f.write("How's it going\n");
for i in range(bin_num1):
	masstot_arr.append(massenc_arr[i] + mass_sink)
	if nh[i] > 0:
		f.write(str(rad[i]) + ' ' + str(nh[i]) + ' ' + str(vrada[i]) + ' ' +  str(vrota[i]) + ' ' + str(masstot_arr[i]) + '\n')
f.close()

################################################################################################

#narr = 200
#dens = []
#nh = []
#rad_nh = []
#for i in range(narr - 1):
#        rad_out = 1.5e13*rad_arr[i+1]
#	rad_in = 1.5e13*rad_arr[i]
#        rad_nh.append(rad_out/1.5e13)
#	volume = (4./3.)*3.14159*(rad_out**3 - rad_in**3) 
#	dens.append(2.e33*(massenc_arr[i+1] - massenc_arr[i])/volume)
#	nh.append(dens[i-1]/fac)


#print 'radA =', radA
#print 'vradA = ', vradA
#print 'nhA =', nhA
#print 'tempA =', tempA
#print 'rad =', rad
#print 'nh =', nh

rmin = 1.e2
rmax = 1.e8
nmin = 1.e-2
nmax = 1.e12

pl.subplot(221)
#pl.plot(rad_nh, nh, 'k')
pl.plot(rad, nh, 'k')
pl.plot(radA, nhA)
pl.plot(radB, nhB, 'r--')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Radius [AU]', fontsize=9)
ax.set_ylabel('number density [cm^{-3}]', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((rmin, rmax, nmin, nmax))

if setup== 3:
        pl.subplot(222)
        pl.plot(rad, temp,'k')
        pl.plot(radA, tempA)
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
	pl.plot(nh_nhA, temp_nhA)
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
        pl.plot(radb, vrad, 'k')
        pl.plot(rmencA, vradA)
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
	pl.plot(radb, vrot, 'k')
	pl.plot(rmencA, vrotA)
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
	pl.plot(radA, presA)
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
	pl.plot(rmencA, mencA)
	pl.plot(rmencB, mencB, 'r--')
	ax = pl.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel('Radius [AU]', fontsize=9)
	ax.set_ylabel('Enclosed Mass [M_{\odot}', fontsize=9)
	pl.xticks(fontsize=9)
	pl.yticks(fontsize=9)
	pl.axis((1e1, 1e8, 1.1, 2e6))

pl.savefig('radprof_com.pdf')

M_i = data["CellMassMsun"].sum()
print "total mass = ", M_i
