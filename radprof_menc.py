from yt.mods import *
import numpy as na
import pylab as pl
from string import rstrip
import fields
#import fields_iso

#import matplotlib
#matplot.use('pdf')
#import matplot.pyplot as pl
#import matplotlib pyplot at pl

#pf = load("/nobackupp7/astacy/orion/chemtest_iso/data.0021.3d.hdf5")
pf = load("/nobackupp7/astacy/orion/data.0000.3d.hdf5")
value, location = pf.h.find_max("Density")
data = pf.h.sphere(location, 3.0/pf['pc'])
print 'location =', location
pc = PlotCollection(pf, center = location)

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


rad_min = min(data['Radius'])
rad_max = max(data['Radius'])
rad_arr = []
massenc_arr = []
print 'rad_arr length = ', len(rad_arr)
narr = 100
narr_doub = 100.0
i_doub = 1.0
for i in range(narr):
	rad_val = 0.5*rad_min + rad_min*((rad_max/rad_min)**(i_doub/narr_doub))
	rad_arr.append(rad_val)    
	i_doub = i_doub + 1.0
	sub_sphere =  pf.h.sphere(location, 1.5e13*rad_arr[i]/pf['cm']) 
	massenc_arr.append(sub_sphere["CellMassMsun"].sum())
	print 'rad_arr =', rad_arr[i]
	print 'massenc =', massenc_arr[i]


#def _massenc(field, data):
#    '''
#    Mass enclosed within given distance. In this problem the center is at density peak.
#    '''
#    rad =  na.sqrt((data['x'] )**2 + (data['y'])**2 + (data['z'])**2)
#    #print "rad[0][0][0] =", rad[0][0][0]
#    leni = len(rad)
#    lenj = len(rad[0])
#    lenk = len(rad[0][0])
#    massenc = na.zeros(leni*lenj*lenk).reshape((leni, lenj, lenk))
#    #print "leni = ", leni, "lenj = ", lenj, "lenk =", lenk 
#    for i in range(leni):
#	for j in range(lenj):
#		for k in range(lenk):	
#			for l in range(narr):
#    				if(rad[i][j][k] < rad_arr[l]):
#					massenc[i][j][k] = massenc_arr[l]
#
#    return massenc
#add_field("enclosed-mass", function=_massenc, take_log=False,
#          units=r'\rm{M_{\odot}}')


print 'rad_min=', rad_min, 'rad_max=', rad_max
#print 'pres_min=', min(data['Pressure']), 'pres_max=', max(data['Pressure'])
print 'temp_min=', min(data['Temperature']), 'temp_max=', max(data['Temperature'])
print 'vrad_min=', min(data['radial-velocity']), 'vrad_max=', max(data['radial-velocity'])
print 'dens_min=', min(data['Density']), 'dens_max=', max(data['Density'])
#print 'enclosed-mass-min', min(data['enclosed-mass']), 'enclosed-mass-max', max(data['enclosed-mass'])

p1 = pc.add_profile_sphere(2., "pc", ['Radius', 'number-density'], weight='CellVolume')
p2 = pc.add_profile_sphere(2., "pc", ['Radius', 'radial-velocity'], weight='CellVolume')

p3 = pc.add_profile_sphere(2., "pc", ['number-density', 'Temperature'], weight='CellVolume')
#p4 = pc.add_profile_sphere(2., "pc", ['Radius', 'enclosed-mass'], weight='CellVolume')

pc.save('radprof')

################read in other data files#######################################
gasdatA = []
turbdatA = []
gasdatB = []
turbdatB = []

with open('/home1/astacy/research_programs/bin_HR10_map_gas.dat', "r") as f:
	gasdatA = [map(float, line.split()) for line in f]

with open('/home1/astacy/research_programs/bin_HR10_map_turb.dat', "r") as f:
        turbdatA = [map(float, line.split()) for line in f]

with open('/home1/astacy/research_programs/bin_HR10_gas.dat', "r") as f:
        gasdatB = [map(float, line.split()) for line in f]

with open('/home1/astacy/research_programs/bin_HR10_turb.dat', "r") as f:
        turbdatB = [map(float, line.split()) for line in f]

radA = []
nhA = []
tempA = []
rmencA = []
vradA = []
mencA = []
presA = []

radB = []
nhB = []
tempB = []
rmencB = []
vradB = []
mencB = []
presB = []

narr = 200
line_arr = []
for i in range(narr):
        line_arr = gasdatA[i]
	radA.append(line_arr[0])
	nhA.append(line_arr[1])
	tempA.append(line_arr[2])

for i in range(narr):
        line_arr = turbdatA[i]
	rmencA.append(line_arr[0])
#        if line_arr[2] < 0:
#        	line_arr[2] = line_arr[2]*(-1)
	vradA.append(line_arr[2])
	mencA.append(line_arr[5])


for i in range(narr):
	line_arr = gasdatB[i]
	radB.append(line_arr[0])
	nhB.append(line_arr[1])
	tempB.append(line_arr[2])

for i in range(narr):
	line_arr = turbdatB[i]
	rmencB.append(line_arr[0])
#        if line_arr[2] < 0:
#                line_arr[2] = line_arr[2]*(-1)
	vradB.append(line_arr[2])
	mencB.append(line_arr[5])

for i in range(narr):
	presA.append(1.36e-18*nhA[i]*tempA[i])
	presB.append(1.36e-18*nhB[i]*tempB[i])      

print 'radA =', radA
print 'vradA = ', vradA
##############################################################################

rad = p1.data['Radius']
nh = p1.data['number-density'] 
rad2 = p2.data['Radius']
vrad = p2.data['radial-velocity'] 
nh3 = p3.data['number-density']
temp = p3.data['Temperature']


print 'radA =', radA
print 'vradA = ', vradA
print 'rad =', rad
print 'nh =', nh
print 'rad2 =', rad2
print 'vrad = ', vrad

pl.subplot(221)
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
pl.axis((1e1, 1e6, 1e2, 1e12))

pl.subplot(222)
pl.plot(nh3, temp,'k')
pl.plot(nhA, tempA)
pl.plot(nhB, tempB, 'r--')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('number density [cm^{-3}]', fontsize=9)
ax.set_ylabel('Temperature [K]', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((1e2, 1e12, 1e2, 1e4))

pl.subplot(223)
pl.plot(rad2, vrad, 'k')
pl.plot(rmencA, vradA)
pl.plot(rmencB, vradB, 'r--')
ax = pl.gca()
ax.set_xscale('log')
ax.set_xlabel('Radius [AU]', fontsize=9)
ax.set_ylabel('radial-velocity [km s^{-1}]', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((1e1, 1e6, -6, 6))

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
pl.axis((1e1, 1e6, 1.1, 3e3))

pl.savefig('radprof_com.png')

M_i = data["CellMassMsun"].sum()
print "total mass = ", M_i
