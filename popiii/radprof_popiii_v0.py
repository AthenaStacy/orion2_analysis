import matplotlib
matplotlib.use('ps')

from yt.mods import *
import numpy as na
import pylab as pl
from string import rstrip
import fields
#import matplotlib.colorbar as cb
#import fields_iso

MRH = 0.76

#datanum = '0023'
datanum = '0024'

pf = load("/global/scratch/minerva/popiii_Bscope/" + 'data.' + datanum + ".3d.hdf5")
value, location = pf.h.find_max("Density")
data = pf.h.sphere(location, 3.0/pf['pc'])
print 'location =', location
pc = PlotCollection(pf, center = location)

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
    rad =  na.sqrt((data['x'] - x)**2 + (data['y'] - y)**2 + (data['z'] - z)**2)
    rad = rad/1.5e13
    return rad
add_field("Radius", function=_Radius, take_log=False,
          units=r'\rm{AU}')

def _mu(field,data):
   h2frac = data["tracer6"]
   hefrac = data["tracer15"]  #should be approximately 0.24
   hfrac = 1. - h2frac - hefrac
   #mu_inv = (1. - h2frac)*0.76 + h2frac*0.76/2.0 + 0.24/4.0 #here h2frac ranges from 0-1
   mu_inv = hfrac + h2frac/2. + hefrac/4.
   mu = 1. / mu_inv
   return  (mu )
add_field("mu",function=_mu,units=r"\rm{Kelvin}",take_log=True)

def _Temp_alt(field,data):
   return  (data["ThermalEnergy"] * (data['tracer1']-1.0) / data["density"] * data['mu'] * 1.67e-24 / 1.38e-16)
add_field("Temp_alt",function=_Temp_alt,units=r"\rm{Kelvin}",take_log=True)

def _nh_alt(field,data):
   return  (data["density"] / data['mu'] / 1.67e-24 )
add_field("nh_alt",function=_nh_alt,units=r"cm^{-3}",take_log=True)
#######################################################################################################

rad_min = min(data['Radius'])
rad_max = max(data['Radius'])
rad_arr = []
print 'rad_arr length = ', len(rad_arr)



print 'rad_min=', rad_min, 'rad_max=', rad_max
print 'mu_min=', min(data['mu']), 'mu_max=', max(data['mu'])
print 'h2_min=', min(data['tracer6']), 'h2_max=', max(data['tracer6'])
print 'temp_min=', min(data['Temp_alt']), 'temp_max=', max(data['Temp_alt'])
print 'vrad_min=', min(data['radial-velocity']), 'vrad_max=', max(data['radial-velocity'])
print 'dens_min=', min(data['Density']), 'dens_max=', max(data['Density'])

bin_num1 = 5000
bin_num2 = 200
rmin = 0.
rmin_menc = 0.
rmax = 3.e5 
#nmin = 1.e4
#nmax = 6.e11
nmin = 1.e3
nmax = 1.e12

#location = [0,0,0]

p1 = pc.add_profile_sphere(2., "pc", ['Radius', 'nh_alt'], weight='CellMassMsun',x_bins=bin_num1, x_log=False, x_bounds = [rmin,rmax], center = location)
p2 = pc.add_profile_sphere(2., "pc", ['Radius', 'radial-velocity'], weight='CellMassMsun', x_bins=bin_num1, x_log=False, x_bounds = [rmin,rmax], center = location)
p3 = pc.add_profile_sphere(2., "pc", ['nh_alt', 'Temp_alt', 'tracer6'], weight='CellMassMsun', x_bins = bin_num2, center = location, x_bounds = [nmin, nmax])
p4 = pc.add_profile_sphere(2., "pc", ['Radius', 'CellMassMsun'], weight=None, accumulation=True, x_bins=bin_num1, x_log=False, x_bounds = [rmin_menc,rmax], center = location)
# pc.add_profile_sphere(1., "pc", ['Radius', 'CellMassMsun'], weight=None, accumulation=True, x_bins=bin_num1, x_log=False, x_bounds = [rmin,rmax], center = location)

pc.save('radprof')

################read in other data files#######################################
dir = '/global/home/users/minerva/gadget_runs/'
basename = 'bin_zoom10_new_cut_ref3'
datanum1 = '7053'
print 'snapnum = ', datanum1

gasdatA = []
gasdatB = []
nhdatA = []
nhdatB = []

with open(dir + basename + '_gas_' + datanum1 + '.dat', "r") as f:
	gasdatA = [map(float, line.split()) for line in f]

with open(dir + basename + '_gas_' + datanum1 + '.dat', "r") as f:
        gasdatB = [map(float, line.split()) for line in f]

with open(dir + basename + '_nhprof_' + datanum1 + '.dat', "r") as f:
        nhdatA = [map(float, line.split()) for line in f]

with open(dir + basename + '_nhprof_' + datanum1 + '.dat', "r") as f:
        nhdatB = [map(float, line.split()) for line in f]


radA = []
nhA = []
h2A = []
tempA = []
rmencA = []
vradA = []
vrotA = []
mencA = []
presA = []
nhA2 = []
tempA2 = []
h2A2 = []

radB = []
nhB = []
h2B = []
tempB = []
rmencB = []
vradB = []
vrotB = []
mencB = []
presB = []
nhB2 = []
tempB2 = []
h2B2 = []

narr = 5000
narr2 = 200
line_arr = []
for i in range(narr):
        line_arr = gasdatA[i]
	radA.append(line_arr[0])
	mencA.append(line_arr[1])
	nhA.append(line_arr[2])
        h2A.append(line_arr[3])
        tempA.append(line_arr[4])
        vradA.append(line_arr[5])
        vrotA.append(line_arr[6])

for i in range(narr):
	line_arr = gasdatB[i]
	radB.append(line_arr[0])
        mencB.append(line_arr[1])
        nhB.append(line_arr[2])
        h2B.append(line_arr[3])
        tempB.append(line_arr[4])
        vradB.append(line_arr[5])
        vrotB.append(line_arr[6])

for i in range(narr2):
        line_arr = nhdatA[i]
        nhA2.append(line_arr[0])
        h2A2.append(line_arr[1])
        tempA2.append(line_arr[6])

for i in range(narr2):
        line_arr = nhdatB[i]
        nhB2.append(line_arr[0])
        h2B2.append(line_arr[1])
        tempB2.append(line_arr[6])

for i in range(narr2):
	presA.append(1.36e-18*nhA[i]*tempA[i])
	presB.append(1.36e-18*nhB[i]*tempB[i])      
	h2A2[i] = h2A2[i]*2.*MRH 
        h2B2[i] = h2B2[i]*2.*MRH 
##############################################################################

rad = p1.data['Radius']
#nh = p1.data['number-density'] 
nh = p1.data['nh_alt']
rad2 = p2.data['Radius']
vrad = p2.data['radial-velocity'] 
#nh3 = p3.data['number-density']
nh3 = p3.data['nh_alt']
temp = p3.data['Temp_alt']
h2 = p3.data['tracer6']
massenc_arr = p4.data['CellMassMsun']
rmenc = p4.data['Radius']

print 'radA =', radA
print 'rad =', rad

print 'nhA = ', nhA
print 'nh =', nh

print 'tempA = ', tempA
print 'temp =', temp

pl.subplot(221)
pl.plot(rad, nh, 'k')
pl.plot(radA, nhA)
#pl.plot(radB, nhB, 'r--')
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
pl.plot(nhA2, tempA2)
#pl.plot(nhB2, tempB2, 'r--')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('number density [cm^{-3}]', fontsize=9)
ax.set_ylabel('Temperature [K]', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((1e2, 1e12, 1e2, 1e4))

pl.subplot(223)
pl.plot(nh3, h2,'k')
pl.plot(nhA2, h2A2)
#pl.plot(nhB2, h2B2, 'r--')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('number density [cm^{-3}]', fontsize=9)
ax.set_ylabel('H2 ', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((1e2, 1e12, 1e-6, 1e0))

#pl.subplot(223)
#pl.plot(rad2, vrad, 'k')
#pl.plot(radA, vradA)
#pl.plot(radB, vradB, 'r--')
#ax = pl.gca()
#ax.set_xscale('log')
#ax.set_xlabel('Radius [AU]', fontsize=9)
#ax.set_ylabel('radial-velocity [km s^{-1}]', fontsize=9)
#pl.xticks(fontsize=10)
#pl.yticks(fontsize=10)
#pl.axis((1e1, 1e6, -6, 6))

pl.subplot(224)
pl.plot(rmenc, massenc_arr,'k')
pl.plot(radA, mencA)
#pl.plot(radB, mencB, 'r--')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Radius [AU]', fontsize=9)
ax.set_ylabel('Enclosed Mass [M_{\odot}', fontsize=9)
pl.xticks(fontsize=9)
pl.yticks(fontsize=9)
pl.axis((1e1, 1e6, 1.1, 3e3))

pl.savefig('radprof.ps')

M_i = data["CellMassMsun"].sum()
print "total mass = ", M_i
