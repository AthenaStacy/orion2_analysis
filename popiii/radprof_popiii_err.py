import matplotlib
matplotlib.use('ps')

from yt.mods import *
import numpy as na
import pylab as pl
from string import rstrip
import fields
#import matplotlib.colorbar as cb
#import fields_iso

igamma = 'tracer1'
iHP =   'tracer2'
iH =   'tracer3'
iHM =  'tracer4'
iH2P =   'tracer5'
iH2 =   'tracer6'
iDP =   'tracer7'
iD =   'tracer8'
iDM =  'tracer9'
iHDP   =  'tracer10'
iHD   =  'tracer11'
iD2P  =  'tracer12'
iD2 =  'tracer13'
iHEP =  'tracer14'
iHE  =  'tracer15'
iHEPP  =  'tracer16'
iELEC = 'tracer17'

MRH = 0.76
h2conv = (2.0/1.3158)
heconv = (4.0/1.3158)

datanum = '0134'
#datanum = '0023'
#datanum = '0024'

pf = load("/global/scratch/minerva/popiii_Bscope5/" + 'data.' + datanum + ".3d.hdf5")
#pf = load("/global/scratch/minerva/popiii_Bscope4/" + 'data.' + datanum + ".3d.hdf5")
#pf = load("/global/scratch/minerva/popiii_Bscope3/" + 'data.' + datanum + ".3d.hdf5")
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
   h2frac = data[iH2] * h2conv #/data['density']
   hefrac = data[iHE] * heconv #/data['density']  #should be approximately 0.24
   hfrac = 1. - h2frac - hefrac
   #mu_inv = (1. - h2frac)*0.76 + h2frac*0.76/2.0 + 0.24/4.0 #here h2frac ranges from 0-1
   mu_inv = hfrac + h2frac/2. + hefrac/4.
   mu = 1. / mu_inv
   return  (mu )
add_field("mu",function=_mu,units=r"\rm{Kelvin}",take_log=True)

def _Temp_alt(field,data):
   return  (data["ThermalEnergy"] * (data[igamma]-1.0) / data["density"] * data['mu'] * 1.67e-24 / 1.38e-16)
add_field("Temp_alt",function=_Temp_alt,units=r"\rm{Kelvin}",take_log=True)

def _nh_alt(field,data):
   return  (data["density"] * MRH / 1.67e-24 )
add_field("nh_alt",function=_nh_alt,units=r"cm^{-3}",take_log=True)

#######################################################################################################

rad_min = min(data['Radius'])
rad_max = max(data['Radius'])
rad_arr = []
print 'rad_arr length = ', len(rad_arr)


print 'rad_min=', rad_min, 'rad_max=', rad_max
print 'mu_min=', min(data['mu']), 'mu_max=', max(data['mu'])
print 'gam_min=', min(data[igamma]), 'gam_max=', max(data[igamma])
print 'hI_min=', min(data[iH]/data['density']), 'hI_max=', max(data[iH]/data['density'])
print 'h2_min=', min(data[iH2]/data['density']), 'h2_max=', max(data[iH2]/data['density'])
print 'hm_min=', min(data[iHM]/data['density']), 'hm_max=', max(data[iHM]/data['density'])
print 'h2p_min=', min(data[iH2P]/data['density']), 'h2p_max=', max(data[iH2P]/data['density'])
print 'dm_min=', min(data[iDM]/data['density']), 'dm_max=', max(data[iDM]/data['density'])
print 'd2_min=', min(data[iD2]/data['density']), 'd2_max=', max(data[iD2]/data['density'])
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
#nmin = 1.e3
#nmax = 1.e12
nmin = 1.e6
nmax = 1.e14

#location = [0,0,0]

p3 = pc.add_profile_sphere(2., "pc", ['nh_alt', 'density', 'Temp_alt', iH2, iHP, iDP, iHD, iHEP], weight='CellMassMsun', x_bins = bin_num2, center = location, x_bounds = [nmin, nmax])

pc.save('radprof')

################read in other data files#######################################
dir = '/global/home/users/minerva/gadget_runs/'
basename = 'bin_zoom10_new_cut_ref3'
datanum1 = '7054'
#datanum1 = '7052'
#datanum1 = '7054'
#datanum1 = '7055'
#datanum1 = '7057'
datanum2 = '7053'
print 'snapnum = ', datanum1

gasdatA = []
gasdatB = []
nhdatA = []
nhdatB = []

with open(dir + basename + '_gas_' + datanum1 + '.dat', "r") as f:
	gasdatA = [map(float, line.split()) for line in f]

with open(dir + basename + '_gas_' + datanum2 + '.dat', "r") as f:
        gasdatB = [map(float, line.split()) for line in f]

with open(dir + basename + '_nhprof_' + datanum1 + '.dat', "r") as f:
        nhdatA = [map(float, line.split()) for line in f]

#with open(dir + basename + '_nhprof_' + datanum2 + '.dat', "r") as f:
with open('/global/home/users/minerva/StandAlone/output.dat', "r") as f:
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
diiA2 = []
hiiA2 = []
hdA2 = []
heiiA2 = []

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
diiB2 = []
hiiB2 = []
hdB2 = []
heiiB2 = []

narr = 5000
narr2 = 200
narr3 = 1379
#narr3 = 200
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
        hdA2.append(line_arr[2])
        diiA2.append(line_arr[3])
        heiiA2.append(line_arr[4])
        hiiA2.append(line_arr[5])
        tempA2.append(line_arr[10])

for i in range(narr3):
        line_arr = nhdatB[i]
        nhB2.append(line_arr[0])
        h2B2.append(line_arr[1])
        hdB2.append(line_arr[2])
        diiB2.append(line_arr[3])
        heiiB2.append(line_arr[4])
        hiiB2.append(line_arr[5])
        tempB2.append(line_arr[10])

print 'tempB2 =', tempB2

for i in range(narr2):
	presA.append(1.36e-18*nhA[i]*tempA[i])
	presB.append(1.36e-18*nhB[i]*tempB[i])      
# 	h2A2[i] = h2A2[i] * 2. / 1.3158 
#        hdA2[i] = hdA2[i] * 3. / 1.3158
#        hiiA2[i] = hiiA2[i] / 1.3158

##############################################################################

#rad = p1.data['Radius']
#nh = p1.data['number-density'] 
#nh = p1.data['nh_alt']
#rad2 = p2.data['Radius']
#vrad = p2.data['radial-velocity'] 
#nh3 = p3.data['number-density']
den = p3.data['density']
nh3 = p3.data['nh_alt']
temp = p3.data['Temp_alt']

#dfac = den
dfac = 1.0

h2 = p3.data[iH2] / dfac
hii = p3.data[iHP] /dfac
dii = p3.data[iDP] /dfac
hd = p3.data[iHD] /dfac
heii = p3.data[iHEP] /dfac

print 'nh3 =', nh3
print 'nhA2 =', nhA2
###############################################################################
#calculate error

temp_err = []
h2_err = []
hd_err = []
dii_err = []
hii_err = []
temp_errB = []
h2_errB = []
hd_errB = []
dii_errB = []
hii_errB = []
dii_ratio = []

for i in range(narr2):
	temp_err.append(np.fabs((temp[i] - tempA2[i])/(temp[i] + tempA2[i])*2))
	h2_err.append(np.fabs((h2[i] - h2A2[i])/(h2[i] + h2A2[i])*2))
        hd_err.append(np.fabs((hd[i] - hdA2[i])/(hd[i] + hdA2[i])*2))
        dii_err.append(np.fabs((dii[i] - diiA2[i])/(dii[i] + diiA2[i])*2))
        hii_err.append(np.fabs((hii[i] - hiiA2[i])/(hii[i] + hiiA2[i])*2))

#for i in range(narr3):
#        temp_errB.append(np.fabs((temp[i] - tempB2[i])/(temp[i] + tempB2[i])*2))
#        h2_errB.append(np.fabs((h2[i] - h2B2[i])/(h2[i] + h2B2[i])*2))
#        hd_errB.append(np.fabs((hd[i] - hdB2[i])/(hd[i] + hdB2[i])*2))
#        dii_errB.append(np.fabs((dii[i] - diiB2[i])/(dii[i] + diiB2[i])*2))
#        hii_errB.append(np.fabs((hii[i] - hiiB2[i])/(hii[i] + hiiB2[i])*2))
#        dii_ratio.append(diiB2[i]/hiiB2[i])

##############################################################################
print 'heii =', heii
#print 'heiiB2 =', heiiB2
print 'dii =', dii
#print 'diiB2 =', diiB2
print 'hii_err =', hii_err

nmin = 1.e6
nmax = 1.e14
emin = 1.e-2
emax = 2.e0

pl.subplot(421)
pl.plot(nh3, temp,'k')
pl.plot(nhA2, tempA2)
pl.plot(nhB2, tempB2, 'b:')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('number density [cm^{-3}]', fontsize=9)
ax.set_ylabel('Temperature [K]', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((nmin, nmax, 1e2, 1e4))

pl.subplot(422)
pl.plot(nhA2, temp_err)
#pl.plot(nhB2, temp_errB, 'b:')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('number density [cm^{-3}]', fontsize=9)
ax.set_ylabel('Temperature error', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((nmin, nmax, emin, emax))

pl.subplot(423)
pl.plot(nh3, h2,'k')
pl.plot(nhA2, h2A2)
pl.plot(nhB2, h2B2, 'b:')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('number density [cm^{-3}]', fontsize=9)
ax.set_ylabel('H2 ', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((nmin, nmax, 1e-4, 1e0))


pl.subplot(424)
pl.plot(nhA2, h2_err)
#pl.plot(nhB2, h2_errB, 'b:')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('number density [cm^{-3}]', fontsize=9)
ax.set_ylabel('H2 error', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((nmin, nmax, emin, emax))

pl.subplot(425)
pl.plot(nh3, hd,'k')
pl.plot(nhA2, hdA2)
pl.plot(nhB2, hdB2, 'b:')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('number density [cm^{-3}]', fontsize=9)
ax.set_ylabel('HD ', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((nmin, nmax, 1e-8, 1e-3))


pl.subplot(426)
pl.plot(nhA2, hd_err)
#pl.plot(nhB2, hd_errB, 'b:')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('number density [cm^{-3}]', fontsize=9)
ax.set_ylabel('HD error', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((nmin, nmax, emin, emax))

#pl.subplot(527)
#pl.plot(nh3, dii,'k')
#pl.plot(nhA2, diiA2)
#pl.plot(nhB2, diiB2, 'b:')
#ax = pl.gca()
#ax.set_xscale('log')
#ax.set_yscale('log')
#ax.set_xlabel('number density [cm^{-3}]', fontsize=9)
#ax.set_ylabel('DII ', fontsize=9)
#pl.xticks(fontsize=10)
#pl.yticks(fontsize=10)
#pl.axis((nmin, nmax, 1e-10, 1e-5))


#pl.subplot(528)
#pl.plot(nhA2, dii_err)
#pl.plot(nhB2, dii_errB, 'b:')
#ax = pl.gca()
#ax.set_xscale('log')
#ax.set_yscale('log')
#ax.set_xlabel('number density [cm^{-3}]', fontsize=9)
#ax.set_ylabel('DII error', fontsize=9)
#pl.xticks(fontsize=10)
#pl.yticks(fontsize=10)
#pl.axis((nmin, nmax, emin, emax))


pl.subplot(427)
pl.plot(nh3, hii,'k')
pl.plot(nhA2, hiiA2)
pl.plot(nhB2, hiiB2, 'b:')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('number density [cm^{-3}]', fontsize=9)
ax.set_ylabel('HII ', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((nmin, nmax, 1e-12, 1e-6))


pl.subplot(4,2,8)
pl.plot(nhA2, hii_err)
#pl.plot(nhB2, hii_errB, 'b:')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('number density [cm^{-3}]', fontsize=9)
ax.set_ylabel('HII error', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((nmin, nmax, emin, emax))

pl.show()
pl.savefig('radprof_err.ps')

M_i = data["CellMassMsun"].sum()
print "total mass = ", M_i
