import matplotlib
matplotlib.use('ps')
import yt
from yt.units import dimensions
from astropy import units as u
#from yt.mods import *
import numpy as na
import pylab as pl
from string import rstrip
#import fields
import fields_bfield
#import matplotlib.colorbar as cb
#import fields_iso
from tracer_def import *

MRH = 0.76

#datanum = '0137'
#datanum = '0023'
datanum = '0131'
#datanum = '0240'

pf_pc = 3e-18

pf = load("/nobackupp7/astacy/popiii_chemtest/" + 'data.' + datanum + ".3d.hdf5")
#pf = load("/nobackupp7/astacy/popiii_bfieldA0/" + 'data.' + datanum + ".3d.hdf5")
value, location = pf.h.find_max("density")

data = pf.h.sphere(location, 3.0/pf_pc)
print 'location =', location
#pc = PlotCollection(pf, center = location)

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


def _nh_alt(field,data):
   return  (data["density"] * MRH / 1.67e-24 )
add_field("nh_alt",function=_nh_alt,units=r"cm^{-3}",take_log=True)

#######################################################################################################
'''
print 'rad_min=', min(data['Radius']), 'rad_max=', max(data['Radius'])
print 'mu_min=', min(data['mu']), 'mu_max=', max(data['mu'])
print 'gam_min=', min(data[igamma]), 'gam_max=', max(data[igamma])
print 'hI_min=', min(data[iH]/data['density']), 'hI_max=', max(data[iH]/data['density'])
print 'h2_min=', min(data[iH2]/data['density']), 'h2_max=', max(data[iH2]/data['density'])
print 'hm_min=', min(data[iHM]/data['density']), 'hm_max=', max(data[iHM]/data['density'])
print 'h2p_min=', min(data[iH2P]/data['density']), 'h2p_max=', max(data[iH2P]/data['density'])
print 'dm_min=', min(data[iDM]/data['density']), 'dm_max=', max(data[iDM]/data['density'])
print 'd2_min=', min(data[iD2]/data['density']), 'd2_max=', max(data[iD2]/data['density'])
print 'temp_min=', min(data['Temp']), 'temp_max=', max(data['Temp'])
print 'dens_min=', min(data['density']), 'dens_max=', max(data['density'])
'''

bin_num1 = 5000
bin_num2 = 200
rmin = 0.
rmin_menc = 0.
rmax = 3.e5 
#nmin = 1.e4
#nmax = 1.e12
nmin = 1.e6
nmax = 1.e14

extrema = {'nh_alt': (nmin, nmax), 'Temp': (100,2000)}
data = data.cut_region("obj['Temp'] < 2e3")

#location = [0,0,0]

#p3 = ProfilePlot(data, "nh_alt", ["density", "Temp", iH2, iHP,  iHD, iHEP], weight_field="cell_mass", n_bins=bin_num2)

p3 = create_profile(data, "nh_alt", ["density", "Temp", iH2, iHP,  iHD, iHEP], weight_field="cell_mass", n_bins=bin_num2, extrema=extrema)

#pc.save('radprof')

################read in other data files#######################################
dir = '/home6/astacy/gadget_runs/'
basename = 'bin_zoom10_new_cut_ref3'
datanum1 = '7054'
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

with open(dir + basename + '_nhprof_' + datanum2 + '.dat', "r") as f:
#with open('/home1/00863/minerva/StandAlone/output.dat', "r") as f:
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
#narr3 = 1379
narr3 = 200
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

profile = p3 #.profiles[0]
nh3 = profile.x
den = profile['density']
temp = profile['Temp']

dfac = den
#dfac = 1.0

h2 = profile[iH2]/dfac 
hii = profile[iHP] /dfac
hd = profile[iHD] /dfac
heii = profile[iHEP] /dfac

print 'nh3 =', nh3
print 'nhA2 =', nhA2
print 'h2 =', h2
print 'hd =', hd
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
hii_errB = []

for i in range(narr2):
	temp_err.append(np.fabs((temp[i].value - tempA2[i])/(temp[i].value  + tempA2[i])*2))
	h2_err.append(np.fabs((h2[i].value  - h2A2[i])/(h2[i].value  + h2A2[i])*2))
        hd_err.append(np.fabs((hd[i].value  - hdA2[i])/(hd[i].value  + hdA2[i])*2))
        hii_err.append(np.fabs((hii[i].value  - hiiA2[i])/(hii[i].value  + hiiA2[i])*2))

#for i in range(narr3):
#        temp_errB.append(np.fabs((temp[i] - tempB2[i])/(temp[i] + tempB2[i])*2))
#        h2_errB.append(np.fabs((h2[i] - h2B2[i])/(h2[i] + h2B2[i])*2))
#        hd_errB.append(np.fabs((hd[i] - hdB2[i])/(hd[i] + hdB2[i])*2))
#        dii_errB.append(np.fabs((dii[i] - diiB2[i])/(dii[i] + diiB2[i])*2))
#        hii_errB.append(np.fabs((hii[i] - hiiB2[i])/(hii[i] + hiiB2[i])*2))
#        dii_ratio.append(diiB2[i]/hiiB2[i])

##############################################################################
print 'heii =', heii
print 'hii_err =', hii_err

nmin = 1.e6
nmax = 1.e13
emin = 1.e-2
emax = 2.e0

pl.subplot(321)
pl.plot(nh3, temp,'k')
pl.plot(nhA2, tempA2, 'g')
pl.plot(nhB2, tempB2, 'b:')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
#ax.set_xlabel('number density [cm^{-3}]', fontsize=9)
ax.set_ylabel('Temperature [K]', fontsize=9)
#pl.xticks(fontsize=6)
ax.xaxis.set_ticklabels([])
pl.yticks(fontsize=10)
pl.axis((nmin, nmax, 1e2, 1e4))

pl.subplot(322)
pl.plot(nhA2, temp_err, 'k')
#pl.plot(nhB2, temp_errB, 'b:')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
#ax.set_xlabel('number density [cm^{-3}]', fontsize=9)
ax.set_ylabel('Temperature error', fontsize=9)
#pl.xticks(fontsize=6)
ax.xaxis.set_ticklabels([])
pl.yticks(fontsize=10)
pl.axis((nmin, nmax, emin, emax))

pl.subplot(323)
pl.plot(nh3, h2,'k')
pl.plot(nhA2, h2A2, 'g')
pl.plot(nhB2, h2B2, 'b:')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
#ax.set_xlabel('number density [cm^{-3}]', fontsize=9)
ax.set_ylabel('H2 ', fontsize=9)
#pl.xticks(fontsize=6)
ax.xaxis.set_ticklabels([])
pl.yticks(fontsize=10)
pl.axis((nmin, nmax, 1e-4, 1e0))


pl.subplot(324)
pl.plot(nhA2, h2_err, 'k')
#pl.plot(nhB2, h2_errB, 'b:')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
#ax.set_xlabel('number density [cm^{-3}]', fontsize=9)
ax.set_ylabel('H2 error', fontsize=9)
#pl.xticks(fontsize=6)
ax.xaxis.set_ticklabels([])
pl.yticks(fontsize=10)
pl.axis((nmin, nmax, emin, emax))

'''
pl.subplot(425)
pl.plot(nh3, hd,'k')
pl.plot(nhA2, hdA2, 'g')
pl.plot(nhB2, hdB2, 'b:')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
#ax.set_xlabel('number density [cm^{-3}]', fontsize=9)
ax.set_ylabel('HD ', fontsize=9)
#pl.xticks(fontsize=6)
ax.xaxis.set_ticklabels([])
pl.yticks(fontsize=10)
pl.axis((nmin, nmax, 1e-8, 1e-3))


pl.subplot(426)
pl.plot(nhA2, hd_err, 'k')
#pl.plot(nhB2, hd_errB, 'b:')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
#ax.set_xlabel('number density [cm^{-3}]', fontsize=9)
ax.set_ylabel('HD error', fontsize=9)
#pl.xticks(fontsize=6)
ax.xaxis.set_ticklabels([])
pl.yticks(fontsize=10)
pl.axis((nmin, nmax, emin, emax))
'''

pl.subplot(325)
pl.plot(nh3, hii,'k')
pl.plot(nhA2, hiiA2, 'g')
pl.plot(nhB2, hiiB2, 'b:')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('number density [cm$^{-3}$]', fontsize=9)
ax.set_ylabel('HII ', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((nmin, nmax, 1e-12, 1e-6))


pl.subplot(326)
pl.plot(nhA2, hii_err, 'k')
#pl.plot(nhB2, hii_errB, 'b:')
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('number density [cm$^{-3}$]', fontsize=9)
ax.set_ylabel('HII error', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((nmin, nmax, emin, emax))

pl.subplots_adjust(top=0.92, bottom=0.1, left=0.10, right=0.95, hspace=0.25,
                    wspace=0.35)

pl.show()
pl.savefig('radprof_err.eps')

