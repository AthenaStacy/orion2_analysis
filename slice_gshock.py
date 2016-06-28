import matplotlib
matplotlib.use('pdf')

from yt.mods import *
import numpy as na
from string import rstrip
#import fields_iso
import fields
import fields_bfield
import fields_wpot

datanum = '0000'

pf = load("/work/00863/minerva/gshock/" + 'data.' + datanum + '.3d.hdf5')

print  pf.h.field_list
print 'dimensions =', pf.dimensionality

########################################################get B-field logs###############
def _Bx_log(field,data):
    return np.log10(np.abs(data["X-magnfield"]))
add_field("Bx-log",function=_Bx_log, take_log=False,
          units=r'Gauss')

def _By_log(field,data):
    return np.log10(np.abs(data["Y-magnfield"]))
add_field("By-log",function=_By_log, take_log=False,
          units=r'Gauss')

def _Bz_log(field,data):
    return np.log10(np.abs(data["Z-magnfield"]))
add_field("Bz-log",function=_Bz_log, take_log=False,
          units=r'Gauss')

def _Bmag_log(field,data):
    return np.log10(np.abs(data["Bmag"]))
add_field("Bmag-log",function=_Bmag_log, take_log=False,
          units=r'Gauss')

def _DivB_log(field,data):
    return np.log10((data["absDivB"]))
add_field("DivB-log",function=_DivB_log, take_log=False,
          units=r'Gauss / cm')


dims = list(pf.domain_dimensions)
DeltaX = pf.domain_right_edge[0] / dims[0]
print 'DeltaX =', DeltaX
def _errDivB(field,data):
    return np.log10(data["absDivB"] * DeltaX / data["Bmag"])
add_field("errDivB",function=_errDivB, take_log=False,
          units=r' ')

def _Temp_alt(field,data):
   #h2frac = data["tracer6"]
   h2frac = 1.e-3
   hefrac = 0.24  #should be approximately 0.24
   hfrac = 1. - h2frac - hefrac
   #mu_inv = (1. - h2frac)*0.76 + h2frac*0.76/2.0 + 0.24/4.0 #here h2frac ranges from 0-1
   mu_inv = hfrac + h2frac/2. + hefrac/4.
   mu = 1. / mu_inv
   #return  (data["ThermalEnergy"] * (data['tracer1']-1.0) / data["density"] * mu * 1.67e-24 / 1.38e-16)
   return  (data["ThermalEnergy"] * (5./3.-1.0) / data["density"] * mu * 1.67e-24 / 1.38e-16)
add_field("Temp_alt",function=_Temp_alt,units=r"\rm{Kelvin}",take_log=True)
######################################################################################


length = 0.1 #length of box we want in pc

#value, location = pf.h.find_max("Density")
#data = pf.h.sphere(location, length/pf['pc'])
#print 'max density location = ', location

value, location = pf.h.find_min("Temp_alt")
data = pf.h.sphere(location, length/pf['pc'])
print 'min Temp location = ', location


#SlicePlot(pf, 'x', "z-velocity", width = (1, 'pc')).save()

#define        iHP   0  (tracer2)
#define        iH    1  (tracer3)
#define        iHM   2  (tracer4)
#define        iH2P  3  (tracer5)
#define        iH2   4  (tracer6)
#define        iDP   5  (tracer7)
#define        iD    6  (tracer8)
#define        iDM   7  (tracer9)
#define        iHDP  8  (tracer10)
#define        iHD   9  (tracer11)
#define        iD2P  10 (tracer12)
#define        iD2   11 (tracer13)
#define        iHEP  12 (tracer14)
#define        iHE   13 (tracer15)
#define        iHEPP 14 (tracer16)
#define        iELEC 15  (tracer17)


pc=PlotCollection(pf, center=[0,0,0])
#pc=PlotCollection(pf, center=location)

dims = list(pf.domain_dimensions)
ngrid_x = dims[0]
ngrid_y = dims[1]
#proj = pf.h.slice(2, 0.0, center=pf.domain_center)
proj = pf.h.slice(2, location[2], center=location)
w = (pf.h.domain_left_edge[0], pf.h.domain_right_edge[0], pf.h.domain_left_edge[1], pf.h.domain_right_edge[1])
frb1 = FixedResolutionBuffer(proj, w, (ngrid_x, ngrid_x), periodic=False)
bx  = frb1['X-magnfield']
by  = frb1['Y-magnfield']
bz  = frb1['Z-magnfield']
bx_face  = frb1['Bx_face']
by_face  = frb1['By_face']
bz_face  = frb1['Bz_face']
rho = frb1['Density']
velx = frb1['x-velocity']
x1  = frb1['x']
y1  = frb1['y']
z1  = frb1['z']


print 'min density=', min(data['density']), 'max density=', max(data['density'])
print 'min Bfield-X=', min(data['X-magnfield']), 'max Bfield-X=', max(data['X-magnfield'])
print 'min Bfield-Y=', min(data['Y-magnfield']), 'max Bfield-Y=', max(data['Y-magnfield'])
print 'min Bfield-Z=', min(data['Z-magnfield']), 'max Bfield-Z=', max(data['Z-magnfield'])
print 'min log-Bfield-X=', min(data['Bx-log']), 'max log-Bfield-X=', max(data['Bx-log'])
print 'min log-Bfield-Y=', min(data['By-log']), 'max log-Bfield-Y=', max(data['By-log'])
print 'min log-Bfield-Z=', min(data['Bz-log']), 'max log-Bfield-Z=', max(data['Bz-log'])

bmin = -1.e-11
bmax =  1.e-11

#p=pc.add_slice('X-magnfield', 'x')
#p.set_zlim(bmin, bmax)
#p=pc.add_slice('Y-magnfield', 'x')
#p.set_zlim(bmin, bmax)
#p=pc.add_slice('Z-magnfield', 'x')
#p.set_zlim(bmin, bmax)

#p=pc.add_slice('X-magnfield', 'y')
#p.set_zlim(bmin, bmax)
#p=pc.add_slice('Y-magnfield', 'y')
#p.set_zlim(bmin, bmax)
#p=pc.add_slice('Z-magnfield', 'y')
#p.set_zlim(bmin, bmax)

#p=pc.add_slice('X-magnfield', 'z')
#p.set_zlim(bmin, bmax)
#p=pc.add_slice('Y-magnfield', 'z')
#p.set_zlim(bmin, bmax)
#p=pc.add_slice('Z-magnfield', 'z')
#p.set_zlim(bmin, bmax)

bmin = -14
bmax = -10
#bmin = -15
#bmax = -12
#bmin = -20
#bmax = -13

#p=pc.add_slice('Bx-log', 'z')
#p.set_zlim(bmin, bmax)
#p.set_cmap("Rainbow18")


emin = 1e-18
emax = 1e-16
#p=pc.add_slice('MagneticEnergy', 'x')
#p.set_zlim(emin, emax)

#p=pc.add_slice('VorticityX', 'x')
#p=pc.add_slice('VorticityMagnitude', 'x')

#p=pc.add_slice('DivB', 'x')
#p=pc.add_slice('DivB', 'y')
#p=pc.add_slice('DivB', 'z')
#p=pc.add_slice('DivB-log', 'x')

tmin = 100.
tmax = 2000.
p=pc.add_slice('Temp_alt', 'x')
#p.set_width(length,'pc')
#p.set_zlim(tmin, tmax)
p.set_cmap("Rainbow18")
p=pc.add_slice('Temp_alt', 'y')
#p.set_width(length,'pc')
#p.set_zlim(tmin, tmax)
p.set_cmap("Rainbow18")
p=pc.add_slice('Temp_alt', 'z')
#p.set_width(length,'pc')
#p.set_zlim(tmin, tmax)
p.set_cmap("Rainbow18")

gmin = 1.65
gmax = 1.67
p=pc.add_slice('tracer1', 'x')
#p.set_zlim(gmin, gmax)
p=pc.add_slice('tracer1', 'y')
#p.set_zlim(gmin, gmax)
p=pc.add_slice('tracer1', 'z')
#p.set_zlim(gmin, gmax)

p=pc.add_slice('tracer6', 'x')
p=pc.add_slice('tracer6', 'y')
p=pc.add_slice('tracer6', 'z')

emin = -3
emax = 1

dmin = 1.e-27
dmax = 1.e-26
p=pc.add_slice('density', 'x')
#p.set_width(length,'pc')
#pc.plots[-1].modify["grids"]()
p=pc.add_slice('density', 'y')
#p.set_width(length,'pc')
#pc.plots[-1].modify["grids"]()
p=pc.add_slice('density', 'z')
#p.set_width(length,'pc')
#pc.plots[-1].modify["grids"]()

p=pc.add_slice('x-velocity', 'x')
#p.set_width(length,'pc')
p=pc.add_slice('x-velocity', 'y')
#p.set_width(length,'pc')
p=pc.add_slice('x-velocity', 'z')
#p.set_width(length,'pc')

p=pc.add_slice('y-velocity', 'x')
#p.set_width(length,'pc')
p=pc.add_slice('y-velocity', 'y')
#p.set_width(length,'pc')
p=pc.add_slice('y-velocity', 'z')
#p.set_width(length,'pc')

p=pc.add_slice('z-velocity', 'x')
#p.set_width(length,'pc')
p=pc.add_slice('z-velocity', 'y')
#p.set_width(length,'pc')
p=pc.add_slice('z-velocity', 'z')
#p.set_width(length,'pc')


p=pc.add_slice('velocity', 'x')
#p.set_width(length,'pc')
p=pc.add_slice('velocity', 'y')
#p.set_width(length,'pc')
p=pc.add_slice('velocity', 'z')
#p.set_width(length,'pc')


print 'total divB = ', data['DivB'].sum()
print 'min density=', min(data['density']), 'max density=', max(data['density'])
print 'min tracer1=', min(data['tracer1']), 'max tracer1=', max(data['tracer1'])
print 'min tracer3=', min(data['tracer3']), 'max tracer3=', max(data['tracer3'])
print 'min tracer6=', min(data['tracer6']), 'max tracer6=', max(data['tracer6'])
#print 'min gravpot=', min(data['gravitational-potential']), 'max gravpot=', max(data['gravitational-potential'])
#print 'min Temp=', min(data['Temperature']), 'max Temp=', max(data['Temperature'])
print 'min Temp=', min(data['Temp_alt']), 'max Temp=', max(data['Temp_alt'])

print 'B average = ', np.mean(data['Bmag'])
print 'dens average =', np.mean(data['density']) 

pc.save('data_'+ datanum)


