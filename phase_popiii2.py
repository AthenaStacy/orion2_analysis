from yt.mods import *
import numpy as na
from string import rstrip
import fields
import fields_bfield
import matplotlib.colorbar as cb

datanum = '0143'

#pf = load("/work/00863/minerva/uni_vargam/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/popiii_Bscope/" + 'data.' + datanum + '.3d.hdf5')
pf = load("/work/00863/minerva/popiii_Bscope2/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/popiii_Bscope3/" + 'data.' + datanum + '.3d.hdf5')

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

def _Temp_alt(field,data):
   h2frac = data[iH2]/data['density']
   hefrac = data[iHE]/data['density']  #should be approximately 0.24
   hfrac = 1. - h2frac - hefrac
   #mu_inv = (1. - h2frac)*0.76 + h2frac*0.76/2.0 + 0.24/4.0 #here h2frac ranges from 0-1
   mu_inv = hfrac + h2frac/2. + hefrac/4.
   mu = 1. / mu_inv
   return  (data["ThermalEnergy"] * (data['tracer1']-1.0) / data["density"] * mu * 1.67e-24 / 1.38e-16)
add_field("Temp_alt",function=_Temp_alt,units=r"\rm{Kelvin}",take_log=True)

def _Radius(field, data):
    '''
    Distance from central point. In this problem the center is at density peak.
    '''
    rad =  na.sqrt((data['x'] )**2 + (data['y'])**2 + (data['z'])**2)
    rad = rad/1.5e13
    return rad
add_field("Radius", function=_Radius, take_log=False,
          units=r'\rm{AU}')

def _H2_abund(field, data):
    abund = data[iH2] / data['density'] 
    return abund
add_field("H2_abund", function=_H2_abund, take_log=False, units='')

def _HP_abund(field, data):
    abund = data[iHP] / data['density']
    return abund
add_field("HP_abund", function=_HP_abund, take_log=False, units='')

######################################################################################

length = 10.0

my_sphere = pf.h.sphere("c", length/pf['pc'])

pc = PlotCollection(pf, "c")

pc.add_phase_object(my_sphere, ["Density", "Temp_alt", "CellMassMsun"])

pc.add_phase_object(my_sphere, ["Density", igamma, "CellMassMsun"])

pc.add_phase_object(my_sphere, ["Density", "H2_abund", "CellMassMsun"])

pc.add_phase_object(my_sphere, ["Density", "HP_abund", "CellMassMsun"])


pc.save("phase")


