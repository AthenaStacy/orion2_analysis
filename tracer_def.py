from yt.mods import *
import numpy as na
from string import rstrip

########################################################################################
#define some stuff######################################################################

h2conv = (2.0/1.3158)
heconv = (4.0/1.3158)
gamma = 1.66667

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

def _dfac(field,data):
    return data["density"]#/data["density"]
add_field("dfac",function=_dfac, take_log=False, units=' ')

def _mu(field,data):
   h2frac = data[iH2] * h2conv /data['dfac']
   hefrac = data[iHE] * heconv /data['dfac']  #should be approximately 0.24
   hfrac = 1. - h2frac - hefrac
   #mu_inv = (1. - h2frac)*0.76 + h2frac*0.76/2.0 + 0.24/4.0 #here h2frac ranges from 0-1
   mu_inv = hfrac + h2frac/2. + hefrac/4.
   mu = 1. / mu_inv
   return  (mu )
add_field("mu",function=_mu,units=r"\rm{Kelvin}",take_log=True)

def _xVelocity(field, data):
    """generate x-velocity from x-momentum and density
    
    """
    return data["X-momentum"]/data["density"]
add_field("x-velocity",function=_xVelocity, take_log=False,
          units=r'\rm{cm}/\rm{s}')

def _yVelocity(field,data):
    """generate y-velocity from y-momentum and density
    
    """
    return data["Y-momentum"]/data["density"]
add_field("y-velocity",function=_yVelocity, take_log=False,
          units=r'\rm{cm}/\rm{s}')

def _zVelocity(field,data):
    """generate z-velocity from z-momentum and density
    
    """
    return data["Z-momentum"]/data["density"]
add_field("z-velocity",function=_zVelocity, take_log=False,
          units=r'\rm{cm}/\rm{s}')


#ThermalEnergy = (specific-internal-energy)*rho
def _ThermalEnergy(field, data):
    """generate thermal (gas energy). Dual Energy Formalism was
    implemented by Stella, but this isn't how it's called, so I'll
    leave that commented out for now.
    """
    #if data.pf["DualEnergyFormalism"]:
    #    return data["Gas_Energy"]
    #else:
    return data["energy-density"]  - 0.5 * data["density"] * (
        data["x-velocity"]**2.0
        + data["y-velocity"]**2.0
        + data["z-velocity"]**2.0 )
add_field("ThermalEnergy", function=_ThermalEnergy,
          units=r"\rm{ergs}/\rm{cm^3}")


def _Temp_alt(field,data):
   h2frac = data[iH2] * h2conv /data['dfac']
   hefrac = data[iHE] * heconv /data['dfac']  #should be approximately 0.24
   hfrac = 1. - h2frac - hefrac
   #mu_inv = (1. - h2frac)*0.76 + h2frac*0.76/2.0 + 0.24/4.0 #here h2frac ranges from 0-1
   mu_inv = hfrac + h2frac/2. + hefrac/4.
   mu = 1. / mu_inv
   #return  (data["ThermalEnergy"] * (data['tracer1']-1.0) / data["density"] * mu * 1.67e-24 / 1.38e-16)
   return  (data["ThermalEnergy"] * (gamma-1.0) / data["density"] * mu * 1.67e-24 / 1.38e-16)
add_field("Temp.",function=_Temp_alt,units=r"\rm{Kelvin}",take_log=True)
########################################################################################
