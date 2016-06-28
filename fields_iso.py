from yt.mods import *
import numpy as na
from string import rstrip

gamma = 1.0001

def _MagneticEnergy(field,data):
    return (data["X-magnfield"]**2 +
            data["Y-magnfield"]**2 +
            data["Z-magnfield"]**2)/2.
add_field("MagneticEnergy", function=_MagneticEnergy, take_log=True,
          units=r"",display_name=r"B^2/8\pi")
ChomboFieldInfo["MagneticEnergy"]._projected_units=r""

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

def _Velocity(field,data):
    """generate velocity from z-momentum and density
    
    """
    return na.sqrt(data["x-velocity"]*data["x-velocity"] + data["y-velocity"]*data["y-velocity"] + data["z-velocity"]*data["z-velocity"])
add_field("velocity",function=_Velocity, take_log=False,
          units=r'\rm{cm}/\rm{s}')

def _Vrad(field,data):
       # Radial velocity: dot product of r_hat with v
       out = np.zeros(data['Density'].shape, dtype='float64')
       vx = data['X-momentum']/data['Density']
       vy = data['Y-momentum']/data['Density']
       vz = data['Z-momentum']/data['Density']
       # Positions
       rmag = np.sqrt( 0.00001 + (data['x'])**2 +(data['y'])**2 + (data['z'])**2 )
       rx = data['x']/rmag
       ry = data['y']/rmag
       rz = data['z']/rmag
       out = (rx*vx + ry*vy + rz*vz)  #with a negative sign out front, infalling is positive
       #out = np.sqrt(rx*rx+ry*ry+rz*rz) #test, should be radial distance
       return out
add_field('Vrad',function=_Vrad, validators=[ValidateSpatial(1,['Density'])], take_log=False)

def _Vmag(field,data):
       # Radial velocity: dot product of r_hat with v
       out = np.zeros(data['Density'].shape, dtype='float64')
       vx = data['X-momentum']/data['Density']
       vy = data['Y-momentum']/data['Density']
       vz = data['Z-momentum']/data['Density']
       # Positions
       rmag = np.sqrt( 0.00001 + (data['x'])**2 +(data['y'])**2 + (data['z'])**2 )
       rx = data['x']/rmag
       ry = data['y']/rmag
       rz = data['z']/rmag
       out = (rx*vx + ry*vy + rz*vz)  #with a negative sign out front, infalling is positive
       #out = np.sqrt(rx*rx+ry*ry+rz*rz) #test, should be radial distance
       return out
add_field('Vmag',function=_Vmag, validators=[ValidateSpatial(1,['Density'])], take_log=False)

# def _rVelocity(field,data):
#     """ generate the radial velocity from the other components

#     """
#     return na.sqrt(data["x-velocity"]**2 + data["y-velocity"]**2 + data['z-velocity']**2)
# add_field("radial-velocity", function=_rVelocity, take_log=False,
#           units=r'\rm{cm}/\rm{s}')


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


def _cs(fields,data):
#    return 132588*1.0001*data["density"]/data["density"]
    return 180000*1.0001*data["density"]/data["density"]
add_field("cs", function=_cs, units=r"\rm{cm}/\rm{s}")


#Temp = (gamma-1) * mu* m_H / k_B  * ThermalEnergy/density
def _Temperature(field,data):
#The 1.650882e12 factor assumes gamma=1.0001, which is NOT CORRECT for ARS case!!!
#    return (data["ThermalEnergy"]/(1.650882e12*data["density"]))
    return  (pow(data['cs'],2) * 1.67e-24 / 1.38e-16) 
add_field("Temperature",function=_Temperature,units=r"\rm{Kelvin}",take_log=False)

def _TotalEnergyDensity(field,data):
    """Total non-magnetic energy, including radiation.
    """
    return data["energy-density"] + data['radiation-energy-density']
add_field("total-energy-density",function=_TotalEnergyDensity,units=r"\rm{ergs}/\rm{cm^3}")

def _NumDens(field,data):
    """Number Density, as derived from mass density
    """
    return (data["density"]/(1.22*1.67e-24))
add_field("number-density",function=_NumDens,units=r"\rm{cm^{-3}}")

#Pressure = (gamma-1)*(specific-internal-energy)*rho
#         =(gamma-1)*u*rho in Gadget
def _Pressure(field,data):
    """M{(Gamma-1.0)*e, where e is thermal energy density
    NB: this will need to be modified for radiation
    """
    return 1.38e-16*data["number-density"]*data["Temperature"]
add_field("Pressure", function=_Pressure, units=r"\rm{dyne}/\rm{cm}^{2}")

def _mach(fields,data):
    return na.fabs(data["Vrad"])/na.sqrt(gamma*data["Pressure"]/data["density"])
add_field("mach", function=_mach, units=r"[ ]", take_log=False)

