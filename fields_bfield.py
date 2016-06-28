from yt.mods import *
import numpy as na
from string import rstrip

gamma = 1.1

no_bfield = 1

cell_volume = (1.23427103e+18/64.) * (1.23427103e+18/64.) * (1.23427103e+18/64.)

if no_bfield == 1:
	def _xbfield(field, data):
		return  (data['density']*1.e-20)
	add_field("X-magnfield", function=_xbfield, take_log=False, units=r'\rm{G}')

        def _ybfield(field, data):
                return  (data['density']*1.e-20)
        add_field("Y-magnfield", function=_ybfield, take_log=False, units=r'\rm{G}')

        def _zbfield(field, data):
                return (data['density']*1.e-20)
        add_field("Z-magnfield", function=_zbfield, take_log=False, units=r'\rm{G}')

def unlimited(right,center,left):
   return 0.5*(right-left)

def _MagneticEnergy(field,data):
    return (data["X-magnfield"]**2 +
            data["Y-magnfield"]**2 +
            data["Z-magnfield"]**2)/2.
add_field("MagneticEnergy", function=_MagneticEnergy, take_log=True,
          units=r"",display_name=r"B^2/8\pi")
ChomboFieldInfo["MagneticEnergy"]._projected_units=r""

def _ufield(field, data):
    UF =  ((data['X-magnfield'] )**2 + (data['Y-magnfield'])**2 + (data['Z-magnfield'])**2)/ (8.0*3.14159)
    return UF
add_field("ufield", function=_ufield, take_log=False,
          units=r'\rm{G}')

def _Bmag(field,data):
    return (np.sqrt(data["X-magnfield"]**2 +
            data["Y-magnfield"]**2 +
            data["Z-magnfield"]**2))
add_field("Bmag", function=_Bmag, take_log=False,
          units=r"",display_name=r"Gauss")

def _vortx(field, data):
   sl_left = slice(None,-2,None)
   sl_right = slice(2,None,None)

   dsx = data['dx'].flat[0]
   dsy = data['dy'].flat[0]
   dsz = data['dz'].flat[0]
   out = na.zeros(data['velocity'].shape, dtype='float64')
   out[1:-1,1:-1,1:-1] =  unlimited(data['z-velocity'][sl_right,1:-1,1:-1], data['z-velocity'][1:-1,1:-1,1:-1], data['z-velocity'][sl_left,1:-1,1:-1])/dsy                                      - unlimited(data['y-velocity'][sl_right,1:-1,1:-1], data['y-velocity'][1:-1,1:-1,1:-1], data['y-velocity'][sl_left,1:-1,1:-1])/dsz
   return out
add_field('vortx', function=_vortx, validators=[ValidateSpatial(1,['velocity'])], take_log=False)

def _vorty(field, data):
   sl_left = slice(None,-2,None)
   sl_right = slice(2,None,None)

   dsx = data['dx'].flat[0]
   dsy = data['dy'].flat[0]
   dsz = data['dz'].flat[0]
   out = na.zeros(data['velocity'].shape, dtype='float64')   
   out[1:-1,1:-1,1:-1] =  unlimited(data['x-velocity'][sl_right,1:-1,1:-1], data['x-velocity'][1:-1,1:-1,1:-1], data['x-velocity'][sl_left,1:-1,1:-1])/dsz                                      - unlimited(data['z-velocity'][sl_right,1:-1,1:-1], data['z-velocity'][1:-1,1:-1,1:-1], data['z-velocity'][sl_left,1:-1,1:-1])/dsx
   return out
add_field('vorty', function=_vorty, validators=[ValidateSpatial(1,['velocity'])], take_log=False)

def _vortz(field, data):
   sl_left = slice(None,-2,None)
   sl_right = slice(2,None,None)

   dsx = data['dx'].flat[0]
   dsy = data['dy'].flat[0]
   dsz = data['dz'].flat[0]
   out = na.zeros(data['velocity'].shape, dtype='float64')   
   out[1:-1,1:-1,1:-1] =  unlimited(data['y-velocity'][sl_right,1:-1,1:-1], data['y-velocity'][1:-1,1:-1,1:-1], data['y-velocity'][sl_left,1:-1,1:-1])/dsx                                      - unlimited(data['x-velocity'][sl_right,1:-1,1:-1], data['x-velocity'][1:-1,1:-1,1:-1], data['x-velocity'][sl_left,1:-1,1:-1])/dsy
   return out
add_field('vortz', function=_vortz, validators=[ValidateSpatial(1,['velocity'])], take_log=False)

#care of YT!
#def _VorticityX(field, data):
#    # We need to set up stencils
#    if data.pf["HydroMethod"] == 2:
#        sl_left = slice(None,-2,None)
#        sl_right = slice(1,-1,None)
#        div_fac = 1.0
#    else:
#        sl_left = slice(None,-2,None)
#        sl_right = slice(2,None,None)
#        div_fac = 2.0
#    new_field = np.zeros(data["z-velocity"].shape, dtype='float64')
#    new_field[1:-1,1:-1,1:-1] = (data["z-velocity"][1:-1,sl_right,1:-1] -
#                                 data["z-velocity"][1:-1,sl_left,1:-1]) \
#                                 / (div_fac*data["dy"].flat[0])
#    new_field[1:-1,1:-1,1:-1] -= (data["y-velocity"][1:-1,1:-1,sl_right] -
#                                  data["y-velocity"][1:-1,1:-1,sl_left]) \
#                                  / (div_fac*data["dz"].flat[0])
#    return new_field


def _Bx_face(field,data):
    sl_right = slice(1,None,None)
    sl_left = slice(None,-1,None)
    f  = data["X-magnfield"][sl_right,0:-1,0:-1]
    f += data["X-magnfield"][sl_left,0:-1,0:-1]
    new_field = np.zeros(data["X-magnfield"].shape, dtype='float64')
    new_field[0:-1,0:-1,0:-1] = f/2.0;
    return new_field
add_field("Bx_face", function=_Bx_face, 
	validators=[ValidateSpatial(ghost_zones=1,
                       fields=["X-magnfield","Y-magnfield","Z-magnfield"])],
	take_log=False, units=r"",display_name=r"Gauss")

def _By_face(field,data):
    sl_right = slice(1,None,None)
    sl_left = slice(None,-1,None)
    f  = data["Y-magnfield"][0:-1,sl_right,0:-1]
    f += data["Y-magnfield"][0:-1,sl_left,0:-1]
    new_field = np.zeros(data["Y-magnfield"].shape, dtype='float64')
    new_field[0:-1,0:-1,0:-1] = f/2.0;
    return new_field
add_field("By_face", function=_By_face, 
	validators=[ValidateSpatial(ghost_zones=1,
                       fields=["X-magnfield","Y-magnfield","Z-magnfield"])],
	take_log=False, units=r"",display_name=r"Gauss")

def _Bz_face(field,data):
    sl_right = slice(1,None,None)
    sl_left = slice(None,-1,None)
    f  = data["Z-magnfield"][0:-1,0:-1,sl_right]
    f += data["Z-magnfield"][0:-1,0:-1,sl_left]
    new_field = np.zeros(data["Z-magnfield"].shape, dtype='float64')
    new_field[0:-1,0:-1,0:-1] = f/2.0;
    return new_field
add_field("Bz_face", function=_Bz_face, 
	validators=[ValidateSpatial(ghost_zones=1,
                       fields=["X-magnfield","Y-magnfield","Z-magnfield"])],
	take_log=False, units=r"",display_name=r"Gauss")

def _DivB(field, data):
    # We need to set up stencils
#    if data.pf["HydroMethod"] == 2:
#        sl_left = slice(None,-2,None)
#        sl_right = slice(1,-1,None)
#        div_fac = 1.0
#    else:
#        sl_left = slice(None,-2,None)
#        sl_right = slice(2,None,None)
#        div_fac = 2.0
    sl_right = slice(1,-1,None)
    sl_left = slice(None,-2,None)
    div_fac = 1.0
    ds = div_fac * data['dx'].flat[0]
    f  = data["X-magnfield"][sl_right,1:-1,1:-1]/ds
    f -= data["X-magnfield"][sl_left ,1:-1,1:-1]/ds
#    f  = data["Bx_face"][sl_right,1:-1,1:-1]/ds
#    f -= data["Bx_face"][sl_left,1:-1,1:-1]/ds
#    if data.pf.dimensionality > 1:
    ds = div_fac * data['dy'].flat[0]
    f += data["Y-magnfield"][1:-1,sl_right,1:-1]/ds
    f -= data["Y-magnfield"][1:-1,sl_left ,1:-1]/ds
#    f += data["By_face"][1:-1,sl_right,1:-1]/ds
#    f -= data["By_face"][1:-1,sl_left,1:-1]/ds
#    if data.pf.dimensionality > 2:
    ds = div_fac * data['dz'].flat[0]
    f += data["Z-magnfield"][1:-1,1:-1,sl_right]/ds
    f -= data["Z-magnfield"][1:-1,1:-1,sl_left ]/ds
#    f += data["Bz_face"][1:-1,1:-1,sl_right]/ds
#    f -= data["Bz_face"][1:-1,1:-1,sl_left ]/ds
    new_field = np.zeros(data["X-magnfield"].shape, dtype='float64')
    new_field[1:-1,1:-1,1:-1] = f
    #new_field[1:None,1:None,1:None] = f
    return new_field
def _convertDivB(data):
    return data.convert("cm")**-1.0
add_field("DivB", function=_DivB,
           validators=[ValidateSpatial(ghost_zones=1,
                       fields=["X-magnfield","Y-magnfield","Z-magnfield"])],
          units=r"\rm{Gauss / cm}", take_log=False,
          convert_function=_convertDivB)

def _absDivB(field, data):
    absdivb =  np.abs(data["DivB"])
    return absdivb
add_field("absDivB", function=_absDivB, take_log=False,
          units=r'\rm{Gauss / cm}')

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

