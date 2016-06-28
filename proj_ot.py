from yt.mods import *
import numpy as na
from string import rstrip
import fields_bfield
import matplotlib.colorbar as cb

datanum = '0050'

pf = load("/work/00863/minerva/orion_Otest/lowBhighRes/" + 'data.' + datanum + '.3d.hdf5')

#prj = ProjectionPlot(pf,2,'z-velocity',weight_field=None,max_level=2)
#prj.save()

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
######################################################################################

dmin = 1.e-1
dmax = 5.e-1
prj = ProjectionPlot(pf, 'z', 'density', weight_field = "Density", width=(1.0,'cm'))
prj.set_zlim('density', dmin, dmax)
#prj.set_cmap("density", "Rainbow + white")
prj.set_cmap("density", "PRISM")
prj.save('data.'+ datanum)

#bmin = -21
#bmax = -18
#bmin = -4
#bmax = 0

prj = ProjectionPlot(pf, 'z', 'Bmag-log', weight_field = "Bmag-log", width=(1.0,'cm'))
#prj.set_zlim('Bmag-log', bmin, bmax)
prj.set_cmap("Bmag-log", "Rainbow18")
prj.save('data.'+ datanum)

prj = ProjectionPlot(pf, 'z', 'Bx-log', weight_field = "Bx-log", width=(1.0,'cm'))
#prj.set_zlim('Bx-log', bmin, bmax)
prj.set_cmap("Bx-log", "Rainbow18")
prj.save('data.'+ datanum)

prj = ProjectionPlot(pf, 'z', 'By-log', weight_field = "By-log", width=(1.0,'cm'))
#prj.set_zlim('By-log', bmin, bmax)
prj.set_cmap("By-log", "Rainbow18")
prj.save('data.'+ datanum)

prj = ProjectionPlot(pf, 'z', 'Bz-log', weight_field = "Bz-log", width=(1.0,'cm'))
#prj.set_zlim('Bz-log', bmin, bmax)
prj.set_cmap("Bz-log", "Rainbow18")
prj.save('data.'+ datanum)











