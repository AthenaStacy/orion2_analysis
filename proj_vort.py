from yt.mods import *
import numpy as na
from string import rstrip
import fields
import fields_bfield
import matplotlib.colorbar as cb

datanum = '0002'

#pf = load("/work/00863/minerva/orion_Otest/lowBhighRes/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion_Otest/lowB/" + 'data.' + datanum + '.3d.hdf5')
pf = load("/work/00863/minerva/orion_Vtest/" + 'data.' + datanum + '.3d.hdf5')


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
#prj.set_zlim('density', dmin, dmax)
#prj.set_cmap("density", "Rainbow + white")
#prj.set_cmap("density", "PRISM")
prj.save('data.'+ datanum)

bmin = -21
bmax = -18
#bmin = -4
#bmax = 0
prj = ProjectionPlot(pf, 'z', 'Bmag-log', weight_field = "Bmag-log", width=(1.0,'cm'))
prj.set_zlim('Bmag-log', bmin, bmax)
prj.set_cmap("Bmag-log", "Rainbow18")
prj.save('data.'+ datanum)

prj = ProjectionPlot(pf, 'z', 'Bx-log', weight_field = "Bx-log", width=(1.0,'cm'))
prj.set_zlim('Bx-log', bmin, bmax)
prj.set_cmap("Bx-log", "Rainbow18")
prj.save('data.'+ datanum)

prj = ProjectionPlot(pf, 'z', 'By-log', weight_field = "By-log", width=(1.0,'cm'))
prj.set_zlim('By-log', bmin, bmax)
prj.set_cmap("By-log", "Rainbow18")
prj.save('data.'+ datanum)

prj = ProjectionPlot(pf, 'z', 'Bz-log', weight_field = "Bz-log", width=(1.0,'cm'))
prj.set_zlim('Bz-log', bmin, bmax)
prj.set_cmap("Bz-log", "Rainbow18")
prj.save('data.'+ datanum)

prj = ProjectionPlot(pf, 'z', 'DivB-log', weight_field = "DivB-log", width=(1.0,'cm'))
#prj.set_zlim('Bz-log', bmin, bmax)
prj.set_cmap("DivB-log", "Rainbow18")
prj.save('data.'+ datanum)

location = [0,0,0]
data = pf.h.sphere(location, 1.e3/pf['pc'])
print 'min vx=', min(data['x-velocity']), 'max vx=', max(data['x-velocity'])
print 'min vy=', min(data['y-velocity']), 'max vy=', max(data['y-velocity'])
print 'min vz=', min(data['z-velocity']), 'max vz=', max(data['z-velocity'])
print 'min vlog=', min(data['vlog']), 'max vlog=', max(data['vlog'])

prj = ProjectionPlot(pf, 'z', 'x-velocity', weight_field = "x-velocity", width=(1.0,'cm'))
prj.save('data.'+ datanum)

prj = ProjectionPlot(pf, 'z', 'y-velocity', weight_field = "y-velocity", width=(1.0,'cm'))
prj.save('data.'+ datanum)

#vmin = np.log10(1.e-1)
vmin = np.log10(1.e-3)
vmax = np.log10(1.5e0)
prj = ProjectionPlot(pf, 'z', 'vlog', weight_field = "vlog", width=(1.0,'cm'))
#prj.set_zlim('vlog', vmin, vmax)
#prj.set_cmap("vlog", "Rainbow18")
prj.save('data.'+ datanum)










