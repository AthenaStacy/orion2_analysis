from yt.mods import *
import numpy as na
from string import rstrip
import fields
import fields_bfield
import matplotlib.colorbar as cb

datanum = '0137'

pf = load("/work/00863/minerva/popiii_Bscope5/" + 'data.' + datanum + '.3d.hdf5')

#prj = ProjectionPlot(pf,2,'z-velocity',weight_field=None,max_level=2)
#prj.save()

#define        iHP   0  (tracer2)
#define        iH    1  (tracer3)
#define        iHM   2  (tracer4)
#define        iH2P  3  (tracer5)
#define        iH2   4  (tracer6)
#define        iDP   5  (tracer7)
#define        iD    6
#define        iDM   7
#define        iHDP  8
#define        iHD   9
#define        iD2P  10
#define        iD2   11
#define        iHEP  12
#define        iHE   13
#define        iHEPP 14
#define        iELEC 15  (tracer17)

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

length = 0.1

dmin = 1.e-1
dmax = 5.e-1
prj = ProjectionPlot(pf, 'z', 'density', weight_field = "Density", width=(length,'pc'))
prj.set_zlim('density', dmin, dmax)
#prj.set_cmap("density", "Rainbow + white")
prj.set_cmap("density", "PRISM")
prj.save('data.'+ datanum)

prj = ProjectionPlot(pf, 'z', 'tracer1', weight_field = "tracer1", width=(length,'pc'))
prj.save('data.'+ datanum)

prj = ProjectionPlot(pf, 'z', 'tracer6', weight_field = "tracer1", width=(length,'pc'))
prj.save('data.'+ datanum)

#bmin = -21
#bmax = -18

prj = ProjectionPlot(pf, 'z', 'Bmag-log', weight_field = "Bmag-log", width=(length,'pc'))
#prj.set_zlim('Bmag-log', bmin, bmax)
prj.set_cmap("Bmag-log", "Rainbow18")
prj.save('data.'+ datanum)

prj = ProjectionPlot(pf, 'z', 'Bx-log', weight_field = "Bx-log", width=(length,'pc'))
#prj.set_zlim('Bx-log', bmin, bmax)
prj.set_cmap("Bx-log", "Rainbow18")
prj.save('data.'+ datanum)

prj = ProjectionPlot(pf, 'z', 'By-log', weight_field = "By-log", width=(length,'pc'))
#prj.set_zlim('By-log', bmin, bmax)
prj.set_cmap("By-log", "Rainbow18")
prj.save('data.'+ datanum)

prj = ProjectionPlot(pf, 'z', 'Bz-log', weight_field = "Bz-log", width=(length,'pc'))
#prj.set_zlim('Bz-log', bmin, bmax)
prj.set_cmap("Bz-log", "Rainbow18")
prj.save('data.'+ datanum)







