from yt.mods import *
import numpy as na
from string import rstrip
import fields
import fields_bfield
import matplotlib.colorbar as cb
import fields_bfield
from tracer_def import *

datanum = '0830'

pf = load("/nobackupp7/astacy//popiii_Bscope17/" + 'data.' + datanum + '.3d.hdf5')


######################################################################################
######################################################################################

#length = 1.0
length = 0.025

#dmin = 1.e-1
#dmax = 5.e-1
prj = ProjectionPlot(pf, 'z', 'density', weight_field = "Density", width=(length,'pc'))
#prj.set_zlim('density', dmin, dmax)
#prj.set_cmap("density", "Rainbow + white")
#prj.set_cmap("density", "PRISM")
prj.annotate_sphere([0,0,0], 1.e-5, {'fill':True})
prj.save('data.'+ datanum)

prj = ProjectionPlot(pf, 'z', 'Temperature', weight_field = "Temperature", width=(length,'pc'))
prj.save('data.'+ datanum)

prj = ProjectionPlot(pf, 'z', 'tracer1', weight_field = "tracer1", width=(length,'pc'))
prj.save('data.'+ datanum)

prj = ProjectionPlot(pf, 'z', 'tracer2', weight_field = "tracer2", width=(length,'pc'))
prj.save('data.'+ datanum)

#bmin = -21
#bmax = -18
bmin = -15
bmax = -13.5

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

prj = ProjectionPlot(pf, 'z', 'X-magnfield', weight_field = "X-magnfield", width=(length,'pc'))
#prj.set_zlim('Bx-log', bmin, bmax)
prj.set_cmap("X-magnfield", "Rainbow18")
prj.save('data.'+ datanum)

prj = ProjectionPlot(pf, 'z', 'Y-magnfield', weight_field = "Y-magnfield", width=(length,'pc'))
#prj.set_zlim('By-log', bmin, bmax)
prj.set_cmap("Y-magnfield", "Rainbow18")
prj.save('data.'+ datanum)

prj = ProjectionPlot(pf, 'z', 'Bmag-log', weight_field = "Bmag-log", width=(length,'pc'))
#prj.set_zlim('Bmag-log', bmin, bmax)
#prj.set_cmap("Y-magnfield", "Rainbow18")
prj.save('data.'+ datanum)




