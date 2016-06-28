from yt.mods import *
import numpy as na
from string import rstrip

datanum = '0050'

pf = load("/work/00863/minerva/orion_Otest/" + 'data.' + datanum + '.3d.hdf5')

#prj = ProjectionPlot(pf,2,'z-velocity',weight_field=None,max_level=2)
#prj.save()

dmin = 1.e-1
dmax = 5.e-1
prj = ProjectionPlot(pf, 'z', 'density', weight_field = "Density", width=(1.0,'cm')).save()

