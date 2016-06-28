from yt.mods import *
import numpy as na
import fields

pf = load("/work/00863/minerva/orion/data.0004.3d.hdf5")
pc = PlotCollection(pf, 'c')
pc.add_profile_sphere(100.0, 'pc', ["number-density", "Temperature"])
#pc.add_profile_sphere(100.0, 'pc', ["number-density", "ThermalEnergy"])
pc.save()
