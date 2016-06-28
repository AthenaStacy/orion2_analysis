from yt.mods import *
import numpy as na
from string import rstrip
import fields

#pf = load("/work/00863/minerva/orion/data.30825.3d.hdf5")
#pf = load("/work/00863/minerva/orion/data.0000.tsc")
pf = load("/work/00863/minerva/orion/data.0000.3d.hdf5")
value, location = pf.h.find_max("Density")
data = pf.h.sphere(location, 30.0/pf['pc'])
#xval = 0.0
#yval = -1.e18
#zval = 1.e18
#data = pf.h.region([-xval, -xval, -xval], [yval, yval, yval], [zval, zval, zval])

#M_i = data["CellMassGrams"].sum()

M_i = data["CellMassGrams"].sum() / 1.98892e33 

print "total mass = ", M_i
