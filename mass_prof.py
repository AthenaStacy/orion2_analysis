from yt.mods import *
import numpy as na
from string import rstrip
import fields

pf = load("/work/00863/minerva/orion/data.0000.3d.hdf5")
value, location = pf.h.find_max("Density")
data = pf.h.sphere(location, 3.0/pf['pc'])

den0  = data['density']
x     = data['x']
temp0 = data['Temperature']
vx0   = data['x-velocity']
time  = pf.current_time

M_i = data["CellMassMsun"].sum()
print "total mass = " 
print M_i
