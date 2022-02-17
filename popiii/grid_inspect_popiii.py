
import yt
from yt.units import dimensions
import matplotlib
matplotlib.use('ps')
from yt.mods import *
import numpy as na
import pylab as pl
import fields_bfield
import tracer_def

dir_name = '/work/00863/minerva/stampede2/popiii_bfieldA/'
datanum = '0245'

ds = load(dir_name + 'data.' + datanum + ".3d.hdf5")
#gs = ds.index.select_grids(ds.index.max_level)
#gs = ds.index.select_grids(0)

level_0  = ds.covering_grid(level=0, left_edge=[0,0.0,0.0],
                                      dims=ds.domain_dimensions)

print(level_0['X-momentum'].shape)

'''
print('min X-momentum = ', np.min(level_0['X-momentum']))
print('min Y-momentum = ', np.min(level_0['Y-momentum']))
print('min Z-momentum = ', np.min(level_0['Z-momentum']))

print('min X-magnfield = ', np.min(level_0['X-magnfield']))
print('min Y-magnfield = ', np.min(level_0['Y-magnfield']))
print('min Z-magnfield = ', np.min(level_0['Z-magnfield']))

print('max X-momentum = ', np.max(level_0['X-momentum']))
print('max Y-momentum = ', np.max(level_0['Y-momentum']))
print('max Z-momentum = ', np.max(level_0['Z-momentum']))

print('max X-magnfield = ', np.max(level_0['X-magnfield']))
print('max Y-magnfield = ', np.max(level_0['Y-magnfield']))
print('max Z-magnfield = ', np.max(level_0['Z-magnfield']))
'''

U_B = (level_0["X-magnfield"]**2 + level_0["Y-magnfield"]**2 + level_0["Z-magnfield"]**2)/2

UB_tot_yt1 = np.sum(level_0['UB_PS'])
UB_tot_yt2 = np.sum(level_0['UB_PS2'])
UB_tot_python = np.sum(U_B)
UB_tot_athena = np.sum(level_0['MagEnergy'])

UK_tot = np.sum(level_0['UK_PS'])
UK_tot_yt = np.sum(level_0['KinEnergy'])


print('UB_tot_yt1 = ', UB_tot_yt1)
print('UB_tot_yt2 = ', UB_tot_yt2)
print('UB_tot_python = ', UB_tot_python)
print('UB_tot_athena = ', UB_tot_athena)


print('UK_tot = ', UK_tot)
print('UK_tot_yt = ', UK_tot_yt)

B_mag_here = np.sqrt(level_0["X-magnfield"]**2 + level_0["Y-magnfield"]**2 + level_0["Z-magnfield"]**2) * np.sqrt(4*3.14159)

B_mag_here2 = (level_0["X-magnfield"]**2 + level_0["Y-magnfield"]**2 + level_0["Z-magnfield"]**2) * 4*3.14159

B_tot_python = np.sum(B_mag_here)
B_tot_python2 = np.sum(B_mag_here2)
B_tot_yt = np.sum(level_0['Bmag'])
B_tot_yt2 = np.sum(level_0['Bmag2'])

print('B_tot_yt = ', B_tot_yt)
print('B_tot_yt2 = ', B_tot_yt2)
print('B_tot_python = ', B_tot_python)
print('B_tot_python2 = ', B_tot_python2)
