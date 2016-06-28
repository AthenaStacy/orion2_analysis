import matplotlib
matplotlib.use('pdf')

from yt.mods import *
import numpy as na
from string import rstrip

#datanum = '0000'
datanum = '0002'

#pf = load("/work/00863/minerva/orion/bfield_ref_maptest_simpB2/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion/bfield_ref3_maptest1/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion/bfield_ref3_maptest1_2lev/" + 'data.' + datanum + '.3d.hdf5')
pf = load("/work/00863/minerva/orion_Vtest/" + 'data.' + datanum + '.3d.hdf5')

print pf.h.field_list

#value, location = pf.h.find_max("Density")
location = [0,0,0]

data = pf.h.sphere(location, 2.0/pf['pc'])

#L = [1.e-5,0,0]
#L = [0,1.e-5,0]
#L = [0,0,1.e-5]

L = [1.e-5,1.e-5,0]

radi = 0.5 
#radi = 5

sl = OffAxisSlicePlot(pf,L,'density',center=location,width=(2.0*radi,'cm'))
#sl = SlicePlot(pf, 'z', 'density', width = (5, 'pc'))
sl.annotate_velocity(factor=16)  # This will put black velocity lines on the plot 
sl.annotate_streamlines('CuttingPlaneBx', 'CuttingPlaneBy',plot_args={'color':'w'})  # This will put white b-field lines on the plot
sl.plots['density'].axes.set_xlabel('AU')
sl.plots['density'].axes.set_ylabel('AU')
sl.set_window_size(6.0)
#ln = 'FACE '+str(int(data['particle_id'][i]))
#ln2 = ' Mass %.3f'%(float(data['particle_mass'][i]/(1.989e33)))
#sl.annotate_title(ln+ln2)

sl.save('data_'+ datanum)
