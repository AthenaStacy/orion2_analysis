from yt.mods import *
import numpy as na
from string import rstrip
import fields
import fields_bfield
import matplotlib.colorbar as cb
from tracer_def import *

#rdir = "/work/00863/minerva/popiii_Bscope6/"
rdir = "/work/00863/minerva/popiii_Bscope10/"

#datanum = '0209'
datanum = '0188'
read_sink = 1

mmin = 1.e-7
mmax = 1.e1

h2min = 1.e-4
h2max = 1.e0

hpmin = 1.e-14
hpmax = 1.e-6

tmin = 1.e2
tmax = 1.e4

#pf = load("/work/00863/minerva/uni_vargam/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/popiii_Bscope/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/popiii_Bscope2/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/popiii_Bscope4/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/popiii_Bscope5/" + 'data.' + datanum + '.3d.hdf5')
pf = load(rdir + 'data.' + datanum + '.3d.hdf5')

#######################################################################
##Read in sink information

res = 16
nlev = 10
pc_to_cm = 3.08567758e18
rad  = 0.5*(pf.domain_right_edge[0]-pf.domain_left_edge[0])*pf['pc']
sink_size =  pf.domain_right_edge[0] / res / (2**(nlev-1))
sink_size = sink_size * 5 / pc_to_cm
print 'sink_size =', sink_size

mass_sink = []
xpos_sink = []
ypos_sink = []
zpos_sink = []
xplot_sink = []
yplot_sink = []
zplot_sink = []
xmom_sink = []
ymom_sink = []
zmom_sink = []
xangmom_sink = []
yangmom_sink = []
zangmom_sink = []
id_sink = []
masstot = []

if read_sink == 1:

	with open(rdir + 'data.' + datanum + '.3d.sink', "r") as f:
        	sinkdat = [map(float, line.split()) for line in f]

	line_arr = sinkdat[0]
	num_sink = int(line_arr[0])
	nmerge = int(line_arr[1])

	for i in range(0,num_sink):
        	line_arr = sinkdat[i+1]
        	mass_sink.append(line_arr[0])
        	xpos_sink.append(line_arr[1])
        	ypos_sink.append(line_arr[2])
        	zpos_sink.append(line_arr[3])
        	xplot_sink.append(line_arr[1]*res/pc_to_cm/(2.*rad) + res*0.5)
        	yplot_sink.append(line_arr[2]*res/pc_to_cm/(2.*rad) + res*0.5)
	        zplot_sink.append(line_arr[3]*res/pc_to_cm/(2.*rad) + res*0.5)
        	xmom_sink.append(line_arr[4])
        	ymom_sink.append(line_arr[5])
        	zmom_sink.append(line_arr[6])
        	xangmom_sink.append(line_arr[7])
        	yangmom_sink.append(line_arr[8])
        	zangmom_sink.append(line_arr[9])
        	id_sink.append(line_arr[10])
        	if(i == 0):
                	masstot.append(mass_sink[0])
        	if(i > 0):
                	masstot[0] = masstot[0] + mass_sink[i]

	print 'num_sink = ', num_sink
	print 'mass_sink = ', mass_sink

########################################################get B-field logs###############

def _Radius(field, data):
    '''
    Distance from central point. In this problem the center is at density peak.
    '''
    rad =  na.sqrt((data['x'] )**2 + (data['y'])**2 + (data['z'])**2)
    rad = rad/1.5e13
    return rad
add_field("Radius", function=_Radius, take_log=False,
          units=r'\rm{AU}')

def _H2_abund(field, data):
    abund = data[iH2] / data['dfac'] 
    return abund
add_field("H_2", function=_H2_abund, take_log=False, units='')

def _HP_abund(field, data):
    abund = data[iHP] / data['dfac']
    return abund
add_field("H^+", function=_HP_abund, take_log=False, units='')

######################################################################################

length = 10.0

my_sphere = pf.h.sphere("c", length/pf['pc'])


if read_sink == 1:
	center1 = [xpos_sink[0], ypos_sink[0], zpos_sink[0]]
	sink1 = pf.h.sphere(center1, sink_size/pf['pc'])
	my_sphere = pf.h.boolean([my_sphere, "NOT", sink1])

pc = PlotCollection(pf, 'c')

p=pc.add_phase_object(my_sphere, ["Density", "Temp.", "CellMassMsun"])
p.set_zlim(mmin, mmax)
p.set_ylim(tmin, tmax)

p=pc.add_phase_object(my_sphere, ["Density", igamma, "CellMassMsun"])
p.set_zlim(mmin, mmax)

p=pc.add_phase_object(my_sphere, ["Density", "H_2", "CellMassMsun"])
p.set_zlim(mmin, mmax)
p.set_ylim(h2min, h2max)

p=pc.add_phase_object(my_sphere, ["Density", "H^+", "CellMassMsun"])
p.set_zlim(mmin, mmax)
p.set_ylim(hpmin, hpmax)

pc.save("phase")


