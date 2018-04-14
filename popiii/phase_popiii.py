import yt
from yt.units import dimensions
#from yt.mods import *
from astropy import units as u
import numpy as na
from string import rstrip
import fields_bfield
import matplotlib.colorbar as cb
from tracer_def import *


rdir = "/nobackupp7/astacy/popiii_bfieldA/"
#rdir = "/nobackupp7/astacy/popiii_bfieldA0/"
#rdir = "/nobackupp7/astacy/popiii_chemtest/"

datanum = '0500'
#datanum = '0250' #right around when first sink forms
#datanum = '0407'
read_sink = 0

dmin = 1.e-21
dmax = 1.e-9

mmin = 1.e-7
mmax = 1.e1

h2min = 1.e-4
h2max = 1.e0

hdmin = 1.e-11
hdmax = 1.e-4

hpmin = 1.e-11
hpmax = 1.e-5

gmin = 1.35
gmax = 1.7

tmin = 1.e2
tmax = 1.e4

bmin = 1.e-20
bmax = 1.e-5

vmin = 1
vmax = 50

tq_min = 1.e8
tq_max = 1.e12

#pf = load("/work/00863/minerva/uni_vargam/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/popiii_Bscope/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/popiii_Bscope2/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/popiii_Bscope4/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/popiii_Bscope5/" + 'data.' + datanum + '.3d.hdf5')
pf = load(rdir + 'data.' + datanum + '.3d.hdf5')

#for i in sorted(pf.field_list):
#  print(i)
#for i in sorted(pf.derived_field_list):
#  print(i)

ad = pf.all_data()
#print 'quantities!', ad.quantities.keys()

pot_data = ad.quantities.min_location('gravitational-potential')
min_pot = pot_data[0] 
x = pot_data[1] 
y = pot_data[2] 
z = pot_data[3]
print pot_data

string_to_write = str(min_pot.value) + ' ' + str(x.value) + ' ' + str(y.value) + ' ' + str(z.value)
f = open('gpot.dat', 'w')
f.write(string_to_write)
f.close()

#######################################################################
##Read in sink information

res = 16
nlev = 10
pc_to_cm = 3.08567758e18
rad  = 0.5*(pf.domain_right_edge[0]-pf.domain_left_edge[0])/pc_to_cm
rad = float(rad)
sink_size =  pf.domain_right_edge[0] / res / (2**(nlev-1))
sink_size = sink_size * 5 / pc_to_cm
print ('sink_size =', sink_size)

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

	print ('num_sink = ', num_sink)
	print ('mass_sink = ', mass_sink)


######################################################################################

length = 10.0

#my_sphere = pf.h.sphere("c", length/pf['pc'])
my_sphere = pf.all_data()

#print ("min gamma = ", min(my_sphere['gamma']), "max gamma= ", max(my_sphere['gamma'])) 
#print ("min h2 = ", min(my_sphere['H_2']), "max h2= ", max(my_sphere['H_2']))
#print ("min temp = ", min(my_sphere['Temp']), "max temp= ", max(my_sphere['Temp']))
#print ("min vel = ", min(my_sphere['velocity']), "max vel= ", max(my_sphere['velocity']))

#print ("min Bmag =", min(my_sphere['Bmag']))
#print ("max Bmag =", max(my_sphere['Bmag']))
#print ("min Bmag-log =", min(my_sphere['Bmag-log']))
#print ("max Bmag-log =", max(my_sphere['Bmag-log']))


print ("min grav_acc = ", min(my_sphere['grav_acc']), "max grav_acc = ", max(my_sphere['grav_acc'])) 
print ("min torque_time = ", min(my_sphere['torque_time']), "max torque_time = ", max(my_sphere['torque_time']))

if read_sink == 1:
	center1 = [xpos_sink[0], ypos_sink[0], zpos_sink[0]]
	sink1 = pf.h.sphere(center1, sink_size/pf['pc'])
	#my_sphere = pf.h.boolean([my_sphere, "NOT", sink1])

#pc = PlotCollection(pf, 'c')

fsize = 20

p = PhasePlot(my_sphere, "density", "torque_time", ["cell_mass"], weight_field=None)
p.set_ylim(tq_min, tq_max)
p.set_xlim(dmin, dmax)
p.set_font({'size': fsize})
p.save('phase_ttime_t'+datanum+'.png')

#p=pc.add_phase_object(my_sphere, ["density", "Temp", "CellMassMsun"])
p = PhasePlot(my_sphere, "density", "Temp", ["cell_mass"], weight_field=None)
p.set_ylim(tmin, tmax)
p.set_xlim(dmin, dmax)
p.set_font({'size': fsize})
p.save('phase_temp_t'+datanum+'.png')

p = PhasePlot(my_sphere, "density", "HD", ["cell_mass"], weight_field=None)
p.set_ylim(hdmin, hdmax)
p.set_xlim(dmin, dmax)
p.set_font({'size': fsize})
p.save('phase_hd_t'+datanum+'.png')

p = PhasePlot(my_sphere, "density", "H_2", ["cell_mass"], weight_field=None)
p.set_ylim(h2min, h2max)
p.set_xlim(dmin, dmax)
p.set_font({'size': fsize})
p.save('phase_h2_t'+datanum+'.png')

p = PhasePlot(my_sphere, "density", "H^+", ["cell_mass"], weight_field=None)
p.set_ylim(hpmin, hpmax)
p.set_xlim(dmin, dmax)
p.set_font({'size': fsize})
p.save('phase_hp_t'+datanum+'.png')

p = PhasePlot(my_sphere, "density", "gamma", ["cell_mass"], weight_field=None)
p.set_ylim(gmin, gmax)
p.set_xlim(dmin, dmax)
p.set_log("gamma", False)
p.set_font({'size': fsize})
p.save('phase_gamma_t'+datanum+'.png')

p = PhasePlot(my_sphere, "density", "velocity", ["cell_mass"], weight_field=None)
p.set_ylim(vmin, vmax)
p.set_xlim(dmin, dmax)
p.set_font({'size': fsize})
p.save('phase_vel_t'+datanum+'.png')

fac= 1e9
p = PhasePlot(my_sphere, "density", "Bmag", ["cell_mass"], weight_field=None)
#p.set_zlim(mmin, mmax)
p.set_ylim(bmin*fac, bmax*fac)
p.set_xlim(dmin, dmax)
x_init = 1.e-20
x = []
y = []
y2 = []
for i in range (0,12):
        x_here = x_init * 10**i
        x.append(x_here)        
        y.append(1.e-16 * (x_here/x_init) ** 0.6666)
        y2.append(1.e-16 * (x_here/x_init) ** 0.9)

print ('xline =', x)
print ('yline =', y)
print ('y2line =', y2)
#p.annotate_line((x[0], y[0]), (x[8], y[8]), coord_system='axis')
#p.annotate_line((x[0], y2[0]), (x[8], y2[8]), coord_system='axis')
p.set_font({'size': fsize})
p.save('phase_bmag_t'+datanum+'.png')

