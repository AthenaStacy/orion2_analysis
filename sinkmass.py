import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as pl

from yt.mods import *
import numpy as na
from string import rstrip
import fields
import fields_wpot

#rdir = "/work/00863/minerva/orion/bfield_ref2_maptest1_2lev/"
#rdir = "/work/00863/minerva/orion/bfield_ref2_maptest1_2lev_ideal/"
rdir = "/work/00863/minerva/orion/bfield_ref3_maptest1/"

filenum = 20

pc_to_cm = 3.08567758e18

# Give me some output and titles
foutname = 'sinkmass.eps'
plttitle = 'Sink Mass'

masstot = []
tacc = []
time = []

#######################################################################
##Read in time and sink information

indnum = 0
for f in range(filenum):

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

	datanum = str(f)

	prefix = '000'
	if(f > 9):
		prefix = '00'
	if(f > 99):
		prefix = '0'
	if(f > 999):
        	prefix = ''
	
	#Give me a file name
	pf = load(rdir + 'data.' + prefix + datanum + '.3d.hdf5')

	# Give me a resolution
	res = pf.domain_dimensions[0]
	res_doub = float(res)

	rad  = 0.5*(pf.domain_right_edge[0]-pf.domain_left_edge[0])*pf['pc']

	with open(rdir + 'data.' + prefix + datanum + '.3d.sink', "r") as f:
        	sinkdat = [map(float, line.split()) for line in f]

	line_arr = sinkdat[0]
	num_sink = int(line_arr[0])
	nmerge = int(line_arr[1])

	for i in range(1,num_sink):
        	line_arr = sinkdat[i]
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
		if(i == 1):
			masstot.append(mass_sink[0])
			time.append(pf.current_time)
			tacc.append(time[indnum] - time[0])
			indnum = indnum + 1
 		if(i > 1):
			masstot[indnum-1] = masstot[indnum-1] + mass_sink[i-1]

	#print "mass_sink =", mass_sink
	#print "xpos_sink =", xpos_sink
	#print "ypos_sink =", ypos_sink
	#print "xplot_sink =", xplot_sink
	#print "yplot_sink =", yplot_sink
	print "tacc =", tacc
	print "num_sink =", num_sink
	print "total sink mass =", masstot
######################################################################    

tacc_yr = []
masstot_msol = []
for i in range(0,indnum-1):
	tacc_yr.append(tacc[i] / 3.14e7)
	masstot_msol.append(masstot[i] / 1.989e33)

mmin = 0
mmax = 1.e3

tmin = 0
tmax= 6.e4

pl.plot(tacc_yr, masstot_msol, 'k')
ax = pl.gca()
ax.set_xlabel('time [yr]', fontsize=9)
ax.set_ylabel('Mass [M_sun]', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((tmin, tmax, mmin, mmax))

pl.savefig('sinkmass.eps')
