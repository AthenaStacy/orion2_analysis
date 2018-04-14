import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as pl

from yt.mods import *
import numpy as na
from string import rstrip
import fields
import fields_wpot

def read_sink_file(rdir, filestart, fileEnd, filenum):
	#pass
	pc_to_cm = 3.08567758e18

	##################################################
	#define some lists
	masstot = []
	num_sink = []
	tacc = []
	time = []

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
	id_sink_final = []
	tacc_sink = []
	id_sink = []

	###################################################
	##Read in time and sink information

	indnum = 0
	for fnum in range(filestart, filestart+filenum):

		mass_sink_here = []
		id_sink_here = []

		prefix = '000'
		if(fnum > 9):
			prefix = '00'
		if(fnum > 99):
			prefix = '0'
		if(fnum > 999):
        		prefix = ''

		datanum = str(fnum)
		datanum2 = str(fnum - fnum%10 + 10)

		if fnum == filestart:
			datanum1 = str(fnum - fnum%10)
			pf = load(rdir + 'data.' + prefix + datanum1 + '.3d.hdf5')
			pf2 = load(rdir + 'data.' + prefix + datanum2 + '.3d.hdf5')
			t1 = pf.current_time
			t2 = pf2.current_time
		
                if fnum%10 == 0 and fnum != filestart:
                        pf = load(rdir + 'data.' + prefix + datanum2 + '.3d.hdf5')
                        t1 = t2
                        t2 = pf.current_time

		if fnum%10 == 0:
			t_current = t1
		else:  
			incrememt = (t2-t1)/10.0
			increment_num = 10.0 - float(fnum%10)
			t_current = t2 - increment_num * (t2-t1) / 10.0


		# Give me a resolution
		res = pf.domain_dimensions[0]
		res_doub = float(res)

		rad  = 0.5*(pf.domain_right_edge[0]-pf.domain_left_edge[0])/pc_to_cm
		rad = float(rad)

		with open(rdir + 'data.' + prefix + datanum + '.3d.sink', "r") as f:
        		sinkdat = [map(float, line.split()) for line in f]

		num_sink_current = int(sinkdat[0][0])
		num_sink.append(num_sink_current)
		nmerge = int(sinkdat[0][1])

		for i in range(1,num_sink_current+1):
        		line_arr = sinkdat[i]
        		mass_sink.append(line_arr[0])
			mass_sink_here.append(line_arr[0])
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
			if line_arr[10] == 3:
				line_arr[10] = 1
        		id_sink.append(line_arr[10])
                	id_sink_here.append(line_arr[10])
			if(i == 1):
				masstot.append(mass_sink_here[0])
				#time.append(pf.current_time)
				time.append(t_current)
				tacc.append(time[indnum] - time[0])
				indnum = indnum + 1
			tacc_sink.append((time[indnum-1] - time[0]))
 			if(i > 1):
				masstot[indnum-1] = masstot[indnum-1] + mass_sink_here[i-1]
			if(fnum == filestart+filenum-1 and i >= 1):
				id_sink_final.append(id_sink_here[i-1])
				print ("id_sink_final = ", id_sink_final)

	########################################################    

	final_sink_num = len(id_sink_final)
	mdict = {}
	tdict = {}
	vxdict = {}
	vydict = {}
	vzdict = {}
	mass_fin = []
	id_fin = []

	print ("len mass_sink = ", len(mass_sink))
	print ("len tacc_sink = ", len(tacc_sink))
	print ('len final_sink_num =', final_sink_num)

        for i in range(len(mass_sink)-1, len(mass_sink)-final_sink_num-1, -1):
                mass_fin.append(mass_sink[i]/1.989e33)
                id_fin.append((id_sink[i], mass_sink[i]))

	for i in range(final_sink_num):
		id_here = id_sink_final[i]
		for j in range(len(mass_sink)):
			if id_here == id_sink[j]:
				if id_here not in mdict.keys():
					mdict[id_here] = []
					tdict[id_here] = []
					vxdict[id_here] = []
					vydict[id_here] = []
					vzdict[id_here] = []
				mdict[id_here].append(mass_sink[j] / 1.989e33)
				tdict[id_here].append(tacc_sink[j] / 3.14e7)
				vxdict[id_here].append(xmom_sink[j] / mass_sink[j])
				vydict[id_here].append(ymom_sink[j] / mass_sink[j])
				vzdict[id_here].append(zmom_sink[j] / mass_sink[j]) 	


	tacc_yr = []
	masstot_msol = []
	for i in range(len(masstot)):
		tacc_yr.append(tacc[i] / 3.14e7)
		masstot_msol.append(masstot[i] / 1.989e33)

	return tacc_yr, id_fin, tdict, mdict, vxdict, vydict, vzdict

def vel_calc(vxdict, vydict, vzdict, id_fin, n):
	length1 = len(tdict[id_fin[0][0]])
	length2 = len(tdict[id_fin[n][0]])
	vx1 = []
	vx2  = []
	vy1 = []
	vy2  = []
	vz1 = []
	vz2 = []
	v = []
	for i in range(length2):
        	vx1.append(vxdict[id_fin[0][0]][i + length1-length2])
        	vx2.append(vxdict[id_fin[n][0]][i])
        	vy1.append(vydict[id_fin[0][0]][i + length1-length2])
        	vy2.append(vydict[id_fin[n][0]][i])
        	vz1.append(vzdict[id_fin[0][0]][i + length1-length2])
        	vz2.append(vzdict[id_fin[n][0]][i])
        	vel = (vx1[i]-vx2[i])*(vx1[i]-vx2[i]) + (vy1[i]-vy2[i])*(vy1[i]-vy2[i]) + (vz1[i]-vz2[i])*(vz1[i]-vz2[i])
        	v.append(vel**0.5 / 1e5)
	return v




type = 'hydro'
#type = 'mhd'
if type == 'hydro':
	rdir = "/nobackupp7/astacy/popiii_bfieldA0/"
	filestart = 244
	fileEnd = 326
	ptitle = 'Hydro'
else:
	rdir = "/nobackupp7/astacy/popiii_bfieldA/"
	filestart = 240
	fileEnd = 460
	ptitle = 'MHD'

filenum = fileEnd - filestart

tacc_yr, id_fin, tdict, mdict, vxdict, vydict, vzdict  = read_sink_file(rdir, filestart, fileEnd, filenum)

id_fin.sort(key=lambda x: x[1], reverse=True)

v00 = vel_calc(vxdict, vydict, vzdict, id_fin, 0)
v01 = vel_calc(vxdict, vydict, vzdict, id_fin, 1)
v02 = vel_calc(vxdict, vydict, vzdict, id_fin, 2)

#now just add initial distances of all those other sinks
v_other = []
t_other = []
for i in range(3,len(id_fin)):
	v = vel_calc(xdict, ydict, zdict, id_fin, i)
	v_other.append(r[0])
	t_other.append(tdict[id_fin[i][0]][0])

print 'v_other = ', v_other
print 't_other = ', t_other
print 'id_fin = ', id_fin

vmin = 0
vmax = 1.1*max([max(v01),max(v02)])

tmin = 0
tmax= 1.1*max(tacc_yr)

fsize = 14

pl.plot(tdict[id_fin[0][0]], v00)
pl.plot(tdict[id_fin[1][0]], v01,'b', label='2nd-most massive sink')
pl.plot(tdict[id_fin[2][0]], v02,'r', label='3rd-most massive sink')
pl.plot(t_other, v_other, 'm.', label='other sinks') 
ax = pl.gca()
ax.set_xlabel('time [yr]', fontsize=fsize)
ax.set_ylabel('velocity relative to main sink [km / s]', fontsize=fsize)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((tmin, tmax, vmin, vmax))
pl.legend(fontsize=10, loc=2)
pl.title(ptitle)
pl.savefig('sinkvel_.eps')
pl.clf()


