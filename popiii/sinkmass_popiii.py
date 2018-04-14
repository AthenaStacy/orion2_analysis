#import matplotlib as pl
#matplotlib.use('ps')
import matplotlib.pyplot as pl
pl.style.use('seaborn-deep')

from yt.mods import *
import numpy as na
from string import rstrip
import fields
import fields_wpot
from sinkmass_g2 import *

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
			print 'field_list', pf.field_list	
			print 'derived_field_list', pf.derived_field_list	
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
				if id_here in mdict.keys():
					mdict[id_here].append(mass_sink[j] / 1.989e33)
					tdict[id_here].append(tacc_sink[j] / 3.14e7) 	
				else:
					mdict[id_here] = []
					mdict[id_here].append(mass_sink[j] / 1.989e33)
					tdict[id_here] = []
					tdict[id_here].append(tacc_sink[j] / 3.14e7)	


	tacc_yr = []
	masstot_msol = []
	for i in range(len(masstot)):
		tacc_yr.append(tacc[i] / 3.14e7)
		masstot_msol.append(masstot[i] / 1.989e33)

	return tacc_yr, masstot_msol, num_sink, mass_fin, id_fin, id_sink, tdict, mdict


rdir = "/nobackupp7/astacy/popiii_bfieldA0/"
filestart = 241
fileEnd = 620
filenum = fileEnd - filestart
tacc_yr, masstot_msol, num_sink, mass_fin, id_fin, id_sink, tdict, mdict = read_sink_file(rdir, filestart, fileEnd, filenum)
#sort (id, mass) array in order of descending mass
id_fin.sort(key=lambda x: x[1], reverse=True)
print ('id_fin =', id_fin)

rdir = "/nobackupp7/astacy/popiii_bfieldA/"
filestart = 246
fileEnd = 570
filenum = fileEnd - filestart
tacc_yr2, masstot_msol2, num_sink2, mass_fin2, id_fin2, id_sink2, tdict2, mdict2 = read_sink_file(rdir, filestart, fileEnd, filenum)
id_fin2.sort(key=lambda x: x[1], reverse=True)

mdot = []
mdot2 = []
tdot = []
tdot2 = []
for i in range(len(masstot_msol)-1):
	mdot_here = (masstot_msol[i+1] - masstot_msol[i]) / (tacc_yr[i+1] - tacc_yr[i])
	if (mdot_here > 0):
		mdot.append(mdot_here)
		tdot.append( (tacc_yr[i] + tacc_yr[i+1])/2)	

for i in range(len(masstot_msol2)-1):
        mdot_here = (masstot_msol2[i+1] - masstot_msol2[i]) / (tacc_yr2[i+1] - tacc_yr2[i])
	if (mdot_here > 0):
	        mdot2.append(mdot_here)
        	tdot2.append( (tacc_yr2[i] + tacc_yr2[i+1])/2)

msink_g2, mtot_g2, nsink_g2, t_msink, t_nsink = read_sink()

mdot_min = 1.e-3
mdot_max = 2.e-1

mmin = 0
mmax = 1.1*max(masstot_msol)

nmin = 0
nmax = 1.1*max(num_sink)

tmin = 0
tmax= 1.1*max(max(tacc_yr),max(tacc_yr2))

fsize = 14

pl.plot(tacc_yr, masstot_msol, 'k', label='total sink mass')
pl.plot(tdict[id_fin[0][0]], mdict[id_fin[0][0]], label='initial sink')
#pl.plot(tdict[id_fin[1][0]], mdict[id_fin[1][0]], 'r', label='2nd-most massive sink')
#pl.plot(tdict[id_fin[2][0]], mdict[id_fin[2][0]],'g', label='3rd-most massive sink')
pl.plot(t_msink, msink_g2, 'b:', label='largest sink from Gadget2')
pl.plot(t_nsink, mtot_g2, 'k--', label='total sink mass from Gadget2')
ax = pl.gca()
ax.set_xlabel('time [yr]', fontsize=fsize)
ax.set_ylabel(r'mass [M$_{\odot}$]', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((tmin, tmax, mmin, mmax))
pl.legend(fontsize=10, loc=2)
pl.title('Hydro')
pl.savefig('sinkmass_noB.eps')
pl.clf()

pl.plot(tacc_yr2, masstot_msol2, 'k', label='total sink mass')
pl.plot(tdict2[id_fin2[0][0]], mdict2[id_fin2[0][0]], label='initial sink')
#pl.plot(tdict2[id_fin2[1][0]], mdict2[id_fin2[1][0]], 'r', label='2nd-most massive sink')
#pl.plot(tdict2[id_fin2[2][0]], mdict2[id_fin2[2][0]],'g', label='3rd-most massive sink')
ax = pl.gca()
ax.set_xlabel('time [yr]', fontsize=fsize)
ax.set_ylabel(r'mass [M$_{\odot}$]', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((tmin, tmax, mmin, mmax))
pl.legend(fontsize=10, loc=2)
pl.title('MHD')
pl.savefig('sinkmass.eps')
pl.clf()

pl.plot(tacc_yr, num_sink, label='hydro')
pl.plot(tacc_yr2, num_sink2, 'k', label='MHD')
pl.plot(t_nsink, nsink_g2, 'k:', label='Gadget2')
ax = pl.gca()
ax.set_xlabel('time [yr]', fontsize=fsize)
ax.set_ylabel('number of sinks', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((tmin, tmax, nmin, nmax))
pl.legend(fontsize=10, loc=2)
pl.savefig('sinknum.eps')
pl.clf()

pl.plot(tdot, mdot, label='hydro')
pl.plot(tdot2, mdot2, 'k', label='MHD')
ax = pl.gca()
ax.set_yscale('log')
ax.set_xlabel('time [yr]', fontsize=fsize)
ax.set_ylabel('dM$_*$/dt [M$_{\odot}$ yr$^{-1}$]', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((tmin, tmax, mdot_min, mdot_max))
pl.legend(fontsize=10, loc=1)
pl.savefig('sinkmdot.eps')
pl.clf()

print ('mass_fin =', mass_fin)
print ('mass_fin2 =', mass_fin2)
bins = numpy.linspace(0, 30, 50)
pl.rcParams["patch.force_edgecolor"] = True
pl.hist(mass_fin,  bins, alpha=0.2, label='hydro', linestyle='dashed', linewidth=2.0)
pl.hist(mass_fin2, bins, alpha=0.2, label='MHD', linestyle='solid')
ax = pl.gca()
pl.xlabel(r'mass [M$_{\odot}$]', fontsize=fsize)
pl.ylabel('number of sinks', fontsize=fsize)
ax.set_ylim([0,7])
pl.legend(fontsize=10, loc=0)
pl.savefig('sinkhist.eps')
pl.clf()

'''
pl.hist(mass_fin2)
ax = pl.gca()
pl.xlabel(r'mass [M$_{\odot}$]', fontsize=fsize)
pl.ylabel('number of sinks', fontsize=fsize)
ax.set_ylim([0,10])
pl.savefig('sinkhist.eps')
pl.clf()
'''


