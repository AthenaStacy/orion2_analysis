import matplotlib
matplotlib.use('ps')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

#from yt.mods import *
import numpy as na
import pylab as pl
#from string import rstrip
#import fields
import fields_bfield
from tracer_def import *

my_rho = YTQuantity(1, 'g / cm**3')
MRH = 0.76
pf_pc = 3e-18

###########################################################################
#read in sink data

rdir = '/work/00863/minerva/stampede2/popiii_bfieldA0/'
filenum = 30

num_sink = []
for f in range(filenum):

	datanum = 240 + f*10
	datanum = str(datanum)
	prefix = '0'
		
	with open(rdir + 'data.' + prefix + datanum + '.3d.sink', "r") as f:
		sinkdat = [line.split() for line in f]
		line_arr = sinkdat[0]
		num_sink_here = int(line_arr[0])
		num_sink.append(num_sink_here)

time_sink = []
with open('time_sink.txt', "r") as f:
	time_dat = [line.split(' ') for line in f]
for i in range(len(time_dat)):
	line = time_dat[i]
	time_sink.append(float(line[1]))

time_sink = np.asarray(time_sink)
time_sink = time_sink/3.14e7 - 8800.0

for i in range(len(time_dat)):
	print('time = ', time_sink[i], ' num = ', num_sink[i])	

sink_form_time = []
sink_y_val = []
for i in range(len(time_dat)-1):
	num_sink_1 = num_sink[i]
	num_sink_2 = num_sink[i+1]
	if num_sink_2 > num_sink_1:
		sink_form_time.append(time_sink[i])
		sink_y_val.append(1.0)
sink_form_time = np.asarray(sink_form_time)
sink_y_val = np.asarray(sink_y_val)

###########################################################################
# read in "max of max" data
time = []
mrat_max = []
mcrit_max = []
mphi_max = []
rho_crit = []
rho_max = []
temp_max = []
beta_max = []
with open('mrat.txt', "r") as f:
	mhd_dat = [line.split(' ') for line in f]
for i in range(len(mhd_dat)):
	line = mhd_dat[i]
	time.append(float(line[1]))
	mrat_max.append(float(line[2]))
	mcrit_max.append(float(line[3]))
	mphi_max.append(float(line[4]))
	rho_crit.append(float(line[5]))
	rho_max.append(float(line[6]))
	temp_max.append(float(line[7]))
	beta_max.append(float(line[8]))
time = np.asarray(time)
mrat_max = np.asarray(mrat_max)
mcrit_max = np.asarray(mcrit_max)
mphi_max = np.asarray(mphi_max)
rho_max = np.asarray(rho_max)
temp_max = np.asarray(temp_max)
time = time/3.14e7 - 9045.0 

time0 = []
mrat_max0 = []
mcrit_max0 = []
mphi_max0 = []
mphi_max0 = []
rho_crit0 = []
rho_max0 = []
temp_max0 = []
beta_max0 = []
with open('mrat0.txt', "r") as f:
        hydro_dat = [line.split(' ') for line in f]
for i in range(len(hydro_dat)):
	print('i = ', i)
	line = hydro_dat[i]
	time0.append(float(line[1]))
	mrat_max0.append(float(line[2]))
	mcrit_max0.append(float(line[3]))
	mphi_max0.append(float(line[4]))
	rho_crit0.append(float(line[5]))
	rho_max0.append(float(line[6]))
	temp_max0.append(float(line[7]))
	beta_max0.append(float(line[8]))
time0 = np.asarray(time0)
mrat_max0 = np.asarray(mrat_max0)
mcrit_max0 = np.asarray(mcrit_max0)
rho_max0 = np.asarray(rho_max0)
temp_max0 = np.asarray(temp_max0) 
time0 = time0/3.14e7 - 8795.0

#for i in range(len(temp_max)):
#	print('temp_max = ', temp_max[i], 'temp_max0 = ', temp_max0[i], 'ratio = ', temp_max[i]/temp_max0[i])



xline = [-8000, -4000, -2000, -1000, 0, 500, 1e3, 1500, 2000]
xline = np.asarray(xline)
yline = [1, 1, 1, 1, 1, 1, 1, 1, 1]
yline = np.asarray(yline)

fac = 1.0

tmin = -1000
tmax = 2000
mmin = 0.0
mmax = 3.

pl.gca()
pl.clf()
fsize = 14
pl.plot(time0, fac*mrat_max0,'b:', linewidth=2.5, label='(M/M$_\mathrm{crit}$)$_\mathrm{max}$ (hydro)')
pl.plot(time, fac*mrat_max,  'k', linewidth=2.0, label='(M/M$_\mathrm{BE}$)$_\mathrm{max}$ (MHD)')
pl.plot(time, fac*mcrit_max,'r--', linewidth=2.0, label='(M/M$_\mathrm{crit}$)$_\mathrm{max}$ (MHD)')
#pl.plot(time, fac*mphi_max,'r--', linewidth=2.0, label='M/M$_\mathrm{\Phi}$ (MHD)')
pl.plot(xline, yline,'g:', linewidth=2.0)
pl.plot(sink_form_time, sink_y_val, 'mo', label='sink formation (Hydro)')
ax = pl.gca()
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('time [yr]', fontsize=fsize)
ax.set_ylabel('(M/M$_\mathrm{crit}$)$_\mathrm{max}$', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.legend(loc='upper left', fontsize=10)
pl.yticks(fontsize=fsize)
pl.axis((tmin, tmax, mmin, mmax))

pl.savefig('mrat_vs_t.eps')

'''
pl.gca()
pl.clf()
fsize = 14
pl.plot(time, rho_max/rho_max0,'k', linewidth=2.0, label=r'max $\rho_\mathrm{MHD}$/$\rho_\mathrm{Hydro}$')
pl.plot(time, temp_max/temp_max0,  'b:', linewidth=2.0, label=r'max $T_\mathrm{MHD}$/$T_\mathrm{Hydro}$')
ax = pl.gca()
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('time [yr]', fontsize=fsize)
ax.set_ylabel(r'max T and $\rho$ ratios', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.legend(loc='upper left', fontsize=10)
pl.yticks(fontsize=fsize)
pl.axis((tmin, tmax, mmin, mmax))

pl.savefig('trat_vs_t.eps')
'''

mmin = 1.e-13
mmax = 1.e-10
pl.gca()
pl.clf()
fsize = 14
pl.plot(time, rho_max,'k', linewidth=2.0, label=r'MHD')
pl.plot(time0, rho_max0,  'b:', linewidth=2.0, label=r'Hydro')
for i in range(len(sink_form_time)):
	xline_here = [sink_form_time[i], sink_form_time[i], sink_form_time[i]]
	yline_here = [mmin, mmax/2, mmax]
	#pl.plot(xline_here, yline_here, 'm:')
ax = pl.gca()
#ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('time [yr]', fontsize=fsize)
ax.set_ylabel(r'$\rho_{\rm max}$ [g cm$^{-3}$]', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.legend(loc='upper left', fontsize=10)
pl.yticks(fontsize=fsize)
pl.axis((tmin, tmax, mmin, mmax))

pl.savefig('rho_max_vs_t.eps')

mmin = 50
mmax = 2000
pl.gca()
pl.clf()
fsize = 14
pl.plot(time, temp_max,'k', linewidth=2.0, label=r'MHD')
pl.plot(time0, temp_max0,  'b:', linewidth=2.0, label=r'Hydro')
ax = pl.gca()
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('time [yr]', fontsize=fsize)
ax.set_ylabel(r'max T', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.legend(loc='upper left', fontsize=10)
pl.yticks(fontsize=fsize)
pl.axis((tmin, tmax, mmin, mmax))

pl.savefig('temp_max_vs_t.eps')

mmin = 1.e-25
mmax = 1.e-19
pl.gca()
pl.clf()
fsize = 14
pl.plot(time, np.power(rho_max, 1.5) * np.power(temp_max,-1.5),'k', linewidth=2.0, label=r'MHD')
pl.plot(time0, np.power(rho_max0, 1.5) * np.power(temp_max0,-1.5),  'b:', linewidth=2.0, label=r'Hydro')
ax = pl.gca()
#ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('time [yr]', fontsize=fsize)
ax.set_ylabel(r'T$^{-1.5} \rho^{1.5}$', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.legend(loc='upper left', fontsize=10)
pl.yticks(fontsize=fsize)
pl.axis((tmin, tmax, mmin, mmax))

pl.savefig('check_max_vs_t.eps')

mmin = 0
mmax = 3
pl.gca()
pl.clf()
fsize = 14
pl.plot(time, rho_max/rho_crit,'k', linewidth=2.0, label=r'MHD')
pl.plot(time0, rho_max0/rho_crit0,  'b:', linewidth=2.0, label=r'Hydro')
pl.plot(sink_form_time, sink_y_val, 'mo', label='sink formation (Hydro)')
pl.plot(xline, yline,'g:', linewidth=2.0)
ax = pl.gca()
#ax.set_yscale('log')
ax.set_xlabel('time [yr]', fontsize=fsize)
ax.set_ylabel(r'max $\rho /\rho_{crit}$', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.legend(loc='upper left', fontsize=10)
pl.yticks(fontsize=fsize)
pl.axis((tmin, tmax, mmin, mmax))

pl.savefig('rho_rhomax_vs_t.eps')

mmin = 1.e-1
mmax = 1.e3
pl.gca()
pl.clf()
fsize = 14
pl.plot(time, beta_max,'k', linewidth=2.0, label=r'MHD')
ax = pl.gca()
ax.set_yscale('log')
ax.set_xlabel('time [yr]', fontsize=fsize)
ax.set_ylabel(r'$\beta_{\rm max}$', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.legend(loc='upper left', fontsize=10)
pl.yticks(fontsize=fsize)
pl.axis((tmin, tmax, mmin, mmax))

pl.savefig('betamax_vs_t.eps')

