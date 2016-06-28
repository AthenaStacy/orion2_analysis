import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as plt

from yt.mods import *
import numpy as na
from string import rstrip
import fields
import fields_wpot

#rdir = "/work/00863/minerva/orion/bfield_ref2_maptest1_2lev/"
#rdir = "/work/00863/minerva/orion/bfield_ref2_maptest1_2lev_ideal/"
rdir = "/work/00863/minerva/orion/bfield_ref3_maptest1/"
datanum = '0230'

pc_to_cm = 3.08567758e18

#Give me a file name
pf = load(rdir + 'data.' + datanum + '.3d.hdf5')


# Give me some output and titles
foutname = 'sinkplot.ps'
plttitle = 'Density + sinks'
#plttitle = 'log(dGmag_Gadget)'

# What direction should I make the slice? (0=x,1=y,2=z)
slicedirection = 2

# Give me a resolution
res = 128
res_doub = 128.0

# What field should I plot
plotme = 'density'
#plotme = 'dGmag_G2'

# Log?
logit = 1


if slicedirection==0:
	xdir = 1
	ydir = 2
elif slicedirection==1:
	xdir = 0
	ydir = 2
else:
	xdir=0
	ydir=1
	slicedirection=2
	


# 0 is the center of the box!!!
c = 0.0*pf.domain_right_edge[slicedirection]
sl = pf.h.slice(slicedirection,c) ##Get the Slice
w = [pf.domain_left_edge[xdir],pf.domain_right_edge[xdir],pf.domain_left_edge[ydir],pf.domain_right_edge[ydir]] #Choose the width, here I choose the entire domain
frb1 = FixedResolutionBuffer(sl,w,(res,res))  #Create FixedResolution Buffer
#frb1[plotme] = frb1[plotme] - frb1[plotme].min()


#Try ticks
ticks = np.arange(0,res+1,res/5.0)
rad  = 0.5*(pf.domain_right_edge[ydir]-pf.domain_left_edge[ydir])*pf['pc']
radt = np.linspace(-1,1,num=len(ticks))
radt = radt * (rad)
radt = np.array(radt,dtype='float')
for i in range(len(radt)):
    radt[i] = '%.2f'%(radt[i])
xticks = []
yticks = []

for i in range(len(radt)):
    xticks.append(radt[i])
    yticks.append(radt[len(radt)-i-1])

#######################################################################
##Read in sink information
with open(rdir + 'data.' + datanum + '.3d.sink', "r") as f:
        sinkdat = [map(float, line.split()) for line in f]

line_arr = sinkdat[0]
num_sink = int(line_arr[0])
nmerge = int(line_arr[1])

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
 	if(i > 1):
		masstot[0] = masstot[0] + mass_sink[i-1]

#print "mass_sink =", mass_sink
#print "xpos_sink =", xpos_sink
#print "ypos_sink =", ypos_sink
#print "xplot_sink =", xplot_sink
#print "yplot_sink =", yplot_sink
print "total sink mass =", masstot[0]
######################################################################    

##Plot each FRB and the difference, include a colorbar
fig=plt.figure()
plt.clf()
ax=plt.gca()

if(logit): image = plt.imshow(np.log10(frb1[plotme]))
else: image = plt.imshow(frb1[plotme])
plt.colorbar(image)
ax.autoscale(False)
ax.plot(xplot_sink, yplot_sink,'k.',)
#ax.plot([64,64],[64,64], 'g')
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xlabel('pc')
ax.set_ylabel('pc')
ax.set(xticklabels=xticks)
ax.set(yticklabels=yticks)
plt.title(plttitle)
plt.savefig(foutname)

value, location = pf.h.find_max("Density")
data = pf.h.sphere(location, 1.e0/pf['pc'])
pc=PlotCollection(pf, center=[0,0,0])
p=pc.add_slice('density', 'z')
plt.plot(xplot_sink, yplot_sink,'k.')
ax.plot(range(11))
pc.save('data')

#fig.show()

