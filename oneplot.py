import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as plt

from yt.mods import *
import numpy as na
from string import rstrip
import fields
import fields_wpot

datanum = '0001'

#Give me a file name
#pf  = load("/work/00863/minerva/orion/" + 'data.' + datanum + '.3d.hdf5.out.0.2pc')
pf = load("/work/00863/minerva/orion/gravpot_small/" + 'data.' + datanum + '.3d.hdf5')
#pf = load("/work/00863/minerva/orion/gravpot_small/" + 'aaron.' + datanum + '.3d.hdf5.out.0.2pc')


# Give me some output and titles
foutname = 'ScriptDensityTest.ps'
#plttitle = 'Density Using Script'
plttitle = 'log(dGmag_Gadget)'

# What direction should I make the slice? (0=x,1=y,2=z)
slicedirection = 0

# Give me a resolution
res = 128

# What field should I plot
#plotme = 'density'
plotme = 'dGmag_G2'

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

##THIS IS NOT ANOTHER FRB. IT IS ACTUALLY A NUMPY ARRAY WITH THE DENSITY DIFFERENCES
#frb1[plotme] = frb1[plotme] - frb1[plotme].min()
#frb2[plotme] = frb2[plotme] - frb2[plotme].min()
#frb3 = np.fabs(frb1[plotme]-frb2[plotme])
#frb3 = 2.0*np.fabs(frb1[plotme]-frb2[plotme])/np.fabs(frb1[plotme]+frb2[plotme])

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
    
##Plot each FRB and the difference, include a colorbar
fig=plt.figure()
plt.clf()
ax=plt.gca()

if(logit): image = plt.imshow(np.log10(frb1[plotme]))
else: image = plt.imshow(frb1[plotme])
plt.colorbar(image)
#mpl.colors.Normalize(vmin=-0.5, vmax=2.5)

ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xlabel('pc')
ax.set_ylabel('pc')
ax.set(xticklabels=xticks)
ax.set(yticklabels=yticks)
plt.title(plttitle)
plt.savefig(foutname)

