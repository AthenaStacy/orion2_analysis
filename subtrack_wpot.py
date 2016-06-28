import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as plt

from yt.mods import *
import numpy as na
from string import rstrip
import fields
import fields_wpot

datanum = '0000'

#Give me a file name
#pf  = load("/work/00863/minerva/orion/gravpot_0.8pc/" + 'data.' + datanum + '.3d.hdf5')
pf   = load("/work/00863/minerva/orion/gravpot_0.4pc/" + 'data.' + datanum + '.3d.hdf5')
pf2  = load("/work/00863/minerva/orion/gravpot_0.2pc/" + 'data.' + datanum + '.3d.hdf5')

#pf  = load("/work/00863/minerva/orion/" + 'data.' + datanum + '.3d.hdf5')
#pf2  = load("/work/00863/minerva/orion/" + 'data.' + datanum + '.3d.hdf5')

# Give me a resolution
res = 128

##########################################################################################################

# What field should I plot

potent = 2

#plotme1 = 'dGdx'
#plotme2 = 'dGdx_G2'

#plotme1 = 'tracer1'
#plotme2 = 'tracer1'

#plotme1 = 'gravitational-potential'
#plotme2 = 'gravitational-potential'

#plotme1 = 'dGmag_G2'
#plotme2 = 'dGmag_G2'

plotme1 = 'dGmag'
plotme2 = 'dGmag_G2'

half_boxsize = 3.08567758e17

value, location = pf.h.find_max("Density")
data = pf.h.sphere(location, half_boxsize/pf['pc'])
gpmin1 = min(data['gravitational-potential'])

value, location = pf2.h.find_max("Density")
data = pf2.h.sphere(location, half_boxsize/pf['pc'])
gpmin2 = min(data['gravitational-potential'])

# What direction should I make the slice? (0=x,1=y,2=z)
slicedirection = 0


if slicedirection==0:
	xdir = 1
	ydir = 2
elif slicedirection==1:
	xdir = 0
	ydir = 2
else:
	xdir=0
	ydir=1


c = 0.0*pf.domain_right_edge[slicedirection]
sl = pf.h.slice(slicedirection,c) ##Get the Slice
gpmin1 = min(sl[plotme1])
w = [-1*half_boxsize, half_boxsize, -1* half_boxsize, half_boxsize]
#w = [pf.domain_left_edge[xdir],pf.domain_right_edge[xdir],pf.domain_left_edge[ydir],pf.domain_right_edge[ydir]] #Choose the width, here I choose the entire domain
frb1 = FixedResolutionBuffer(sl,w,(res,res))  #Create FixedResolution Buffer

##Do the same (THIS CAN BE A DIFFERENT DATASET)
c = 0.0*pf2.domain_right_edge[slicedirection]
sl = pf2.h.slice(slicedirection,c)
gpmin2 = min(sl[plotme2])
w = [-1*half_boxsize, half_boxsize, -1* half_boxsize, half_boxsize]
#w = [pf2.domain_left_edge[xdir],pf2.domain_right_edge[xdir],pf2.domain_left_edge[ydir],pf2.domain_right_edge[ydir]]
frb2 = FixedResolutionBuffer(sl,w,(res,res))

#if we need to swith signs for the gravitational potential...
if(potent):
        frb1[plotme1] = frb1[plotme1] - gpmin1
	frb2[plotme2] = frb2[plotme2] - gpmin2
	print 'gpmin1 =', gpmin1, 'gpmin2 =', gpmin2
	
#print 'min =', min(frb1[plotme1])[0]

##THIS IS NOT ANOTHER FRB. IT IS ACTUALLY A NUMPY ARRAY WITH THE DENSITY DIFFERENCES
#frb3 = frb1[plotme1]-frb2[plotme2]
frb3 = (frb1[plotme1]-frb2[plotme2])/((frb1[plotme1]+frb2[plotme2])/2)

#Try ticks
ticks = np.arange(0,res+1,res/5.0)
#rad  = (pf.domain_right_edge[ydir]-pf.domain_left_edge[ydir])*pf['pc']
rad  = (2*half_boxsize)*pf['pc']
radt = np.linspace(-0.5,0.5,num=len(ticks))
radt = radt * (rad)
radt = np.array(radt,dtype='float')
for i in range(len(radt)):
    radt[i] = '%.2f'%(radt[i])
xticks = []
yticks = []

for i in range(len(radt)):
    xticks.append(radt[i])
    yticks.append(radt[len(radt)-i-1])
   
ymin = 0
#ymax = 5.e-4
ymax = 1.54e12 

##Plot each FRB and the difference, include a colorbar
fig=plt.figure()
plt.clf()
ax=plt.gca()
image = plt.imshow(frb1[plotme1])
plt.colorbar(image)
#plt.clim(ymin,ymax)  #color bar limits
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xlabel('pc')
ax.set_ylabel('pc')
ax.set(xticklabels=xticks)
ax.set(yticklabels=yticks)
plt.savefig('pf1_data.ps')

plt.clf()
ax=plt.gca()
image = plt.imshow(frb2[plotme2])
plt.colorbar(image)
#plt.clim(ymin,ymax)  #color bar limits
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xlabel('pc')
ax.set_ylabel('pc')
ax.set(xticklabels=xticks)
ax.set(yticklabels=yticks)
plt.savefig('pf2_data.ps')

plt.clf()
ax = plt.gca()
image = plt.imshow(frb3)
plt.colorbar(image)
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xlabel('pc')
ax.set_ylabel('pc')
ax.set(xticklabels=xticks)
ax.set(yticklabels=yticks)
plt.savefig('subtract.ps')
