import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as plt

from yt.mods import *
import numpy as na
from string import rstrip
import fields

datanum1 = '0021'
datanum2 = '0000'

#Give me a file name
pf  = load("/nobackupp7/astacy/orion/bfield_comp1/" + 'data.' + datanum1 + '.3d.hdf5')
pf2 = load("/nobackupp7/astacy/orion/bfield_comp2/" + 'data.' + datanum2 + '.3d.hdf5')


# Give me a resolution
res = 128

# What field should I plot
plotme = 'vortx'

logit = 1

# What direction should I make the slice? (0=x,1=y,2=z)
slicedirection = 1


if slicedirection==0:
	xdir = 1
	ydir = 2
elif slicedirection==1:
	xdir = 0
	ydir = 2
else:
	xdir=0
	ydir=1
	

dd = pf.h.all_data()
dd2 = pf.h.all_data()




c = 0.0*pf.domain_right_edge[slicedirection]
sl = pf.h.slice(slicedirection,c) ##Get the Slice
w = [pf.domain_left_edge[xdir],pf.domain_right_edge[xdir],pf.domain_left_edge[ydir],pf.domain_right_edge[ydir]] #Choose the width, here I choose the entire domain
frb1 = FixedResolutionBuffer(sl,w,(res,res))  #Create FixedResolution Buffer

##Do the same (THIS CAN BE A DIFFERENT DATASET)
c = 0.0*pf2.domain_right_edge[slicedirection]
sl = pf2.h.slice(slicedirection,c)
w = [pf2.domain_left_edge[xdir],pf2.domain_right_edge[xdir],pf2.domain_left_edge[ydir],pf2.domain_right_edge[ydir]]
frb2 = FixedResolutionBuffer(sl,w,(res,res))

##THIS IS NOT ANOTHER FRB. IT IS ACTUALLY A NUMPY ARRAY WITH THE DENSITY DIFFERENCES
frb3 = (frb1[plotme]-frb2[plotme])/((frb1[plotme]+frb2[plotme])/2.0)
#frb3 = np.fabs((frb1[plotme] - frb2[plotme])/frb1[plotme])

print 'frb1 =', frb1[plotme]
print 'frb2 =', frb2[plotme]

frb3sq = frb3*frb3
print 'frb3 =', frb3
print 'length =', len(frb3)
print 'frb3sq =', frb3sq
print 'length =', len(frb3sq)

#Calculate rms error
ms = np.mean(frb3sq)
print 'mean = ', ms
print 'rms error =', np.sqrt(ms) 

#Calculate mass-weighted median value
val1 = frb1[plotme] 
val2 = frb2[plotme]
high = []
low = []
all = []
median = np.median(frb1[plotme])
print 'median val1 =', median
print 'length val1 =', len(val1)

for i in range(len(val1)):
    for j in range(len(val1)):
        all.append((val1[i,j] - val2[i,j])/((val1[i,j]+val2[i,j])/2))
        if val1[i,j] < median:
            low.append((val1[i,j] - val2[i,j])/((val1[i,j]+val2[i,j])/2))
        if val1[i,j] > median:
            high.append((val1[i,j] - val2[i,j])/((val1[i,j]+val2[i,j])/2))

all = np.array(all)
print 'len =', len(all)
error_all_sq = all*all
ms = np.mean(error_all_sq)
print 'mean error all = ', ms
print 'rms error all =', np.sqrt(ms)

high = np.array(high)
print 'len =', len(high)
error_high_sq = high*high
ms = np.mean(error_high_sq)
print 'mean error high = ', ms
print 'rms error high =', np.sqrt(ms)

low = np.array(low)
print 'len =', len(low)
error_low_sq = low*low
ms = np.mean(error_low_sq)
print 'mean error low = ', ms
print 'rms error low =', np.sqrt(ms)

#Try ticks
ticks = np.arange(0,res+1,res/5.0)
rad  = (pf.domain_right_edge[ydir]-pf.domain_left_edge[ydir])*pf['pc']
print 'rad=', rad
radt = np.linspace(-0.5,0.5,num=len(ticks))
print 'radt=', radt
radt = radt * (rad)

radt = np.array(radt,dtype='float')
print 'radt=', radt
for i in range(len(radt)):
    radt[i] = '%.2f'%(radt[i])
xticks = []
yticks = []

for i in range(len(radt)):
    xticks.append(radt[i])
    yticks.append(radt[len(radt)-i-1])

print 'radt=', radt
print 'xticks=', xticks
    
##Plot each FRB and the difference, include a colorbar
fig=plt.figure()
plt.clf()
ax=plt.gca()
if(logit): image = plt.imshow(np.log10(frb1[plotme]))
else: image = plt.imshow(frb1[plotme])
plt.colorbar(image)
#plt.clim(-2,0)  #color bar limits
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xlabel('pc')
ax.set_ylabel('pc')
ax.set(xticklabels=xticks)
ax.set(yticklabels=yticks)
plt.savefig('pf1_data.ps')

plt.clf()
ax=plt.gca()
if(logit): image = plt.imshow(np.log10(frb2[plotme]))
else: image = plt.imshow(frb2[plotme])
plt.colorbar(image)
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xlabel('pc')
ax.set_ylabel('pc')
ax.set(xticklabels=xticks)
ax.set(yticklabels=yticks)
plt.savefig('pf2_data.ps')

err_low  = -.2
err_high = .2
plt.clf()
ax = plt.gca()
image = plt.imshow(frb3)
plt.colorbar(image)
plt.clim(err_low,err_high)  #color bar limits
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xlabel('pc')
ax.set_ylabel('pc')
ax.set(xticklabels=xticks)
ax.set(yticklabels=yticks)
plt.savefig('subtract.ps')
