import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as plt

from yt.mods import *
import numpy as np
from string import rstrip
import fields
import fields_wpot
import sys

datanum = '0000'

fname = "/work/00863/minerva/orion/" + 'data.' + datanum + '.3d.hdf5'


plt2=1

if(plt2):
   fname2 = "/work/00863/minerva/orion/" + 'data.' + datanum + '.3d.hdf5'

sname = 'TestLine'


#Load data
pf = load(fname)

if(plt2): pf2= load(fname2)

# Chose center of data (selecting highest density point)
c = pf.h.find_max('Density')[1]
c2 = pf2.h.find_max('Density')[1]

gpmin = pf.h.find_min('gravitational-potential')[0]
gpmax = pf.h.find_max('gravitational-potential')[0]
if(plt2):
   gpmin2 = pf.h.find_min('gravpot_G2')[0]
   gpmax2 = pf.h.find_max('gravpot_G2')[0]

################################################################################################################################

# Get rays through this point
rayx = pf.h.ortho_ray(0,(c[1],c[2]))
rayy = pf.h.ortho_ray(1,(c[0],c[2]))
rayz = pf.h.ortho_ray(2,(c[0],c[1]))

if(plt2):
   rayx2 = pf2.h.ortho_ray(0,(c2[1],c2[2]))
   rayy2 = pf2.h.ortho_ray(1,(c2[0],c2[2]))
   rayz2 = pf2.h.ortho_ray(2,(c2[0],c2[1]))

# Save the distance axes since they are used often
unit = 'pc'
axX = rayx['x']*pf[unit]
axY = rayy['y']*pf[unit]
axZ = rayz['z']*pf[unit]

if(plt2):
   axX2 = rayx2['x']*pf[unit]
   axY2 = rayy2['y']*pf[unit]
   axZ2 = rayz2['z']*pf[unit]


## Custom fields  ##

#dGmag error
if(plt2):
    dGm_errX = 2.0*np.fabs(rayx['dGmag']-rayx2['dGmag_G2'])/(np.fabs(rayx['dGmag']) + np.fabs(rayx2['dGmag_G2']))
    dGm_errY = 2.0*np.fabs(rayy['dGmag']-rayy2['dGmag_G2'])/(np.fabs(rayy['dGmag']) + np.fabs(rayy2['dGmag_G2']))
    dGm_errZ = 2.0*np.fabs(rayz['dGmag']-rayz2['dGmag_G2'])/(np.fabs(rayz['dGmag']) + np.fabs(rayz2['dGmag_G2']))
    #dGm_errX = np.fabs(rayx['dGmag']-rayx2['dGmag_G2'])/(np.fabs(rayx2['dGmag_G2']))
    #dGm_errY = np.fabs(rayy['dGmag']-rayy2['dGmag_G2'])/(np.fabs(rayy2['dGmag_G2']))
    #dGm_errZ = np.fabs(rayz['dGmag']-rayz2['dGmag_G2'])/(np.fabs(rayz2['dGmag_G2']))

# Make some plots subplot(numrows,numcol,fignum=0,1,...,numrow*numcol-1)
plt.subplot(221)
plt.title('Black=x, Green=y, Red=z')
plotme = 'Density'
plt.semilogy(axX,rayx[plotme],color='black')
plt.semilogy(axY,rayy[plotme],color='green')
plt.semilogy(axZ,rayz[plotme],color='red')
if(plt2):
   plt.semilogy(axX2,rayx2[plotme],'--',color='black')
   plt.semilogy(axY2,rayy2[plotme],'--',color='green')
   plt.semilogy(axZ2,rayz2[plotme],'--',color='red')

plt.ylabel(plotme)

plt.subplot(222)
plt.title('Solid=out, Dash=Gadget')
plotme = 'dGmag'
plt.semilogy(axX,rayx[plotme],color='black')
plt.semilogy(axY,rayy[plotme],color='green')
plt.semilogy(axZ,rayz[plotme],color='red')
if(plt2):
   plotme2 = 'dGmag_G2'
   plt.semilogy(axX2,rayx2[plotme2],'--',color='black')
   plt.semilogy(axY2,rayy2[plotme2],'--',color='green')
   plt.semilogy(axZ2,rayz2[plotme2],'--',color='red')
plt.ylabel(plotme)


#Custom field
plt.subplot(223)
#plotme = 'Vmag'
#plt.semilogy(axX,rayx[plotme],color='black')
#plt.semilogy(axY,rayy[plotme],color='green')
#plt.semilogy(axZ,rayz[plotme],color='red')
if(plt2):
   plt.semilogy(axX,dGm_errX,color='black')
   plt.semilogy(axY,dGm_errY,color='green')
   plt.semilogy(axZ,dGm_errZ,color='red')
plt.ylim(-2, 2)  
plt.ylabel('dGmag_error')
plt.xlabel('Distance (' + unit + ')')

plt.subplot(224)
plotme = 'gravitational-potential'
print 'gpmin = ', gpmin
plt.semilogy(axX,rayx[plotme]-gpmin,color='black')
plt.semilogy(axY,rayy[plotme]-gpmin,color='green')
plt.semilogy(axZ,rayz[plotme]-gpmin,color='red')
if(plt2):
   print 'gpmin2 = ', gpmin2
   print 'gpmin2 = ', gpmin2
   plotme2 = 'gravpot_G2'
   plt.semilogy(axX2,rayx2[plotme2]-gpmin2 ,'--',color='black')
   plt.semilogy(axY2,rayy2[plotme2]-gpmin2 ,'--',color='green')
   plt.semilogy(axZ2,rayz2[plotme2]-gpmin2 ,'--',color='red')
plt.ylabel(plotme)
plt.xlabel('Distance (' + unit + ')')


plt.tight_layout()
plt.savefig(sname + ".ps")

