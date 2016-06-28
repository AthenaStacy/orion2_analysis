import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as plt

from yt.mods import *

import numpy as np
#import sys
import fields_wpot
#import scipy

triad1 = [239./256,0./256,42./256]
triad2 = [3./256,137./256,156./256]
triad3 = [168./256,186./256,47./256]
schemename = 'Red=x, Blue=y, Grellow=z'


CompareNum = 2

fname  = '/work/00863/minerva/orion/gravpot_0.4pc_064/data.0000.3d.hdf5'
fname2 = '/work/00863/minerva/orion/gravpot_0.8pc_128/data.0000.3d.hdf5'
fname3 = '/work/00863/minerva/orion/gravpot_1.6pc_256/data.0000.3d.hdf5'
fname4 = 'data.0001.512.1.5pc'

sname = 'TestComp'
DeltaX = 9.64274244e15
DeltaX = 1.0

#Load data
pf = load(fname)

if(CompareNum==2):
    pf2= load(fname2)
elif(CompareNum==3):
    pf2= load(fname3)
else:
    pf2= load(fname4)

# Chose center of data (selecting highest density point)
c = pf.h.find_max('Density')[1]
c2 = pf2.h.find_max('Density')[1]




# Get rays through this point
rayx = pf.h.ortho_ray(0,(c[1],c[2]))
rayy = pf.h.ortho_ray(1,(c[0],c[2]))
rayz = pf.h.ortho_ray(2,(c[0],c[1]))

rayx2 = pf2.h.ortho_ray(0,(c2[1],c2[2]))
rayy2 = pf2.h.ortho_ray(1,(c2[0],c2[2]))
rayz2 = pf2.h.ortho_ray(2,(c2[0],c2[1]))

# Save the distance axes since they are used often
unit = 'pc'
axX = rayx['x']*pf[unit]
axY = rayy['y']*pf[unit]
axZ = rayz['z']*pf[unit]

axX2 = rayx2['x']*pf[unit]
axY2 = rayy2['y']*pf[unit]
axZ2 = rayz2['z']*pf[unit]



## Since ax(XYZ) is 128 and ax(XYZ)2 is only 64, cut out the edges for error analysis
## also, the 64 seems to duplicate the middle cell... and the right side is shifted right...
## so get rid of entry 32 for the 64 array and shrink the 128 
# new axes
#Xforerr = axX[0:63]
#Yforerr = axY[0:63]
#Zforerr = axZ[0:63]
addme = 0
extra = 0
Xforerr = axX[0+addme:63+extra+addme]
Yforerr = axY[0+addme:63+extra+addme]
Zforerr = axZ[0+addme:63+extra+addme]

# new 64 data (omit entry 32)
omitme = 0
#dGmag1x = np.delete(rayx['dGmag'],omitme)
#dGmag1y = np.delete(rayy['dGmag'],omitme)
#dGmag1z = np.delete(rayz['dGmag'],omitme)
dGmag1x = rayx['dGmag'][0:64]
dGmag1y = rayx['dGmag'][0:64]
dGmag1z = rayx['dGmag'][0:64]

if(CompareNum==2):
    # shrink the 128 data to have 64 entries
    dGmag2x = rayx2['dGmag'][32:96+extra]
    dGmag2y = rayy2['dGmag'][32:96+extra]
    dGmag2z = rayz2['dGmag'][32:96+extra]
elif(CompareNum==3):
    # shrink the 256 data to have 64 entries
    dGmag2x = rayx2['dGmag'][97:160+extra]
    dGmag2y = rayy2['dGmag'][97:160+extra]
    dGmag2z = rayz2['dGmag'][97:160+extra]
else:
    # shrink the 512 data to have 64 entries
    dGmag2x = rayx2['dGmag'][225:288+extra]
    dGmag2y = rayy2['dGmag'][225:288+extra]
    dGmag2z = rayz2['dGmag'][225:288+extra]


## Custom fields  ##

#dGmag error
#dGm_errX = 2.0*np.fabs(rayx['dGmag']-rayx2['dGmag'])/(rayx['dGmag']+rayx2['dGmag'])
#dGm_errY = 2.0*np.fabs(rayy['dGmag']-rayy2['dGmag'])/(rayy['dGmag']+rayy2['dGmag'])
#dGm_errZ = 2.0*np.fabs(rayz['dGmag']-rayz2['dGmag'])/(rayz['dGmag']+rayz2['dGmag'])
dGm_errX = 2.0*np.fabs((dGmag2x-dGmag1x)/(dGmag2x+dGmag1x))
dGm_errY = 2.0*np.fabs((dGmag2y-dGmag1y)/(dGmag2y+dGmag1y))
dGm_errZ = 2.0*np.fabs((dGmag2z-dGmag1z)/(dGmag2z+dGmag1z))

# Make some plots subplot(numrows,numcol,fignum=0,1,...,numrow*numcol-1)
ax1 = plt.subplot(221)
plt.title(schemename)
plotme = 'gravitational-potential'
gpmin = pf.h.find_min('gravitational-potential')[0]
gpmin2 = pf2.h.find_min('gravitational-potential')[0]

plt.semilogy(axX,rayx[plotme]-gpmin+4e10,color=triad1)
plt.semilogy(axY,rayy[plotme]-gpmin+4e10,color=triad2)
plt.semilogy(axZ,rayz[plotme]-gpmin+4e10,color=triad3)

plt.semilogy(axX2,rayx2[plotme]-gpmin2+4e10,'--',color=triad1)
plt.semilogy(axY2,rayy2[plotme]-gpmin2+4e10,'--',color=triad2)
plt.semilogy(axZ2,rayz2[plotme]-gpmin2+4e10,'--',color=triad3)

ax1.set_xlim([-0.1,0.1])
plt.ylabel(plotme)



ax2=plt.subplot(222)
plt.title('Solid=64, Dash=128')
plotme = 'dGmag'
plt.semilogy(Xforerr,dGmag1x,color=triad1)
plt.semilogy(Yforerr,10*dGmag1y,color=triad2)
plt.semilogy(Zforerr,0.1*dGmag1z,color=triad3)
#plt.semilogy(axX2,rayx2[plotme],'--',color='black')
#plt.semilogy(axY2,rayy2[plotme],'--',color='green')
#plt.semilogy(axZ2,rayz2[plotme],'--',color='red')
plt.semilogy(Xforerr,dGmag2x,'--',color=triad1)
plt.semilogy(Yforerr,10*dGmag2y,'--',color=triad2)
plt.semilogy(Zforerr,0.1*dGmag2z,'--',color=triad3)
plt.ylabel(plotme)





#Custom field
ax3=plt.subplot(223)
#plotme = 'Vmag'
#plt.semilogy(axX,rayx[plotme],color='black')
#plt.semilogy(axY,rayy[plotme],color='green')
#plt.semilogy(axZ,rayz[plotme],color='red')
plt.semilogy(Xforerr,dGm_errX,color=triad1)
plt.semilogy(Yforerr,dGm_errY,color=triad2)
plt.semilogy(Zforerr,dGm_errZ,color=triad3)

plt.ylabel('dGmag_error')
plt.xlabel('Distance (' + unit + ')')







ax4=plt.subplot(224)
plotme = 'gravitational-potential'
gpmin = pf.h.find_min('gravitational-potential')[0]
gpmin2 = pf2.h.find_min('gravitational-potential')[0]

addme = 0
extra = 1
shift = 0
Xgp1 = axX[0+addme:63+extra+addme]
Gp1x = (rayx[plotme]-1*gpmin+1)/DeltaX
Gp1y = (rayy[plotme]-1*gpmin+1)/DeltaX
Gp1z = (rayz[plotme]-1*gpmin+1)/DeltaX
Gp2x = (rayx2[plotme]-1*gpmin2+1)/DeltaX
Gp2y = (rayy2[plotme]-1*gpmin2+1)/DeltaX
Gp2z = (rayz2[plotme]-1*gpmin2+1)/DeltaX

#Gp1x = (Gp1x[2:64]-Gp1x[0:62])/(2*DeltaX)
#Gp1y = (Gp1y[2:64]-Gp1y[0:62])/(2*DeltaX)
#Gp1z = (Gp1z[2:64]-Gp1z[0:62])/(2*DeltaX)

#if(CompareNum==2):
    # shrink the 128 data to have 64-1 entries
    #Gp2x = (Gp2x[2:128]-Gp2x[0:126])/(2*DeltaX)
    #Gp2y = (Gp2y[2:128]-Gp2y[0:126])/(2*DeltaX)
    #Gp2z = (Gp2z[2:128]-Gp2z[0:126])/(2*DeltaX)

#Gpdiffx = np.fabs(Gp1x - Gp2x[1+32:63+32])/(Gp1x + Gp2x[1+32:63+32])
#Gpdiffy = np.fabs(Gp1y - Gp2y[1+32:63+32])/(Gp1y + Gp2y[1+32:63+32])
#Gpdiffz = np.fabs(Gp1z - Gp2z[1+32:63+32])/(Gp1z + Gp2z[1+32:63+32])
Gpdiffx = 2.0*np.fabs(Gp1x - Gp2x[0+32:64+32])/(Gp1x + Gp2x[0+32:64+32])
Gpdiffy = 2.0*np.fabs(Gp1y - Gp2y[0+32:64+32])/(Gp1y + Gp2y[0+32:64+32])
Gpdiffz = 2.0*np.fabs(Gp1z - Gp2z[0+32:64+32])/(Gp1z + Gp2z[0+32:64+32])

shift=0
Gpdiffx = 1.0*np.fabs(Gp1x - Gp2x[shift+0+32:shift+64+32])
Gpdiffy = 1.0*np.fabs(Gp1y - Gp2y[shift+0+32:shift+64+32])
Gpdiffz = 1.0*np.fabs(Gp1z - Gp2z[shift+0+32:shift+64+32])

#plt.plot(axX[1:63],Gp1x,color='black')
#plt.plot(axY[1:63],Gp1y,color='green')
#plt.plot(axZ[1:63],Gp1z,color='red')
#plt.plot(axX2[1:127],Gp2x,'--',color='black')
#plt.plot(axY2[1:127],Gp2y,'--',color='green')
#plt.plot(axZ2[1:127],Gp2z,'--',color='red')


expand=1
plt.semilogy(axX[1-expand:63+expand],Gpdiffx,'.',color=triad1,alpha=0.6)
plt.semilogy(axY[1-expand:63+expand],Gpdiffy,'.',color=triad2,alpha=0.6)
plt.semilogy(axZ[1-expand:63+expand],Gpdiffz,'.',color=triad3,alpha=0.6)


plt.ylabel(r'$\Phi$/$\Delta$x err')
plt.xlabel('Distance (' + unit + ')')


plt.tight_layout()
plt.savefig(sname + ".pdf")

