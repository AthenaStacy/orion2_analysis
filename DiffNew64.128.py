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

#fname  = '/work/00863/minerva/orion/gravpot_0.8pc_128/data.0000.3d.hdf5'
#fname2 = '/work/00863/minerva/orion/gravpot_0.8pc_128/data.0000.3d.hdf5'

fname =  '/work/00863/minerva/orion/gravpot_0.5pc_256_ref3/data.0000.3d.hdf5'
fname2 = '/work/00863/minerva/orion/gravpot_1pc_512_ref3//data.0000.3d.hdf5'

cen = 127
cen2 = 255
half_grid_num = cen + 1

sname = 'TestNew'


#Load data
pf = load(fname)
pf2= load(fname2)

# Chose center of data (selecting highest density point)
c = pf.h.find_max('Density')[1]
c2 = pf2.h.find_max('Density')[1]

dims = list(pf.domain_dimensions)
dims2 = list(pf2.domain_dimensions)

cube = pf.h.covering_grid(level=0,left_edge=pf.h.domain_left_edge,dims=dims)
cube2 = pf2.h.covering_grid(level=0,left_edge=pf2.h.domain_left_edge,dims=dims2)


# Create the axes
axis = cube['dx'][1,1,:]
for i in range(len(axis)): axis[i] = i*axis[i]
for i in range(len(axis)): axis[i] = axis[i] - axis[-1]/2.0
axis2 = cube2['dx'][1,1,:]
for i in range(len(axis2)): axis2[i] = i*axis2[i]
for i in range(len(axis2)): axis2[i] = axis2[i] - axis2[-1]/2.0
unit = 'pc'
axis = axis*pf[unit]
axis2= axis2*pf[unit]

DeltaX = cube['dx'][1,1,1]



d = cube['Density']
d2 = cube2['Density']

g  = cube['gravitational-potential']
g2 = cube2['gravitational-potential']
#g2 = cube2['gravpot_G2']


def CenDiff(g,i,j,k,idx):
    if idx==1:
        r = g[i+1,j,k]
        l = g[i-1,j,k]
    elif idx==2:
        r = g[i,j+1,k]
        l = g[i,j-1,k]
    else:
        r = g[i,j,k+1]
        l = g[i,j,k-1]
    return (r-l)/2.0



## Custom fields  ##
accXx=[]
accX2x=[]
accYx=[]
accY2x=[]
accZx=[]
accZ2x=[]

accXy=[]
accX2y=[]
accYy=[]
accY2y=[]
accZy=[]
accZ2y=[]

accXz=[]
accX2z=[]
accYz=[]
accY2z=[]
accZz=[]
accZ2z=[]


# Acceleration
for i in range(1,len(axis)-1,1): accXx.append( CenDiff(g,i,cen,cen,1)/DeltaX )
for i in range(1,len(axis2)-1,1): accX2x.append( CenDiff(g2,i,cen2,cen2,1)/DeltaX )
for i in range(1,len(axis)-1,1): accYx.append( CenDiff(g,i,cen,cen,2)/DeltaX )
for i in range(1,len(axis2)-1,1): accY2x.append( CenDiff(g2,i,cen2,cen2,2)/DeltaX )
for i in range(1,len(axis)-1,1): accZx.append( CenDiff(g,i,cen,cen,3)/DeltaX )
for i in range(1,len(axis2)-1,1): accZ2x.append( CenDiff(g2,i,cen2,cen2,3)/DeltaX )

for i in range(1,len(axis)-1,1): accXy.append( CenDiff(g,cen,i,cen,1)/DeltaX )
for i in range(1,len(axis2)-1,1): accX2y.append( CenDiff(g2,cen2,i,cen2,1)/DeltaX )
for i in range(1,len(axis)-1,1): accYy.append( CenDiff(g,cen,i,cen,2)/DeltaX )
for i in range(1,len(axis2)-1,1): accY2y.append( CenDiff(g2,cen2,i,cen2,2)/DeltaX )
for i in range(1,len(axis)-1,1): accZy.append( CenDiff(g,cen,i,cen,3)/DeltaX )
for i in range(1,len(axis2)-1,1): accZ2y.append( CenDiff(g2,cen2,i,cen2,3)/DeltaX )

for i in range(1,len(axis)-1,1): accXz.append( CenDiff(g,cen,cen,i,1)/DeltaX )
for i in range(1,len(axis2)-1,1): accX2z.append( CenDiff(g2,cen2,cen2,i,1)/DeltaX )
for i in range(1,len(axis)-1,1): accYz.append( CenDiff(g,cen,cen,i,2)/DeltaX )
for i in range(1,len(axis2)-1,1): accY2z.append( CenDiff(g2,cen2,cen2,i,2)/DeltaX )
for i in range(1,len(axis)-1,1): accZz.append( CenDiff(g,cen,cen,i,3)/DeltaX )
for i in range(1,len(axis2)-1,1): accZ2z.append( CenDiff(g2,cen2,cen2,i,3)/DeltaX )




accMagX = []
for i in range(len(accXx)): accMagX.append( np.sqrt(accXx[i]**2 + accYx[i]**2 + accZx[i]**2 ))
accMagX2 = []
for i in range(len(accX2x)): accMagX2.append( np.sqrt(accX2x[i]**2 + accY2x[i]**2 + accZ2x[i]**2) )
accMagY = []
for i in range(len(accXy)): accMagY.append( np.sqrt(accXy[i]**2 + accYy[i]**2 + accZy[i]**2 ))
accMagY2 = []
for i in range(len(accX2y)): accMagY2.append( np.sqrt(accX2y[i]**2 + accY2y[i]**2 + accZ2y[i]**2) )
accMagZ = []
for i in range(len(accXz)): accMagZ.append( np.sqrt(accXz[i]**2 + accYz[i]**2 + accZz[i]**2 ))
accMagZ2 = []
for i in range(len(accX2z)): accMagZ2.append( np.sqrt(accX2z[i]**2 + accY2z[i]**2 + accZ2z[i]**2) )


Saxis = axis[1:dims[1]-1]
Saxis2 = axis2[1:dims2[1]-1]



# Time for plots!

# Make some plots subplot(numrows,numcol,fignum=0,1,...,numrow*numcol-1)
ax1=plt.subplot(221)
plt.title(schemename)
plotme = 'gravitational-potential'
gpmin = g.min()
gpmin2 = g2.min()

plt.semilogy(axis,g[:,cen,cen]-gpmin+4e10,color=triad1)
plt.semilogy(axis,g[cen,:,cen]-gpmin+4e10,color=triad2)
plt.semilogy(axis,g[cen,cen,:]-gpmin+4e10,color=triad3)

plt.semilogy(axis2,g2[:,cen2,cen2]-gpmin2+4e10,'--',color=triad1)
plt.semilogy(axis2,g2[cen2,:,cen2]-gpmin2+4e10,'--',color=triad2)
plt.semilogy(axis2,g2[cen2,cen2,:]-gpmin2+4e10,'--',color=triad3)

ax1.set_xlim([-0.15,0.15])
plt.ylabel(plotme)



ax2=plt.subplot(222)
plt.title('Solid=64, Dash=128')
plotme = '1D Acc in colored dir'
plt.plot(Saxis,accXx,color=triad1)
plt.plot(Saxis,accYy,color=triad2)
plt.plot(Saxis,accZz,color=triad3)
#plt.semilogy(Yforerr,10*dGmag1y,color=triad2)
#plt.semilogy(Zforerr,0.1*dGmag1z,color=triad3)
plt.plot(Saxis2,accX2x,'--',color=triad1)
plt.plot(Saxis2,accY2y,'--',color=triad2)
plt.plot(Saxis2,accZ2z,'--',color=triad3)
#plt.semilogy(Yforerr,10*dGmag2y,'--',color=triad2)
#plt.semilogy(Zforerr,0.1*dGmag2z,'--',color=triad3)
plt.ylabel(plotme)





#Custom field
ax3=plt.subplot(223)
#plotme = 'Vmag'

plt.semilogy(Saxis,accMagX,color=triad1)
plt.semilogy(Saxis,accMagY,color=triad2)
plt.semilogy(Saxis,accMagZ,color=triad3)


plt.semilogy(Saxis2,accMagX2,'--',color=triad1)
plt.semilogy(Saxis2,accMagY2,'--',color=triad2)
plt.semilogy(Saxis2,accMagZ2,'--',color=triad3)

#plt.semilogy(Yforerr,dGm_errY,color=triad2)
#plt.semilogy(Zforerr,dGm_errZ,color=triad3)

plt.ylabel('Mag of Acc')
plt.xlabel('Distance (' + unit + ')')




ax4=plt.subplot(224)

accErrorX = []
for i in range(len(accMagX)): accErrorX.append( 2.0*np.fabs(accMagX2[i+half_grid_num]-accMagX[i])/( accMagX2[i+half_grid_num]+accMagX[i]) )
accErrorY = []
for i in range(len(accMagY)): accErrorY.append( 2.0*np.fabs(accMagY2[i+half_grid_num]-accMagY[i])/( accMagY2[i+half_grid_num]+accMagY[i]) )
accErrorZ = []
for i in range(len(accMagZ)): accErrorZ.append( 2.0*np.fabs(accMagZ2[i+half_grid_num]-accMagZ[i])/( accMagZ2[i+half_grid_num]+accMagZ[i]) )


plt.semilogy(Saxis,accErrorX,color=triad1)
plt.semilogy(Saxis,accErrorY,color=triad2)
plt.semilogy(Saxis,accErrorZ,color=triad3)

plt.ylabel(r'Acc err')
plt.xlabel('Distance (' + unit + ')')


plt.tight_layout()
plt.savefig(sname + ".pdf")

