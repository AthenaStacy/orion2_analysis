import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as plt

from yt.mods import *

import numpy as np
#import sys
import fields_wpot
#import scipy

au_cm = 1.49597871e13
au_pc = 4.84813681e-6  
pc_cm = 3.08567758e18
Msun = 1.98892e33

triad1 = [239./256,0./256,42./256]
triad2 = [3./256,137./256,156./256]
triad3 = [168./256,186./256,47./256]
schemename = 'Red=x, Blue=y, Grellow=z'

rdir = '/home1/00863/minerva/research_programs/'
file_prefix = 'bin_zoom10_new_cut_ref3_wpot'
#file_prefix = 'bin_HR10_ref3_wpot'

snapnum = '7131'

with open(rdir + file_prefix + '_gas_' + snapnum +'.dat', "r") as f:
        gasdatB = [map(float, line.split()) for line in f]

with open(rdir + file_prefix + '_gpot_' + snapnum +'.dat', "r") as f:
        potdatB = [map(float, line.split()) for line in f]

densx = []
densy = []
densz = []
radpot_pos = []
radpot_neg = []
gpotx = []
gpoty = []
gpotz = []
gaccxLF = []
gaccxRF = []
gaccyLF = []
gaccyRF = []
gacczLF = []
gacczRF = []
gaccx = []
gaccy = []
gaccz = []

narr = 256

fac = 1.2195*1.67e-24
for i in range(narr):
        line_arr = gasdatB[i]
        densx.append(line_arr[2]*fac)
        densy.append(line_arr[3]*fac)
        densz.append(line_arr[4]*fac)
for i in range(narr):
        line_arr = potdatB[i]
        radpot_pos.append(line_arr[0])
	radpot_pos[i] = (radpot_pos[i] ) / pc_cm 
	radpot_neg.append(-radpot_pos[i])
        gpotx.append(line_arr[1])
        gpoty.append(line_arr[2])
        gpotz.append(line_arr[3])
        gaccx.append(line_arr[4])
        gaccy.append(line_arr[5])
        gaccz.append(line_arr[6])

for i in range(narr):
        gaccxLF.append( np.abs(gpotx[0] - gpotx[narr-1]) / (radpot_pos[narr-1] ) )
        gaccxRF.append( np.abs(gpotx[0] - gpotx[narr-1]) / (radpot_pos[narr-1] ) )
        gaccyLF.append( np.abs(gpoty[0] - gpoty[narr-1]) / (radpot_pos[narr-1] ) )
        gaccyRF.append( np.abs(gpoty[0] - gpoty[narr-1]) / (radpot_pos[narr-1] ) )
        gacczLF.append( np.abs(gpotz[0] - gpotz[narr-1]) / (radpot_pos[narr-1] ) )
        gacczRF.append( np.abs(gpotz[0] - gpotz[narr-1]) / (radpot_pos[narr-1] ) )

radpot = []
gpotXx = []
gpotYy = []
gpotZz = []
gaccXx = []
gaccYy = []
gaccZz = []
DeltaX = radpot_pos[narr-1] / float(narr) 
print 'DeltaX =', DeltaX
for i in range(narr):
	radpot.append(radpot_pos[i] - DeltaX)
        gpotXx.append(gpotx[i])
        gpotYy.append(gpoty[i])
        gpotZz.append(gpotz[i])
	gaccXx.append(gaccx[i])
	gaccYy.append(gaccy[i])
	gaccZz.append(gaccz[i])

print 'radpot =', radpot
print 'gpotZz =', gpotZz
print 'gaccZz =', gaccZz
print 'gaccxLF =', gaccxLF

fname = '/work/00863/minerva/orion/gadget2pot_ref3_20pc_256'
g2 = np.fromfile(fname, dtype=float, count=-1)
print 'g2[10] =', g2[10]
print 'g2[10000] =', g2[10000]
g2 = numpy.reshape(g2, [256,256,256])
print 'g2[10] =', g2[255,255,255]
print 'g2[10000] =', g2[0,0,0]

fname  = '/work/00863/minerva/orion/gravpot_20pc_256_ref3/data.0000.3d.hdf5'
fname2 = '/work/00863/minerva/orion/gravpot_20pc_256_ref3/data.0000.3d.hdf5'
cen =  (narr / 2) - 1
cen2 = (narr / 2) - 1

sname = 'TestNew'

#Load data
pf = load(fname)
pf2= load(fname2)

# Chose center of data (selecting highest density point) - but actually c and c2 are not used at all later in the code
c = pf.h.find_max('Density')[1]
print 'max density c =', c
c = pf.h.find_min('gravitational-potential')[1]
print 'min potential c =', c

dims = list(pf.domain_dimensions)
dims2 = list(pf2.domain_dimensions)

cube = pf.h.covering_grid(level=0,left_edge=pf.h.domain_left_edge,dims=dims)
cube2 = pf2.h.covering_grid(level=0,left_edge=pf2.h.domain_left_edge,dims=dims2)

# Create the axes
axis = cube['dx'][1,1,:]
length = len(axis)
print 'length =', length
for i in range(len(axis)): axis[i] = i*axis[i]
for i in range(len(axis)): axis[i] = axis[i] - axis[-1]/2.0
axis2 = cube2['dx'][1,1,:]
print 'length2 =', len(axis2)
for i in range(len(axis2)): axis2[i] = i*axis2[i]
for i in range(len(axis2)): axis2[i] = axis2[i] - axis2[-1]/2.0
unit = 'pc'
axis = axis*pf[unit]
axis2= axis2*pf[unit]

DeltaX = cube['dx'][0,0,0]
#DeltaX = np.abs(DeltaX) * 2. / float(length)
print 'DeltaX = ', DeltaX

i = 127
print 'x[0] =', cube['dx'][1,1,i]
print 'y[0] =', cube['dy'][2,2,i]
print 'z[0] =', cube['dz'][1,1,i]

d = cube['Density']
#d2 = cube2['Density']
#print 'central density = ', d[cen,cen,cen]

g  = cube['gravitational-potential']
#g2 = cube2['gravitational-potential']
#g2 = cube2['gravpot_G2']

print 'gravpot dims =', cube['gravitational-potential'].shape
print 'gpot[0,0,0,] =', g[0,0,0]
print 'gpot[255,255,255] =', g[255,255,255]
print 'gpot[cen,cen,cen] =', g[cen,cen,cen]

def CenDiff(g,i,j,k,idx):
    if idx==1:
        r = g[i+1,j,k]
        l = g[i,j,k]
    elif idx==2:
        r = g[i,j+1,k]
        l = g[i,j,k]
    else:
        r = g[i,j,k+1]
        l = g[i,j,k]
    return np.fabs(r-l)/2.0



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


print 'length1 = ', len(accXx)
print 'length2 = ', len(accX2x)

eps = 1.e-15
accMagX = []
for i in range(len(accXx)): accMagX.append( np.sqrt(accXx[i]**2 + accYx[i]**2 + accZx[i]**2 ) + eps)
accMagX2 = []
for i in range(len(accX2x)): accMagX2.append( np.sqrt(accX2x[i]**2 + accY2x[i]**2 + accZ2x[i]**2) )
accMagY = []
for i in range(len(accXy)): accMagY.append( np.sqrt(accXy[i]**2 + accYy[i]**2 + accZy[i]**2 ) )
accMagY2 = []
for i in range(len(accX2y)): accMagY2.append( np.sqrt(accX2y[i]**2 + accY2y[i]**2 + accZ2y[i]**2) + eps)
accMagZ = []
for i in range(len(accXz)): accMagZ.append( np.sqrt(accXz[i]**2 + accYz[i]**2 + accZz[i]**2 ))
accMagZ2 = []
for i in range(len(accX2z)): accMagZ2.append( np.sqrt(accX2z[i]**2 + accY2z[i]**2 + accZ2z[i]**2) + eps)


Saxis = axis[1:dims[1]-1]
Saxis2 = axis2[1:dims2[1]-1]

print 'axis =', axis[0:dims[1]]

# Time for plots!

# Make some plots subplot(numrows,numcol,fignum=0,1,...,numrow*numcol-1)
ax1=plt.subplot(221)
plt.title(schemename)
plotme = 'gravitational-potential'
gpmin = g.min()
gpmin2 = g2.min()

gpotXx_min = 1.e30
gpotYy_min = 1.e30
gpotZz_min = 1.e30

for i in range(narr):
        if gpotXx[i] < gpotXx_min:
                gpotXx_min = gpotXx[i]
        if gpotYy[i] < gpotYy_min:
                gpotYy_min = gpotYy[i]
        if gpotZz[i] < gpotZz_min:
                gpotZz_min = gpotZz[i]

for i in range(narr):
        gpotXx[i] = gpotXx[i] - gpotXx_min
        gpotYy[i] = gpotYy[i] - gpotYy_min
        gpotZz[i] = gpotZz[i] - gpotZz_min

gp_shift = 0

plt.semilogy(axis,g[:,cen,cen]-gpmin+gp_shift,color=triad1)
plt.semilogy(axis,g[cen,:,cen]-gpmin+gp_shift,color=triad2)
plt.semilogy(axis,g[cen,cen,:]-gpmin+gp_shift,color=triad3)

plt.semilogy(axis2,g2[:,cen2,cen2]-gpmin2+gp_shift,'--',color=triad1)
plt.semilogy(axis2,g2[cen2,:,cen2]-gpmin2+gp_shift,'--',color=triad2)
plt.semilogy(axis2,g2[cen2,cen2,:]-gpmin2+gp_shift,'--',color=triad3)

#plt.plot(radpot, gpotXx, '--', color=triad1)
#plt.plot(radpot, gpotYy, '--', color=triad2)
#plt.plot(radpot, gpotZz, '--', color=triad3)

#ax1.set_xlim([-0.15,0.15])
plt.ylabel(plotme)



ax2=plt.subplot(222)
plt.title('Solid=Orion, Dash=Gadget')
plotme = '1D Acc in colored dir'
plt.plot(Saxis,accXx,color=triad1)
plt.plot(Saxis,accYy,color=triad2)
plt.plot(Saxis,accZz,color=triad3)

plt.plot(Saxis2,accX2x,'--',color=triad1)
plt.plot(Saxis2,accY2y,'--',color=triad2)
plt.plot(Saxis2,accZ2z,'--',color=triad3)

#plt.plot(radpot, gaccXx, '--', color=triad1)
#plt.plot(radpot, gaccYy, '--', color=triad2)
#plt.plot(radpot, gaccZz, '--', color=triad3)

ax2.set_yscale('log')
ax2.set_ylim([1.e-9,1.e-4])



#Custom field
ax3=plt.subplot(223)

grid_shift = 0
densx_arr = []
densy_arr = []
densz_arr = []
for i in range(len(axis)): 
	densx_arr.append(d[i, cen, cen]) 
	densy_arr.append(d[cen, i, cen])
	densz_arr.append(d[cen, cen, i])

plt.semilogy(axis,densx_arr,color=triad1)
plt.semilogy(axis,densy_arr,color=triad2)
plt.semilogy(axis,densz_arr,color=triad3)

plt.semilogy(radpot,densx,'--',color=triad1)
plt.semilogy(radpot,densy,'--',color=triad2)
plt.semilogy(radpot,densz,'--',color=triad3)

ax3.set_xlim([-1,1])
plt.ylabel('Density [g cm^-3]')
plt.xlabel('Distance (' + unit + ')')


print 'length1 =', len(accMagX)
print 'length2 =', len(accXx)

ax4=plt.subplot(224)

accErrorX = []
for i in range(len(accMagX)): accErrorX.append( 2.0*np.fabs(accMagX2[i+grid_shift]-accMagX[i])/( np.fabs(accMagX2[i+grid_shift]) + np.fabs(accMagX[i]) ) )
#for i in range(len(accXx)): accErrorX.append( 2.0*np.fabs(gaccXx[i+grid_shift]-accXx[i])/( np.fabs(gaccXx[i+grid_shift]) + np.fabs(accXx[i]) ) )
accErrorY = []
for i in range(len(accMagY)): accErrorY.append( 2.0*np.fabs(accMagY2[i+grid_shift]-accMagY[i])/( np.fabs(accMagY2[i+grid_shift]) + np.fabs(accMagY[i])) )
#for i in range(len(accYy)): accErrorY.append( 2.0*np.fabs(gaccYy[i+grid_shift]-accYy[i])/( np.fabs(gaccYy[i+grid_shift]) + np.fabs(accYy[i]) ) )
accErrorZ = []
for i in range(len(accMagZ)): accErrorZ.append( 2.0*np.fabs(accMagZ2[i+grid_shift]-accMagZ[i])/( np.fabs(accMagZ2[i+grid_shift]) + np.fabs(accMagZ[i])) )
#for i in range(len(accZz)): accErrorZ.append( 2.0*np.fabs(gaccZz[i+grid_shift]-accZz[i])/( np.fabs(gaccZz[i+grid_shift]) + np.fabs(accZz[i]) ) )

plt.semilogy(Saxis,accErrorX,color=triad1)
plt.semilogy(Saxis,accErrorY,color=triad2)
plt.semilogy(Saxis,accErrorZ,color=triad3)
plt.axis((-1.0, 1.0, 1.e-3, 1.e1))

plt.ylabel(r'Acc err')
plt.xlabel('Distance (' + unit + ')')

for i in range(len(accXx)): print 'i=', i, 'radpot =', radpot[i+grid_shift], 'saxis =', Saxis[i], 'gaccXx = ', gaccXx[i+grid_shift], 'accXx =', accXx[i]

plt.tight_layout()
plt.savefig(sname + ".pdf")

