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
#file_prefix = 'bin_HR10_ref3_wpot'

file_prefix = 'bin_HR10_wpot'

snapnum = '0002'

#file_prefix = 'bin_zoom10_new_cut_ref3_wpot'
#file_prefix = 'bin_zoom10_new_10pc_ref3_wpot'
#snapnum = '7131'
#snapnum = '6901'

#file_prefix = 'bin_zoom9_new_cut_ref3_wpot'
#snapnum = '9451'

with open(rdir + file_prefix + '_gas_' + snapnum +'.dat', "r") as f:
        gasdatB = [map(float, line.split()) for line in f]

with open(rdir + file_prefix + '_gpot_' + snapnum +'.dat', "r") as f:
        potdatB = [map(float, line.split()) for line in f]

with open(rdir + file_prefix + '_gas2_' + snapnum +'.dat', "r") as f:
        gasdatB2 = [map(float, line.split()) for line in f]

with open(rdir + file_prefix + '_gpot2_' + snapnum +'.dat', "r") as f:
        potdatB2 = [map(float, line.split()) for line in f]

with open(rdir + file_prefix + '_gas3_' + snapnum +'.dat', "r") as f:
        gasdatB3 = [map(float, line.split()) for line in f]

with open(rdir + file_prefix + '_gpot3_' + snapnum +'.dat', "r") as f:
        potdatB3 = [map(float, line.split()) for line in f]


radpot_pos = []
densx = []
densy = []
densz = []
gpotXx = []
gpotYy = []
gpotZz = []
gaccXx = []
gaccYy = []
gaccZz = []

radpot_pos2 = []
densx2 = []
densy2 = []
densz2 = []
gpotXx2 = []
gpotYy2 = []
gpotZz2 = []
gaccXx2 = []
gaccYy2 = []
gaccZz2 = []

radpot_pos3 = []
densx3 = []
densy3 = []
densz3 = []
gpotXx3 = []
gpotYy3 = []
gpotZz3 = []
gaccXx3 = []
gaccYy3 = []
gaccZz3 = []


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
        gpotXx.append(line_arr[1])
        gpotYy.append(line_arr[2])
        gpotZz.append(line_arr[3])
        gaccXx.append(line_arr[4])
        gaccYy.append(line_arr[5])
        gaccZz.append(line_arr[6])

for i in range(narr):
        line_arr = gasdatB2[i]
        densx2.append(line_arr[2]*fac)
        densy2.append(line_arr[3]*fac)
        densz2.append(line_arr[4]*fac)
for i in range(narr):
        line_arr = potdatB2[i]
        radpot_pos2.append(line_arr[0])
        radpot_pos2[i] = (radpot_pos2[i] ) / pc_cm
        gpotXx2.append(line_arr[1])
        gpotYy2.append(line_arr[2])
        gpotZz2.append(line_arr[3])
        gaccXx2.append(line_arr[4])
        gaccYy2.append(line_arr[5])
        gaccZz2.append(line_arr[6])

for i in range(narr):
        line_arr = gasdatB3[i]
        densx3.append(line_arr[2]*fac)
        densy3.append(line_arr[3]*fac)
        densz3.append(line_arr[4]*fac)
for i in range(narr):
        line_arr = potdatB3[i]
        radpot_pos3.append(line_arr[0])
        radpot_pos3[i] = (radpot_pos3[i] ) / pc_cm
        gpotXx3.append(line_arr[1])
        gpotYy3.append(line_arr[2])
        gpotZz3.append(line_arr[3])
        gaccXx3.append(line_arr[4])
        gaccYy3.append(line_arr[5])
        gaccZz3.append(line_arr[6])


radpot = []
DeltaX = radpot_pos[narr-1] / float(narr) 
print 'DeltaX =', DeltaX
for i in range(narr):
	radpot.append(radpot_pos[i] - DeltaX)

radpot2 = []
DeltaX = radpot_pos2[narr-1] / float(narr)
print 'DeltaX =', DeltaX
for i in range(narr):
        radpot2.append(radpot_pos2[i] - DeltaX)

radpot3 = []
DeltaX = radpot_pos3[narr-1] / float(narr)
print 'DeltaX =', DeltaX
for i in range(narr):
        radpot3.append(radpot_pos3[i] - DeltaX)

print 'radpot =', radpot
print 'gpotZz =', gpotZz
print 'gaccZz =', gaccZz

fname  = '/work/00863/minerva/orion/data.0000.3d.hdf5'

#fname  = '/work/00863/minerva/orion/gravpot_2pc_256Nest_HR10/data.0000.3d.hdf5'
#fname2  = '/work/00863/minerva/orion/gravpot_2pc_256Nest_HR10/data.0000.3d.hdf5'

sname = 'TestNew'

#Load data
pf = load(fname)

# Chose center of data (selecting highest density point) - but actually c and c2 are not used at all later in the code
c = pf.h.find_max('Density')[1]
print 'max density c =', c
c = pf.h.find_min('gravitational-potential')[1]
print 'min potential c =', c

length = pf.h.domain_right_edge[1] - pf.h.domain_left_edge[1]
i_pmin = pf.domain_dimensions[1] * c[0] / length + pf.domain_dimensions[1]/2
j_pmin = pf.domain_dimensions[1] * c[1] / length + pf.domain_dimensions[1]/2
k_pmin = pf.domain_dimensions[1] * c[2] / length + pf.domain_dimensions[1]/2
 
print 'i_pmin =', i_pmin, 'j_pmin =', j_pmin, 'k_pmin =', k_pmin

dims = list(pf.domain_dimensions)

nlev = 2

cube = pf.h.covering_grid(level=0, left_edge=1.0*pf.h.domain_left_edge, right_edge=1.0*pf.h.domain_right_edge, dims=pf.domain_dimensions)
cube2 = pf.h.covering_grid(level=nlev-1, left_edge=0.5*pf.h.domain_left_edge, right_edge=0.5*pf.h.domain_right_edge, dims=pf.domain_dimensions)
cube3 = pf.h.covering_grid(level=nlev,   left_edge=0.25*pf.h.domain_left_edge, right_edge=0.25*pf.h.domain_right_edge, dims=pf.domain_dimensions)

# Create the axes
axis = cube['dx'][1,1,:]
length = len(axis)
print 'length_axis =', length
for i in range(len(axis)): axis[i] = i*axis[i]
for i in range(len(axis)): axis[i] = axis[i] - axis[-1]/2.0

axis2 = cube2['dx'][1,1,:]
print 'length_axis2 =', len(axis2)
for i in range(len(axis2)): axis2[i] = i*axis2[i]
for i in range(len(axis2)): axis2[i] = axis2[i] - axis2[-1]/2.0
unit = 'pc'

axis3 = cube3['dx'][1,1,:]
print 'length_axis3 =', len(axis3)
for i in range(len(axis3)): axis3[i] = i*axis3[i]
for i in range(len(axis3)): axis3[i] = axis3[i] - axis3[-1]/2.0

unit = 'pc'
axis = axis*pf[unit]
axis2= axis2*pf[unit]
axis3= axis3*pf[unit]

DeltaX = cube['dx'][0,0,0]
print 'DeltaX = ', DeltaX

DeltaX2 = cube2['dx'][0,0,0]
print 'DeltaX2 = ', DeltaX2

DeltaX3 = cube3['dx'][0,0,0]
print 'DeltaX3 = ', DeltaX3

cen =  len(axis)  / 2  - 1
cen2 = len(axis2) / 2 - 1
cen3 = len(axis3) / 2 - 1

d = cube['Density']
d2 = cube2['Density']
d3 = cube3['Density']
print 'central density = ', d2[cen2,cen2,cen2]


g  = cube['gravitational-potential']
g2 = cube2['gravitational-potential']
g3 = cube3['gravitational-potential']

print 'gpot[0,0,0,] =', g[0,0,0]
print 'gpot[cen,cen,cen] =', g[cen,cen,cen]

#def CenDiff(g,i,j,k,idx):
#    if idx==1:
#        r = g[i+1,j,k]
#        l = g[i-1,j,k]
#    elif idx==2:
#        r = g[i,j+1,k]
#        l = g[i,j-1,k]
#    else:
#        r = g[i,j,k+1]
#        l = g[i,j,k-1]
#    return np.fabs(r-l)/2.0

def CenDiff(g,i,j,k,idx):
    if idx==1:
	if i < i_pmin:
        	r = g[i+1,j,k]   
        	l = g[i-1,j,k]
	if i >= i_pmin:
                r = g[i+1,j,k]
                l = g[i-1,j,k]
    elif idx==2:
	if j < j_pmin:
        	r = g[i,j+1,k] 
        	l = g[i,j-1,k]
	if j >= j_pmin:
                r = g[i,j+1,k]
                l = g[i,j-1,k]
    else:
	if k < k_pmin:
        	r = g[i,j,k+1] 
        	l = g[i,j,k-1]
        if k >= k_pmin:
                r = g[i,j,k+1]
                l = g[i,j,k-1]
    return np.fabs(r-l)/2.0


## Custom fields  ##
accXx=[]
accYx=[]
accZx=[]

accXy=[]
accYy=[]
accZy=[]

accXz=[]
accYz=[]
accZz=[]

accXx2=[]
accYx2=[]
accZx2=[]

accXy2=[]
accYy2=[]
accZy2=[]

accXz2=[]
accYz2=[]
accZz2=[]

accXx3=[]
accYx3=[]
accZx3=[]

accXy3=[]
accYy3=[]
accZy3=[]

accXz3=[]
accYz3=[]
accZz3=[]

# Acceleration
gshift = 1
for i in range(1,len(axis)-gshift,1): accXx.append( CenDiff(g,i,cen,cen,1)/DeltaX )
for i in range(1,len(axis)-gshift,1): accYx.append( CenDiff(g,i,cen,cen,2)/DeltaX )
for i in range(1,len(axis)-gshift,1): accZx.append( CenDiff(g,i,cen,cen,3)/DeltaX )

for i in range(1,len(axis)-gshift,1): accXy.append( CenDiff(g,cen,i,cen,1)/DeltaX )
for i in range(1,len(axis)-gshift,1): accYy.append( CenDiff(g,cen,i,cen,2)/DeltaX )
for i in range(1,len(axis)-gshift,1): accZy.append( CenDiff(g,cen,i,cen,3)/DeltaX )

for i in range(1,len(axis)-gshift,1): accXz.append( CenDiff(g,cen,cen,i,1)/DeltaX )
for i in range(1,len(axis)-gshift,1): accYz.append( CenDiff(g,cen,cen,i,2)/DeltaX )
for i in range(1,len(axis)-gshift,1): accZz.append( CenDiff(g,cen,cen,i,3)/DeltaX )

for i in range(1,len(axis2)-gshift,1): accXx2.append( CenDiff(g2,i,cen2,cen2,1)/DeltaX2 )
for i in range(1,len(axis2)-gshift,1): accYx2.append( CenDiff(g2,i,cen2,cen2,2)/DeltaX2 )
for i in range(1,len(axis2)-gshift,1): accZx2.append( CenDiff(g2,i,cen2,cen2,3)/DeltaX2 )

for i in range(1,len(axis2)-gshift,1): accXy2.append( CenDiff(g2,cen2,i,cen2,1)/DeltaX2 )
for i in range(1,len(axis2)-gshift,1): accYy2.append( CenDiff(g2,cen2,i,cen2,2)/DeltaX2 )
for i in range(1,len(axis2)-gshift,1): accZy2.append( CenDiff(g2,cen2,i,cen2,3)/DeltaX2 )

for i in range(1,len(axis2)-gshift,1): accXz2.append( CenDiff(g2,cen2,cen2,i,1)/DeltaX2 )
for i in range(1,len(axis2)-gshift,1): accYz2.append( CenDiff(g2,cen2,cen2,i,2)/DeltaX2 )
for i in range(1,len(axis2)-gshift,1): accZz2.append( CenDiff(g2,cen2,cen2,i,3)/DeltaX2 )

for i in range(1,len(axis3)-gshift,1): accXx3.append( CenDiff(g3,i,cen3,cen3,1)/DeltaX3 )
for i in range(1,len(axis3)-gshift,1): accYx3.append( CenDiff(g3,i,cen3,cen3,2)/DeltaX3 )
for i in range(1,len(axis3)-gshift,1): accZx3.append( CenDiff(g3,i,cen3,cen3,3)/DeltaX3 )

for i in range(1,len(axis3)-gshift,1): accXy3.append( CenDiff(g3,cen3,i,cen3,1)/DeltaX3 )
for i in range(1,len(axis3)-gshift,1): accYy3.append( CenDiff(g3,cen3,i,cen3,2)/DeltaX3 )
for i in range(1,len(axis3)-gshift,1): accZy3.append( CenDiff(g3,cen3,i,cen3,3)/DeltaX3 )

for i in range(1,len(axis3)-gshift,1): accXz3.append( CenDiff(g3,cen3,cen3,i,1)/DeltaX3 )
for i in range(1,len(axis3)-gshift,1): accYz3.append( CenDiff(g3,cen3,cen3,i,2)/DeltaX3 )
for i in range(1,len(axis3)-gshift,1): accZz3.append( CenDiff(g3,cen3,cen3,i,3)/DeltaX3 )

print 'length_accXx = ', len(accXx)

accMagX = []
for i in range(len(accXx)): accMagX.append( np.sqrt(accXx[i]**2 + accYx[i]**2 + accZx[i]**2 ) )
accMagY = []
for i in range(len(accXy)): accMagY.append( np.sqrt(accXy[i]**2 + accYy[i]**2 + accZy[i]**2 ) )
accMagZ = []
for i in range(len(accXz)): accMagZ.append( np.sqrt(accXz[i]**2 + accYz[i]**2 + accZz[i]**2 ))

accMagX2 = []
for i in range(len(accXx2)): accMagX2.append( np.sqrt(accXx2[i]**2 + accYx2[i]**2 + accZx2[i]**2 ) )
accMagY2 = []
for i in range(len(accXy2)): accMagY2.append( np.sqrt(accXy2[i]**2 + accYy2[i]**2 + accZy2[i]**2 ) )
accMagZ2 = []
for i in range(len(accXz2)): accMagZ2.append( np.sqrt(accXz2[i]**2 + accYz2[i]**2 + accZz2[i]**2 ))

accMagX3 = []
for i in range(len(accXx3)): accMagX3.append( np.sqrt(accXx3[i]**2 + accYx3[i]**2 + accZx3[i]**2 ) )
accMagY3 = []
for i in range(len(accXy3)): accMagY3.append( np.sqrt(accXy3[i]**2 + accYy3[i]**2 + accZy3[i]**2 ) )
accMagZ3 = []
for i in range(len(accXz3)): accMagZ3.append( np.sqrt(accXz3[i]**2 + accYz3[i]**2 + accZz3[i]**2 ))

Saxis = axis[1:(len(axis))-gshift]
Saxis2 = axis2[1:(len(axis2))-gshift]
Saxis3 = axis3[1:(len(axis3))-gshift]

print 'length_Saxis =', len(Saxis)

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

plt.plot(radpot, gpotXx, '--', color=triad1)
plt.plot(radpot, gpotYy, '--', color=triad2)
plt.plot(radpot, gpotZz, '--', color=triad3)

#ax1.set_xlim([-0.15,0.15])
plt.ylabel(plotme)
plt.xlabel('Distance (' + unit + ')')


ax2=plt.subplot(222)
plt.title('Solid=Orion, Dash=Gadget')
plotme = '1D Acc in colored dir'


length = len(Saxis)
for i in range(length/2 - length/4 + 2, length/2 + length/4 - 2):
	accXx[i] = 0 
	accYy[i] = 0
	accZz[i] = 0
length = len(Saxis2)
for i in range(length/2 - length/4 + 2, length/2 + length/4 - 2):
        accXx2[i] = 0
        accYy2[i] = 0
        accZz2[i] = 0

length = len(radpot)
for i in range(length/2 - length/4 + 2, length/2 + length/4 - 2):
        gaccXx[i] = 0
        gaccYy[i] = 0
        gaccZz[i] = 0
length = len(radpot2)
for i in range(length/2 - length/4 + 2, length/2 + length/4 - 2):
        gaccXx2[i] = 0
        gaccYy2[i] = 0
        gaccZz2[i] = 0

plt.plot(Saxis,accXx,color=triad1)
plt.plot(Saxis,accYy,color=triad2)
plt.plot(Saxis,accZz,color=triad3)

plt.plot(Saxis2,accXx2,color=triad1)
plt.plot(Saxis2,accYy2,color=triad2)
plt.plot(Saxis2,accZz2,color=triad3)

plt.plot(Saxis3,accXx3,color=triad1)
plt.plot(Saxis3,accYy3,color=triad2)
plt.plot(Saxis3,accZz3,color=triad3)

plt.plot(radpot, gaccXx, '--', color=triad1)
plt.plot(radpot, gaccYy, '--', color=triad2)
plt.plot(radpot, gaccZz, '--', color=triad3)

plt.plot(radpot2, gaccXx2, '--', color=triad1)
plt.plot(radpot2, gaccYy2, '--', color=triad2)
plt.plot(radpot2, gaccZz2, '--', color=triad3)

plt.plot(radpot3, gaccXx3, '--', color=triad1)
plt.plot(radpot3, gaccYy3, '--', color=triad2)
plt.plot(radpot3, gaccZz3, '--', color=triad3)

ax2.set_yscale('log')
ax2.set_ylim([1.e-9,1.e-4])

plt.ylabel(r'Acc [cm s^-2]')
plt.xlabel('Distance (' + unit + ')')


#Custom field
ax3=plt.subplot(223)

grid_shift = 1

densx_arr = []
densy_arr = []
densz_arr = []

densx2_arr = []
densy2_arr = []
densz2_arr = []

densx3_arr = []
densy3_arr = []
densz3_arr = []

for i in range(len(axis)):
        densx_arr.append(d[i, cen, cen])
        densy_arr.append(d[cen, i, cen])
        densz_arr.append(d[cen, cen, i])

for i in range(len(axis2)): 
	densx2_arr.append(d2[i, cen2, cen2]) 
	densy2_arr.append(d2[cen2, i, cen2])
	densz2_arr.append(d2[cen2, cen2, i])

for i in range(len(axis3)):
        densx3_arr.append(d3[i, cen3, cen3])
        densy3_arr.append(d3[cen3, i, cen3])
        densz3_arr.append(d3[cen3, cen3, i])


length = len(axis2)
for i in range(length/2 - length/4, length/2 + length/4):
        densx2_arr[i] = 0
        densy2_arr[i] = 0
        densz2_arr[i] = 0
length = len(axis3)
for i in range(length/2 - length/4, length/2 + length/4):
        densx3_arr[i] = 0
        densy3_arr[i] = 0
        densz3_arr[i] = 0

length = len(radpot2)
for i in range(length/2 - length/4, length/2 + length/4):
        densx2[i] = 0
        densy2[i] = 0
        densz2[i] = 0
length = len(radpot3)
for i in range(length/2 - length/4, length/2 + length/4):
        densx3[i] = 0
        densy3[i] = 0
        densz3[i] = 0

plt.semilogy(axis,densx_arr,color=triad1)
plt.semilogy(axis,densy_arr,color=triad2)
plt.semilogy(axis,densz_arr,color=triad3)

plt.semilogy(axis2,densx2_arr,color=triad1)
plt.semilogy(axis2,densy2_arr,color=triad2)
plt.semilogy(axis2,densz2_arr,color=triad3)

plt.semilogy(axis3,densx3_arr,color=triad1)
plt.semilogy(axis3,densy3_arr,color=triad2)
plt.semilogy(axis3,densz3_arr,color=triad3)

plt.semilogy(radpot,densx,'--',color=triad1)
plt.semilogy(radpot,densy,'--',color=triad2)
plt.semilogy(radpot,densz,'--',color=triad3)

plt.semilogy(radpot2,densx2,'--',color=triad1)
plt.semilogy(radpot2,densy2,'--',color=triad2)
plt.semilogy(radpot2,densz2,'--',color=triad3)

plt.semilogy(radpot3,densx3,'--',color=triad1)
plt.semilogy(radpot3,densy3,'--',color=triad2)
plt.semilogy(radpot3,densz3,'--',color=triad3)

#ax3.set_xlim([-0.3,0.3])
plt.ylabel('Density [g cm^-3]')
plt.xlabel('Distance (' + unit + ')')


print 'length1 =', len(accMagX)
print 'length2 =', len(accXx)

ax4=plt.subplot(224)

accErrorX = []
for i in range(len(accXx)): accErrorX.append( 2.0*np.fabs(gaccXx[i+grid_shift]-np.fabs(accMagX[i]))/( np.fabs(gaccXx[i+grid_shift]) + np.fabs(accMagX[i]) ) )
accErrorY = []
for i in range(len(accYy)): accErrorY.append( 2.0*np.fabs(gaccYy[i+grid_shift]-np.fabs(accMagY[i]))/( np.fabs(gaccYy[i+grid_shift]) + np.fabs(accMagY[i]) ) )
accErrorZ = []
for i in range(len(accZz)): accErrorZ.append( 2.0*np.fabs(gaccZz[i+grid_shift]-np.fabs(accMagZ[i]))/( np.fabs(gaccZz[i+grid_shift]) + np.fabs(accMagZ[i]) ) )
accErrorX2 = []
for i in range(len(accXx2)): accErrorX2.append( 2.0*np.fabs(gaccXx2[i+grid_shift]-np.fabs(accMagX2[i]))/( np.fabs(gaccXx2[i+grid_shift]) + np.fabs(accMagX2[i]) ) )
accErrorY2 = []
for i in range(len(accYy2)): accErrorY2.append( 2.0*np.fabs(gaccYy2[i+grid_shift]-np.fabs(accMagY2[i]))/( np.fabs(gaccYy2[i+grid_shift]) + np.fabs(accMagY2[i]) ) )
accErrorZ2 = []
for i in range(len(accZz2)): accErrorZ2.append( 2.0*np.fabs(gaccZz2[i+grid_shift]-np.fabs(accMagZ2[i]))/( np.fabs(gaccZz2[i+grid_shift]) + np.fabs(accMagZ2[i]) ) )
accErrorX3 = []
for i in range(len(accXx3)): accErrorX3.append( 2.0*np.fabs(gaccXx3[i+grid_shift]-np.fabs(accMagX3[i]))/( np.fabs(gaccXx3[i+grid_shift]) + np.fabs(accMagX3[i]) ) )
accErrorY3 = []
for i in range(len(accYy3)): accErrorY3.append( 2.0*np.fabs(gaccYy3[i+grid_shift]-np.fabs(accMagY3[i]))/( np.fabs(gaccYy3[i+grid_shift]) + np.fabs(accMagY3[i]) ) )
accErrorZ3 = []
for i in range(len(accZz3)): accErrorZ3.append( 2.0*np.fabs(gaccZz3[i+grid_shift]-np.fabs(accMagZ3[i]))/( np.fabs(gaccZz3[i+grid_shift]) + np.fabs(accMagZ3[i]) ) )

length = len(Saxis)
for i in range(length/2 - length/4, length/2 + length/4):
        accErrorX[i] = 0
        accErrorY[i] = 0
        accErrorZ[i] = 0

length = len(Saxis2)
for i in range(length/2 - length/4, length/2 + length/4):
        accErrorX2[i] = 0
        accErrorY2[i] = 0
        accErrorZ2[i] = 0

length = len(Saxis3)
for i in range(length/2 - 3, length/2 + 2):
        accErrorX3[i] = 0
        accErrorY3[i] = 0
        accErrorZ3[i] = 0


plt.semilogy(Saxis,accErrorX,color=triad1)
plt.semilogy(Saxis,accErrorY,color=triad2)
plt.semilogy(Saxis,accErrorZ,color=triad3)

plt.semilogy(Saxis2,accErrorX2,color=triad1)
plt.semilogy(Saxis2,accErrorY2,color=triad2)
plt.semilogy(Saxis2,accErrorZ2,color=triad3)

plt.semilogy(Saxis3,accErrorX3,color=triad1)
plt.semilogy(Saxis3,accErrorY3,color=triad2)
plt.semilogy(Saxis3,accErrorZ3,color=triad3)

#ax4.set_xlim(-0.6, 0.6)
ax4.set_ylim(1.e-3, 1.e1)

plt.ylabel(r'Acc err')
plt.xlabel('Distance (' + unit + ')')

for i in range(len(accXx)): print 'i=', i, 'radpot =', radpot[i+grid_shift], 'saxis =', Saxis[i], 'gaccXx = ', gaccXx[i+grid_shift], 'accXx =', accXx[i]

plt.tight_layout()
plt.savefig(sname + ".pdf")

