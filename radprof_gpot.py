import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as pl
#import pylab as pl

from yt.mods import *
import numpy as na
from string import rstrip
import fields
#import fields_iso

au_cm = 1.49597871e13
Msun = 1.98892e33

#wdir = '/work/00863/minerva/orion//maptest_iso_ref_newvel/'
wdir = '/work/00863/minerva/orion/gravpot_6.4pc_256_ref3b/'
#wdir = '/work/00863/minerva/orion/gravpot_6.4pc_256/'
#wdir = '/work/00863/minerva/orion/gravpot_3.2pc_256/'
#wdir = '/work/00863/minerva/orion/gravpot_1.6pc_256/'
#wdir = '/work/00863/minerva/orion/'
rdir = '/home1/00863/minerva/research_programs/'

#file_prefix = 'bin_HR10'
#snapnum = '0007'

file_prefix = 'bin_HR10_ref3_wpot'
snapnum = '7131'

#pf = load("/nobackupp7/astacy/orion/chemtest_iso/data.0021.3d.hdf5")
pf = load(wdir+"data.0000.3d.hdf5")
value, location = pf.h.find_max("Density")
data = pf.h.sphere(location, 2.0/pf['pc'])
pc = PlotCollection(pf, center = location)

print pf.h.field_list

bin_num_square = 128
bin_doub = 128.0
DeltaX = pf.domain_right_edge[0] / bin_doub

x=location[0]
y=location[1]
z=location[2]

massenc_com = []
vx = []
vy = []
vz = []

vx.append(0)
vy.append(0)
vz.append(0)

print 'location =', location, 'x=', x, 'y=', y, 'z=', z, 'vx =', vx, 'vy =', vy, 'vz =', vz

def _rVelocity(field, data):
    '''
    The infall velocity. In this problem the center is at density peak.
    '''
    vr = (data['x-velocity'] - vx[0])*(data['x'] - x) + (data['y-velocity'] - vy[0])*(data['y'] - y) + (data['z-velocity'] - vz[0])*(data['z'] - z)
    vr = (vr) / na.sqrt((data['x'] - x)**2 + (data['y'] - y)**2 + (data['z'] - z)**2)
    vr = vr/1.e5
    return vr
add_field("radial-velocity", function=_rVelocity, take_log=False,
          units=r'\rm{km}/\rm{s}')

def _Radius(field, data):
    '''
    Distance from central point. In this problem the center is at density peak.
    '''
    rad =  na.sqrt((data['x'] - x)**2 + (data['y'] - y)**2 + (data['z'] - z)**2)
    rad = rad/1.5e13
    return rad
add_field("Radius", function=_Radius, take_log=False,
          units=r'\rm{AU}')

def _rotVelocity(field, data):
    '''
    Rotational velocity. In this problem the center is at density peak.
    '''
    vrotx = (data['z-velocity'] - vz[0])*(data['y'] - y) - (data['y-velocity'] - vy[0])*(data['z'] - z)
    vroty = (data['x-velocity'] - vx[0])*(data['z'] - z) - (data['z-velocity'] - vz[0])*(data['x'] - x)
    vrotz = (data['y-velocity'] - vy[0])*(data['x'] - x) - (data['x-velocity'] - vx[0])*(data['y'] - y)
    vrot = na.sqrt(vrotx**2 + vroty**2 + vrotz**2)
    vrot = (vrot) / na.sqrt((data['x'] - x)**2 + (data['y'] - y)**2 + (data['z'] - z)**2)
    vrot = vrot/1.e5
    return vrot
add_field("rot-velocity", function=_rotVelocity, take_log=False,
          units=r'\rm{km}/\rm{s}')

######x vel's
def _xVelocity(field, data):
    velx = (data['x-velocity'] - vx[0])
    return velx/1.e5
add_field("velx", function=_xVelocity, take_log=False,
          units=r'\rm{km}/\rm{s}')

def _xvrad(field, data):
    vradx = data['radial-velocity']*(data['x'] - x)/(au_cm*data['Radius'])
    return vradx
add_field("vradx", function=_xvrad, take_log=False,
          units=r'\rm{km}/\rm{s}')

def _xvrot(field, data):
    vrotx = (data['z-velocity'] - vz[0])*(data['y'] - y) - (data['y-velocity'] - vy[0])*(data['z'] - z)
    return vrotx/1.e5/(au_cm*data['Radius'])
add_field("vrotx", function=_xvrot, take_log=False,
          units=r'\rm{km}/\rm{s}')

def _xmomentum(field, data):
    return (data["X-momentum"]*data["CellMassMsun"]/data["density"])
add_field("xmomentum", function=_xmomentum, take_log=False, units=r'\rm{km}/\rm{s}')

#######y vel's
def _yVelocity(field, data):
    vely = (data['y-velocity'] - vy[0])
    return vely/1.e5
add_field("vely", function=_yVelocity, take_log=False,
          units=r'\rm{km}/\rm{s}')

def _yvrad(field, data):
    vrady = data['radial-velocity']*(data['y'] - y)/(au_cm*data['Radius'])
    return vrady
add_field("vrady", function=_yvrad, take_log=False,
          units=r'\rm{km}/\rm{s}')

def _yvrot(field, data):
    vroty = (data['x-velocity'] - vx[0])*(data['z'] - z) - (data['z-velocity'] - vz[0])*(data['x'] - x)
    return vroty/1.e5/(au_cm*data['Radius'])
add_field("vroty", function=_yvrot, take_log=False,
          units=r'\rm{km}/\rm{s}')

def _ymomentum(field, data):
    return (data["Y-momentum"]*data["CellMassMsun"]/data["density"])
add_field("ymomentum", function=_ymomentum, take_log=False, units=r'\rm{km}/\rm{s}')

######z vel's
def _zVelocity(field, data):
    velz = (data['z-velocity'] - vz[0])
    return velz/1.e5
add_field("velz", function=_zVelocity, take_log=False,
          units=r'\rm{km}/\rm{s}')

def _zvrad(field, data):
    vradz = data['radial-velocity']*(data['z'] - z)/(au_cm*data['Radius'])
    return vradz
add_field("vradz", function=_zvrad, take_log=False,
          units=r'\rm{km}/\rm{s}')

def _zvrot(field, data):
    vrotz = (data['y-velocity'] - vy[0])*(data['x'] - x) - (data['x-velocity'] - vx[0])*(data['y'] - y)
    return vrotz/1.e5/(au_cm*data['Radius'])
add_field("vrotz", function=_zvrot, take_log=False,
          units=r'\rm{km}/\rm{s}')

def _rotVelocity(field, data):
    '''
    Rotational velocity. In this problem the center is at density peak.
    '''
    vel_rot = na.power(data['vrotx']*data['vrotx'] + data['vroty']*data['vroty'] + data['vrotz']*data['vrotz'], 0.5)
    return vel_rot
add_field("rot-velocity", function=_rotVelocity, take_log=False,
          units=r'\rm{km}/\rm{s}')

def _zmomentum(field, data):
    return (data["Z-momentum"]*data["CellMassMsun"]/data["density"])
add_field("zmomentum", function=_zmomentum, take_log=False, units=r'\rm{km}/\rm{s}')

################read in other data files#######################################

with open(rdir + file_prefix + '_gas_' + snapnum +'.dat', "r") as f:
        gasdatB = [map(float, line.split()) for line in f]

with open(rdir + file_prefix + '_energy_' + snapnum +'.dat', "r") as f:
        edatB = [map(float, line.split()) for line in f]

with open(rdir + file_prefix + '_gpot_' + snapnum +'.dat', "r") as f:
        potdatB = [map(float, line.split()) for line in f]

radB = []
nhB = []
densB = []
tempB = []
h2B = []
elecB = []
hdB = []
rmencB = []
csB = []
vradB = []
vrotB = []
turbB = []
mencB = []
presB = []
nh_nhB = []
temp_nhB = []
xmomB = []
ymomB = []
zmomB = []
eintB = []
gpotB = []
gaccB = []

radpot = []
gpotxL = []
gpotxR = []
gpotyL = []
gpotyR = []
gpotzL = []
gpotzR = []
gaccxLF = []
gaccxRF = []
gaccyLF = []
gaccyRF = []
gacczLF = []
gacczRF = []
#gaccxL = np.zeros(bin_num_square)
#gaccxR = np.zeros(bin_num_square)
#gaccyL = np.zeros(bin_num_square)
#gaccyR = np.zeros(bin_num_square)
#gacczL = np.zeros(bin_num_square)
#gacczR = np.zeros(bin_num_square)
gaccxL = []
gaccxR = []
gaccyL = []
gaccyR = []
gacczL = []
gacczR = []


narr1 = 200
line_arr = []
fac = 1.2195*1.67e-24
narr_vel = 200

narr = 128
for i in range(narr):
	line_arr = gasdatB[i]
	radB.append(line_arr[0])
	nhB.append(line_arr[1])
        densB.append(line_arr[1]*fac)
	tempB.append(line_arr[2])
        h2B.append(line_arr[3])
        elecB.append(line_arr[4])
        hdB.append(line_arr[5])
	mencB.append(line_arr[6])
	presB.append(1.38e-16*nhB[i]*tempB[i])

for i in range(narr):
        line_arr = edatB[i]
        eintB.append(line_arr[0])
        xmomB.append(line_arr[1])
        ymomB.append(line_arr[2])
        zmomB.append(line_arr[3])
        gpotB.append(line_arr[4])
for i in range(narr):
        if i < 1:
		gaccB.append( (gpotB[1] - gpotB[0]) / DeltaX ) 
	if i >= 1:
		gaccB.append( (gpotB[i] - gpotB[i-1]) / DeltaX ) 

for i in range(narr):
        line_arr = potdatB[i]
	radpot.append(line_arr[0])
        gpotxL.append(line_arr[1])
        gpotxR.append(line_arr[2])
        gpotyL.append(line_arr[3])
        gpotyR.append(line_arr[4])
        gpotzL.append(line_arr[5])
        gpotzR.append(line_arr[6])
        gaccxL.append(line_arr[7])
        gaccxR.append(line_arr[8])
        gaccyL.append(line_arr[9])
        gaccyR.append(line_arr[10])
        gacczL.append(line_arr[11])
        gacczR.append(line_arr[12])

for i in range(narr):
        gaccxLF.append( np.abs(gpotxL[narr-1] - gpotxR[narr-1]) / (radpot[narr-1] ) )
        gaccxRF.append( np.abs(gpotxR[narr-1] - gpotxL[narr-1]) / (radpot[narr-1] ) )
        gaccyLF.append( np.abs(gpotyL[narr-1] - gpotyR[narr-1]) / (radpot[narr-1] ) )
        gaccyRF.append( np.abs(gpotyR[narr-1] - gpotyL[narr-1]) / (radpot[narr-1] ) )
        gacczLF.append( np.abs(gpotzL[narr-1] - gpotzR[narr-1]) / (radpot[narr-1] ) )
        gacczRF.append( np.abs(gpotzR[narr-1] - gpotzL[narr-1]) / (radpot[narr-1] ) )

print 'gpotxL =', gpotxL
print 'gaccxL =', gaccxL
print 'gaccxLF =', gaccxLF

##############################################################################

rmin = 100.
rmax = 3.e5
rmax_arr = pf.domain_right_edge[0]
nmin = 1e2
nmax = 1e16
#bin_num1 = 7500
bin_num1 = 200
bin_num2 = 200

i_doub = 1.0

rad_arr = []
rad_neg = []
gpot_arr = []
gacc_arr = []
massenc_arr = []

gpotxL_arr = []
gaccxL_arr = []
gpotxR_arr = []
gaccxR_arr = []
gpotyL_arr = []
gaccyL_arr = []
gpotyR_arr = []
gaccyR_arr = []
gpotzL_arr = []
gacczL_arr = []
gpotzR_arr = []
gacczR_arr = []

dfac = 0.99
location = [0.0,0.0,0.0]

for i in range(bin_num_square):
        if i > 0:
		Lval_old = Lval
		Rval_old = Rval
        Lval = -rmax_arr  * (i_doub/bin_doub)
        Rval =  rmax_arr  * (i_doub/bin_doub)

	if i > 0:
                data_in = data

        data = pf.h.region(location, [dfac*Lval, dfac*Lval, dfac*Lval], [dfac*Rval, dfac*Rval, dfac*Rval])
	#massenc_arr.append(data["CellMassGrams"].sum() / 1.98892e33 ) #WARNING - must modify fields.py appropriately for CellMassGrams to be correct	 
	massenc_arr.append(data["CellMassMsun"].sum() )
	if i == 0:
		gpot_arr.append( (data['gravitational-potential'] * data['CellMassMsun']).sum() / (massenc_arr[i]) )
                gpotxL_arr.append( (data['gravitational-potential'] * data['CellMassMsun']).sum() / (massenc_arr[i]) )
                gpotxR_arr.append( (data['gravitational-potential'] * data['CellMassMsun']).sum() / (massenc_arr[i]) )
                gpotyL_arr.append( (data['gravitational-potential'] * data['CellMassMsun']).sum() / (massenc_arr[i]) )
                gpotyR_arr.append( (data['gravitational-potential'] * data['CellMassMsun']).sum() / (massenc_arr[i]) )
                gpotzL_arr.append( (data['gravitational-potential'] * data['CellMassMsun']).sum() / (massenc_arr[i]) )
                gpotzR_arr.append( (data['gravitational-potential'] * data['CellMassMsun']).sum() / (massenc_arr[i]) )		

        if i > 0:
                data_ann = pf.h.boolean([data, "NOT", data_in])
		gpot_arr.append( (data_ann['gravitational-potential'] * data_ann['CellMassMsun']).sum() / data_ann['CellMassMsun'].sum() )

		data_subtract = pf.h.region(location, [dfac*Lval_old, dfac*Lval, dfac*Lval], [dfac*Rval, dfac*Rval, dfac*Rval])
                dataxL  = pf.h.boolean([data, "NOT", data_subtract])
		if(dataxL['CellMassMsun'].sum() > 0):
			gpotxL_arr.append( (dataxL['gravitational-potential'] * dataxL['CellMassMsun']).sum() / dataxL['CellMassMsun'].sum() )
		else:
			gpotxL_arr.append(0)

		data_subtract = pf.h.region(location, [dfac*Lval, dfac*Lval, dfac*Lval], [dfac*Rval_old, dfac*Rval, dfac*Rval])
                dataxR  = pf.h.boolean([data, "NOT", data_subtract])
		if(dataxR['CellMassMsun'].sum() > 0):
                	gpotxR_arr.append( (dataxR['gravitational-potential'] * dataxR['CellMassMsun']).sum() / dataxR['CellMassMsun'].sum() )
		else: 
			gpotxR_arr.append(0)

		data_subtract = pf.h.region(location, [dfac*Lval, dfac*Lval_old, dfac*Lval], [dfac*Rval, dfac*Rval, dfac*Rval])
                datayL  = pf.h.boolean([data, "NOT", data_subtract])
		if(datayL['CellMassMsun'].sum() > 0):
                	gpotyL_arr.append( (datayL['gravitational-potential'] * datayL['CellMassMsun']).sum() / datayL['CellMassMsun'].sum() )
		else:
			gpotyL_arr.append(0)

		data_subtract = pf.h.region(location, [dfac*Lval, dfac*Lval, dfac*Lval], [dfac*Rval, dfac*Rval_old, dfac*Rval])
                datayR  = pf.h.boolean([data, "NOT", data_subtract])
		if(datayR['CellMassMsun'].sum() > 0):
                	gpotyR_arr.append( (datayR['gravitational-potential'] * datayR['CellMassMsun']).sum() / datayR['CellMassMsun'].sum() )
		else:
			gpotyR_arr.append(0)

		data_subtract = pf.h.region(location, [dfac*Lval, dfac*Lval, dfac*Lval_old], [dfac*Rval, dfac*Rval, dfac*Rval])
                datazL  = pf.h.boolean([data, "NOT", data_subtract])
		if(datazL['CellMassMsun'].sum() > 0):
                	gpotzL_arr.append( (datazL['gravitational-potential'] * datazL['CellMassMsun']).sum() / datazL['CellMassMsun'].sum() )
		else:
			gpotzL_arr.append(0)	

		data_subtract = pf.h.region(location, [dfac*Lval, dfac*Lval, dfac*Lval], [dfac*Rval, dfac*Rval, dfac*Rval_old])
                datazR  = pf.h.boolean([data, "NOT", data_subtract])
		if(datazL['CellMassMsun'].sum() > 0):
                	gpotzR_arr.append( (datazR['gravitational-potential'] * datazR['CellMassMsun']).sum() / datazR['CellMassMsun'].sum() )
		else:
			gpotzR_arr.append(0)
        rad_arr.append(Rval / au_cm )
        rad_neg.append(Lval / au_cm )
        i_doub = i_doub + 1.0


for i in range(bin_num_square):
        if i == 0:
                gacc_arr.append(np.abs(gpot_arr[1] - gpot_arr[0]) / DeltaX)
                gaccxL_arr.append(np.abs(gpotxR_arr[1] - gpotxL_arr[0]) / DeltaX )
                gaccxR_arr.append(np.abs(gpotxL_arr[1] - gpotxR_arr[0]) / DeltaX)
                gaccyL_arr.append(np.abs(gpotyR_arr[1] - gpotyL_arr[0]) / DeltaX)
                gaccyR_arr.append(np.abs(gpotyL_arr[1] - gpotyR_arr[0]) / DeltaX)
                gacczL_arr.append(np.abs(gpotzR_arr[1] - gpotzL_arr[0]) / DeltaX)
                gacczR_arr.append(np.abs(gpotzL_arr[1] - gpotzR_arr[0]) / DeltaX)
        if i > 0:
                gacc_arr.append( np.abs(gpot_arr[i] - gpot_arr[i-1]) / DeltaX )
                gaccxL_arr.append( np.abs(gpotxL_arr[i] - gpotxL_arr[i-1]) / DeltaX )
                gaccxR_arr.append( np.abs(gpotxR_arr[i] - gpotxR_arr[i-1]) / DeltaX )
                gaccyL_arr.append( np.abs(gpotyL_arr[i] - gpotyL_arr[i-1]) / DeltaX )
                gaccyR_arr.append( np.abs(gpotyR_arr[i] - gpotyR_arr[i-1]) / DeltaX )
                gacczL_arr.append( np.abs(gpotzL_arr[i] - gpotzL_arr[i-1]) / DeltaX )
                gacczR_arr.append( np.abs(gpotzR_arr[i] - gpotzR_arr[i-1]) / DeltaX )


#####################################################################################

print 'gpotxL_arr =', gpotxL_arr

gpotB_min = 1.e30
gpot_arr_min = 1.e30
gpotxL_arr_min = 1.e30
gpotxR_arr_min = 1.e30
gpotyL_arr_min = 1.e30
gpotyR_arr_min = 1.e30
gpotzL_arr_min = 1.e30
gpotzR_arr_min = 1.e30
gpotxL_min = 1.e30
gpotxR_min = 1.e30
gpotyL_min = 1.e30
gpotyR_min = 1.e30
gpotzL_min = 1.e30
gpotzR_min = 1.e30

for i in range(bin_num_square):
        if gpot_arr[i] < gpot_arr_min:
		gpot_arr_min = gpot_arr[i]
        if gpotB[i] < gpotB_min:
		gpotB_min = gpotB[i]
        if gpotxL_arr[i] < gpotxL_arr_min and np.abs(gpotxL_arr[i]) > 0.1:
                gpotxL_arr_min = gpotxL_arr[i]
        if gpotxR_arr[i] < gpotxR_arr_min and np.abs(gpotxR_arr[i]) > 0.1:
                gpotxR_arr_min = gpotxR_arr[i]
        if gpotyL_arr[i] < gpotyL_arr_min and np.abs(gpotyL_arr[i]) > 0.1:
                gpotyL_arr_min = gpotyL_arr[i]
        if gpotyR_arr[i] < gpotyR_arr_min and np.abs(gpotyR_arr[i]) > 0.1:
                gpotyR_arr_min = gpotyR_arr[i]
        if gpotzL_arr[i] < gpotzL_arr_min and np.abs(gpotzL_arr[i]) > 0.1:
                gpotzL_arr_min = gpotzL_arr[i]
        if gpotzR_arr[i] < gpotzR_arr_min and np.abs(gpotzR_arr[i]) > 0.1:
                gpotzR_arr_min = gpotzR_arr[i]

for i in range(bin_num_square):
        if gpotxL[i] < gpotxL_min:
                gpotxL_min = gpotxL[i]
        if gpotxR[i] < gpotxR_min:
                gpotxR_min = gpotxR[i]
        if gpotyL[i] < gpotyL_min:
                gpotyL_min = gpotyL[i]
        if gpotyR[i] < gpotyR_min:
                gpotyR_min = gpotyR[i]
        if gpotzL[i] < gpotzL_min:
                gpotzL_min = gpotzL[i]
        if gpotzR[i] < gpotzR_min:
                gpotzR_min = gpotzR[i]


for i in range(bin_num_square):
	gpot_arr[i] = gpot_arr[i] - gpot_arr_min
        gpotxL_arr[i] = gpotxL_arr[i] - gpotxL_arr_min
        gpotxR_arr[i] = gpotxR_arr[i] - gpotxR_arr_min
        gpotyL_arr[i] = gpotyL_arr[i] - gpotyL_arr_min
        gpotyR_arr[i] = gpotyR_arr[i] - gpotyR_arr_min
        gpotzL_arr[i] = gpotzL_arr[i] - gpotzL_arr_min
        gpotzR_arr[i] = gpotzR_arr[i] - gpotzR_arr_min

	gpotB[i] = gpotB[i] - gpotB_min
	gpotxL[i] = gpotxL[i] - gpotxL_min
	gpotxR[i] = gpotxR[i] - gpotxR_min
	gpotyL[i] = gpotyL[i] - gpotyL_min
	gpotyR[i] = gpotyR[i] - gpotyR_min
	gpotzL[i] = gpotzL[i] - gpotzL_min
	gpotzR[i] = gpotzR[i] - gpotzR_min

print 'gpotxR_arr_min', gpotxR_arr_min, 'gpotyR_arr_min', gpotyR_arr_min, 'gpotzR_arr_min', gpotyL_arr_min
print 'gpotxR_min', gpotxR_min, 'gpotyR_min', gpotyR_min, 'gpotzR_min', gpotyL_min
print 'gpot_arr =', gpot_arr
print 'gacc_arr =', gacc_arr
print 'gpotxL_arr =', gpotxL_arr
print 'gpotxL_arr =', gpotxL_arr
print 'gacxL_arr =', gaccxL_arr
print 'gpotyR_arr =', gpotyR_arr
print 'gaccyR_arr =', gaccyR_arr
#####################################################################################
#calculate errors
gpot_err = []
gacc_err = []
gaccxL_err = []
gaccxR_err = []
gaccyL_err = []
gaccyR_err = []
gacczL_err = []
gacczR_err = []

gaccxL_err_pn = []
gaccxR_err_pn = []
gaccyL_err_pn = []
gaccyR_err_pn = []
gacczL_err_pn = []
gacczR_err_pn = []

for i in range(bin_num_square):
#        if i == 0:
#                gpot_err.append(1.e-10)
#                gacc_err.append(1.e-10)
#                gaccxL_err.append(1.e-10)
#                gaccxR_err.append(1.e-10)
#                gaccyL_err.append(1.e-10)
#                gaccyR_err.append(1.e-10)
#                gacczL_err.append(1.e-10)
#                gacczR_err.append(1.e-10)
	gpot_err.append(np.abs(gpot_arr[i] - gpotB[i])/gpotB[i])
        #gacc_err.append(np.abs(gacc_arr[i] - gaccB[i])/gaccB[i])
	gacc_err.append(np.abs(gacc_arr[i] - gaccB[i]))
	gaccxL_err.append(np.abs(gaccxL_arr[i] - gaccxL[i]))
        gaccxR_err.append(np.abs(gaccxR_arr[i] - gaccxR[i]))
        gaccyL_err.append(np.abs(gaccyL_arr[i] - gaccyL[i]))
        gaccyR_err.append(np.abs(gaccyR_arr[i] - gaccyR[i]))
        gacczL_err.append(np.abs(gacczL_arr[i] - gacczL[i]))
        gacczR_err.append(np.abs(gacczR_arr[i] - gacczR[i]))
#####################################################################################

rmin = -2.e5
rmax = 2.e5

pl.subplot(321)
pl.plot(rad_neg, gpotxL, 'r--')
pl.plot(radB, gpotxR, 'r--')
pl.plot(rad_neg, gpotxL_arr, 'k')
pl.plot(radB, gpotxR_arr, 'k')
ax = pl.gca()
ax.set_yscale('log')
ax.set_xlabel('Radius [AU]', fontsize=9)
ax.set_ylabel('X-face grav. pot.', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((rmin, rmax, 1.e10, 1.e13))

pl.subplot(322)
pl.plot(rad_neg, gaccxL_arr, 'k')
pl.plot(rad_arr, gaccxR_arr, 'k')
pl.plot(rad_neg, gaccxL,'r--')
pl.plot(rad_arr, gaccxR,'r--')
pl.plot(rad_neg, gaccxL_err,'b:')
pl.plot(rad_arr, gaccxR_err,'b:')
pl.plot(rad_neg, gaccxLF,'g--')
pl.plot(rad_arr, gaccxRF,'g--')
ax = pl.gca()
ax.set_yscale('log')
ax.set_xlabel('Radius [AU]', fontsize=9)
ax.set_ylabel('X-face grav. acc.', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((rmin, rmax, 1.e-9, 1.e-4))

pl.subplot(323)
pl.plot(rad_neg, gpotyL,'r--')
pl.plot(radB, gpotyR, 'r--')
pl.plot(rad_neg, gpotyL_arr, 'k')
pl.plot(radB, gpotyR_arr, 'k')
ax = pl.gca()
ax.set_yscale('log')
ax.set_xlabel('Radius [AU]', fontsize=9)
ax.set_ylabel('Y-face grav. pot.', fontsize=9)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
pl.axis((rmin, rmax, 1.e10, 1.e13))

pl.subplot(324)
pl.plot(rad_neg, gaccyL_arr, 'k')
pl.plot(rad_arr, gaccyR_arr, 'k')
pl.plot(rad_neg, gaccyL,'r--')
pl.plot(rad_arr, gaccyR,'r--')
pl.plot(rad_neg, gaccyL_err,'b:')
pl.plot(rad_arr, gaccyR_err,'b:')
pl.plot(rad_neg, gaccyLF,'g--')
pl.plot(rad_arr, gaccyRF,'g--')
ax = pl.gca()
ax.set_yscale('log')
ax.set_xlabel('Radius [AU]', fontsize=9)
ax.set_ylabel('Y-face grav. acc.', fontsize=9)
pl.xticks(fontsize=9)
pl.yticks(fontsize=9)
pl.axis((rmin, rmax, 1.e-9, 1.e-4))

pl.subplot(325)
pl.plot(rad_neg, gpotzL,'r--')
pl.plot(radB, gpotzR, 'r--')
pl.plot(rad_neg, gpotzL_arr, 'k')
pl.plot(radB, gpotzR_arr, 'k')
ax = pl.gca()
ax.set_yscale('log')
ax.set_xlabel('Radius [AU]', fontsize=9)
ax.set_ylabel('Z-face grav. pot.', fontsize=9)
pl.xticks(fontsize=9)
pl.yticks(fontsize=9)
pl.axis((rmin, rmax, 1.e10, 1.e13))

pl.subplot(326)
pl.plot(rad_neg, gacczL_arr, 'k')
pl.plot(rad_arr, gacczR_arr, 'k')
pl.plot(rad_neg, gacczL,'r--')
pl.plot(rad_arr, gacczR,'r--')
pl.plot(rad_neg, gacczL_err,'b:')
pl.plot(rad_arr, gacczR_err,'b:')
pl.plot(rad_neg, gacczLF,'g--')
pl.plot(rad_arr, gacczRF,'g--')
ax = pl.gca()
ax.set_yscale('log')
ax.set_xlabel('Radius [AU]', fontsize=9)
ax.set_ylabel('Z--face grav. pot.', fontsize=9)
pl.xticks(fontsize=9)
pl.yticks(fontsize=9)
pl.axis((rmin, rmax, 1.e-9, 1.e-4))


pl.savefig('plot.ps')

M_i = data["CellMassMsun"].sum()
print "total mass = ", M_i
