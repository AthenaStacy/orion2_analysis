import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as pl
#import pylab as pl

from yt.mods import *
import numpy as na
from string import rstrip
import fields
import fields_bfield

pc_cm = 3.08567758e18
au_cm = 1.49597871e13
Msun = 1.98892e33
rho_nh = 1./(1.2195*1.6726e-24)
kB = 1.3806e-16

#shrink_fac = 0.125
#my_level = 3

shrink_fac = 0.0625
my_level = 4

shrink_fac2 = 0.015625
my_level2 = 6

#datanum1 = '0000'
#datanum2 = '0000'
#datanum1 = '0020'   #--corresponding times and data steps for orion128
#datanum2 = '0746'
#datanum1 = '0023'
#datanum2 = '0827'  #--

#datanum1 = '0025'
#datanum2 = '0200'
datanum1 = '0029'  #--corresponding times and data steps for orion256
datanum2 = '0679' 

odir1 = '/work/00863/minerva/orion/bfield_comp_orion256/'
odir2 = '/work/00863/minerva/orion/'
odir3 = '/work/00863/minerva/orion//bfield_comp_vpot/'
rdir = '/home1/00863/minerva/gadget_runs/'

file_level = 'lev4_'

file1 = odir1 + 'data.' + datanum1 + '.3d.hdf5'
#file1 = "/work/00863/minerva/orion_Btest/bfield_homolog/" + 'data.' + datanum1 + '.3d.hdf5'
#file1 = "/work/00863/minerva/orion_Btest/" + 'data.' + datanum1 + '.3d.hdf5'
#file1 = "/work/00863/minerva/orion_Gtest/" + 'data.' + datanum1 + '.3d.hdf5'

file2 = rdir+"bfield_comp_"+datanum2

file_dens = odir2+"gadget2dens_ideal_"+file_level+datanum2
file_vx = odir2+"gadget2vx_ideal_"+file_level+datanum2
file_vy = odir2+"gadget2vy_ideal_"+file_level+datanum2
file_vz = odir2+"gadget2vz_ideal_"+file_level+datanum2
file_egy = odir2+"gadget2egy_ideal_"+file_level+datanum2
file_bx = odir2+"gadget2bx_ideal_"+file_level+datanum2
file_by = odir2+"gadget2by_ideal_"+file_level+datanum2
file_bz = odir2+"gadget2bz_ideal_"+file_level+datanum2

file_bx_vp = odir3+"gadget2bx_ideal_"+file_level+datanum2
file_by_vp = odir3+"gadget2by_ideal_"+file_level+datanum2
file_bz_vp = odir3+"gadget2bz_ideal_"+file_level+datanum2

pf = load(file1)

print pf.h.field_list

dims = list(pf.domain_dimensions)
bin_num = dims[0] * 1
bin_doub = float(bin_num)

##############################################################
#read from gadget mapped-onto-grid files
gref = bin_num
print 'gref =', gref

dens_gadget = np.fromfile(file_dens, dtype=float, count=-1)
dens_gadget = numpy.reshape(dens_gadget, [gref,gref,gref])

vx_gadget = np.fromfile(file_vx, dtype=float, count=-1)
vx_gadget = numpy.reshape(vx_gadget, [gref,gref,gref])

vy_gadget = np.fromfile(file_vy, dtype=float, count=-1)
vy_gadget = numpy.reshape(vy_gadget, [gref,gref,gref])

vz_gadget = np.fromfile(file_vz, dtype=float, count=-1)
vz_gadget = numpy.reshape(vz_gadget, [gref,gref,gref])

egy_gadget = np.fromfile(file_egy, dtype=float, count=-1)
egy_gadget = numpy.reshape(egy_gadget, [gref,gref,gref])

Bx_gadget = np.fromfile(file_bx, dtype=float, count=-1)
Bx_gadget = numpy.reshape(Bx_gadget, [gref,gref,gref])

By_gadget = np.fromfile(file_by, dtype=float, count=-1)
By_gadget = numpy.reshape(By_gadget, [gref,gref,gref])

Bz_gadget = np.fromfile(file_bz, dtype=float, count=-1)
Bz_gadget = numpy.reshape(Bz_gadget, [gref,gref,gref])

B_gadget = np.sqrt(Bx_gadget*Bx_gadget + By_gadget*By_gadget + Bz_gadget*Bz_gadget)

Bx_vp = np.fromfile(file_bx_vp, dtype=float, count=-1)
Bx_vp = numpy.reshape(Bx_vp, [gref,gref,gref])

By_vp = np.fromfile(file_by_vp, dtype=float, count=-1)
By_vp = numpy.reshape(By_vp, [gref,gref,gref])

Bz_vp = np.fromfile(file_bz, dtype=float, count=-1)
Bz_vp = numpy.reshape(Bz_vp, [gref,gref,gref])

B_vp = np.sqrt(Bx_vp*Bx_vp + By_vp*By_vp + Bz_vp*Bz_vp)

gamma = 1.1
n_gadget = dens_gadget * rho_nh
temp_gadget = egy_gadget * (gamma - 1) / rho_nh / kB

max_arr = np.argmax(dens_gadget)
max_arr = np.unravel_index(max_arr, [gref,gref,gref])
print 'max_arr =', max_arr
#################################################################



DeltaX = 2. * pf.domain_right_edge[0] / bin_doub * shrink_fac
rmax = pf.domain_right_edge[0] / pc_cm * shrink_fac

value, location = pf.h.find_max("Density")
data = pf.h.sphere(location, rmax/pf['pc'])
print 'max density location = ', location

pc = PlotCollection(pf, center = location)

print pf.h.field_list


x=location[0]
y=location[1]
z=location[2]

def _rVelocity(field, data):
    '''
    The infall velocity. In this problem the center is at density peak.
    '''
    vr = data['x-velocity']*data['x'] + data['y-velocity']*data['y'] + data['z-velocity']*data['z']
    vr = (vr) / na.sqrt(data['x']**2 + data['y']**2 + data['z']**2)
    vr = vr/1.e5
    return vr
add_field("radial-velocity", function=_rVelocity, take_log=False,
          units=r'\rm{km}/\rm{s}')

def _Radius(field, data):
    '''
    Distance from central point. In this problem the center is at density peak.
    '''
    rad =  na.sqrt((data['x'] )**2 + (data['y'])**2 + (data['z'])**2)
    rad = rad/pc_cm
    return rad
add_field("Radius", function=_Radius, take_log=False,
          units=r'\rm{AU}')

take_ray = 1
take_profile = 0

print 'temp_min=', min(data['Temperature']), 'temp_max=', max(data['Temperature'])
print 'vrad_min=', min(data['radial-velocity']), 'vrad_max=', max(data['radial-velocity'])
print 'dens_min=', min(data['Density']), 'dens_max=', max(data['Density'])

if take_profile == 1:
	p1 = pc.add_profile_sphere(rmax, "pc", ['Radius', 'number-density'], weight='CellVolume', x_bins=bin_num)
	p2 = pc.add_profile_sphere(rmax, "pc", ['Radius', 'radial-velocity', 'x-velocity', 'y-velocity', 'z-velocity', 'Temperature'], weight='CellVolume', x_log=False, x_bins=bin_num)
	p3 = pc.add_profile_sphere(rmax, "pc", ['Radius', 'density', 'Temperature', 'Bmag', 'ufield', 'MagneticEnergy', 'DivB', 'absDivB'], weight='CellVolume', x_log=False, x_bins=bin_num)

	rad_arr = p3.data['Radius']

# Chose center of data (selecting highest density point)
#rayx = pf.h.ortho_ray(0,(0,0))

# Save the distance axes since they are used often
unit = 'pc'

oref = gref
cube = pf.h.covering_grid(level=my_level,left_edge=pf.h.domain_left_edge*shrink_fac,dims=[oref,oref,oref])
cube2 = pf.h.covering_grid(level=my_level2,left_edge=pf.h.domain_left_edge*shrink_fac2,dims=[oref,oref,oref])

if take_ray == 1:
	axX = cube['x']*pf[unit]
	density = cube['density']

        axX_ref = cube2['x']*pf[unit]
        density_ref = cube2['density']

	max_arr_orion = np.argmax(density)
	max_arr_orion = np.unravel_index(max_arr_orion, [oref,oref,oref])
	print 'max_arr_orion =', max_arr_orion

        max_arr_orion_ref = np.argmax(density_ref)
        max_arr_orion_ref = np.unravel_index(max_arr_orion_ref, [oref,oref,oref])
        print 'max_arr_orion =', max_arr_orion

	bfield =  cube['Bmag']
	ufield =  cube['ufield'] 
	#absDivB = cube['absDivB']
        vx =      cube['x-velocity'] / 1.e5
        vy =      cube['y-velocity'] / 1.e5
        vz =      cube['z-velocity'] / 1.e5
        temp =    cube['Temperature']

	bfield_ref = cube2['Bmag']	

	axX = axX[:, max_arr_orion[1], max_arr_orion[2]]
	density = density[:, max_arr_orion[1], max_arr_orion[2]]
	bfield =  bfield[:, max_arr_orion[1], max_arr_orion[2]]
	ufield = ufield[:, max_arr_orion[1], max_arr_orion[2]]
	#absDivB = absDivB[:, max_arr_orion[1], max_arr_orion[2]]
	vx = vx[:, max_arr_orion[1], max_arr_orion[2]]
	vy = vy[:, max_arr_orion[1], max_arr_orion[2]]
	vz = vz[:, max_arr_orion[1], max_arr_orion[2]]
	temp = temp[:, max_arr_orion[1], max_arr_orion[2]]
	
	density_ref = density_ref[:, max_arr_orion_ref[1], max_arr_orion_ref[2]]
        bfield_ref =  bfield_ref[:, max_arr_orion_ref[1], max_arr_orion_ref[2]]

	rmin = -rmax
	
	n_gadget = n_gadget[:, max_arr[1], max_arr[2]]
	vx_gadget = vx_gadget[:, max_arr[1], max_arr[2]] / 1.e5
        vy_gadget = vy_gadget[:, max_arr[1], max_arr[2]] / 1.e5
        vz_gadget = vz_gadget[:, max_arr[1], max_arr[2]] / 1.e5
        temp_gadget = temp_gadget[:, max_arr[1], max_arr[2]]
	B_gadget = B_gadget[:, max_arr[1], max_arr[2]]
	B_vp = B_vp[:, max_arr[1], max_arr[2]]

center_expect = max_arr[0]
cen_shift = max_arr_orion[0] - center_expect

x_gadget = np.zeros(gref)
gref_doub = float(gref)
for i in range(gref):
	i_doub = float(i)
	#x_gadget[i] = (rmin + (rmax-rmin) * i_doub/gref_doub)
	#x_gadget[i] = (rmin*pc_cm + DeltaX/2. + i_doub*DeltaX) / pc_cm
	x_gadget[i] = (rmin*pc_cm  + (2*cen_shift+1)*DeltaX/2 + i_doub*DeltaX) / pc_cm

if take_profile == 1:
	axX =     p3.data['Radius']
	density = p3.data['density']
        bfield =  p3.data['Bmag']
        ufield =  p3.data['ufield']
        absDivB = p3.data['absDivB']	
	vx =      p2.data['x-velocity'] / 1.e5
        vy =      p2.data['y-velocity'] / 1.e5
        vz =      p2.data['z-velocity'] / 1.e5
	temp =    p2.data['Temperature'] 
	rmin = 0

nh = density * rho_nh
nh_ref = density_ref * rho_nh
ufield_norm = ufield / (density**1.33333)
#errDivB = absDivB * DeltaX / bfield

print 'bfield =', bfield
print 'ufield_norm =', ufield_norm
#print 'error divB =', errDivB
print 'total divB = ', data['DivB'].sum()
print 'velx = ', vx
print 'velz = ', vz

###################################################################
read_sph = 0

narr1 = bin_num
if take_profile == 1: 
	narr1 = narr1*2

if read_sph == 1:
	with open(file2) as f:
        	Bdat = [map(float, line.split()) for line in f]

	line_arr = []
	x_gadget = []
	n_gadget = []
	B_gadget = []
	vx_gadget = []
	vy_gadget = []
	vz_gadget = []
	temp_gadget = []

	for i in range(narr1):
        	line_arr = Bdat[i]
        	x_gadget.append(line_arr[0])
        	n_gadget.append(line_arr[1])
        	B_gadget.append(line_arr[5])
		vx_gadget.append(line_arr[6])
		vy_gadget.append(line_arr[7])
		vz_gadget.append(line_arr[8])
		temp_gadget.append(line_arr[9])
#####################################################################

dens_min = 1.e-28
dens_max = 1.e-26

nmin = 1.e3
nmax = 1.e13

bmin = 1.e-22
bmax = 1.e-16

umin = 1.e8
umax = 1.e15

DivBmin = -1.e-22
DivBmax = 1.e-22

tmin = 0
tmax = 2000

vmin = -10
vmax = 10

plot_B = 1
plot_vel = 0

size = 8.5

##############################################################################
#divide data into left and right sides...
x_gadget_left = []
axX_left = []
nh_left = []
n_gadget_left = []
B_gadget_left = []
B_vp_left = []
bfield_left = []

x_gadget_right = []
axX_right = []
nh_right = []
n_gadget_right = [] 
B_gadget_right = []
B_vp_right = []
bfield_right = []

for i in range(gref):
	if x_gadget[i] <= 0:
		x_gadget_left.append(x_gadget[i])
		n_gadget_left.append(n_gadget[i])
		B_gadget_left.append(B_gadget[i])
		B_vp_left.append(B_vp[i])
        if x_gadget[i] > 0:
                x_gadget_right.append(x_gadget[i])
                n_gadget_right.append(n_gadget[i])
		B_gadget_right.append(B_gadget[i])
		B_vp_right.append(B_vp[i])
for i in range(oref):
        if axX[i] <= 0:
                axX_left.append(axX[i])
                nh_left.append(nh[i])
                bfield_left.append(bfield[i])
        if axX[i] > 0:
                axX_right.append(axX[i])
                nh_right.append(nh[i])
                bfield_right.append(bfield[i])
for i in range(oref):
	print 'i =', i, 'x =', axX[i], 'nh =', nh[i], 'bfield =', bfield[i]
########################################################################################

print 'x_gadget =', x_gadget
print 'n_gadget =', n_gadget

print 'axX.shape =', axX.shape
print 'x_gadget.shape =', x_gadget.shape
print 'n_gadget.shape =', n_gadget.shape

if plot_vel == 1:

        pl.subplot(321)
        #pl.plot(axX, density,'k')
        pl.plot(axX, nh,'k')
        pl.plot(x_gadget, n_gadget)
        ax = pl.gca()
        ax.set_yscale('log')
        ax.set_xlabel('Radius [pc]', fontsize=size)
        ax.set_ylabel('nh [cm^-3]', fontsize=size)
        pl.xticks(fontsize=size)
        pl.yticks(fontsize=size)
        #pl.axis((rmin, rmax, dens_min, dens_max))
        pl.axis((rmin, rmax, nmin, nmax))

        pl.subplot(322)
        pl.plot(axX, temp, 'k')
        pl.plot(x_gadget, temp_gadget)
        ax = pl.gca()
        ax.set_ylabel('Temp [K]', fontsize=size)
        ax.set_xlabel('Radius [pc]', fontsize=size)
        pl.xticks(fontsize=size)
        pl.yticks(fontsize=size)
        pl.axis((rmin, rmax, tmin, tmax))

        pl.subplot(323)
        pl.plot(axX, vx, 'k')
        pl.plot(x_gadget, vx_gadget)
        ax = pl.gca()
        ax.set_ylabel('x-velocity [km / s]', fontsize=size)
        ax.set_xlabel('Radius [pc]', fontsize=size)
        pl.xticks(fontsize=size)
        pl.yticks(fontsize=size)
        pl.axis((rmin, rmax, vmin, vmax))

        pl.subplot(324)
        pl.plot(axX, vy, 'k')
        pl.plot(x_gadget, vy_gadget)
        ax = pl.gca()
        ax.set_ylabel('y-velocity [km / s]', fontsize=size)
        ax.set_xlabel('Radius [pc]', fontsize=size)
        pl.xticks(fontsize=size)
        pl.yticks(fontsize=size)
        pl.axis((rmin, rmax, vmin, vmax))

        pl.subplot(325)
        pl.plot(axX, vz, 'k')
        pl.plot(x_gadget, vz_gadget)
        ax = pl.gca()
        ax.set_ylabel('z-velocity [km / s]', fontsize=size)
        ax.set_xlabel('Radius [pc]', fontsize=size)
        pl.xticks(fontsize=size)
        pl.yticks(fontsize=size)
        pl.axis((rmin, rmax, vmin, vmax))



if plot_B == 1:

	pl.subplot(221)
	pl.title('Black=Orion, Blue=Roberts, Red=Vec. Pot.', fontsize = 10)
	#pl.plot(axX, density,'k')
	pl.plot(axX, nh,'k')
	pl.plot(x_gadget, n_gadget)
	ax = pl.gca()
	#ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel('X [pc]', fontsize=9)
	ax.set_ylabel('nh [cm^-3]', fontsize=9)
	pl.xticks(fontsize=10)
	pl.yticks(fontsize=10)
	#pl.axis((rmin, rmax, dens_min, dens_max))
	pl.axis((rmin, rmax, nmin, nmax))

	pl.subplot(222)
	pl.plot(axX, bfield, 'k')
	pl.plot(x_gadget, B_gadget)
	pl.plot(x_gadget, B_vp, 'r')
	ax = pl.gca()
	#ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_ylabel('B-field [G]', fontsize=9)
	ax.set_xlabel('X [pc]', fontsize=9)
	pl.xticks(fontsize=10)
	pl.yticks(fontsize=10)
	pl.axis((rmin, rmax, bmin, bmax))

#pl.subplot(222)
#pl.plot(axX, ufield_norm,'k')
#ax = pl.gca()
##ax.set_xscale('log')
#ax.set_yscale('log')
#ax.set_xlabel('axX [pc]', fontsize=9)
#ax.set_ylabel('U / rho^(4/3) [cgs]', fontsize=9)
#pl.xticks(fontsize=10)
#pl.yticks(fontsize=10)
#pl.axis((rmin, rmax, umin, umax))

	pl.subplot(223)
	pl.plot(nh, bfield, 'k')
	pl.plot(n_gadget, B_gadget, 'b')
	pl.plot(n_gadget, B_vp, 'r')
        #pl.plot(nh_right, bfield_right, 'k')
        #pl.plot(n_gadget_right, B_gadget_right, 'b')
        #pl.plot(n_gadget_right, B_vp_right, 'r')
	ax = pl.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel('nh [cm^-3]', fontsize=9)
	ax.set_ylabel('B-field [G]', fontsize=9)
	pl.xticks(fontsize=10)
	pl.yticks(fontsize=10)
	pl.axis((nmin, nmax, bmin, bmax))

        pl.subplot(224)
        pl.plot(nh_ref, bfield_ref, 'k')
        ax = pl.gca()
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('nh [cm^-3]', fontsize=9)
        ax.set_ylabel('B-field [G], refined', fontsize=9)
        pl.xticks(fontsize=10)
        pl.yticks(fontsize=10)
        pl.axis((nmin, 10*nmax, bmin, bmax))

	#pl.subplot(224)
	#pl.plot(axX, errDivB,'k')
	#ax = pl.gca()
	#ax.set_xscale('log')
	#ax.set_yscale('log')
	#ax.set_xlabel('X [pc]', fontsize=9)
	#ax.set_ylabel('error DivB (divB * cell_size / B)', fontsize=9)
	#pl.xticks(fontsize=10)
	#pl.yticks(fontsize=10)
	#pl.axis((rmin, rmax, 1.e-10, 1.e0))

pl.savefig('bfield.eps')

M_i = data["CellMassMsun"].sum()
print "total mass = ", M_i
