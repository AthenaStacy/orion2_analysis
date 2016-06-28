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

shrink_fac = 1.0
my_level = 0

shrink_fac2 = 1.0
my_level2 = 0

datanum1 = '0000'  
datanum2 = '5979' 

odir1 = '/work/00863/minerva/orion_Btest/'
odir2 = '/work/00863/minerva/orion/bfield_comp_vpot_fwards/'
odir3 = '/work/00863/minerva/orion/bfield_comp_vpot_fwards/'
rdir = '/home1/00863/minerva/gadget_runs/'

file_level = ''

file1 = odir1 + 'data.' + datanum1 + '.3d.hdf5'
#file1 = "/work/00863/minerva/orion_Btest/bfield_homolog/" + 'data.' + datanum1 + '.3d.hdf5'
#file1 = "/work/00863/minerva/orion_Btest/" + 'data.' + datanum1 + '.3d.hdf5'
#file1 = "/work/00863/minerva/orion_Gtest/" + 'data.' + datanum1 + '.3d.hdf5'

file2 = rdir+"bfield_comp_vpot_"+datanum2

file_dens = odir2+"gadget2dens_ideal_"+file_level+datanum2
file_vx = odir2+"gadget2vx_ideal_"+file_level+datanum2
file_vy = odir2+"gadget2vy_ideal_"+file_level+datanum2
file_vz = odir2+"gadget2vz_ideal_"+file_level+datanum2
file_egy = odir2+"gadget2egy_ideal_"+file_level+datanum2
file_bx = odir2+"gadget2bx_ideal_"+file_level+datanum2
file_by = odir2+"gadget2by_ideal_"+file_level+datanum2
file_bz = odir2+"gadget2bz_ideal_"+file_level+datanum2

#file_bx_vp = odir3+"gadget2bx_ideal_"+file_level+datanum2
#file_by_vp = odir3+"gadget2by_ideal_"+file_level+datanum2
#file_bz_vp = odir3+"gadget2bz_ideal_"+file_level+datanum2

file_bx_vp = odir3+"gadget2bx_conv_"+file_level+datanum2
file_by_vp = odir3+"gadget2by_conv_"+file_level+datanum2
file_bz_vp = odir3+"gadget2bz_conv_"+file_level+datanum2

#file_bx_vp = odir3+"gadget2ax_ideal_"+file_level+datanum2
#file_by_vp = odir3+"gadget2ay_ideal_"+file_level+datanum2
#file_bz_vp = odir3+"gadget2az_ideal_"+file_level+datanum2

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
vx_gadget = vx_gadget/1.e5

vy_gadget = np.fromfile(file_vy, dtype=float, count=-1)
vy_gadget = numpy.reshape(vy_gadget, [gref,gref,gref])
vy_gadget = vy_gadget/1.e5

vz_gadget = np.fromfile(file_vz, dtype=float, count=-1)
vz_gadget = numpy.reshape(vz_gadget, [gref,gref,gref])
vz_gadget = vz_gadget/1.e5

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
temp_gadget = egy_gadget * (gamma - 1) / rho_nh / kB

max_arr = np.argmax(dens_gadget)
max_arr = np.unravel_index(max_arr, [gref,gref,gref])
print 'max_arr =', max_arr
#################################################################



DeltaX = 2. * pf.domain_right_edge[0] / bin_doub * shrink_fac
rmax = pf.domain_right_edge[0] / pc_cm * shrink_fac

value, location = pf.h.find_max("Density")
data = pf.h.sphere(location, 2*rmax/pf['pc'])
print 'max density location = ', location

location = [0,0,0]
pc = PlotCollection(pf, center = location)

print pf.h.field_list

x=location[0]
y=location[1]
z=location[2]

################################################################################################
#add fields
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
####################################################################################################

take_ray = 1
take_profile = 0

print 'temp_min=', min(data['Temperature']), 'temp_max=', max(data['Temperature'])
print 'vrad_min=', min(data['radial-velocity']), 'vrad_max=', max(data['radial-velocity'])
print 'dens_min=', min(data['Density']), 'dens_max=', max(data['Density'])


# Chose center of data (selecting highest density point)
#rayx = pf.h.ortho_ray(0,(0,0))

# Save the distance axes since they are used often
unit = 'pc'

oref = gref
max_arr_orion = [oref/2, oref/2, oref/2]
rmin = -rmax

cube = pf.h.covering_grid(level=my_level,left_edge=pf.h.domain_left_edge*shrink_fac,dims=[oref,oref,oref])
cube2 = pf.h.covering_grid(level=my_level2,left_edge=pf.h.domain_left_edge*shrink_fac2,dims=[oref,oref,oref])

axX = cube['x']*pf[unit]
density = cube['density']

axX_ref = cube2['x']*pf[unit]
density_ref = cube2['density']

max_arr_orion = np.argmax(density)
max_arr_orion = np.unravel_index(max_arr_orion, [oref,oref,oref])
print 'max_arr_orion =', max_arr_orion

max_arr_orion_ref = np.argmax(density_ref)
max_arr_orion_ref = np.unravel_index(max_arr_orion_ref, [oref,oref,oref])
print 'max_arr_orion_ref =', max_arr_orion_ref

rad =  cube['Radius']
bfield =  cube['Bmag']
vx =      cube['x-velocity'] / 1.e5
vy =      cube['y-velocity'] / 1.e5
vz =      cube['z-velocity'] / 1.e5
temp =    cube['Temperature']

bfield_ref = cube2['Bmag']

if take_ray == 1:
	axX = axX[:, max_arr_orion[1], max_arr_orion[2]]
	density = density[:, max_arr_orion[1], max_arr_orion[2]]
	bfield =  bfield[:, max_arr_orion[1], max_arr_orion[2]]
	vx = vx[:, max_arr_orion[1], max_arr_orion[2]]
	vy = vy[:, max_arr_orion[1], max_arr_orion[2]]
	vz = vz[:, max_arr_orion[1], max_arr_orion[2]]
	temp = temp[:, max_arr_orion[1], max_arr_orion[2]]
	
	density_ref = density_ref[:, max_arr_orion_ref[1], max_arr_orion_ref[2]]
        bfield_ref =  bfield_ref[:, max_arr_orion_ref[1], max_arr_orion_ref[2]]
	
	dens_gadget = dens_gadget[:, max_arr[1], max_arr[2]]
	vx_gadget = vx_gadget[:, max_arr[1], max_arr[2]] 
        vy_gadget = vy_gadget[:, max_arr[1], max_arr[2]] 
        vz_gadget = vz_gadget[:, max_arr[1], max_arr[2]] 
        temp_gadget = temp_gadget[:, max_arr[1], max_arr[2]]
	B_gadget = B_gadget[:, max_arr[1], max_arr[2]]
	B_vp = B_vp[:, max_arr[1], max_arr[2]]

center_expect = max_arr[0]
cen_shift = max_arr_orion[0] - center_expect

x_gadget = np.zeros(gref)
gref_doub = float(gref)
for i in range(gref):
	i_doub = float(i)
	x_gadget[i] = (rmin*pc_cm  + (2*cen_shift+1)*DeltaX/2 + i_doub*DeltaX) / pc_cm

if take_profile == 1:
	bin_num_doub = float(bin_num)
	array_min = 0
	array_max = rmax

	array_bin = []
	for i in range(bin_num):
        	i_doub = float(i)
        	array_bin.append(array_min + (array_max-array_min)*i_doub/bin_num_doub)

	lmax = oref 
	rad_gadget = np.zeros((gref,gref,gref))
	x1_gadget  = np.zeros(gref)
	x2_gadget  = np.zeros(gref)
	x3_gadget  = np.zeros(gref)

        for i in range(gref):
		i_doub = float(i)
                print 'i=', i
		x1_gadget[i] = -rmax + (i_doub) * DeltaX / pc_cm + DeltaX/2/pc_cm 
                x2_gadget[i] = -rmax + (i_doub) * DeltaX / pc_cm + DeltaX/2/pc_cm
                x3_gadget[i] = -rmax + (i_doub) * DeltaX / pc_cm + DeltaX/2/pc_cm

        for i in range(gref):
                for j in range(gref):
                        for k in range(gref):
				rad_gadget[i,j,k] = x1_gadget[i]*x1_gadget[i] + x2_gadget[j]*x2_gadget[j] +x3_gadget[k]*x3_gadget[k] 
				rad_gadget[i,j,k] = np.sqrt(rad_gadget[i,j,k])

	num        =  np.zeros(bin_num)
	dens_sum   = np.zeros(bin_num)
	vx_sum     = np.zeros(bin_num)
	vy_sum     = np.zeros(bin_num)
        vz_sum     = np.zeros(bin_num)
        temp_sum   = np.zeros(bin_num)
        bfield_sum = np.zeros(bin_num)

        num_gadget        =  np.zeros(bin_num)
        dens_sum_gadget   = np.zeros(bin_num)
        vx_sum_gadget     = np.zeros(bin_num)
        vy_sum_gadget     = np.zeros(bin_num)
        vz_sum_gadget     = np.zeros(bin_num)
        temp_sum_gadget   = np.zeros(bin_num)
        B_sum_gadget = np.zeros(bin_num)
        B_vp_sum_gadget = np.zeros(bin_num)

	for i in range(lmax):
        	for j in range(lmax):
                	for k in range(lmax):
                        	for l in range(bin_num-1):
                                        if rad[i,j,k] > array_bin[l] and rad[i,j,k] < array_bin[l+1]:
                                                num[l] = num[l] +  1.0
                                        	dens_sum[l] = dens_sum[l] + density[i,j,k]
                                        	vx_sum[l] = vx_sum[l] + vx[i,j,k]
                                        	vy_sum[l] = vy_sum[l] + vy[i,j,k]
                                                vz_sum[l] = vz_sum[l] + vz[i,j,k]
                                                temp_sum[l] = temp_sum[l] + temp[i,j,k]
                                                bfield_sum[l] = bfield_sum[l] + bfield[i,j,k]
                                        if rad_gadget[i,j,k] > array_bin[l] and rad_gadget[i,j,k] < array_bin[l+1]:
                                                num_gadget[l] = num_gadget[l] +  1.0
                                                dens_sum_gadget[l] = dens_sum_gadget[l] + dens_gadget[i,j,k]
                                                vx_sum_gadget[l] = vx_sum_gadget[l] + vx_gadget[i,j,k]
                                                vy_sum_gadget[l] = vy_sum_gadget[l] + vy_gadget[i,j,k]
                                                vz_sum_gadget[l] = vz_sum_gadget[l] + vz_gadget[i,j,k]
                                                temp_sum_gadget[l] = temp_sum_gadget[l] + temp_gadget[i,j,k]
                                                B_sum_gadget[l] = B_sum_gadget[l] + B_gadget[i,j,k]
						B_vp_sum_gadget[l] = B_vp_sum_gadget[l] + B_vp[i,j,k]

	axX    = array_bin
	x_gadget= array_bin
	rmin = 0

        density= np.zeros(bin_num)
        vx     = np.zeros(bin_num)
        vy     = np.zeros(bin_num)
        vz     = np.zeros(bin_num)
        temp   = np.zeros(bin_num)
        bfield = np.zeros(bin_num)
        dens_gadget   = np.zeros(bin_num)
        vx_gadget     = np.zeros(bin_num)
        vy_gadget     = np.zeros(bin_num)
        vz_gadget     = np.zeros(bin_num)
        temp_gadget   = np.zeros(bin_num)
        B_gadget      = np.zeros(bin_num)
        B_vp          = np.zeros(bin_num)

        for i in range(lmax):
		density[i] = dens_sum[i]   / num[i]
        	bfield[i]  = bfield_sum[i] / num[i]
		vx[i]      = vx_sum[i]     / num[i]      
        	vy[i]      = vy_sum[i]     / num[i]     
        	vz[i]      = vz_sum[i]     / num[i]    
		temp[i]    = temp_sum[i]   / num[i]

                dens_gadget[i] = dens_sum_gadget[i] / num_gadget[i]
                B_gadget[i]    = B_sum_gadget[i]    / num_gadget[i]
                B_vp[i]        = B_vp_sum_gadget[i] / num_gadget[i]
                vx_gadget[i]   = vx_sum_gadget[i]   / num_gadget[i]
                vy_gadget[i]   = vy_sum_gadget[i]   / num_gadget[i]
                vz_gadget[i]   = vz_sum_gadget[i]   / num_gadget[i]
                temp_gadget[i] = temp_sum_gadget[i] / num_gadget[i]

nh = density * rho_nh
nh_gadget = dens_gadget * rho_nh
if take_ray == 1:
	nh_ref = density_ref * rho_nh

print 'bfield =', bfield
print 'velx = ', vx
print 'velz = ', vz

###################################################################
read_sph = 0
narr1 = bin_num/2

if read_sph == 1:
	with open(file2) as f:
        	Bdat = [map(float, line.split()) for line in f]

	line_arr = []
	x_gadget = []
	nh_gadget = []
	B_gadget = []
	B_vp = []
	vx_gadget = []
	vy_gadget = []
	vz_gadget = []
	temp_gadget = []

	for i in range(narr1):
        	line_arr = Bdat[i]
        	x_gadget.append(line_arr[0])
        	nh_gadget.append(line_arr[1])
        	B_gadget.append(line_arr[5])
		vx_gadget.append(line_arr[6])
		vy_gadget.append(line_arr[7])
		vz_gadget.append(line_arr[8])
		temp_gadget.append(line_arr[9])
        for i in range(narr1):
		B_vp.append(line_arr[5])
#####################################################################

dens_min = 1.e-28
dens_max = 1.e-26

nmin = 1.e3
nmax = 1.e13

bmin = 1.e-20
bmax = 1.e-10

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

print 'x_gadget =', x_gadget
print 'nh_gadget =', nh_gadget
if take_profile == 1:
	print 'num_gadget=', num_gadget
	print 'num=', num

if take_ray == 1:
	print 'axX.shape =', axX.shape
	print 'x_gadget.shape =', x_gadget.shape
	print 'nh_gadget.shape =', nh_gadget.shape

if plot_vel == 1:

        pl.subplot(321)
        #pl.plot(axX, density,'k')
        pl.plot(axX, nh,'k')
        pl.plot(x_gadget, nh_gadget)
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
	pl.plot(x_gadget, nh_gadget)
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
        #ax.set_ylabel('Vec. Pot.', fontsize=9)
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
	pl.plot(nh_gadget, B_gadget, 'b')
	pl.plot(nh_gadget, B_vp, 'r')
        #pl.plot(nh_right, bfield_right, 'k')
        #pl.plot(nh_gadget_right, B_gadget_right, 'b')
        #pl.plot(nh_gadget_right, B_vp_right, 'r')
	ax = pl.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel('nh [cm^-3]', fontsize=9)
	ax.set_ylabel('B-field [G]', fontsize=9)
        #ax.set_ylabel('Vec. Pot.', fontsize=9)
	pl.xticks(fontsize=10)
	pl.yticks(fontsize=10)
	pl.axis((nmin, nmax, bmin, bmax))

        pl.subplot(224)
        #pl.plot(nh_ref, bfield_ref, 'k')
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
