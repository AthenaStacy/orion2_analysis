import matplotlib
matplotlib.use('ps')

from yt.mods import *
import numpy as na
import pylab as pl
#from string import rstrip
import fields_bfield
import tracer_def

MRH = 0.76
h2conv = (2.0/1.3158)
heconv = (4.0/1.3158)
pc_to_cm = 3.08567758e18

def load_file(dir, datanum):

	pf = load(dir + 'data.' + datanum + ".3d.hdf5")
	value, location = pf.h.find_max("density")
	data = pf.h.sphere(location, 3.0*pc_to_cm)
	print ('location =', location)

	data = data.cut_region(["obj['density'] > 4.0e-19"])

	bin_num2 = 50
	nmin = 1.e-20
	nmax = 1.e-9

	#p = pc.add_profile_sphere(2., "pc", ['Density', 'Bmag'], weight='CellMassMsun', x_bins = bin_num2, center = location, x_bounds = [nmin, nmax])
	plot = ProfilePlot(data, "n_h", ['Bmag2'], n_bins = bin_num2, weight_field="cell_volume")

	profile = plot.profiles[0]
	nh = profile.x
	bmag = np.sqrt(profile['Bmag2'])
	#bmag =profile['Bmag']

	return nh, bmag

'''
den_gadget = [5.94338E-25, 1.39513E-24, 4.96551E-24, 1.19045E-23, 2.51893E-23, 5.32961E-23, 1.03767E-22, 1.93797E-22, 3.77363E-22, 8.49733E-22, 2.12331E-21, 5.08935E-21, 1.24548E-20, 3.38221E-20, 7.93962E-20, 1.54543E-19, 5.27591E-19, 1.21298E-18, 3.2936E-18, 9.13224E-18, 1.97235E-17,4.34896E-17, 7.95201E-17, 1.45396E-16]
b_gadget = [3.47E-11, 7.13E-11, 2.04E-10, 6.55E-10, 1.83E-09, 4.68E-09, 1.27E-08, 3.35E-08, 1.07E-07, 2.27E-07, 5.51E-07, 1.26E-06, 2.82E-06, 6.85E-06, 1.49E-05, 2.74E-05, 7.02E-05, 1.48E-04, 3.05E-04, 7.61E-04, 0.001360717, 0.002118973, 0.003210206, 0.004600912]
'''

den_gadget = [2.5353638e-26,   4.0743740e-26,   6.5475901e-26,   1.0522091e-25,   1.6909185e-25,
   2.7173356e-25,   4.3668058e-25,   7.0175335e-25,   1.1277300e-24,   1.8122819e-24,
   2.9123683e-24,   4.6802263e-24,   7.5212052e-24,   1.2086709e-23,   1.9423557e-23,
   3.1213992e-23,   5.0161425e-23,   8.0610306e-23,   1.2954216e-22,   2.0817650e-22,
   3.3454324e-22,   5.3761679e-22,   8.6395979e-22,   1.3883983e-21,   2.2311802e-21,
   3.5855452e-21,   5.7620351e-21,   9.2596937e-21,   1.4880483e-20,   2.3913202e-20,
   3.8428943e-20,   6.1755955e-20,   9.9242924e-20,   1.5948515e-19,   2.5629530e-19,
   4.1187115e-19,   6.6188432e-19,   1.0636592e-18,   1.7093182e-18,   2.7469048e-18,
   4.4143251e-18,   7.0938947e-18,   1.1400016e-17,  1.8320029e-17,   2.9440592e-17,
   4.7311558e-17,   7.6030522e-17,   1.2218232e-16,   1.9634917e-16]

den_gadget = np.asarray(den_gadget)
nh_gadget = den_gadget / 1.67e-24 * .76

b_gadget = [4.0405926e-12,   5.1267997e-12,   6.8901578e-12,   9.5083125e-12,   1.3468342e-11,
   1.9389195e-11,   2.8110089e-11,   3.8397422e-11,   5.0160857e-11,   7.8511168e-11,
   1.2891716e-10,   2.3446289e-10,   4.2062647e-10,   8.1297054e-10,   1.7883228e-09,
   2.8125155e-09,   7.2976238e-09,   1.4638656e-08,   3.0624952e-08,   5.7106274e-08,
   9.9304308e-08,   1.6778045e-07,   2.9406524e-07,   4.9065663e-07,   7.1840074e-07,
   1.0299489e-06,   1.6076816e-06,   2.5628763e-06,   3.8764825e-06,   5.9786005e-06,
   8.9687602e-06,   1.3232808e-05,   1.9969269e-05,   2.9737359e-05,   4.2966386e-05,
   6.3547005e-05,   9.3293571e-05,   0.00014333730,   0.00020191851,   0.00030134383,
   0.00045198875,   0.00066851871,    0.0010280378,    0.0015085077,    0.0021590413,
    0.0029591806,    0.0037736749,    0.0049152486,    0.0077241384]

b_gadget = np.asarray(b_gadget)

bfac_g = 1.8 #np.sqrt(4*3.14159)
b_gadget = b_gadget / bfac_g 
dfac_gadget = np.power(nh_gadget, -0.6666667)

for i in range(len(b_gadget)):
	print('nh_gadget = ', nh_gadget[i], 'bfield = ', b_gadget[i])

rho_turk = []
u_turk = []
with open('turk_bfield.csv', "r") as f:
	turkdat = [line.split(',') for line in f]
for i in range(len(turkdat)):
	line = turkdat[i]
	rho_turk.append(float(line[0]))
	u_turk.append(float(line[1]))
rho_turk = np.asarray(rho_turk)
nh_turk = rho_turk / 1.67e-24 * .76

u_turk = u_turk * (np.power(rho_turk, 1.333333))
b_turk = np.power(u_turk, 0.5) * 8.0 * 3.14159 * 1.e3 / bfac_g 
#extra 1e3 factor is to make Turk initial b-field (1e-14) compatible with Gadget initial bfield (1e-11)
dfac_turk = np.power(nh_turk, -0.66667)

nh_lowres = []
b_lowres = []
with open('gadget_lowres_bfield.csv', "r") as f:
        lowres_dat = [line.split(',') for line in f]
for i in range(len(lowres_dat)):
        line = lowres_dat[i]
        nh_lowres.append(float(line[0]))
        b_lowres.append(float(line[1])/bfac_g)
dfac_lowres = np.power(nh_lowres, -0.66667)

print('nh_turk = ', nh_turk)
print('b_turk = ', b_turk)

fac = 1.e10
fac2 = 1.e9

x_init = 1.e-25
x = []
y = []
y2 = []
for i in range (0,16):
        x_here = x_init * 10**i
        x.append(x_here)
        y.append(fac*5.e-20 * (x_here/x_init) ** 0.66667)
        y2.append(fac*5.e-20 * (x_here/x_init) ** 0.9)
x = np.asarray(x)
dfac_line = np.power(x, -0.6666667)

print ('xline =', x)
print ('yline =', y*dfac_line)
print ('y2line =', y2*dfac_line)

dir_name =  '/nobackupp12/astacy/popiii_Bscope12/'
datanum = '0010'
nh0, bmag0 = load_file(dir_name, datanum)
dfac0 = np.power(nh0, -0.6666667)

#dir_name = '/nobackupp7/astacy/popiii_bfac3/'
dir_name = '/nobackupp12/astacy/popiii_bfieldA/'

datanum  = '0100'
nh, bmag = load_file(dir_name, datanum)
dfac = np.power(nh, -0.6666667)

datanum = '0200'
nh2, bmag2 = load_file(dir_name, datanum)
dfac2 = np.power(nh2, -0.6666667)

datanum = '0245'  #first sink forms
nh3, bmag3 = load_file(dir_name, datanum)
dfac3 = np.power(nh3, -0.666667)

datanum = '0250'
nh4, bmag4 = load_file(dir_name, datanum)
dfac4 = np.power(nh4, -0.6666667)

datanum = '0350'
nh5, bmag5 = load_file(dir_name, datanum)
dfac5 = np.power(nh5, -0.6666667)

print('nh0 = ', nh0)
print('bmag0 = ', bmag0)

dfac_line2 = np.power(1.e-10, -0.66667)
dfac_line3 = np.power(1.e-7, -0.66667)
dfac_line4 = np.power(1.e-1, -0.66667)

nmin = 1.e-1
nmax = 1.e14
#bmin = 0.1*min(bmag5*dfac5)
bmin = min(b_gadget*dfac_gadget)
bmax = 1.e4*max(bmag5*dfac5)

bmax = bmax.value

fsize = 12

#pl.plot(nh0, bmag0*fac2*dfac0,'r.', linewidth=2.0, label='initial B-field, no div. clean')
#line1, = pl.plot(nh, bmag*dfac,'k', linewidth=2.0, label="initial B-field")
#pl.plot(nh2, bmag2*dfac2,'g--', linewidth=2.0, label='1000 yr before first sink forms')
pl.plot(nh_gadget, b_gadget*dfac_gadget, 'k--', linewidth=2.0, label='-9000 yr')
line2, = pl.plot(nh3, bmag3*dfac3,'c', linewidth=2.0, label='0 yr, first sink forms')
pl.plot(nh4, bmag4*dfac4,'-.', linewidth=2.0, label='100 yr')
line3, = pl.plot(nh5, bmag5*dfac5,'m--', linewidth=2.0, label='1000 yr')
#pl.plot(x, y*dfac_line, 'b:')
#pl.plot(x, y2*dfac_line, 'b:')
pl.plot(nh_turk, b_turk*dfac_turk, 'b:', linewidth=2.0, label=r'Turk ea 2011 ($\times 10^{3}$)')
#pl.plot(nh_lowres, b_lowres*dfac_lowres, 'g:', linewidth=2.0, label='Gadget low-res')

ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'n$_{\rm H}$ [cm$^3$]', fontsize=fsize)
ax.set_ylabel(r'B-field magnitude / n$_{\rm H}^{2/3}$ [G / cm$^{-3}$]', fontsize=fsize)
#ax.annotate(r'$B\propto\rho^{2/3}$', xy=(1e-10,1e-10*fac*dfac_line2),xytext=(1.e-13, 1e7), fontsize=12)
#ax.annotate(r'$B\propto\rho^{0.9}$', xy=(1e-10,1e-10*fac*dfac_line2),xytext=(1.e-13, 1.e10), fontsize=12)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((nmin, nmax, bmin, bmax))
pl.tight_layout()
pl.legend(loc=2, fontsize=10)

pl.savefig('nprof_bfield.eps')

