import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as plt
import numpy as np
from yt.mods import *

import fields
import fields_bfield
import fft

pc_cm = 3.08567758e18
au_cm = 1.49597871e13
Msun = 1.98892e33
rho_nh = 1./(1.2195*1.6726e-24)

fname1 = 'bmin_vort'
fname2 = 'bmax_vort'
fname3 = 'brms_vort'

###################################################################
narr1 = 9

num_rot = []
bmin_orion_midres = []
bmin_orion_hires = []
bmin_orion = []
bmin_gadget_midres = []
bmin_gadget_hires = []
bmin_anal = []

bmax_orion_midres = []
bmax_orion_hires = []
bmax_orion = []
bmax_gadget_midres = []
bmax_gadget_hires = []
bmax_anal = []

brms_orion_midres = []
brms_orion_hires = []
brms_orion = []
brms_gadget_midres = []
brms_gadget_hires = []
brms_anal = []

with open(fname1, "r") as f:
        Bdat = [map(float, line.split()) for line in f]

for i in range(narr1):
        line_arr = Bdat[i]
        num_rot.append(line_arr[0])
        bmin_orion_midres.append(line_arr[1])
        bmin_orion_hires.append(line_arr[2])
        bmin_gadget_midres.append(line_arr[3])
        bmin_gadget_hires.append(line_arr[4])
        bmin_anal.append(line_arr[5])

f.close

with open(fname2, "r") as f:
        Bdat = [map(float, line.split()) for line in f]

for i in range(narr1):
        line_arr = Bdat[i]
        bmax_orion_midres.append(line_arr[1])
        bmax_orion_hires.append(line_arr[2])
        bmax_gadget_midres.append(line_arr[3])
        bmax_gadget_hires.append(line_arr[4])
        bmax_anal.append(line_arr[5])

f.close

with open(fname3, "r") as f:
        Bdat = [map(float, line.split()) for line in f]

for i in range(narr1):
        line_arr = Bdat[i]
        brms_orion_midres.append(line_arr[1])
        brms_orion_hires.append(line_arr[2])
        brms_gadget_midres.append(line_arr[3])
        brms_gadget_hires.append(line_arr[4])
        brms_anal.append(line_arr[5])

f.close

#####################################################################

do_bfield = 1

print 'num_rot =', num_rot
print 'brms_orion_hires', brms_orion_hires

xmin = 1
xmax = 16
ymin = 1.e-20
ymax = 1.e-16

plt.title('Black = Orion low res, Magenta = Orion high res, Blue = SPH low res, Red = SPH high res, Green = analytic', fontsize = 10)
plt.plot(num_rot, brms_orion_midres,  'k')
plt.plot(num_rot, brms_orion_hires, 'm')
plt.plot(num_rot, brms_gadget_midres)
plt.plot(num_rot, brms_gadget_hires, 'r')
plt.plot(num_rot, brms_anal, 'g--', linewidth=3.0)
ax = plt.gca()
ax.set_yscale('log')
ax.set_xlabel('number of rotations', fontsize=9)
ax.set_ylabel('Bmag', fontsize=9)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.axis((xmin, xmax, ymin, ymax))

plt.savefig("b_evol.eps")


