from yt.mods import *
import sys
import pylab as plt
import numpy as np
from tabular import *
from math import pi

G = 6.6726e-8
c_iso = 1.8e4
t0 = 1.29302897416e+12 # s
msol = 1.98892e33

# look for sink particle file dumps
pfiles=os.popen("ls data.*.sink").read().strip(' \n').split("\n")
mass = []
massana = []
time = []
for file in pfiles:
    f=open(file)
    f.readline()
    time.append(load(file[0:len(file)-5]+'.hdf5').current_time)
    mass.append(float(f.readline().split(' ')[0])/msol)
    massana.append(DELTA**2*alpha[0]*(DELTA+v[0])*c_iso**3*(time[-1]+t0)/G/msol)

mass = np.array(mass)
massana = np.array(massana)

if os.path.exists("data.0021.3d.sink") or True:
    L2Error = np.sum(np.sqrt(((mass-massana)/massana)**2))/len(mass)
    LInfError = np.max(np.abs((mass-massana)/massana))
else:
    L2Error = 1e100
    LInfError = 1e100

# as of 10/10/2011 (the first cut at the code for this problem), 
# the L2 and Linf errors are:
# L2Error = 0.0367239649748
# LInfError = 0.0828887270188

# on 3/8/2011, ATM raised the tolerances to reflect changes in the
# accretion model. The new errors are:
#L2Error = 0.0725675475759
#LInfError = 0.108279486431

print L2Error
print LInfError
if L2Error < 0.045 and LInfError < 0.065:
    print 'Pass!'
    writer=open('test_pass', 'w')
    writer.close()
else:
    print 'FAIL!'

