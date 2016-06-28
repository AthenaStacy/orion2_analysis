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

# look for the last dump
fn=os.popen("ls data.*.hdf5").read().strip(' \n').split("\n")[-1]

pf = load(fn)
cen = (pf.domain_left_edge + pf.domain_right_edge) / 2
ray = pf.h.ray([pf.domain_left_edge[0],cen[1],cen[2]], [pf.domain_right_edge[0],cen[1],cen[2]])

x = ray['x']
den = ray['density']
vx = np.abs(ray['x-velocity'])
time = pf.current_time

# compute analytic solution from self-similar model of Shu (1977)
# tabulated shu solution
xana = c_iso*(time+t0)*np.linspace(DELTA, (len(alpha)+1)*DELTA, len(alpha), endpoint=True)
denana = np.array(alpha)/(4*pi*G*(time+t0)**2)
vxana = np.array(v)*c_iso
xana = np.concatenate((-xana[::-1],xana))
denana = np.concatenate((denana[::-1],denana))
vxana = np.concatenate((vxana[::-1],vxana))

kau = 1.49598e16 # cm
plt.clf()
plt.plot(x/kau,den)
plt.plot(xana/kau,denana)
plt.gca().set_xlabel(r'x [10$^3$ AU]')
plt.gca().set_ylabel(r'$\rho$ [g cm$^{-3}$]')
plt.gca().set_yscale('log')
plt.gca().set_xlim(-3.5,3.5)
plt.savefig('density.png')

plt.clf()
plt.plot(x/kau,vx)
plt.plot(xana/kau,vxana)
plt.gca().set_xlabel('x [1000 AU]')
plt.gca().set_ylabel(r'-v$_r$ [cm s$^{-1}$]')
plt.gca().set_yscale('log')
plt.gca().set_xlim(-3.5,3.5)
plt.savefig('v.png')

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

plt.clf()
plt.plot(time,mass)
plt.plot(time,massana)
plt.gca().set_xlabel('t [s]')
plt.gca().set_ylabel(r'M$_{\rm sink}$ [M$_\odot$]')
plt.savefig('m.png')
