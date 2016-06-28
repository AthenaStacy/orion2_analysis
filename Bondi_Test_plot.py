from yt.mods import *
import pylab as pl
import fields
from readBondi import readBondi

def _rVelocity(field, data):
    '''
    The infall velocity. In this problem the center is at 0, 0, 0.
    '''
    vr = data['x-velocity']*data['x'] + data['y-velocity']*data['y'] + data['z-velocity']*data['z']
    vr = abs(vr) / na.sqrt(data['x']**2 + data['y']**2 + data['z']**2)
    return vr
add_field("radial-velocity", function=_rVelocity, take_log=False,
          units=r'\rm{cm}/\rm{s}')

frame = sys.argv[1] # in python, sys.argv[0] is actually the name of the script
fn = 'data.0' + str(frame) + '.3d.hdf5'

pf = load(fn)
pc = PlotCollection(pf, center = [0, 0, 0])

# define some parameters
r_bondi = 2*3.926e12
r_big = 4.0 * r_bondi
r_small = 0.25 * r_bondi
rho_inf = 1e-18
cs = 1.3e7

# make plot and adjust
#p = pc.add_profile_sphere(3.2e13, "cm", ['Radius', 'radial-velocity'], weight='CellVolume', x_bounds = (0.25*8e12, 4*8e12))
#u = p.data['radial-velocity'] / cs
p1 = pc.add_profile_sphere(3.2e13, "cm", ['Radius', 'density'], weight='CellVolume', x_bounds = (0.25*8e12, 8*8e12))
p2 = pc.add_profile_sphere(3.2e13, "cm", ['Radius', 'radial-velocity'], weight='CellVolume', x_bounds = (0.25*8e12, 8*8e12))
a = p1.data['density'] / rho_inf
u = p2.data['radial-velocity'] / cs
r = p1.data['Radius'] / r_bondi

pl.subplot(211)
pl.plot(r, a, 'k.')
xs, us, alphas = readBondi('tabular.h')
pl.plot(xs, alphas, 'r')
ax = pl.gca()
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylabel('Non-dimensional density')
pl.axis((0.25, 4.0, 1.0, 20.0))

pl.subplot(212)
pl.plot(r, u, 'k.')
pl.plot(xs, us, 'r')
ax = pl.gca()
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylabel('Non-dimensional infall velocity')
ax.set_xlabel('Bondi Radii')
pl.axis((0.25, 4.0, 0.01, 2.5))

# save
pl.savefig('prof')
