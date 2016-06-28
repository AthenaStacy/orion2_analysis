from yt.mods import *
import matplotlib.colorbar as cb
import sys
import fields

fn = 'data.0006.3d.hdf5'
field = 'radiation-energy-density' 
dir = 2

pf = load(fn)
c = (pf.domain_left_edge + pf.domain_right_edge) / 2
pc = PlotCollection(pf,center = c)
p = pc.add_slice(field,dir)
p.modify['grids']()
p._redraw_image()
p.save_image('Blast_Result')

