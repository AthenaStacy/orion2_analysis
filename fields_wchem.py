from yt.mods import *
import numpy as na
from string import rstrip


def _H2(field, data):
   out = data['tracer1']
   return out
#add_field('H2', function=_H2, validators=[ValidateSpatial(1,['tracer1'])], take_log=False)
add_field('H2', function=_H2, take_log=False)

def _HD(field, data):
   out = data['tracer2']
   return out
#add_field('HD', function=_HD, validators=[ValidateSpatial(1,['tracer2'])], take_log=False)
add_field('HD', function=_HD, take_log=False)

def _elec(field, data):
   out = data['tracer3']
   return out
#add_field('elec', function=_elec, validators=[ValidateSpatial(1,['tracer3'])], take_log=False)
add_field('elec', function=_elec, take_log=False)
