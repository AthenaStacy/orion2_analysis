from yt.mods import *
import numpy as na
from string import rstrip


def unlimited(right,center,left):
   return 0.5*(right-left)

def _gravpot_G2(field, data):
   out = data['tracer4']
   return out
add_field('gravpot_G2', function=_gravpot_G2, take_log=False)

def _dGdx_G2(field, data):
   sl_left = slice(None,-2,None)
   sl_right = slice(2,None,None)
   ds = data['dx'].flat[0]
   out = na.zeros(data['gravpot_G2'].shape, dtype='float64')
   out[1:-1,1:-1,1:-1] = unlimited(data['gravpot_G2'][sl_right,1:-1,1:-1],data['gravpot_G2'][1:-1,1:-1,1:-1],data['gravpot_G2'][sl_left ,1:-1,1:-1])/ds
   return out
add_field('dGdx_G2', function=_dGdx_G2, validators=[ValidateSpatial(1,['gravpot_G2'])], take_log=False)


def _dGdy_G2(field, data):
   sl_left = slice(None,-2,None)
   sl_right = slice(2,None,None)
   ds = data['dy'].flat[0]
   out = na.zeros(data['gravpot_G2'].shape, dtype='float64')
   out[1:-1,1:-1,1:-1] = unlimited(data['gravpot_G2'][1:-1,sl_right,1:-1],data['gravpot_G2'][1:-1,1:-1,1:-1],data['gravpot_G2'][1:-1,sl_left,1:-1])/ds
   return out
add_field('dGdy_G2', function=_dGdy_G2, validators=[ValidateSpatial(1,['gravpot_G2'])], take_log=False)


def _dGdz_G2(field, data):
   sl_left = slice(None,-2,None)
   sl_right = slice(2,None,None)
   ds = data['dz'].flat[0]
   out = na.zeros(data['gravpot_G2'].shape, dtype='float64')
   out[1:-1,1:-1,1:-1] = unlimited(data['gravpot_G2'][1:-1,1:-1,sl_right],data['gravpot_G2'][1:-1,1:-1,1:-1],data['gravpot_G2'][1:-1,1:-1,sl_left])/ds
   return out
add_field('dGdz_G2', function=_dGdz_G2, validators=[ValidateSpatial(1,['gravpot_G2'])], take_log=False)


def _dGmag_G2(field,data):
       out = na.zeros(data['gravitational-potential'].shape, dtype='float64')
       out[1:-1,1:-1,1:-1] = na.sqrt( data['dGdx_G2'][1:-1,1:-1,1:-1]*data['dGdx_G2'][1:-1,1:-1,1:-1] + data['dGdy_G2'][1:-1,1:-1,1:-1]*data['dGdy_G2'][1:-1,1:-1,1:-1] + data['dGdz_G2'][1:-1,1:-1,1:-1]*data['dGdz_G2'][1:-1,1:-1,1:-1] )
       return out
add_field('dGmag_G2',function=_dGmag_G2, validators=[ValidateSpatial(1,['gravpot_G2'])], take_log=False)

def _dGdx(field, data):
   sl_left = slice(None,-2,None)
   sl_right = slice(2,None,None)
   ds = data['dx'].flat[0]
   out = np.zeros(data['gravitational-potential'].shape, dtype='float64')
   out[1:-1,1:-1,1:-1] = -unlimited(data['gravitational-potential'][sl_right,1:-1,1:-1],data['gravitational-potential'][1:-1,1:-1,1:-1],data['gravitational-potential'][sl_left ,1:-1,1:-1])/ds
   return out
add_field('dGdx', function=_dGdx, validators=[ValidateSpatial(1,['gravitational-potential'])], take_log=False)


def _dGdy(field, data):
   sl_left = slice(None,-2,None)
   sl_right = slice(2,None,None)
   ds = data['dy'].flat[0]
   out = na.zeros(data['gravitational-potential'].shape, dtype='float64')
   out[1:-1,1:-1,1:-1] = -unlimited(data['gravitational-potential'][1:-1,sl_right,1:-1],data['gravitational-potential'][1:-1,1:-1,1:-1],data['gravitational-potential'][1:-1 ,sl_left,1:-1])/ds
   return out
add_field('dGdy', function=_dGdy, validators=[ValidateSpatial(1,['gravitational-potential'])], take_log=False)


def _dGdz(field, data):
   sl_left = slice(None,-2,None)
   sl_right = slice(2,None,None)
   ds = data['dz'].flat[0]
   out = na.zeros(data['gravitational-potential'].shape, dtype='float64')
   out[1:-1,1:-1,1:-1] = -unlimited(data['gravitational-potential'][1:-1,1:-1,sl_right],data['gravitational-potential'][1:-1,1:-1,1:-1],data['gravitational-potential'][1:-1 ,1:-1,sl_left])/ds
   return out
add_field('dGdz', function=_dGdz, validators=[ValidateSpatial(1,['gravitational-potential'])], take_log=False)


def _dGmag(field,data):
       out = na.zeros(data['gravitational-potential'].shape, dtype='float64')
       out = na.sqrt( data['dGdx']*data['dGdx'] + data['dGdy']*data['dGdy'] + data['dGdz']*data['dGdz'])
       return out
add_field('dGmag',function=_dGmag, validators=[ValidateSpatial(1,['gravitational-potential'])], take_log=False)


