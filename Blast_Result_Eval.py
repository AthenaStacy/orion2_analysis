#!/home1/atmyers/yt/bin/python2.6

# This script is called like "./make_slice.py XXX '<field_name>'" .
# where XXX is the frame number and <field_name> is what you want to plot
# for instance, 'python make_slice.py 005 'density' make a density slice
# from the output file data.0005.hdf5

from yt.mods import *
import matplotlib.colorbar as cb
import sys
import fields
import glob
import os
import pylab as plt

# this generates a list of plt files from whatever directory you call the script from
list_of_names = glob.glob('data.*.hdf5')
list_of_names.sort()
N = len(list_of_names)

vals = na.zeros(N)
times = na.zeros(N)

# loop through those files and make your plots
i = 0
for fn in list_of_names:
    print fn
    pf = load(fn)
    times[i] = pf.current_time
    data = pf.h.all_data()
    vals[i] = data.quantities["WeightedAverageQuantity"]("total-energy-density", "CellVolume", lazy_reader=True)
    i = i + 1

N = len(vals)
diff = na.zeros(N-1)

for n in na.arange(N-1):
    diff[n] = abs(vals[n+1]-vals[n])/vals[n]

if any(diff > 2.0e5):
    print "FAIL!"
else:
    print "PASS!"
    writer = open('test_pass', 'w')
    writer.close()
