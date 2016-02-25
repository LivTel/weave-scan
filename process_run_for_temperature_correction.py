import sys
import argparse
import copy
import os

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import pyfits

from functions import read_data_file, correct_stage_positions

F_CALX = "calibration_files/45855915_cal.dat"
F_CALY = "calibration_files/45855915_cal.dat"

parser = argparse.ArgumentParser()
parser.add_argument('--f', help="path to file", action="store", type=str)
parser.add_argument('--de', help="maximum error in distance measurement", action="store", type=float, default=1.0)
parser.add_argument('--o', help="order of fit", action="store", type=int, default=2)
parser.add_argument('--tr', help="CSV range of temperatures over which to fit", action="store", type=str, default='23,25')
parser.add_argument('--tcal', help="temperature at which to correct at (should be mid-range ish)", action="store", type=float, default=24)
parser.add_argument('--p', help="make plots?", action="store_true")
args = parser.parse_args()

args.tr = ([float(T) for T in args.tr.split(',')])

print "processing file " + args.f

# READ INPUT DATA FILE
# --------------------
this_time, x, y, d, d_err, t, h = read_data_file(args.f, args.de, None)

c = np.polyfit([ti for ti in t if ti>args.tr[0] and ti<args.tr[1]], [di for di, ti in zip(d, t) if ti>args.tr[0] and ti<args.tr[1]], args.o)

# we correct for the position by finding how much by calibrating the measured distance to 
# what it would be at a reference temperature of [t_cal]

d_cal = np.polyval(c, args.tcal) 	# the calibration measurement at our reference temperature
offsets = np.polyval(c, t) - d_cal	# the offsets that need to be applied for these measured temperatures to make them equivalent to being taken at the reference temperature

c2 = np.polyfit([ti for ti in t if ti>args.tr[0] and ti<args.tr[1]], [oi for oi, ti in zip(offsets, t) if ti>args.tr[0] and ti<args.tr[1]], args.o)

np.save("t_cor", c2)

if args.p:
    plt.subplot(211)
    plt.plot(t, d, 'kx')
    plt.plot(t, np.polyval(c, t), 'r-')
    plt.ylabel("sensor reading (mm)")
    plt.subplot(212)
    plt.plot(t, offsets, 'kx')
    plt.plot(t, np.polyval(c2, t), 'r-')
    plt.ylabel("offset (mm)")
    plt.xlabel("temperature (deg)")
    plt.show()


