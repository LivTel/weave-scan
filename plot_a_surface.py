import sys
import argparse
import itertools
import copy
import os

from collections import Counter
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate, optimize, signal
from scipy.interpolate import griddata, interp1d
import pyfits

from functions import print_stats, polyfit2d, polyval2d, read_data_file, correct_stage_positions, correct_for_temperature

F_CALX = "calibration_files/45855915_cal.dat"
F_CALY = "calibration_files/45855915_cal.dat"
F_INSTRESP = "calibration_files/instresp.npy"
F_CORT = "calibration_files/t_cor.npy"

parser = argparse.ArgumentParser()
parser.add_argument('d', help="path to data directory", action="store", type=str)
parser.add_argument('--w', help="window to use data points (x1,x2,y1,y2).", action="store", type=str, default='10,290,10,290')
parser.add_argument('--de', help="maximum error in distance measurement", action="store", type=float, default=1.0)
parser.add_argument('--nc', help="don't apply stage calibration files", action="store_false")
parser.add_argument('--nf', help="don't apply flat calibration", action="store_false")
parser.add_argument('--nt', help="don't apply temperature correction file", action="store_false")
parser.add_argument('--p', help="make plots", action="store_true")
args = parser.parse_args()

args.w = ([float(n) for n in args.w.split(',')])

# READ INPUT DATA FILE
# --------------------
print "using directory " + args.d
TIME = []
x = []
y = []
d = []
t = []
h = []
for f in os.listdir(args.d):
    print "processing file " + f
    if f.endswith(".log"):
        print "ignoring file"
        continue
    fi = args.d.rstrip('/') + '/' + f

    # READ INPUT DATA FILE
    # --------------------
    this_run_time, this_run_x, this_run_y, this_run_d, this_run_d_err, this_run_t, this_run_h = read_data_file(fi, args.de, args.w)

    # CORRECT FOR NONLINEAR MOTION OF STAGES
    # --------------------------------------
    if args.nc:
        this_run_x, this_run_y = correct_stage_positions(F_CALX, F_CALY, this_run_x, this_run_y)

    # CORRECT FOR TEMPERATURE EXCURSIONS
    # ----------------------------------
    if args.nt:
        this_run_d = correct_for_temperature(F_CORT, this_run_d, this_run_t)

    TIME.append(np.asarray(this_run_time))
    x.append(np.mean(this_run_x))
    y.append(np.mean(this_run_y))
    d.append(np.mean(this_run_d))
    t.append(np.mean(this_run_t))
    h.append(np.mean(this_run_h))
    print

TIME = np.asarray(TIME)
x = np.asarray(x) # required by poly[fit||val]2d
y = np.asarray(y) # required by poly[fit||val]2d
d = np.asarray(d) # required by poly[fit||val]2d
t = np.asarray(t)
h = np.asarray(h)

# SUBTRACT INSTRUMENT RESPONSE OFF
# --------------------------------
if args.nf:
    print "Subtracting instrument response..."
    it = np.load(F_INSTRESP)
    for idx in range(len(x)): 
        this_x = x[idx]
        this_y = y[idx]
        d[idx] = d[idx] - interpolate.bisplev(this_x, this_y, it)
print_stats(d)

# SUBTRACT MEAN LEVEL OFF (GROSS X/Y TILT)
# ----------------------------------------
print "Subtracting mean plane (x/y tilt)..."
m = polyfit2d(x, y, d, order=1)							# find coeffs
d = d - polyval2d(x, y, m)							# remove mean level
print_stats(d)

if args.p:
    # SHOW DEVIATION
    # --------------
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_trisurf(x, y, d)
    plt.show()



