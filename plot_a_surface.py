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
from scipy.signal import medfilt
import pyfits

from functions import print_stats, polyfit2d, polyval2d, read_data_file, correct_stage_positions, correct_for_temperature

F_CALX = "calibration_files/45855915_cal.dat"
F_CALY = "calibration_files/45855915_cal.dat"
F_INSTRESP = "calibration_files/instresp.npy"
F_CORT = "calibration_files/t_cor.npy"

parser = argparse.ArgumentParser()
parser.add_argument('--f', help="path to file or directory", action="store", type=str)
parser.add_argument('--w', help="window to use data points (x1,x2,y1,y2).", action="store", type=str, default='20,200,248,280')
parser.add_argument('--i', help="increment of interpolated grid as CSV tuple (x,y)", action="store", type=str, default='5,5')
parser.add_argument('--de', help="maximum error in repeatability of distance", action="store", type=float, default=0.5)
parser.add_argument('--int', help="do interpolation of result and smooth", action="store_true")
parser.add_argument('--m', help="kernel size of median filter (if --in set)", action="store", type=float, default=3)
parser.add_argument('--o', help="order of interpolate for low frequency structure (if --in set)", action="store", type=int, default=5)
parser.add_argument('--s', help="smoothing factor for interpolate (if --in set)", action="store", type=int, default=100)
parser.add_argument('--nc', help="don't apply stage calibration files", action="store_false")
parser.add_argument('--nf', help="don't apply flat calibration", action="store_false")
parser.add_argument('--nt', help="don't apply temperature correction file", action="store_false")
parser.add_argument('--p', help="make plots", action="store_true")
args = parser.parse_args()

args.w = ([float(n) for n in args.w.split(',')])
args.i = ([float(n) for n in args.i.split(',')])

TIME = []
x = []
y = []
d = []
t = []
h = []
# READ INPUT DATA FILE
# --------------------

if os.path.isdir(args.f):
    print "using directory " + args.f
    for f in os.listdir(args.f):
        print "processing file " + f
        if f.endswith(".log"):
            print "ignoring file"
            continue
        fi = args.f.rstrip('/') + '/' + f

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
        t.append(np.asarray(this_run_t))
        h.append(np.asarray(this_run_h))
        print
else:
    print "using file " + args.f

    # READ INPUT DATA FILE
    # --------------------
    this_run_time, this_run_x, this_run_y, this_run_d, this_run_d_err, this_run_t, this_run_h = read_data_file(args.f, args.de, args.w)
    # CORRECT FOR NONLINEAR MOTION OF STAGES
    # --------------------------------------
    if args.nc:
        this_run_x, this_run_y = correct_stage_positions(F_CALX, F_CALY, this_run_x, this_run_y)

    # CORRECT FOR TEMPERATURE EXCURSIONS
    # ----------------------------------
    if args.nt:
        this_run_d = correct_for_temperature(F_CORT, this_run_d, this_run_t)

    TIME.append(np.asarray(this_run_time))
    x.append(np.asarray(this_run_x))
    y.append(np.asarray(this_run_y))
    d.append(np.asarray(this_run_d))
    t.append(np.asarray(this_run_t))
    h.append(np.asarray(this_run_h))
    print

TIME = np.asarray(TIME)
x = np.asarray(x) # required by poly[fit||val]2d
y = np.asarray(y) # required by poly[fit||val]2d
d = np.asarray(d) # required by poly[fit||val]2d
t = np.asarray(t)
h = np.asarray(h)

# REMOVE NaNs
# -----------
idx_nonan = [np.where(np.logical_and(np.isnan(x) == False, np.isnan(y) == False, np.isnan(d) == False))][0]
x = x[idx_nonan]
y = y[idx_nonan]
d = d[idx_nonan]

# SUBTRACT MEAN LEVEL OFF (GROSS X/Y TILT)
# ----------------------------------------
print "Subtracting mean plane (x/y tilt)..."
m = polyfit2d(x, y, d, order=1)		# find coeffs
d = d - polyval2d(x, y, m)		# remove mean level
print_stats(d)

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

#FIXME
'''import scipy.signal	
fig = plt.figure()
plt.scatter(x, y, s=5, c=h, marker='x', vmax=np.percentile(h, 95), vmin=np.percentile(h, 5))  
plt.colorbar()
plt.show()
exit(0)'''


if args.p:
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_zlim(-20,20)
    ax.plot_trisurf(x, y, d)
    plt.show()
    
    # INTERPOLATE SURFACE ONTO REGULAR GRID AND SMOOTH
    # ------------------------------------------------
    if args.int:
        print "Smoothing surface..."
        xx_interpolated, yy_interpolated = np.mgrid[args.w[0]:args.w[1]:args.i[0], args.w[2]:args.w[3]:args.i[1]]

        zz_interpolated = griddata((x, y), d, (xx_interpolated, yy_interpolated), method='linear')

        xx_interpolated = xx_interpolated[1:-1, 1:-1]
        yy_interpolated = yy_interpolated[1:-1, 1:-1]
        zz_interpolated = zz_interpolated[1:-1, 1:-1]

        zz_interpolated = medfilt(zz_interpolated, args.m)

        it = interpolate.bisplrep(xx_interpolated, yy_interpolated, zz_interpolated, kx=args.o, ky=args.o, s=args.s)
        fitted_s = copy.deepcopy(zz_interpolated)
        for idx_j, row in enumerate(xx_interpolated): 
            for idx_i, col in enumerate(row):
                this_x = xx_interpolated[idx_j, idx_i]
                this_y = yy_interpolated[idx_j, idx_i]
                fitted_s[idx_j, idx_i] = interpolate.bisplev(this_x, this_y, it)

        print_stats(zz_interpolated)

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_zlim(-20,20)
        ax.plot_surface(xx_interpolated, yy_interpolated, zz_interpolated, cmap=cm.Reds)
        plt.show()

#FIXME
import pyfits
pyfits.writeto("test.fits", zz_interpolated)






