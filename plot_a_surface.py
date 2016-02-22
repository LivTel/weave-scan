# RMB: circle fitting code removed as of 07/12/15

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

from functions import print_stats, polyfit2d, polyval2d, read_data_file, correct_stage_positions

F_CALX = "calibration_files/45855915_cal.dat"
F_CALY = "calibration_files/45855915_cal.dat"
F_FLAT = "calibration_files/instresp.dat"
F_FLATX = "calibration_files/instresp_x.dat"
F_FLATY = "calibration_files/instresp_y.dat"

parser = argparse.ArgumentParser()
parser.add_argument('f', help="path to file", action="store", type=str)
parser.add_argument('--w', help="window to grid over, defaults to min/max of data (x1,x2,y1,y2)", action="store", type=str, default='5, 290, 5, 290')
parser.add_argument('--i', help="increment of interpolated grid as CSV tuple (x,y)", action="store", type=str, default='5, 5')
parser.add_argument('--de', help="maximum error in distance measurement", action="store", type=float, default=1.0)
parser.add_argument('--nc', help="don't apply stage calibration files", action="store_false")
parser.add_argument('--nf', help="don't apply flat calibration", action="store_false")
parser.add_argument('--m', help="median filter size as CSV tuple (x,y)", action="store", type=str, default='3,3')
parser.add_argument('--z', help="z-depth for mean-levelled plots, or rather distance, as CSV tuple (high, low)", action="store", type=str, default='15, -15')
parser.add_argument('--p2D', help="plot 2D", action="store_true")
parser.add_argument('--p3D', help="plot 3D", action="store_true")
parser.add_argument('--rs', help="row stride (3d plot only)", action="store", type=float, default=5)
parser.add_argument('--cs', help="column stride (3d plot only)", action="store", type=float, default=5)
parser.add_argument('--fits', help="make fits file", action="store_true")
args = parser.parse_args()

args.i = ([float(n) for n in args.i.split(',')])
args.w = ([float(n) for n in args.w.split(',')])
args.m = ([int(n) for n in args.m.split(',')])
args.z = ([float(n) for n in args.z.split(',')])

# READ INPUT DATA FILE
# --------------------
x, y, d, d_err = read_data_file(args.f, args.de, args.w)

# CORRECT FOR NONLINEAR MOTION OF STAGES
# --------------------------------------
if args.nc:
    x, y = correct_stage_positions(F_CALX, F_CALY, x, y)

# SHOW DEVIATION FOR NO INTERPOLATION
# -----------------------------------
#plt.scatter(x, y, c=d, s=15, linewidth=0, marker='s', cmap="autumn")
#plt.colorbar()
#plt.show()

# INTERPOLATE DATA ON TO REGULAR GRID 
# -----------------------------------
  
x = np.array(x) # required by poly[fit||val]2d
y = np.array(y) # required by poly[fit||val]2d
d = np.array(d) # required by poly[fit||val]2d

xx_interpolated, yy_interpolated = np.mgrid[args.w[0]:args.w[1]:args.i[0], args.w[2]:args.w[3]:args.i[1]]
zz_interpolated = griddata((x, y), d, (xx_interpolated, yy_interpolated), method='linear')
xx_interpolated = xx_interpolated[1:-1,1:-1]		# trim off border pixels where interpolation will have yielded NaN
yy_interpolated = yy_interpolated[1:-1,1:-1]		# trim off border pixels where interpolation will have yielded NaN
zz_interpolated = zz_interpolated[1:-1,1:-1]		# trim off border pixels where interpolation will have yielded NaN
args.w = [args.w[0]+1, args.w[1]-1, args.w[2]+1, args.w[3]-1]

surfaces = []
extents = []

# SET WINDOW IF NOT SET AT COMMAND LINE
# -------------------------------------
if not any(args.w):     # we can now set the default window to min/max of x,y 		                                             
  args.w = [round(min(x), 2), round(max(x), 2), round(min(y), 2), round(max(y), 2)]

print_stats(zz_interpolated)
surfaces.append(zz_interpolated)
extents.append(args.w)

# MEDIAN FILTER TO SMOOTH NOISE
# -----------------------------
print "Applying median filter..."
zz_interpolated = signal.medfilt(np.nan_to_num(zz_interpolated), kernel_size=[args.m[0],args.m[1]])
xx_interpolated = xx_interpolated[(args.m[1]-1)/2:-(args.m[1]-1)/2,(args.m[0]-1)/2:-(args.m[0]-1)/2]		# trim off border pixels where filter will have yielded 0
yy_interpolated = yy_interpolated[(args.m[1]-1)/2:-(args.m[1]-1)/2,(args.m[0]-1)/2:-(args.m[0]-1)/2]		# trim off border pixels where filter will have yielded 0
zz_interpolated = zz_interpolated[(args.m[1]-1)/2:-(args.m[1]-1)/2,(args.m[0]-1)/2:-(args.m[0]-1)/2]		# trim off border pixels where filter will have yielded 0
args.w = [args.w[0]+(args.m[0]-1)/2, args.w[1]-(args.m[0]-1)/2, args.w[2]+(args.m[1]-1)/2, args.w[3]-(args.m[1]-1)/2]
print_stats(zz_interpolated)
surfaces.append(zz_interpolated)
extents.append(args.w)

# SUBTRACT MEAN LEVEL OFF (GROSS X/Y TILT)
# ----------------------------------------
print "Subtracting mean plane (x/y tilt)..."
xx_interpolated_1d = xx_interpolated.flatten()							# cast to 1d
yy_interpolated_1d = yy_interpolated.flatten()							
zz_interpolated_1d = zz_interpolated.flatten()
m = polyfit2d(xx_interpolated_1d, yy_interpolated_1d, zz_interpolated_1d, order=1)		# find coeffs
zz_interpolated_1d = zz_interpolated_1d - polyval2d(xx_interpolated_1d, yy_interpolated_1d, m)	# remove mean level
zz_interpolated = zz_interpolated_1d.reshape(zz_interpolated.shape)				# reshape
print_stats(zz_interpolated)
surfaces.append(zz_interpolated)
extents.append(args.w)

# SUBTRACT INSTRUMENT RESPONSE OFF
# --------------------------------

if args.nf:
    print "Subtracting instrument response (SURFACE)..."
    s = np.loadtxt(F_FLAT)
    xx_interpolated_1d = xx_interpolated.flatten()							# cast to 1d
    yy_interpolated_1d = yy_interpolated.flatten()							
    zz_interpolated_1d = zz_interpolated.flatten()
    zz_interpolated_1d_surf = copy.deepcopy(zz_interpolated_1d)
    zz_interpolated_1d_surf = zz_interpolated_1d_surf - polyval2d(xx_interpolated_1d, yy_interpolated_1d, s)
    zz_interpolated_surf = zz_interpolated_1d_surf.reshape(zz_interpolated.shape)			# reshape
    print_stats(zz_interpolated_surf)
    surfaces.append(zz_interpolated_surf)
    extents.append(args.w)

    # this is fitting per row/column
    print "Subtracting instrument response in x (1D)..."
    instresp_x = np.loadtxt("calibration_files/instresp_x.dat")
    xx_interpolated = xx_interpolated.transpose()
    yy_interpolated = yy_interpolated.transpose()
    zz_interpolated = copy.deepcopy(zz_interpolated.transpose())
    for idx, (this_xx, this_yy, this_zz) in enumerate(zip(xx_interpolated, yy_interpolated, zz_interpolated)):
        for x in instresp_x:
            y = x[0]
            m = x[1:]
            if y == this_yy[0]: 
                zz_interpolated[idx] = zz_interpolated[idx] - np.polyval(m, xx_interpolated[idx])
    xx_interpolated = xx_interpolated.transpose()	# switch back
    yy_interpolated = yy_interpolated.transpose()
    zz_interpolated = zz_interpolated.transpose()

    print_stats(zz_interpolated)
    surfaces.append(zz_interpolated)
    extents.append(args.w)

    print "Fitting low frequency structure in y (1D)..."
    instresp_y = np.loadtxt("calibration_files/instresp_y.dat")
    for idx, (this_xx, this_yy, this_zz) in enumerate(zip(xx_interpolated, yy_interpolated, zz_interpolated)):
        for y in instresp_y:
            x = y[0]
            m = y[1:]
            if x == this_xx[0]: 
                zz_interpolated[idx] = zz_interpolated[idx] - np.polyval(m, yy_interpolated[idx])

    print_stats(zz_interpolated)
    surfaces.append(zz_interpolated)
    extents.append(args.w)

# PLOTS
# -----
# 2D
if args.p2D:
    plt.figure()
    plt.subplot(321)
    plt.title("data in")
    plt.imshow(surfaces[0].transpose(), interpolation=None, origin='lower', extent=extents[0])
    plt.colorbar()
    plt.subplot(322)
    plt.title("median filtered")
    plt.imshow(surfaces[1].transpose(), interpolation=None, origin='lower', extent=extents[1])
    plt.colorbar()
    plt.subplot(323)
    plt.title("mean tilt removed")
    plt.imshow(surfaces[2].transpose(), vmin=args.z[1], vmax=args.z[0], interpolation=None, origin='lower', extent=extents[2])
    plt.colorbar()
    plt.subplot(324)
    plt.title("surface response removed (2D mode)")
    plt.imshow(surfaces[3].transpose(), vmin=args.z[1], vmax=args.z[0], interpolation=None, origin='lower', extent=extents[3])
    plt.colorbar()
    plt.subplot(325)
    plt.title("residual after x fitting (1D mode)")
    plt.imshow(surfaces[4].transpose(), vmin=args.z[1], vmax=args.z[0], interpolation=None, origin='lower', extent=extents[4])
    plt.colorbar()
    plt.subplot(326)
    plt.title("residual after y fitting (1D mode)")
    plt.imshow(surfaces[5].transpose(), vmin=args.z[1], vmax=args.z[0], interpolation=None, origin='lower', extent=extents[5])
    plt.colorbar()
    plt.tight_layout()
    plt.show()
 
    plt.subplot(211)
    plt.title("mean deviation in x")
    #plt.plot(np.mean(surfaces[3].transpose(), axis=0))	# 2D mode
    plt.plot(np.mean(surfaces[5].transpose(), axis=0))	# 1D mode
    plt.subplot(212)
    plt.title("mean deviation in y")
    #plt.plot(np.mean(surfaces[3].transpose(), axis=1))	# 2D mode
    plt.plot(np.mean(surfaces[5].transpose(), axis=1))	# 1D mode
    plt.show()

# 3D  
if args.p3D:
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(xx_interpolated, yy_interpolated, zz_interpolated, rstride=args.rs, cstride=args.cs, cmap=cm.coolwarm, alpha=0.3)
    plt.show()

# fits
if args.fits:
    hdr = pyfits.core.Header()
    pyfits.writeto(os.path.basename(args.f) + ".fits", zz_interpolated, hdr)
    pyfits.writeto(os.path.basename(args.f) + ".surf.fits", zz_interpolated_surf, hdr)

