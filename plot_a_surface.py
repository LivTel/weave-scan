# RMB: circle fitting code removed as of 07/12/15

import sys
import argparse
import itertools
from collections import Counter

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate, optimize, signal
from scipy.interpolate import griddata, interp1d
import pyfits

from functions import print_stats, polyfit2d, polyval2d, read_data_file, correct_stage_positions

parser = argparse.ArgumentParser()
parser.add_argument('f', help="path to file", action="store", type=str)
parser.add_argument('--w', help="window to grid over, defaults to min/max of data (x1,x2,y1,y2)", action="store", type=str, default='5, 290, 5, 290')
parser.add_argument('--i', help="increment of interpolated grid as CSV tuple (x,y)", action="store", type=str, default='1, 1')
parser.add_argument('--de', help="maximum error in distance measurement", action="store", type=float, default=1.0)
parser.add_argument('--ca', help="don't apply stage calibration files", action="store_false")
parser.add_argument('--cax', help="calibration file for x if stage calibration is to be applied", action="store", type=str, default='calibration_files/45855915_cal.dat')
parser.add_argument('--cay', help="calibration file for y if stage calibration is to be applied", action="store", type=str, default='calibration_files/45855916_cal.dat')
parser.add_argument('--fl', help="don't apply flat calibration", action="store_false")
parser.add_argument('--flf', help="calibration file if flat calibration is to be applied", action="store", type=str, default='calibration_files/instresp.dat')
parser.add_argument('--z', help="z-depth for mean-levelled plots, or rather distance, as CSV tuple (high, low)", action="store", type=str, default='10, -10')
parser.add_argument('--p2D', help="plot 2D", action="store_true")
parser.add_argument('--p3D', help="plot 3D", action="store_true")
parser.add_argument('--rs', help="row stride (3d plot only)", action="store", type=float, default=5)
parser.add_argument('--cs', help="column stride (3d plot only)", action="store", type=float, default=5)
parser.add_argument('--fits', help="make fits file", action="store_true")
args = parser.parse_args()

args.i = ([float(n) for n in args.i.split(',')])
args.w = ([float(n) for n in args.w.split(',')])
args.z = ([float(n) for n in args.z.split(',')])

# READ INPUT DATA FILE
# --------------------
x, y, d, d_err = read_data_file(args.f, args.de, args.w)

# CORRECT FOR NONLINEAR MOTION OF STAGES
# --------------------------------------
if args.ca:
    x, y = correct_stage_positions(args.cax, args.cay, x, y)

## DEBUG TO SHOW DEVIATION FOR NO INTERPOLATION
'''plt.scatter(x, y, c=d, s=15, linewidth=0, marker='s', cmap="autumn")
plt.colorbar()
plt.show()'''

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

if args.fl:
    print "Subtracting instrument response..."
    s = np.loadtxt(args.flf)
    xx_interpolated_1d = xx_interpolated.flatten()							# cast to 1d
    yy_interpolated_1d = yy_interpolated.flatten()							
    zz_interpolated_1d = zz_interpolated.flatten()
    zz_interpolated_1d = zz_interpolated_1d - polyval2d(xx_interpolated_1d, yy_interpolated_1d, s)
    zz_interpolated = zz_interpolated_1d.reshape(zz_interpolated.shape)					# reshape
    print_stats(zz_interpolated)
    surfaces.append(zz_interpolated)
    extents.append(args.w)

# PLOTS
# -----
# 2D
if args.p2D:
    plt.figure()
    plt.subplot(221)
    plt.title("data in")
    plt.imshow(surfaces[0].transpose(), interpolation=None, origin='lower', extent=extents[0])
    plt.colorbar()
    plt.subplot(222)
    plt.title("mean tilt removed")
    plt.imshow(surfaces[1].transpose(), vmin=args.z[1], vmax=args.z[0], interpolation=None, origin='lower', extent=extents[1])
    plt.colorbar()
    plt.subplot(223)
    plt.title("instrument response removed")
    plt.imshow(surfaces[2].transpose(), vmin=args.z[1], vmax=args.z[0], interpolation=None, origin='lower', extent=extents[2])
    plt.colorbar()
    plt.tight_layout()
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
    pyfits.writeto(args.f + ".fits", zz_interpolated, hdr)



