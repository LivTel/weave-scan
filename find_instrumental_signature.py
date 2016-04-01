import sys
import argparse
import copy

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate, optimize, signal
from scipy.interpolate import griddata, interp1d
import pyfits

from functions import print_stats, polyfit2d, polyval2d, read_processed_scan_file, correct_stage_positions

parser = argparse.ArgumentParser()
parser.add_argument('--f', help="path to processed scan output", action="append", type=str)
parser.add_argument('--w', help="window to use data points (x1,x2,y1,y2).", action="store", type=str, default='10,290,10,290')
parser.add_argument('--i', help="increment of interpolated grid as CSV tuple (x,y)", action="store", type=str, default='5,5')
parser.add_argument('--de', help="maximum error in repeatability of distance", action="store", type=float, default=1)
parser.add_argument('--o', help="order of interpolate for low frequency structure", action="store", type=int, default=5)
parser.add_argument('--s', help="smoothing factor for interpolate", action="store", type=int, default=3)
parser.add_argument('--p', help="make plots?", action="store_true")
args = parser.parse_args()

args.w = ([float(n) for n in args.w.split(',')])
args.i = ([float(n) for n in args.i.split(',')])

# READ INPUT PROCESSED SCAN FILES
# -------------------------------

x = []
y = []
d = []
d_err = []
zz_interpolated_all = []

w = args.w
xx_interpolated, yy_interpolated = np.mgrid[w[0]:w[1]:args.i[0], w[2]:w[3]:args.i[1]]
for fi in args.f:
    this_x, this_y, this_d, this_d_err = read_processed_scan_file(fi, args.de)

    # INTERPOLATE DATA ON TO REGULAR GRID 
    # -----------------------------------
    x = np.asarray(this_x) # required by poly[fit||val]2d
    y = np.asarray(this_y) # required by poly[fit||val]2d
    d = np.asarray(this_d) # required by poly[fit||val]2d

    this_zz_interpolated = griddata((x, y), d, (xx_interpolated, yy_interpolated), method='linear')

    # SUBTRACT MEAN LEVEL OFF (GROSS X/Y TILT)
    # ----------------------------------------
    print "Subtracting mean plane (x/y tilt)..."
    this_xx_interpolated = xx_interpolated[1:-1, 1:-1]
    this_yy_interpolated = yy_interpolated[1:-1, 1:-1]
    this_zz_interpolated = this_zz_interpolated[1:-1, 1:-1]

    this_xx_interpolated_1d = this_xx_interpolated.flatten()								# cast to 1d
    this_yy_interpolated_1d = this_yy_interpolated.flatten()							
    this_zz_interpolated_1d = this_zz_interpolated.flatten()
    m = polyfit2d(this_xx_interpolated_1d, this_yy_interpolated_1d, this_zz_interpolated_1d, order=1)			# find coeffs
    this_zz_interpolated_1d = this_zz_interpolated_1d - polyval2d(this_xx_interpolated_1d, this_yy_interpolated_1d, m)	# remove mean level
    this_zz_interpolated = this_zz_interpolated_1d.reshape(this_zz_interpolated.shape)					# reshape

    zz_interpolated_all.append(this_zz_interpolated)

xx_interpolated = xx_interpolated[1:-1, 1:-1]
yy_interpolated = yy_interpolated[1:-1, 1:-1]
zz_interpolated = np.mean(zz_interpolated_all, axis=0)

surfaces = []
extents = []

print_stats(zz_interpolated)
surfaces.append(zz_interpolated)
extents.append(w)

# GENERATE INTERPOLATE FOR LOW FREQUENCY STRUCTURE
# ------------------------------------------------
# produces a single instresp.dat file.
print "Interpolating 2D surface..."
it = interpolate.bisplrep(xx_interpolated, yy_interpolated, zz_interpolated, kx=args.o, ky=args.o, s=args.s)
fitted_s = copy.deepcopy(zz_interpolated)
for idx_j, row in enumerate(xx_interpolated): 
    for idx_i, col in enumerate(row):
        this_x = xx_interpolated[idx_j, idx_i]
        this_y = yy_interpolated[idx_j, idx_i]
        fitted_s[idx_j, idx_i] = interpolate.bisplev(this_x, this_y, it)

surfaces.append(fitted_s)			# add fitted surface to list
extents.append(w)
surfaces.append(fitted_s-zz_interpolated)	# add residual surface to list
extents.append(w)

print_stats(fitted_s-zz_interpolated)
np.save("instresp", it)
print "Wrote file instresp.dat."
print

# PLOTS
# -----
if args.p:
  # 2D
  plt.figure()
  plt.subplot(221)
  plt.title("data in")
  plt.imshow(surfaces[0].transpose(), interpolation=None, origin='lower', extent=extents[0])
  plt.colorbar()
  plt.subplot(222)
  plt.title("fitted surface to be removed")
  plt.imshow(surfaces[1].transpose(), interpolation=None, origin='lower', extent=extents[1])
  plt.colorbar()
  plt.subplot(223)
  plt.title("residual after surface removal")
  plt.imshow(surfaces[2].transpose(), interpolation=None, origin='lower', extent=extents[2])
  plt.colorbar()
  plt.tight_layout()
  plt.show()

  fig = plt.figure()
  ax = fig.gca(projection='3d')
  ax.set_zlim(-20,20)
  ax.plot_surface(xx_interpolated, yy_interpolated, surfaces[1].transpose(), cmap=cm.Reds)
  plt.show()

