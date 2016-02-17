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

from functions import print_stats, polyfit2d, polyval2d, read_data_file, correct_stage_positions

F_CALX = "calibration_files/45855915_cal.dat"
F_CALY = "calibration_files/45855915_cal.dat"

parser = argparse.ArgumentParser()
parser.add_argument('f', help="path to flat", action="store", type=str)
parser.add_argument('--w', help="window, defaults to min/max of data (x1,x2,y1,y2)", action="store", type=str, default='10,290,10,290')
parser.add_argument('--i', help="increment of interpolated grid as CSV tuple (x,y)", action="store", type=str, default='5,5')
parser.add_argument('--de', help="maximum error in distance measurement", action="store", type=float, default=1.0)
parser.add_argument('--nc', help="don't apply stage calibration files", action="store_false")
parser.add_argument('--m', help="median filter size as CSV tuple (x,y)", action="store", type=str, default='3,3')
parser.add_argument('--o', help="order of fit for low frequency structure", action="store", type=int, default='3')
parser.add_argument('--p', help="plot", action="store_true")
args = parser.parse_args()

args.i = ([float(n) for n in args.i.split(',')])
args.m = ([int(n) for n in args.m.split(',')])
args.w = ([float(n) for n in args.w.split(',')])

# READ INPUT DATA FILE
# --------------------
x, y, d, d_err = read_data_file(args.f, args.de, args.w)

# CORRECT FOR NONLINEAR MOTION OF STAGES
# --------------------------------------
if args.nc:
    x, y = correct_stage_positions(F_CALX, F_CALY, x, y)

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
# ------------------------------------------------
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

# FIT LOW FREQUENCY STRUCTURE
# ----------------------------
# this is fitting for surface
# produces a single instresp.dat file.
print "Fitting 2D surface..."
xx_interpolated_1d = xx_interpolated.flatten()							# cast to 1d
yy_interpolated_1d = yy_interpolated.flatten()							
zz_interpolated_1d = zz_interpolated.flatten()
m = polyfit2d(xx_interpolated_1d, yy_interpolated_1d, zz_interpolated_1d, order=args.o)		# find coeffs
s = polyval2d(xx_interpolated_1d, yy_interpolated_1d, m)
zz_interpolated_1d_surf = copy.deepcopy(zz_interpolated_1d)
zz_interpolated_1d_surf = zz_interpolated_1d_surf - s						# remove mean level
zz_interpolated_1d_surf = zz_interpolated_1d_surf.reshape(zz_interpolated.shape)		# reshape

surfaces.append(s.reshape(zz_interpolated.shape))						# add fitted surface to list
extents.append(args.w)
surfaces.append(zz_interpolated_1d_surf)							# add residual surface to list
extents.append(args.w)

print_stats(zz_interpolated_1d_surf)
np.savetxt("instresp.dat", m)
print "Wrote file instresp.dat."
print

# this is fitting per row/column
# produces an instresp_[xy?].dat file for each axis.
## x
print "Fitting low frequency structure in x... (1D)"
xx_interpolated = xx_interpolated.transpose()
yy_interpolated = yy_interpolated.transpose()
zz_interpolated = copy.deepcopy(zz_interpolated.transpose())
instresp_x = []
for idx, (this_xx, this_yy, this_zz) in enumerate(zip(xx_interpolated, yy_interpolated, zz_interpolated)):
    m = np.polyfit(this_xx, this_zz, deg=args.o)
    instresp_x.append(np.insert(m, 0, this_yy[0]))
    zz_interpolated[idx] = zz_interpolated[idx] - np.polyval(m, xx_interpolated[idx])
    if this_yy[0] == 275:
      plt.plot(xx_interpolated[idx], np.polyval(m, xx_interpolated[idx]))
plt.show()
xx_interpolated = xx_interpolated.transpose()	# switch back!
yy_interpolated = yy_interpolated.transpose()
zz_interpolated = zz_interpolated.transpose()

surfaces.append(zz_interpolated)
extents.append(args.w)

print_stats(zz_interpolated)
np.savetxt("instresp_x.dat", instresp_x)
print "Wrote file instresp_x.dat."
print

## y
print "Fitting low frequency structure in y... (1D)"
instresp_y = []
for idx, (this_xx, this_yy, this_zz) in enumerate(zip(xx_interpolated, yy_interpolated, zz_interpolated)):
    m = np.polyfit(this_yy, this_zz, deg=args.o)
    instresp_y.append(np.insert(m, 0, this_xx[0]))
    zz_interpolated[idx] = zz_interpolated[idx] - np.polyval(m, yy_interpolated[idx])
    plt.plot(yy_interpolated[idx], np.polyval(m, yy_interpolated[idx]))
plt.show()
surfaces.append(zz_interpolated)
extents.append(args.w)

print_stats(zz_interpolated)
np.savetxt("instresp_y.dat", instresp_y)
print "Wrote file instresp_y.dat."
print

# PLOTS
# -----
# 2D
if args.p:
  plt.figure()
  plt.subplot(431)
  plt.title("data in")
  plt.imshow(surfaces[0].transpose(), interpolation=None, origin='lower', extent=extents[0])
  plt.colorbar()
  plt.subplot(432)
  plt.title("median filtered")
  plt.imshow(surfaces[1].transpose(), interpolation=None, origin='lower', extent=extents[1])
  plt.colorbar()
  plt.subplot(433)
  plt.title("tilt removed")
  plt.imshow(surfaces[2].transpose(), interpolation=None, origin='lower', extent=extents[2])
  plt.colorbar()
  plt.subplot(434)
  plt.title("fitted surface to be removed")
  plt.imshow(surfaces[3].transpose(), interpolation=None, origin='lower', extent=extents[3])
  plt.colorbar()
  plt.subplot(435)
  plt.title("residual after surface removal")
  plt.imshow(surfaces[4].transpose(), interpolation=None, origin='lower', extent=extents[4])
  plt.colorbar()
  plt.subplot(436)
  plt.title("residual after x fitting")
  plt.imshow(surfaces[5].transpose(), interpolation=None, origin='lower', extent=extents[5])
  plt.colorbar()
  plt.subplot(437)
  plt.title("residual after y fitting")
  plt.imshow(surfaces[6].transpose(), interpolation=None, origin='lower', extent=extents[6])
  plt.colorbar()
  plt.tight_layout()
  plt.show()

