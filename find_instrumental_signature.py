import sys
import itertools
import argparse

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate, optimize, signal
from scipy.interpolate import griddata, interp1d
import pyfits

def print_stats(array):
    print
    print "Percentile (0.5) of data is: " + '\t\t\t' + str(round(np.nanpercentile(array, 0.5), 2)) + "um"
    print "Percentile (99.5) of data is: " + '\t\t\t' + str(round(np.nanpercentile(array, 99.5), 2)) + "um"
    print "Peak-to-peak amplitude of structure is: " + '\t' + str(round(np.nanpercentile(array, 99.5)-np.nanpercentile(array, 0.5), 2)) + "um"
    print "Percentile (99.5) of absolute error is: " + '\t' + str(round(np.nanpercentile(abs(array), 99.5), 2)) + "um"
    print 

def polyfit2d(x, y, z, order=1, linear=False):
    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
        if linear & (i != 0.) & (j != 0.):
            G[:, k] = 0
    m, _, _, _ = np.linalg.lstsq(G, z)
    return m

def polyval2d(x, y, m):
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z

parser = argparse.ArgumentParser()
parser.add_argument('f', help="path to flat", action="store", type=str)
parser.add_argument('--w', help="window, defaults to min/max of data (x1,x2,y1,y2)", action="store", type=str, default='10,290,10,290')
parser.add_argument('--i', help="increment of interpolated grid as CSV tuple (x,y)", action="store", type=str, default='1,1')
parser.add_argument('--de', help="maximum error in distance measurement", action="store", type=float, default=1.0)
parser.add_argument('--ca', help="apply stage calibration files", action="store_true")
parser.add_argument('--cax', help="calibration file for x if stage calibration is to be applied", action="store", type=str, default='calibration_files/45855915_cal.dat')
parser.add_argument('--cay', help="calibration file for y if stage calibration is to be applied", action="store", type=str, default='calibration_files/45855916_cal.dat')
parser.add_argument('--m', help="median filter size as CSV tuple (x,y)", action="store", type=str, default='15,15')
parser.add_argument('--o', help="order of fit for low frequency structure", action="store", type=int, default='3')
parser.add_argument('--p2D', help="plot 2D", action="store_true")
parser.add_argument('--p3D', help="plot 3D", action="store_true")
parser.add_argument('--rs', help="row stride (3d plot only)", action="store", type=float, default=5)
parser.add_argument('--cs', help="column stride (3d plot only)", action="store", type=float, default=5)
parser.add_argument('--fits', help="make fits file", action="store_true")
args = parser.parse_args()

args.i = ([float(n) for n in args.i.split(',')])
args.m = ([int(n) for n in args.m.split(',')])
args.w = ([float(n) for n in args.w.split(',')])

# READ INPUT DATA FILE
# --------------------.
print "Reading input data..."
x       = []
y       = []
d       = []
d_err 	= []
with open(args.f) as data:
  for line in data:
    if line.startswith('#'):
      continue
    res = line.split()
    this_x = float(res[0])
    this_y = float(res[1])
    this_d = float(res[2])
    if any(args.w):				# if we don't have a default for the window
      if any([this_x < args.w[0],
              this_x > args.w[1],
              this_y < args.w[2],
              this_y > args.w[3]]):
        continue

    try:
      this_d_err = float(res[3])    
      if this_d_err>args.de or this_d==0:       # additional quality checks
        continue
      d_err.append(this_d_err)
    except IndexError:
      x.append(this_x)
      y.append(this_y)
      d.append(this_d)
    x.append(this_x)
    y.append(this_y)
    d.append(this_d)

# CORRECT FOR NONLINEAR MOTION OF STAGES
# --------------------------------------
if args.ca:
  print "Applying calibration to stages..."
  stage_x_rd = []
  stage_x_act = []
  with open(args.cax) as cal:
    for line in cal:
      try:
        stage_x_rd.append(float(line.split()[0]))
        stage_x_act.append(float(line.split()[1]))
      except ValueError:
        pass
  fx = interp1d(stage_x_rd, stage_x_act, kind='cubic')

  stage_y_rd = []
  stage_y_act = []
  with open(args.cay) as cal:
    for line in cal:
      try:
        stage_y_rd.append(float(line.split()[0]))
        stage_y_act.append(float(line.split()[1]))
      except ValueError:
        pass
  fy = interp1d(stage_y_rd, stage_y_act, kind='cubic')

  x = fx(x)
  y = fy(y)

x = np.array(x)	# required by poly[fit||val]2d
y = np.array(y) # required by poly[fit||val]2d
d = np.array(d) # required by poly[fit||val]2d

# INTERPOLATE DATA ON TO REGULAR GRID 
# -----------------------------------
xx_interpolated, yy_interpolated = np.mgrid[args.w[0]:args.w[1]:args.i[0], args.w[2]:args.w[3]:args.i[1]]
zz_interpolated = griddata((x, y), d, (xx_interpolated, yy_interpolated), method='linear')
xx_interpolated = xx_interpolated[1:-1,1:-1]		# trim off border pixels where interpolation will have yielded NaN
yy_interpolated = yy_interpolated[1:-1,1:-1]		# trim off border pixels where interpolation will have yielded NaN
zz_interpolated = zz_interpolated[1:-1,1:-1]		# trim off border pixels where interpolation will have yielded NaN
args.w = [args.w[0]+1, args.w[1]-1, args.w[2]+1, args.w[3]-1]

surfaces = []
extents = []

# SET WINDOW AND ZDEPTH IF NOT SET AT COMMAND LINE
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
# this is fitting of mean in each axis independenly
'''print "Fitting low frequency structure in x..."
xx_interpolated_mean = np.mean(xx_interpolated, axis=1)
zz_interpolated_mean = np.mean(zz_interpolated, axis=1)
m = np.polyfit(xx_interpolated_mean, zz_interpolated_mean, deg=args.o)

xx_interpolated = xx_interpolated.T
zz_interpolated = zz_interpolated.T
for idx, xx in enumerate(xx_interpolated):
  zz_interpolated[idx] = zz_interpolated[idx] - np.polyval(m, xx_interpolated[idx])
xx_interpolated = xx_interpolated.T
zz_interpolated = zz_interpolated.T
print_stats(zz_interpolated)

print "Fitting low frequency structure in y..."
yy_interpolated_mean = np.mean(yy_interpolated, axis=0)
zz_interpolated_mean = np.mean(zz_interpolated, axis=0)
m = np.polyfit(yy_interpolated_mean, zz_interpolated_mean, deg=args.o)

for idx, yy in enumerate(yy_interpolated):
  zz_interpolated[idx] = zz_interpolated[idx] - np.polyval(m, yy_interpolated[idx])
print_stats(zz_interpolated)'''

# this is fitting for surface
xx_interpolated_1d = xx_interpolated.flatten()							# cast to 1d
yy_interpolated_1d = yy_interpolated.flatten()							
zz_interpolated_1d = zz_interpolated.flatten()
m = polyfit2d(xx_interpolated_1d, yy_interpolated_1d, zz_interpolated_1d, order=args.o)		# find coeffs
s = polyval2d(xx_interpolated_1d, yy_interpolated_1d, m)
zz_interpolated_1d = zz_interpolated_1d - s							# remove mean level
zz_interpolated = zz_interpolated_1d.reshape(zz_interpolated.shape)				# reshape

surfaces.append(s.reshape(zz_interpolated.shape))
extents.append(args.w)

np.savetxt("instresp.dat", m)

# this is fitting per row/column
'''print "Fitting low frequency structure in x..."
xx_interpolated = xx_interpolated.transpose()
zz_interpolated = zz_interpolated.transpose()
for idx, (this_xx, this_zz) in enumerate(zip(xx_interpolated, zz_interpolated)):
  m = np.polyfit(this_xx, this_zz, deg=args.o)
  zz_interpolated[idx] = zz_interpolated[idx] - np.polyval(m, xx_interpolated[idx])
xx_interpolated = xx_interpolated.transpose()	# switch back
zz_interpolated = zz_interpolated.transpose()
print_stats(zz_interpolated)

print "Fitting low frequency structure in y..."
for idx, (this_yy, this_zz) in enumerate(zip(yy_interpolated, zz_interpolated)):
  m = np.polyfit(this_yy, this_zz, deg=args.o)
  zz_interpolated[idx] = zz_interpolated[idx] - np.polyval(m, yy_interpolated[idx])'''

print_stats(zz_interpolated)
surfaces.append(zz_interpolated)
extents.append(args.w)

# PLOTS
# -----
# 2D
if args.p2D:
  plt.figure()
  plt.subplot(231)
  plt.title("data in")
  plt.imshow(surfaces[0].transpose(), interpolation=None, origin='lower', extent=extents[0])
  plt.colorbar()
  plt.subplot(232)
  plt.title("median filtered")
  plt.imshow(surfaces[1].transpose(), interpolation=None, origin='lower', extent=extents[1])
  plt.colorbar()
  plt.subplot(233)
  plt.title("tilt removed")
  plt.imshow(surfaces[2].transpose(), interpolation=None, origin='lower', extent=extents[2])
  plt.colorbar()
  plt.subplot(234)
  plt.title("fitted surface to be removed")
  plt.imshow(surfaces[3].transpose(), interpolation=None, origin='lower', extent=extents[3])
  plt.colorbar()
  plt.subplot(235)
  plt.title("residual after surface removal")
  plt.imshow(surfaces[4].transpose(), interpolation=None, origin='lower', extent=extents[4])
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
  pyfits.writeto("instresp.fits", zz_interpolated, hdr)




