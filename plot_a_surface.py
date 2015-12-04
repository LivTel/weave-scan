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

def calc_R(x,y, xc, yc):
  return np.sqrt((x-xc)**2 + (y-yc)**2)

def f(c, x, y):
  Ri = calc_R(x, y, *c)
  return Ri - Ri.mean()

def leastsq_circle(x,y,x_m,y_m):
  center_estimate = x_m, y_m
  center, ier = optimize.leastsq(f, center_estimate, args=(x,y))
  xc, yc = center
  Ri     = calc_R(x, y, *center)
  R      = np.median(Ri)
  rms    = np.sqrt(np.mean((Ri - R)**2))
  return xc, yc, R, rms

parser = argparse.ArgumentParser()
parser.add_argument('f', help="path to file", action="store", type=str)
parser.add_argument('--w', help="window to grid over, defaults to min/max of data (x1,x2,y1,y2)", action="store", type=str, default='0, 0, 0, 0')
parser.add_argument('--i', help="increment of interpolated grid as CSV tuple (x,y)", action="store", type=str, default='1, 1')
parser.add_argument('--de', help="maximum error in distance measurement", action="store", type=float, default=1.0)
parser.add_argument('--ca', help="apply stage calibration files", action="store_true")
parser.add_argument('--cax', help="calibration file for x if stage calibration is to be applied", action="store", type=str, default='calibration_files/45855915_cal.dat')
parser.add_argument('--cay', help="calibration file for y if stage calibration is to be applied", action="store", type=str, default='calibration_files/45855916_cal.dat')
parser.add_argument('--fl', help="add window for flattening as CSV tuple (x,y)", action='append')
parser.add_argument('--f1', help="maximum number of points with same x or y before the coordinate is disregarded", action="store", type=str, default=50)
parser.add_argument('--fc', help="attempt a circular fit to surface (surface must be isolated by NaN!)", action="store_true")
parser.add_argument('--fco', help="explictly define coordinates of circle centre best guess as CSV tuple (x, y)", action="store", type=str, default='0, 0')
parser.add_argument('--p2D', help="plot 2D", action="store_true")
parser.add_argument('--z', help="z-depth, or rather distance, as CSV tuple (high, low)", action="store", type=str, default='0, 0')
parser.add_argument('--p3D', help="plot 3D", action="store_true")
parser.add_argument('--rs', help="row stride (3d plot only)", action="store", type=float, default=5)
parser.add_argument('--cs', help="column stride (3d plot only)", action="store", type=float, default=5)
parser.add_argument('--fits', help="make fits file", action="store_true")
args = parser.parse_args()

args.i = ([float(n) for n in args.i.split(',')])
args.w = ([float(n) for n in args.w.split(',')])
args.z = ([float(n) for n in args.z.split(',')])
args.fco = ([float(n) for n in args.fco.split(',')])
if args.fl is not None:
  args.fl = [tuple(w.split(',')) for w in args.fl]

#args.fl = [(150,12), (290, 150), (200, 293), (8,200), (50, 50), (250, 50), (53,264)]		# cryostat
#args.fl = [(70, 115), (70, 205), (235, 115), (235,205)]					# window
#args.fl = [(50,90), (60,225), (248,90), (248,225)]						# retainer plate

# read data file
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
    if any(args.w):                           # if we don't have a default for the window
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
    
if args.ca:
  # load stage calibration files
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

  x_old = list(x)

  # correct for stage positions using calibration file
  x = fx(x)
  y = fy(y)

  ## DEBUG FOR STAGE CALIBRATION CURVE
  '''plt.plot(stage_x_rd, [x1-x2 for x1, x2 in zip(stage_x_act,stage_x_rd)])
  plt.plot(x_old, [x1-x2 for x1, x2 in zip(x,x_old)])
  plt.show()'''

# we can now set the default window to min/max of x,y
if not any(args.w):                                                   
  args.w = [round(min(x), 2), round(max(x), 2), round(min(y), 2), round(max(y), 2)]

## DEBUG TO SHOW DEVIATION FOR NO INTERPOLATION
plt.scatter(x, y, c=d, s=15, linewidth=0, marker='s', cmap="autumn")
plt.colorbar()
plt.show()
  
# create regular grid from data
xx_interpolated, yy_interpolated = np.mgrid[args.w[0]:args.w[1]:args.i[0], args.w[2]:args.w[3]:args.i[1]]
zz_interpolated = griddata((x, y), d, (xx_interpolated, yy_interpolated), method='linear', fill_value=0)

'''# FLATTEN
# -------

if args.fl is not None:
  v = []
  for c in args.fl:
    this_window_centre = np.where(np.logical_and(xx_interpolated==float(c[0]), yy_interpolated==float(c[1])))
    this_window_centre_x = this_window_centre[0][0]
    this_window_centre_y = this_window_centre[1][0]

    v.append(np.nanmedian(zz_interpolated[this_window_centre_x-1:this_window_centre_x+2, this_window_centre_y-1:this_window_centre_y+2]))

  m = polyfit2d(np.array([float(x) for x,y in args.fl]), np.array([float(y) for x,y in args.fl]), np.array(v))

  for jj in range(zz_interpolated.shape[0]):
    for ii in range(zz_interpolated.shape[1]):
      if zz_interpolated[jj, ii] > 0:
        zz_interpolated[jj, ii] = (zz_interpolated[jj, ii]-polyval2d(xx_interpolated[jj, ii], yy_interpolated[jj, ii], m))*-1'''
# -------

print_stats(zz_interpolated)

# MEDIAN FILTER TO SMOOTH NOISE
# -----------------------------
args.m=[3,3]
print "Applying median filter..."
zz_interpolated = signal.medfilt(np.nan_to_num(zz_interpolated), kernel_size=[args.m[0],args.m[1]])
xx_interpolated = xx_interpolated[(args.m[1]-1)/2:-(args.m[1]-1)/2,(args.m[0]-1)/2:-(args.m[0]-1)/2]		# trim off border pixels where filter will have yielded 0
yy_interpolated = yy_interpolated[(args.m[1]-1)/2:-(args.m[1]-1)/2,(args.m[0]-1)/2:-(args.m[0]-1)/2]		# trim off border pixels where filter will have yielded 0
zz_interpolated = zz_interpolated[(args.m[1]-1)/2:-(args.m[1]-1)/2,(args.m[0]-1)/2:-(args.m[0]-1)/2]		# trim off border pixels where filter will have yielded 0
args.w = [args.w[0]+(args.m[0]-1)/2, args.w[1]-(args.m[0]-1)/2, args.w[2]+(args.m[1]-1)/2, args.w[3]-(args.m[1]-1)/2]
print_stats(zz_interpolated)

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

s = np.loadtxt("calibration_files/instresp.dat")
# this is fitting for surface
xx_interpolated_1d = xx_interpolated.flatten()							# cast to 1d
yy_interpolated_1d = yy_interpolated.flatten()							
zz_interpolated_1d = zz_interpolated.flatten()
zz_interpolated_1d = zz_interpolated_1d - polyval2d(xx_interpolated_1d, yy_interpolated_1d, s)	# remove mean level
zz_interpolated = zz_interpolated_1d.reshape(zz_interpolated.shape)				# reshape

print_stats(zz_interpolated)

# if we haven't explictly set a zdepth, use percentiles
if not any(args.z):                                                   
  args.z = (np.nanpercentile(zz_interpolated, 0.5), np.nanpercentile(zz_interpolated, 99.5))

# CIRCLE FITTING
# --------------
if args.fc: 
  # find edges of data where pixels have >=3 adjacent NaNs
  x_edge = []
  y_edge = []
  for jj in range(zz_interpolated.shape[0]):
    for ii in range(zz_interpolated.shape[1]):
      if np.isnan(zz_interpolated[jj, ii]):                             # if this pixel is itself NaN
        continue
      else:
        c_NaN = 0
        for jj_window in range(-1,2):                                   # make a 3x3 pixel window around this pixel
          for ii_window in range(-1,2):
            if jj_window == ii_window == 0:                             # ignoring the current pixel
              continue
            elif jj+jj_window >= zz_interpolated.shape[0] or\
                 jj+jj_window <= 0 or\
                 ii+ii_window >= zz_interpolated.shape[1] or\
                 ii+ii_window <= 0:                                     # making sure it's not out of bounds
              continue
            elif np.isnan(zz_interpolated[jj+jj_window, ii+ii_window]): # ascertaining has many bordering NaNs
              c_NaN = c_NaN + 1
        if c_NaN >= 3:                                                  # and assigning as edge if 3 or more are detected
          x_edge.append(xx_interpolated[jj+jj_window, ii+ii_window])
          y_edge.append(yy_interpolated[jj+jj_window, ii+ii_window])

  # filter for erroneous points where there are more than args.f1 values in the same coordinate
  ## x
  for k, v in Counter(x_edge).iteritems():
    if v>=args.f1:
      for idx, xi in enumerate(x_edge):
        if xi==k:
          x_edge[idx] = np.nan
          y_edge[idx] = np.nan
  ## y
  for k, v in Counter(y_edge).iteritems():
    if v>=args.f1:
      for idx, yi in enumerate(y_edge):
        if yi==k:
          x_edge[idx] = np.nan
          y_edge[idx] = np.nan
        
  # remove NaNs, otherwise we have a problem with fit
  x_edge = np.array(x_edge)
  y_edge = np.array(y_edge)
  x_edge = x_edge[np.logical_not(np.isnan(x_edge))]
  y_edge = y_edge[np.logical_not(np.isnan(y_edge))]

  # set best guess circle centre from input argument, otherwise guess using min/max of zz_interpolated array
  if not any(args.fco):      
    args.fco[0] = xx_interpolated.flatten()[np.nanargmin(zz_interpolated)]
    args.fco[1] = yy_interpolated.flatten()[np.nanargmin(zz_interpolated)]
  circ_fit_vals = leastsq_circle(x_edge, y_edge, args.fco[0], args.fco[1])
  print "Fit (circle) x centroid at " + str(round(circ_fit_vals[0], 3)) + "mm."
  print "Fit (circle) y centroid at " + str(round(circ_fit_vals[1], 3)) + "mm."
  print "Fit (circle) radius is " + str(round(circ_fit_vals[2], 3)) + "mm." 
  print "Fit (circle) RMS of radius is " + str(round(circ_fit_vals[3]*1000, 1)) + "um."
  
# PLOTS
# -----
# 2D
if args.p2D:
  if args.fc:
    plt.axis('equal')
    theta_fit = np.linspace(-np.pi, np.pi, 180)
    x_fit = circ_fit_vals[0] + circ_fit_vals[2]*np.cos(theta_fit)
    y_fit = circ_fit_vals[1] + circ_fit_vals[2]*np.sin(theta_fit)
    plt.plot(x_fit, y_fit, 'b-' , label="fitted circle", lw=2)
    plt.plot([circ_fit_vals[0]], [circ_fit_vals[1]], 'bo', mec='y', mew=1, label="centroid") 
    plt.plot(x_edge, y_edge, 'rx', label='data', mew=1)
    plt.grid() 
    plt.legend(numpoints=1)
  plt.imshow(zz_interpolated.transpose(), interpolation=None, vmin=args.z[0], vmax=args.z[1], origin='lower', extent=args.w)
  plt.colorbar()
  plt.show()

# 3D  
if args.p3D:
  fig = plt.figure()
  ax = fig.gca(projection='3d')
  ax.plot_surface(xx_interpolated, yy_interpolated, zz_interpolated, rstride=args.rs, cstride=args.cs, cmap=cm.coolwarm, alpha=0.3)
  if args.fc:
    ax.plot([circ_fit_vals[0], circ_fit_vals[0]], [circ_fit_vals[1], circ_fit_vals[1]], [args.z[0], args.z[1]], c='blue')
  plt.show()



# fits
if args.fits:
  hdr = pyfits.core.Header()
  pyfits.writeto(args.f + ".fits", zz_interpolated, hdr)



