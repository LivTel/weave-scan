import sys
import argparse

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from mayavi import mlab
from scipy import interpolate
from scipy.interpolate import griddata

NORMALISE = False
INCREMENT = (0.05, 0.05)

parser = argparse.ArgumentParser()
parser.add_argument("-n", help="Normalise", action="store_true")
args = parser.parse_args()

for fi_idx, fi in enumerate(sys.argv[1:]):
  x 		= []
  y 		= []
  d 		= []
  d_err 	= []
  WINDOW_COORDS = []
  with open(fi) as f:
    for line in f:
      if line.startswith('#'):
        if 'WINDOW' in line:
          xl = float(line.split(':')[1].split(',')[0])
          xh = float(line.split(':')[1].split(',')[1])
          yl = float(line.split(':')[1].split(',')[2])
          yh = float(line.split(':')[1].split(',')[3])
          WINDOW_COORDS.append((xl,xh))
          WINDOW_COORDS.append((yl,yh))
      elif line == '':
        pass
      else:
        res = line.split()
        this_x = float(res[0])
        this_y = float(res[1])
        this_d = float(res[2])
        this_d_err = float(res[3])
        if this_x < WINDOW_COORDS[0][0] or\
           this_x > WINDOW_COORDS[0][1] or\
           this_y < WINDOW_COORDS[1][0] or\
           this_y > WINDOW_COORDS[1][1]:
          pass
        else:
          # additional quality criteria
          if this_d_err>2.0 or this_d==0:
            pass
          else:
            x.append(this_x)
            y.append(this_y)
            d.append(this_d*-1)
            d_err.append(this_d_err)
           
  xx_interpolated, yy_interpolated = np.mgrid[WINDOW_COORDS[0][0]:WINDOW_COORDS[0][1]:0.05, WINDOW_COORDS[1][0]:WINDOW_COORDS[1][1]:0.05]
  zz_interpolated = griddata((x, y), d, (xx_interpolated, yy_interpolated), method='cubic')
  
  zdepth = None
  if zdepth is None:
    zdepth = (np.nanpercentile(zz_interpolated, 0.5), np.nanpercentile(zz_interpolated, 99.5))
  
  plt.imshow(zz_interpolated, vmin=zdepth[0], vmax=zdepth[1], origin='lower', aspect='auto', extent=[WINDOW_COORDS[0][0],WINDOW_COORDS[0][1],WINDOW_COORDS[1][0],WINDOW_COORDS[1][1]])
  plt.colorbar()
  plt.show()

  '''if NORMALISE:
    # interpolate surface to smooth and remove Nones
    f_interpolated = interpolate.RectBivariateSpline([y[0] for y in yy], xx[0], dd)
 
    # construct flat surface using corners of window from interpolated (smoothed) surface
    x1_surf = WINDOW_COORDS[0][0]
    x2_surf = WINDOW_COORDS[0][1]
    y1_surf = WINDOW_COORDS[1][0]
    y2_surf = WINDOW_COORDS[1][1]
    d_x1_y1_surf = float(f_interpolated(WINDOW_COORDS[1][0],WINDOW_COORDS[0][0]))
    d_x1_y2_surf = float(f_interpolated(WINDOW_COORDS[1][1],WINDOW_COORDS[0][0]))
    d_x2_y1_surf = float(f_interpolated(WINDOW_COORDS[1][0],WINDOW_COORDS[0][1]))
    d_x2_y2_surf = float(f_interpolated(WINDOW_COORDS[1][1],WINDOW_COORDS[0][1]))
  
    xx_surf, yy_surf = np.meshgrid([x1_surf, x2_surf], [y1_surf, y2_surf])
    dd_surf = [[d_x1_y1_surf ,d_x2_y1_surf],[d_x1_y2_surf ,d_x2_y2_surf]] 
    f_surface = interpolate.interp2d(xx_surf, yy_surf, dd_surf, kind='linear')

    # subtract flat surface from interpolated surface
    dd_subtracted = np.zeros(dd.shape)
    for idx_j in range(xx.shape[0]):
      for idx_i in range(yy.shape[1]):
        this_x = xx[idx_j, idx_i]
        this_y = yy[idx_j, idx_i]
        dd_subtracted[idx_j, idx_i] = f_interpolated(this_y, this_x) - f_surface(this_x, this_y)
    dd = dd_subtracted

  # get rid of offset
  xx_min = np.zeros(xx.shape)
  xx_min.fill(np.min(xx))
  xx_nooffset = xx-xx_min
  yy_min = np.zeros(yy.shape)
  yy_min.fill(np.min(yy))
  yy_nooffset = yy-yy_min


  # 3d
  ax = fig.gca(projection='3d')
  ax.plot_surface(xx_nooffset, yy_nooffset, dd, rstride=1, cstride=1, cmap='Oranges')

  # 3d (lines)
  plt.plot(np.mean(xx_nooffset, axis=0), [0 for i in range(len(xx_nooffset[0]))], np.mean(dd, axis=0), pl_colors[fi_idx] + 'x-', linewidth=3.0, label=fi)
  plt.plot([0 for i in range(len(yy_nooffset.transpose()[0]))], np.mean(yy_nooffset.transpose(), axis=0), np.mean(dd.transpose(), axis=0), pl_colors[fi_idx] + 'x-', linewidth=3.0)
  plt.legend()
  
  
  # 2d
  plt.subplot(211)
  plt.plot(np.mean(xx_nooffset, axis=0), np.median(dd, axis=0), pl_colors[fi_idx] + 'x-', linewidth=2.0, label=fi + ' (x) ')
  plt.legend()
  plt.xlabel('position along shortest axis of plate (mm)')
  plt.ylabel('deformation from flat surface (um)')
  plt.subplot(212)
  plt.plot(np.mean(yy_nooffset.transpose(), axis=0), np.median(dd.transpose(), axis=0), pl_colors[fi_idx] + 'x--', linewidth=2.0, label=fi + ' (y) ')
  plt.xlabel('position along longest axis of plate (mm)')
  plt.ylabel('deformation from flat surface (um)')
  plt.legend(loc='lower right')'''
  
#plt.show()



