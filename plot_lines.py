import sys
import numpy as np
import pylab as plt

for f in sys.argv[1].split(','):
	# read data file
	x       = []
	y       = []  
	d       = []
	d_err 	= []
	with open(f) as data:
	  for line in data:
	    if line.startswith('#'):
	      continue
	    res = line.split()
	    this_x = float(res[0])
	    this_y = float(res[1])
	    this_d = float(res[2])
	    this_d_err = float(res[3])
	    
	    x.append(this_x)
	    y.append(this_y)
	    d.append(this_d)
	    d_err.append(this_d_err)

	c = np.polyfit(x, d, 1)
	d_slope = np.polyval(c, x)
	plt.plot(x, d-d_slope, label="y= " + f)
plt.legend()
plt.ylim([-50,50])
plt.show()
