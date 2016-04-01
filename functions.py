import itertools

import numpy as np
from scipy.interpolate import interp1d

def correct_for_temperature(cat, d, t):
    print "Applying temperature correction..."
    c = np.load(cat)
    offsets = np.polyval(c, t)
    d_cal = d - offsets
    return d_cal

def correct_stage_positions(cax, cay, x, y):
    print "Applying calibration to stages..."
    stage_x_rd = []
    stage_x_act = []
    with open(cax) as cal:
        for line in cal:
            try:
                stage_x_rd.append(float(line.split()[0]))
                stage_x_act.append(float(line.split()[1]))
            except ValueError:
                pass
    fx = interp1d(stage_x_rd, stage_x_act, kind='cubic')

    stage_y_rd = []
    stage_y_act = []
    with open(cay) as cal:
        for line in cal:
            try:
                stage_y_rd.append(float(line.split()[0]))
                stage_y_act.append(float(line.split()[1]))
            except ValueError:
                pass
    fy = interp1d(stage_y_rd, stage_y_act, kind='cubic')

    x = fx(x)
    y = fy(y)

    return x, y

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

def print_stats(array):
    print
    print "5th percentile of data is: " + '\t\t\t' + str(round(np.nanpercentile(array, 5), 2)) + "um"
    print "95th percentile of data is: " + '\t\t\t' + str(round(np.nanpercentile(array, 95), 2)) + "um"
    print "Peak-to-peak amplitude of structure is: " + '\t' + str(round(np.nanpercentile(array, 95)-np.nanpercentile(array, 5), 2)) + "um"
    print "Half peak-to-peak amplitude of structure is: " + '\t' + str(round((np.nanpercentile(array, 95)-np.nanpercentile(array, 5))/2, 2)) + "um"
    print 

def read_data_file(f, de, w):
    print "Reading input data..."
    TIME  = []
    x	  = []
    y	  = []
    d	  = []
    d_err = [] 
    t 	  = []
    h     = []
    with open(f) as data:
        for line in data:
            if line.startswith('#'):
                continue
            res = line.split()
            this_time = float(res[0])
            try:
                this_x = float(res[1])
                this_y = float(res[2])
            except ValueError:
                this_x = -1
                this_y = -1
            this_d = float(res[3])
            this_d_err = float(res[4])   
            this_t = float(res[5])   
            this_h = float(res[6])   
            if w is not None:
                if any(w):				# if we don't have a default for the window
                    if any([this_x < w[0],
                           this_x > w[1],
                           this_y < w[2],
                           this_y > w[3]]):
                         continue
            if this_d_err>de or this_d==0:       	# additional quality checks
                continue
            TIME.append(this_time)
            x.append(this_x)
            y.append(this_y)
            d.append(this_d)
            d_err.append(this_d_err)
            t.append(this_t)
            h.append(this_h)
    
    return TIME, x, y, d, d_err, t, h

def read_processed_scan_file(f, de):
    print "Reading input data..."
    x	  = []
    y	  = []
    d	  = []
    d_err = []
    with open(f) as data:
        for line in data:
            if line.startswith('#'):
                continue
            res = line.split()
            this_x = float(res[0])
            this_y = float(res[1])
            this_d = float(res[2])
            try:
                this_d_err = float(res[3])    
                if this_d_err>de or this_d==0:       	# additional quality checks
                    continue
                d_err.append(this_d_err)
            except IndexError:
                x.append(this_x)
                y.append(this_y)
                d.append(this_d)
            x.append(this_x)
            y.append(this_y)
            d.append(this_d)
    
    return x, y, d, d_err
