import sys
import argparse
import copy
import os

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import pyfits

from functions import read_data_file, correct_stage_positions, correct_for_temperature

F_CALX = "calibration_files/45855915_cal.dat"
F_CALY = "calibration_files/45855915_cal.dat"
F_CORT = "calibration_files/t_cor.npy"

parser = argparse.ArgumentParser()
parser.add_argument('--f', help="path to flat run", action="store", type=str)
parser.add_argument('--w', help="window to use data points (x1,x2,y1,y2).", action="store", type=str, default='10,290,10,290')
parser.add_argument('--de', help="maximum error in distance measurement", action="store", type=float, default=1.0)
parser.add_argument('--nc', help="don't apply stage calibration files", action="store_false")
parser.add_argument('--nt', help="don't apply temperature correction file", action="store_false")
parser.add_argument('--p', help="make plots?", action="store_true")
args = parser.parse_args()

args.w = ([float(n) for n in args.w.split(',')])

run_data_x = []
run_data_y = []
run_data_d = []
run_data_d_err = []
for f in os.listdir(args.f):
    print "processing file " + f
    if f.endswith(".log"):
        print "ignoring file"
        continue
    fi = args.f.rstrip('/') + '/' + f
    # READ INPUT DATA FILE
    # --------------------
    this_time, x, y, d, d_err, t, h = read_data_file(fi, args.de, args.w)

    # CORRECT FOR NONLINEAR MOTION OF STAGES
    # --------------------------------------
    if args.nc:
        x, y = correct_stage_positions(F_CALX, F_CALY, x, y)

    # CORRECT FOR TEMPERATURE EXCURSIONS
    # ----------------------------------
    if args.nt:
        d = correct_for_temperature(F_CORT, d, t)

    run_data_x.append(np.mean(x))
    run_data_y.append(np.mean(y))
    run_data_d.append(np.mean(d))
    run_data_d_err.append(np.mean(d_err))
    print

with open("flat.dat", 'w') as out:
    for i in range(len(run_data_x)):
        if not np.isnan(run_data_x[i]):
            out.write(str(round(run_data_x[i], 2)) + '\t' + str(round(run_data_y[i], 2)) + '\t' + str(round(run_data_d[i], 2)) + '\t' + str(round(run_data_d_err[i], 2)) + '\n')

if args.p:
    plt.scatter(run_data_x, run_data_y, c=run_data_d, s=50, linewidth=0, marker='s', cmap="autumn")
    plt.colorbar()
    plt.show()



