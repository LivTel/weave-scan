import sys
import numpy as np
import pylab as plt
import argparse

from functions import read_data_file, correct_stage_positions

parser = argparse.ArgumentParser()
parser.add_argument('--fl', help="add file to list", action='append')
parser.add_argument('--de', help="maximum error in distance measurement", action="store", type=float, default=1.0)
parser.add_argument('--w', help="window, defaults to min/max of data (x1,x2,y1,y2)", action="store", type=str, default='1,300,1,300')
parser.add_argument('--ca', help="don't apply stage calibration files", action="store_false")
parser.add_argument('--cax', help="calibration file for x if stage calibration is to be applied", action="store", type=str, default='calibration_files/45855915_cal.dat')
parser.add_argument('--cay', help="calibration file for y if stage calibration is to be applied", action="store", type=str, default='calibration_files/45855916_cal.dat')
args = parser.parse_args()

args.w = ([float(n) for n in args.w.split(',')])

for f in args.fl:
    # READ INPUT DATA FILE(S)
    # -----------------------
    x, y, d, de = read_data_file(f, args.de, args.w)

    # CORRECT FOR NONLINEAR MOTION OF STAGES
    # --------------------------------------
    if args.ca:
        x, y = correct_stage_positions(args.cax, args.cay, x, y)

    c = np.polyfit(x, d, 1)
    d_slope = np.polyval(c, x)
    plt.plot(x, d-d_slope, label="y= " + f)
plt.legend()
plt.ylim([-50,50])
plt.show()
