import os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--f', help="file title prefix", action="store", type=str, default="res")
    parser.add_argument('--i', help="increment to use (mm)", action="store", type=float, default=0.1)
    parser.add_argument('--wx', help="add window centres x (mm)", action="append", type=str)
    parser.add_argument('--wy', help="add window centres y (mm)", action="append", type=str)
    parser.add_argument('--b', help="size of box to average measurements over", type=float, default=1.5)
    parser.add_argument('--n', help="number of times to do scan", type=int, default=1)
    args = parser.parse_args()
    
    for ni in range(args.n):
        for idx_wx, wiy in enumerate(args.wy):
            yi = float(wiy.strip())
            for idx_wy, wix in enumerate(args.wx):
                xi = float(wix.strip())
                file_title = args.f + "_" + str(ni) + "_" + str(idx_wx) + "_" + str(idx_wy) + "_" + str(xi) + "_" + str(yi) + "_" + str(args.i)
                lower_xi = xi - args.b/2.
                upper_xi = xi + args.b/2.
                lower_yi = yi - args.b/2.
                upper_yi = yi + args.b/2.
                os.system("python scan.py " + file_title + " --w " + str(lower_xi) + "," + str(upper_xi) + "," + str(lower_yi) + "," + str(upper_yi) + " --sxi " + str(args.i) + " --syi " + str(args.i) + " --m --s")

