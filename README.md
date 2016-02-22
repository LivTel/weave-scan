The procedure is as follows.

We must first generate a "flat" to map the distortion of the translation stages. We have been doing this with granite.

A "run" is defined as a sequence of scans, with each scan containing several apertures taken over the length and width of the flat surface.

Scans are conducted using "multi\_scan\_window.csh". We must be wary of temperature fluctuations, so it's better to schedule the scan using the "at" command, and schedule it for late evening (~8pm) with the aircon off:

$ sudo at 20:00 February 22
$ at> ./multi\_scan\_window.csh [RESULTS\_FOLDER]

Empirically this has been found to give us a repeatability of a few micron, typically 1-2.

We process this scan using "process_scan.py". This produces an output file (typically "flat.0") containing calibrated x/y positions (from the stage calibration files), and corresponding averaged d values with respective errors.

$ python process\_scan.py --d results/1/ --d results/2/ --d results/3/ --p

And then find the instrumental signature of these averaged results

$ python find\_instrumental\_signature.py flat.0 --p

This outputs an interpolate "instresp.dat.npy". Load this by moving the file to calibration_files/. 

This interpolate is loaded in "plot\_a\_surface.py":

$ python plot\_a\_surface.py results/3/ --p




