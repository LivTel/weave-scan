weave-scan
=============

# Overview

This package provides tools to drive the Thorlabs LTS300/M translation stages and STIL confocal sensor both
independently and in unison (as a surface scanner). 

Calibration and visualation scripts are also included.

# Translation Stage Command Set

The serial number, [**serial**], for each stage is printed on the side nearest to the motor drive. You must 
be superuser to execute these commands, or use `sudo`.

* Home the stages (req. after power cycle)

`sudo python home.py [serial]`

* Move a stage to an absolute position

`sudo python goto.py [serial] [position_in_mm]`

* Move a stage to a relative position

`sudo python move.py [serial] [position_in_mm]`

* Get the current position of a stage

`sudo python get_position.py [serial]`

# Sensor and Scanning Command Set

Scanning functionality is provided by the `scan.py` script. You must be superuser to execute these commands, 
or use `sudo`. `scan.py` expects any temperature and humidity data to be sent tab-separated over serial (ttyUSB0, though 
this can be changed with the `--dt` flag). The baud rate can also be changed with the `--dtb` flag (default is 9600). 
Similar flags exist for the stage device name and baud (`--ds` and `--dsb` respectively).

* Print help blurb

`sudo python scan.py --h`

* Initiate continous scan of the sensor head, output data to [**file**]

`sudo python scan.py [file] --r`

* Move sensor head to starting position (`--m`) defined by window (`--w`) parameter, begin scan (`--s`) and 
read data to [**file**]. Scan intervals are set by the `--sxi` and `--syi` flags. Units are in mm. Window 
coordinates are measured as (*x1*,*x2*,*y1*,*y2*), e.g. to scan a window from 0 to 300mm in both x and y, 
with an increment of 1mm:

`sudo python scan.py [file] --m --w 0,300,0,300 --s --sxi 1 --syi 1`

It is also possible to scan in *aperture* mode using `multi_scan.py`. In this mode, a series of apertures with 
x and y coordinates defined by the `--wx` and `--wy` *appendable* flags are scanned over a box of size `--b`, `--n` times 
(default 1) with an increment of `--i`. Seperate files outputted for each, e.g. for four apertures of (150, 150),
 (150, 200), (200, 150), and (200, 200):

`sudo python multi_scan.py --wx 150 --wx 200 --wy 150 --wy 200 --b 1 --i 0.1 --n 1`

## Enabling a Web Interface

To enable a web interface for environment monitoring, symlink the `weave-scan` directory to a directory served by 
apache, e.g. `/var/www/html`. The (`--p`) flag must have been used when invoking scans in order to generate 
depth, temperature and humidity plots.

# Data Calibration

Data recorded by the scanner is subject to instrument errors, both repeatable and random, arising from the 
shape of the translation stages (giving error of *z* ~20micron) and temperature fluctuations (giving error of 
*z* ~5micron). There are a series of calibration files that can be used to compensate for these effects if the 
desired depth accuracy is to be lower than a few tens of micron. A description of the procedure to generate these files 
is given in the following sections.

Additonally, two files are provided as per specification from the manufacturer to map encoder values to actual values 
for both the x and y stages. These are found in `calibration_files/4585591[56?]_cal.dat`. These files will 
automatically be used by any routine with stage x/y position correction as an option unless the `--nc` flag 
is specified.

## Building a Temperature Calibration File

Depending on ambient air temperature, the scanner reading drifts non-linearly to the order of ~3micron/deg. To correct for 
this, the sensor is read out without moving the stages and the air-con turned off (only if the temperature has already 
stabilised at the air-con set temperature). Obviously the temperature sensor must be attached. 

### 1. Generating the Data

Data can be generated by:

`sudo python scan.py test --r`

which produces a tab-separated-variable file (TSV), [**correction\_file**] with lines like:

`1456321257.35	-	-	6471.11937359	0.461643712486	22.94	32.26`

corresponding to time (unix, s), x stage position (blank, not required), y stage position (blank, not required), distance (micron), 
distance error, temperature (degC) and humidity (%).

### 2. Processing

We process this data using the script `process_run_for_temperature_correction.py`:

`python process_run_for_temperature_correction.py --f [correction_file] --p`

which derives a correction for the temperature range defined by the `--tr` flag. The `--p` flag allows visual inspection 
of the resulting fit. Don't use too high an order (`--o`) otherwise it may barrel if it needs to extrapolate. However, if 
your [**correction\_file**] covers the range of temperatures over which the scans will happen, it be more accurate to use 
a higher order.

This outputs a correction array file `t_cor.npy` which can be loaded by moving to the `calibration_files/` directory. The
calibration will automatically be used by any routine with temperature correction as an option unless the `--nt` flag is 
specified. If no temperature file is generated, you must also use the `--nt` flag otherwise it will fault when looking for 
the file.

## Building an Instrumental Flat

We generate an instrumental **flat** (or flats) in order to map the warping of the translation stages, which 
produces a shift in z reading for a flat surface as a function of x and y. This should be done everytime 
the stages are moved; granite can act as an optical flat for this purpose. Typically the instrumental response 
follows a saddle-like shape:

![alt text](https://raw.githubusercontent.com/LivTel/weave-scan/master/img/instsurf.png "Instrumental shape")

### 1. Generate the Data

Unfortunately, although granite is flat to ~1um, the surface quality appears can be fairly rough 
(maybe 1-10um depending on grain size?). We counter this by measuring the average sensor reading in 
a series of **apertures** on the granite. Each **run** has a series of **scans**, with each scan making 
measurements of a single aperture.

Consequently, the input should be a series of files generated using the `multi_scan.py` script. Although 
we try to compensate for temperature, it is best to keep temperature fluctuations to a minimum, so, if 
possible, it's better to schedule the scan using the "at" command for late evening (~8pm) with the aircon off:

`$ sudo at 20:00 February 22`

`$ at> python ./multi_scan.py [parameters]` 

where a description of [parameters] is given above and in the help blurb `--h`.

This produces another TSV file with entries as discussed in **Building a Temperature Calibration File** above.

### 2. Processing

Once the data has been acquired, we process the run using the `process_run_for_flat.py` script. This 
produces a single output file, `flat.dat`, containing (optionally position calibrated) stage x/y positions 
and corresponding (optionally temperature calibrated) averaged measurement values/errors. Adding several 
directories to the command will use the average over all runs, but each run must be exactly the same in terms 
of aperture x/y coordinates for this to make sense. 

`$ python process_run_for_flat.py --d results/1/ --d results/2/ --d results/3/`

Once the average flat is made, we can find the instrumental signature by removing the mean plane and finding a 
suitable interpolant function. This is done by using the `find_instrumental_signature.py` script. The interpolation 
increment (`--i`), interpolant order (`--o`) and smoothing factor (`--s`) may need to be adjusted to get the 
script to output sensibly. The plots flag (`--p`) is particularly useful for eyeballing the result. Be wary 
of overfitting.

`$ python find_instrumental_signature.py flat.dat --i 5,5 --p`

This outputs a correction array file `instresp.npy` which can be loaded by moving to the `calibration_files/` 
directory.  The calibration will automatically be used by any routine with instrumental flag correction as an option 
unless the `--nf` flag is specified. If no instrumental flat file is generated, you must also use the `--nf` flag 
otherwise it will fault when looking for the file.

# Plotting a Surface

To plot a surface, try the `plot_a_surface.py` script. Although the script is fairly trivial to invoke, plotting 
3D surfaces can be heavily nuanced by the variety data being plotted. It is only really provided as a template.

**You will need to alter this script to match what you're trying to plot!**

`python plot_a_surface.py results/granite_5_36_1.5x1.5_0.1x0.1/`







