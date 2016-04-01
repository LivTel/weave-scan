weave-scan
=============

# Overview

This package provides tools to drive the Thorlabs LTS300/M translation stages and STIL confocal sensor both
independently and in unison (as a scan).

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




# To enable web tracking of environment and sensor readings

Symlink the weave-scan directory to a directory served by apache, e.g. /var/www/html.

# To plot a surface

To plot a surface, use the `plot_a_surface.py` script:

`python plot_a_surface.py results/granite_5_36_1.5x1.5_0.1x0.1/`

# To generate a flat

We must generate a **flat** to map the distortion of the translation stages. This should presumably be done everytime we move the stages. We have been doing this with granite. 

Unfortunately, although the granite is flat to ~1um, the surface quality appears to be fairly rough (maybe 1-10um depending on grain size?). We counter this by measuring the average sensor reading in a series of **apertures** on the granite. Each **run** has a series of **scans**, with each scan making measurements of a single aperture.

A series of python scripts are used to generate the flat. Each has a help flag (--h) which will give some insight into parameters, but generally speaking the defaults should be OK for granite.

Scans are conducted using the `multi_scan.py` script. Although we try to compensate for temperature, it is best to keep temperature fluctuations to a minimum, so, if possible, it's better to schedule the scan using the "at" command for late evening (~8pm) with the aircon off. In room 3.02 this has been empirically found to give a fluctuation of 1-2 micron for a 12h run:

`$ sudo at 20:00 February 22`

`$ at> python ./multi_scan.py` 

**TODO: PARAMETERS**

This produces a file for each aperture. The file contains information on time, stage position, measurement, measurement error and some environment readings.

We then process the run using the `process_run_for_flat.py` script. This produces a single output file, "flat.dat", containing (optionally position calibrated) stage x/y positions and corresponding (optionally temperature calibrated) averaged measurement values/errors. Adding several directories to the command will use the average over all runs, but each run must be exactly the same in terms of aperture x/y coordinates for this to make sense. 

`$ python process_run_for_flat.py --d results/1/ --d results/2/ --d results/3/`

Once the flat is made, we can find the instrumental signature after removing the mean plane. This is done by using the `find_instrumental_signature.py` script. The only thing to note here is that we want our interpolation increment (flag --i) to match that of the scan increment.

`$ python find_instrumental_signature.py flat.dat --i 55,55 --p`

This outputs a 2D interpolate to the file "instresp.npy". Load this by moving the file to the calibration_files/ directory. 

# To generate a temperature correction file

We must generate a correction file to compensate for the change in temperature over the scan of the order ~3um/deg, but non-linear. This is done by allowing the sensor to read out without moving the stages and turning the air-con off (presuming the temperature has already stabilised at the air-con set temperature). This will enable a good range of readings.

**TODO: SCRIPT CMD?**

We then process this run using the script `process_run_for_temperature_correction.py`.

`python process_run_for_temperature_correction.py --f input_file --p`

The --p flag allows visual inspection of the fit. Don't use too high an order (--o) otherwise it may barrel if/when it needs to extrapolate. Ideally your `input_file` should cover the range of temperatures over which the scans will happen, in which case a higher order may give you more accuracy. A second order fit looks good to ~0.3 micron.

This outputs a correction file "t\_cor.npy". Load this by moving the file to the calibration_files/ directory. 








