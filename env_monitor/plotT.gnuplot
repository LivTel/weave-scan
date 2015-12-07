#!/usr/bin/gnuplot

set terminal png
set output "/var/www/html/T_full.png"
set xlabel "Time"
set ylabel "Temperature (K)"
set yr [280:290]
set xdata time
set timefmt "%s" 
set format x "%H:%M %m/%d"
set xtic rotate by -90
set offsets 0, 0, 0, 0 
plot "T.dat" using 1:5 title 'Granite' with lines lt 1 lc 1, \
"T.dat" using 1:6 title 'Bench' with lines lt 1 lc 2, \
"T.dat" using 1:7 title 'Plate' with lines lt 1 lc 3

set terminal png
set output "/var/www/html/T.png"
set xlabel "Time"
set ylabel "Temperature (K)"
set autoscale
set xdata time
set timefmt "%s"
set format x "%H:%M"
set xtic rotate by -90
set offsets 0, 0, 10, 10
plot "<(tail -600 T.dat)" using 1:5 title 'Granite' with lines lt 1 lc 1, \
"<(tail -1200 T.dat)" using 1:6 title 'Bench' with lines lt 1 lc 2, \
"<(tail -1200 T.dat)" using 1:7 title 'Plate' with lines lt 1 lc 3

set output "/var/www/html/V_full.png"
set xlabel "Time"
set ylabel "Forward voltage drop (volts)"
set yr [0:1]
set xdata time
set timefmt "%s"
set format x "%H:%M %m/%d"
set xtic rotate by -90
set offsets 0, 0, 0, 0
plot "T.dat" using 1:2 title 'Granite' with lines lt 1 lc 1, \
"T.dat" using 1:3 title 'Bench' with lines lt 1 lc 2, \
"T.dat" using 1:4 title 'Plate' with lines lt 1 lc 3

set output "/var/www/html/V.png"
set xlabel "Time"
set ylabel "Forward voltage drop (volts)"
set autoscale
set xdata  time
set timefmt "%s"
set format x "%H:%M"
set xtic rotate by -90
set offsets 0, 0, 0.02, 0.02
plot "<(tail -600 T.dat)" using 1:2 title 'Granite' with lines lt 1 lc 1, \
"<(tail -1200 T.dat)" using 1:3 title 'Bench' with lines lt 1 lc 2, \
"<(tail -1200 T.dat)" using 1:4 title 'Plate' with lines lt 1 lc 3

set terminal png
set output "/var/www/html/S_full.png"
set xlabel "Time"
set ylabel "Sensor Reading (um)"
set xdata time
set timefmt "%s" 
set format x "%H:%M %m/%d"
set xtic rotate by -90
set offsets 0, 0, 10, 10 
plot "S.dat" using 1:2:3 with yerrorbars

set terminal png
set output "/var/www/html/S.png"
set xlabel "Time"
set ylabel "Sensor Reading (um)"
set autoscale
set xdata time
set timefmt "%s" 
set format x "%H:%M %m/%d"
set xtic rotate by -90
set offsets 0, 0, 0, 0 
plot "<(tail -1200 S.dat)" using 1:2

pause 10; refresh; reread
