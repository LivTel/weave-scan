#/usr/bin/env python

import serial
import sys
import time

if (len(sys.argv) == 4):
  port = sys.argv[1]
  baud = sys.argv[2]
  filename = sys.argv[3]
else:
  port = "/dev/ttyACM0"
  baud = 9600
  filename = "T.dat"

ser = serial.Serial(port, baud)
f = open(filename, 'w', 0) # require 0 so that output is fed directly to file, otherwise file doesn't get written to

for x in range(0, 3):
  ser.readline()	   # flush some crap from input

while True:
  timestamp = str(int(time.time()))	# 3600 for DST
  line = ser.readline()
  #print timestamp + "\t" + line
  f.write(timestamp + "\t" + line)


