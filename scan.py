import serial
import sys
import time
import argparse
import math
import datetime

import pylibftdi
import numpy as np

import pyAPT

READING_DELAY   = 0.0                               # delay in seconds between consecutive sensor readings
N_READINGS      = 100                               # number of sensor readings to take
MAX_READ_WAIT   = (READING_DELAY*N_READINGS) + 15   # maximum period to wait for sensor reads
SETTLE_TIME     = 0.0                               # time to wait after moving stages before taking the first reading
SERIALS         = [45855915, 45855916]              # serial numbers of stages (x, y)

def readSensor(n_readings, reading_delay, max_read_wait):
  ser = serial.Serial(
    port='/dev/ttyAMA0',
    baudrate=115200,
    parity=serial.PARITY_NONE,
    stopbits=serial.STOPBITS_ONE,
    bytesize=serial.EIGHTBITS
  )
  if ser.isOpen():
    ser.close()
  ser.open()    
  lastByte = None
  thisByte = None
  
  # these are byte order bits which tell us which byte is next
  nextIsMSB1 = False
  nextIsMSB2 = False
  nextIsLSB1 = False
  nextIsLSB2 = False
  
  # these are the values for the bytes themselves, which will be 
  # fully populated when both thisByte == lastByte == 255
  MSB1 = None
  MSB2 = None
  LSB1 = None
  LSB2 = None

  finish     = False    # bit to tell us when to finish
  c_readings = 0        # counter for number of readings
  readings   = []
  
  begin           = time.time()
  nextReadingTime = begin + reading_delay
  maxReadingTime  = begin + max_read_wait
  while True:
    while ser.inWaiting() > 0:
      lastByte = thisByte
      thisByte = ord(ser.read(1))   # read one byte
      if nextIsMSB1:
        nextIsMSB1 = False
        nextIsMSB2 = True
        MSB1 = thisByte
      elif nextIsMSB2:
        nextIsMSB2 = False
        nextIsLSB1 = True
        MSB2 = thisByte
      elif nextIsLSB1:
        nextIsLSB1 = False
        nextIsLSB2 = True
        LSB1 = thisByte
      elif nextIsLSB2:
	nextIsLSB2 = False
	LSB2 = thisByte
	MSB = (MSB1*256) + MSB2
	LSB = (LSB1*256) + LSB2
  
        now = time.time()
        if now > nextReadingTime:
          val = ((MSB*pow(2,15) + LSB) * (12000./pow(2,30)))
          readings.append(val)
          c_readings = c_readings + 1
          nextReadingTime = now + reading_delay

      ## is this a new reading? if so, reset our byte order bits.
      if thisByte == 255 and lastByte == 255:
        nextIsMSB1 = True
        nextIsMSB2 = False
        nextIsLSB1 = False
        nextIsLSB2 = False

      ## have we taken enough readings yet?
      if c_readings >= n_readings:
        finish = True
       
      ## or have we ran out of time? 
      now = time.time() 
      if max_read_wait > now:
        ser.close()
        return None, None       # we failed
            
      if finish:
        ser.close()
        return np.mean(readings), np.std(readings)
    
def begin(home=False):
  if home:
    print "homing stages.."
    print "x homed (" + str(homeStage(SERIALS[0], SETTLE_TIME)) + "mm)"
    print "y homed (" + str(homeStage(SERIALS[1], SETTLE_TIME)) + "mm)"
  print "moving to starting position.."
  print "x stage at: " + str(moveStageAbs(SERIALS[0], args.w[0][0], SETTLE_TIME)) + "mm"
  print "y stage at: " + str(moveStageAbs(SERIALS[1], args.w[1][0], SETTLE_TIME)) + "mm"
  print "ready"
  
def raster_scan(outfile):
  X_R = args.w[0][1]-args.w[0][0]
  Y_R = args.w[1][1]-args.w[1][0]  
  serial_x = SERIALS[0]
  serial_y = SERIALS[1]
  x_s = getStagePos(serial_x)
  for idx_y, y in enumerate(np.arange(0, Y_R+args.syi, args.syi)):
    if idx_y == 0:					# start the clock if first scan
      n_scans = 0
      t_start = time.time()
    else:
      moveStageRel(serial_y, args.syi, SETTLE_TIME)	# don't want to increment this the first time!
    moveStageAbs(serial_x, x_s, SETTLE_TIME)
    for x in np.arange(0, X_R+args.sxi, args.sxi):
      with open(outfile, "a") as f:
        mean, stdev = readSensor(N_READINGS, READING_DELAY, MAX_READ_WAIT)
        f.write(str(time.time()) + '\t' + str(getStagePos(serial_x)) + '\t' + str(getStagePos(serial_y)) + '\t' + str(mean) + '\t' + str(stdev) + '\n')
      moveStageRel(serial_x, args.sxi, SETTLE_TIME)

      n_scans = n_scans+1
      t_elapsed = time.time()-t_start
      rate = n_scans/t_elapsed
      n_scans_tot = len(np.arange(0, Y_R+args.syi, args.syi))*len(np.arange(0, X_R+args.sxi, args.sxi))
      n_scans_left = n_scans_tot-n_scans
      print
      print("scan " + str(n_scans) + " of " + str(int(n_scans_tot)))
      print("-> approximate time elapsed " + str(round(t_elapsed,1)) + "s")
      print("-> approximate rate of " + str(round(rate,2)) + " scans/s")
      print("-> approximate time left " + str(round(n_scans_left/rate,1)) + "s")
      print("-> approximate time of completion " + str(datetime.datetime.fromtimestamp(int(time.time()+(n_scans_left/rate))).strftime('%Y-%m-%d %H:%M:%S')))
      print 
   
def moveStageRel(stage_serial, dist, settle_time):
  with pyAPT.MTS50(serial_number=stage_serial) as con:
    con.move(dist)
    time.sleep(settle_time)
    return con.position()

def moveStageAbs(stage_serial, pos, settle_time):
  with pyAPT.MTS50(serial_number=stage_serial) as con:
    con.goto(pos)
    time.sleep(settle_time)
    return con.position()

def getStagePos(stage_serial):
  with pyAPT.MTS50(serial_number=stage_serial) as con:
    return con.position()    

def homeStage(stage_serial, settle_time):
  with pyAPT.MTS50(serial_number=stage_serial) as con:
    con.home()  
    time.sleep(settle_time)
    return con.position()
      
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('f', help="path to output file", action="store", type=str, default="results")
  parser.add_argument('--r', action="store_true", help="read sensor")
  parser.add_argument('--m', action="store_true", help="move stages to starting position")
  parser.add_argument('--s', action="store_true", help="start scan")
  parser.add_argument('--w', help="window to scan over (x1,x2,y1,y2) (mm)", action="store", type=str)
  parser.add_argument('--sxi', help="scan x interval (mm)", action="store", type=float, default=0.01)
  parser.add_argument('--syi', help="scan y interval (mm)", action="store", type=float, default=0.01)
  args = parser.parse_args()

  if args.r:
    while True:
      res = readSensor(N_READINGS, READING_DELAY, MAX_READ_WAIT)
      print str(time.time()), '\t', str(res[0]), '\t', str(res[1])
      with open(args.f, 'a') as f:
        f.write(str(time.time()) + '\t' + '-' + '\t' + '-' + '\t' + str(res[0]) + '\t' + str(res[1]) + '\n')

  if args.m or args.s:
    if args.w is None:
      print("Window not defined.")
      exit(0)
    args.w = ([float(n) for n in args.w.split(',')])
    args.w = [(args.w[0], args.w[1]), (args.w[2], args.w[3])]

  if args.m: 
    begin()
  if args.s:
    raster_scan(args.f)


	

