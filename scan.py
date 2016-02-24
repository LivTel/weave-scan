import serial
import sys
import time
import argparse
import math
import datetime

import pylibftdi
import numpy as np
import pylab as plt

import pyAPT

def readEnvironment(attempts):
    t = -1
    h = -1
    if T_SERIAL is not None:
        T_SERIAL.flushInput()
        T_SERIAL.flushOutput()
        for i in range(0,attempts):
            try:
                line = T_SERIAL.readline()
                t = float(line.split('\t')[0].strip(' \t\r\n\0'))
                h = float(line.split('\t')[1].strip(' \t\r\n\0'))
                assert t < 30
                assert t > 10
                assert h > 0
                assert h < 100
                break
            except AssertionError:
                continue
            except ValueError:
                continue
            except OSError:
                continue
            except IndexError:
                continue
    return t, h

def readSensor(n_readings, reading_delay, max_read_wait):
  if S_SERIAL is not None:
      S_SERIAL.flushInput()
      S_SERIAL.flushOutput()   
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
        while S_SERIAL.inWaiting() > 0:
          lastByte = thisByte
          thisByte = ord(S_SERIAL.read(1))   # read one byte
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
            return None, None       # we failed
            
          if finish:
            return np.mean(readings), np.std(readings)
    
def begin(home=False):
  if home:
    print "homing stages.."
    print "x homed (" + str(homeStage(args.sx, args.wf)) + "mm)"
    print "y homed (" + str(homeStage(args.sy, args.wf)) + "mm)"
  print "moving to starting position.."
  print "x stage at: " + str(moveStageAbs(args.sx, args.w[0][0], args.wf)) + "mm"
  print "y stage at: " + str(moveStageAbs(args.sy, args.w[1][0], args.wf)) + "mm"
  print "ready"
  
def raster_scan():
  X_R = args.w[0][1]-args.w[0][0]
  Y_R = args.w[1][1]-args.w[1][0]  
  serial_x = args.sx
  serial_y = args.sy
  x_s = getStagePos(serial_x)
  for idx_y, y in enumerate(np.arange(0, Y_R+args.syi, args.syi)):
    if idx_y == 0:					# start the clock if first scan
      n_scans = 0
      t_start = time.time()
      TIME = []
      T = []
      H = []
      D = []
      DE = []
    else:
      moveStageRel(serial_y, args.syi, args.wf)		# don't want to increment this the first time!
    moveStageAbs(serial_x, x_s, args.wf)
    for x in np.arange(0, X_R+args.sxi, args.sxi):
      with open(args.f, "a") as f:
        this_time = time.time()
        mean, stdev = readSensor(args.n, args.rd, args.wr)
        t, h = readEnvironment(args.nt)
        f.write(str(this_time) + '\t' + str(getStagePos(serial_x)) + '\t' + str(getStagePos(serial_y)) + '\t' + str(mean) + '\t' + str(stdev) + '\t' + str(t) + '\t' + str(h) + '\n')
        TIME.append(this_time)
        T.append(t)
        H.append(h)
        D.append(mean)
        DE.append(stdev)
        if args.p:
          # temperature
          plt.plot([ti-t_start for ti in TIME], T, 'ko-')
          plt.gca().ticklabel_format(useOffset=False)
          plt.xlabel("Time since start (s)")
          plt.ylabel("Temperature (deg)")
          plt.savefig("T.png")
          plt.clf()
          # humidity
          plt.plot([ti-t_start  for ti in TIME], H, 'ko-')
          plt.gca().ticklabel_format(useOffset=False)
          plt.xlabel("Time since start (s)")
          plt.ylabel("Humidity (%)")
          plt.savefig("H.png")
          plt.clf()
          # distance
          plt.plot([ti-t_start  for ti in TIME], D, 'ko-')
          plt.gca().ticklabel_format(useOffset=False)
          plt.xlabel("Time since start (s)")
          plt.ylabel("Distance measured (micron)")
          plt.savefig("D.png")
          plt.clf()
      moveStageRel(serial_x, args.sxi, args.wf)

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
  parser.add_argument('--rd', help="delay in seconds between consecutive sensor readings (s)", action="store", type=float, default=0)
  parser.add_argument('--n', help="number of sensor readings to take", action="store", type=int, default=100)
  parser.add_argument('--wr', help="maximum period to wait for sensor reads (s)", action="store", type=float, default=10)
  parser.add_argument('--wf', help="time to wait after moving stages before taking the first reading (s)", action="store", type=float, default=0)
  parser.add_argument('--sx', help="serial number of x stage", action="store", type=int, default=45855915)
  parser.add_argument('--sy', help="serial number of y stage", action="store", type=int, default=45855916)
  parser.add_argument('--dt', help="serial environment device location", action="store", type=str, default="/dev/ttyUSB0")
  parser.add_argument('--ds', help="serial stage device location", action="store", type=str, default="/dev/ttyAMA0")
  parser.add_argument('--dtb', help="serial environment device baud", action="store", type=int, default=9600)
  parser.add_argument('--dsb', help="serial stage device baud", action="store", type=int, default=115200)
  parser.add_argument('--nt', help="number of attempts to obtain temperature", action="store", type=int, default=3)
  parser.add_argument('--p', help="make plot?", action="store_true")
  args = parser.parse_args()

  T_SERIAL = serial.Serial(
    port=args.dt,
    baudrate=args.dtb,
    parity=serial.PARITY_NONE,
    stopbits=serial.STOPBITS_ONE,
    bytesize=serial.EIGHTBITS
  )
  if T_SERIAL.isOpen():
    T_SERIAL.close()
  T_SERIAL.open()  

  S_SERIAL = serial.Serial(
    port=args.ds,
    baudrate=args.dsb,
    parity=serial.PARITY_NONE,
    stopbits=serial.STOPBITS_ONE,
    bytesize=serial.EIGHTBITS
  )
  if S_SERIAL.isOpen():
    S_SERIAL.close()
  S_SERIAL.open() 

  if args.r:
    TIME = []
    T = []
    H = []
    D = []
    DE = []
    st_time = time.time()
    while True:
      this_time = time.time()
      t, h = readEnvironment(args.nt)
      mean, stdev = readSensor(args.n, args.rd, args.wr)
      d = float(mean)
      de = float(stdev)
      TIME.append(this_time)
      T.append(t)
      H.append(h)
      D.append(d)
      DE.append(de)
      print str(this_time), '\t', str(d), '\t', str(de), '\t', str(t), '\t', str(h)
      with open(args.f, 'a') as f:
        f.write(str(this_time) + '\t' + '-' + '\t' + '-' + '\t' + str(d) + '\t' + str(de) + '\t' + str(t) + '\t' + str(h) + '\n')
      if args.p:
          # T
          plt.plot([ti-st_time for ti in TIME], T, 'ko-')
          plt.gca().ticklabel_format(useOffset=False)
          plt.xlabel("Time since start (s)")
          plt.ylabel("Temperature (deg)")
          plt.savefig("T.png")
          plt.clf()
          # H
          plt.plot([ti-st_time for ti in TIME], H, 'ko-')
          plt.gca().ticklabel_format(useOffset=False)
          plt.xlabel("Time since start (s)")
          plt.ylabel("Humidity (%)")
          plt.savefig("H.png")
          plt.clf()
          # D
          plt.plot([ti-st_time for ti in TIME], D, 'ko-')
          plt.gca().ticklabel_format(useOffset=False)
          plt.xlabel("Time since start (s)")
          plt.ylabel("Distance measured (micron)")
          plt.savefig("D.png")
          plt.clf()

  if args.m or args.s:
    if args.w is None:
      print("Window not defined.")
      exit(0)
    args.w = ([float(n) for n in args.w.split(',')])
    args.w = [(args.w[0], args.w[1]), (args.w[2], args.w[3])]

  if args.m: 
    begin()
  if args.s:
    raster_scan()


	

