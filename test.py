import serial
import sys
import time

import pylibftdi
import numpy as np

import pyAPT

READING_DELAY   = 0.05                              # delay in seconds between consecutive sensor readings
N_READINGS      = 30                                # number of sensor readings to take
MAX_READ_WAIT   = (READING_DELAY*N_READINGS) + 15   # maximum period to wait for sensor reads
SETTLE_TIME     = 2                                 # time to wait after moving stages before taking the first reading
SERIALS         = [45855916, 45855915]              # serial numbers of stages (x, y)

SCAN_WINDOW     = [(15, 160), (79, 299)]            # [(x_lo, x_hi), (y_lo, y_hi)] scan window
## RASTER SCAN
RASTER_X_INTERVAL   = 10
RASTER_Y_INTERVAL   = 10

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
  print "x stage at: " + str(moveStageAbs(SERIALS[0], SCAN_WINDOW[0][0], SETTLE_TIME)) + "mm"
  print "y stage at: " + str(moveStageAbs(SERIALS[1], SCAN_WINDOW[1][0], SETTLE_TIME)) + "mm"
  print "ready"
  
def raster_scan():
  X_R = SCAN_WINDOW[0][1] - SCAN_WINDOW[0][0]
  Y_R = SCAN_WINDOW[1][1] - SCAN_WINDOW[1][0]  
  serial_x = SERIALS[0]
  serial_y = SERIALS[1]
  x_s = getStagePos(serial_x)
  for y in range(0, Y_R+RASTER_Y_INTERVAL, RASTER_Y_INTERVAL):
    moveStageRel(serial_y, RASTER_Y_INTERVAL, SETTLE_TIME)
    moveStageAbs(serial_x, x_s, SETTLE_TIME)
    for x in range(0, X_R+RASTER_X_INTERVAL, RASTER_X_INTERVAL):
      with open("results", "a") as f:
        mean, stdev = readSensor(N_READINGS, READING_DELAY, MAX_READ_WAIT)
        f.write(str(getStagePos(serial_x)) + '\t' + str(getStagePos(serial_y)) + '\t' + str(mean) + '\t' + str(stdev) + '\n')
      moveStageRel(serial_x, RASTER_X_INTERVAL, SETTLE_TIME)
   
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
  if sys.argv[1] == 'r':
    print readSensor(N_READINGS, READING_DELAY, MAX_READ_WAIT)
  elif sys.argv[1] == 'm': 
    begin()
  elif sys.argv[1] == 's':
    raster_scan()
  elif sys.argv[1] == 'b':
    begin()
    raster_scan()
	

