# Python Modules
import numpy as np
import struct
import math
import numpy
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
import sys
import getopt
import os
import time
from multiprocessing import Pool
from athena_read import athdf as athdf
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
    FileTransferSpeed, FormatLabel, Percentage, \
        ProgressBar, ReverseBar, RotatingMarker, \
            SimpleProgress, Timer
from blessings import Terminal
from termcolor import colored
from DiskPlanet import torque, angles, SimulationSetup

# Simulation setup (put frequently changed parameters here)
inc = 0.0                          # inclination of the planet
timestep = 0.062831853071795864769 # step between outputs

try:
  opts, args = getopt.getopt(sys.argv[1:], "t:a:n:p:",["torque","angles"])
    
except getopt.GetoptError:
  print "Unknown flag"
  sys.exit(2)
    
# Clear terminal window
term = Terminal()
sys.stderr.write("\x1b[2J\x1b[H")

class Writer():
  def __init__(self, location):
    self.location = location
  def write(self, string):
    with term.location(*self.location):
      print(string)

flagT = 0
flagA = 0

# Read arguments
for opt, arg in opts:
  if opt in ("-n"):
    total = int(arg)
  elif opt in ("-p"):
    proc = int(arg)
  elif opt in ("--torque"):
    flagT = 1
    if not os.path.exists("sigma"):
      os.makedirs("sigma")
    if not os.path.exists("Nforce"):
      os.makedirs("Nforce")
    if not os.path.exists("NforcePert"):
      os.makedirs("NforcePert")
    if not os.path.exists("Tforce"):
      os.makedirs("Tforce")
    if not os.path.exists("Rforce"):
      os.makedirs("Rforce")
    if not os.path.exists("RforcePert"):
      os.makedirs("RforcePert")
    if not os.path.exists("angular_momentum"):
      os.makedirs("angular_momentum")
  elif opt in ("--angles"):
    flagA = 1
    if not os.path.exists("inc_angle"):
      os.makedirs("inc_angle")
    if not os.path.exists("prec_angle"):
      os.makedirs("prec_angle")

arguments = np.arange(proc)

def worker(n):
  x_pos = 0
  y_pos = n+1
  location = (x_pos, y_pos)
  writer = Writer(location)
  widgets = [colored("".join(["Processor number ","%02d"%n,": "]),'red'), Percentage(),' ', Bar(marker="#"), ' ', ETA()]
  pbar = ProgressBar(widgets=widgets, maxval = (total+1)//proc,fd=writer)
  pbar.start()
  for i in range(0,(total+1)//proc):
    if n*((total+1)//proc) + i < total+1:
      T = n*((total+1)//proc) + i
      setup = SimulationSetup(GM = 1.0, GM1 = 1.0e-4, R0 = 1.0, Rsoft = 0.025, inclination = inc, time = timestep, number = T, rho_0 = 1.0, zeta = 2.0, rho_floor = 1.0e-9,hor=0.1)
      filename = ".".join(["temp","out1","%05d"%T,"athdf"])
      if flagT == 1:
        torque(filename, setup)
      if flagA == 1:
        angles(filename, setup)
    pbar.update(i+1)
  pbar.finish()

print "Calculating migration torque and surface density..."
print("\033[%d;0HStart..."%(proc+2))
pool = Pool(proc)
pool.map(worker,arguments.tolist())
print("\033[%d;6H...finish!"%(proc+2))
for i in range(proc+4):
	print "\n"
