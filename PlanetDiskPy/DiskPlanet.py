"""
Package for working with the simulation output for disk-planet.cpp problem generator.

Contains the following classses:
SimulationSetup(
    GM          - mass of the central object
    GM1         - mass of the planet
    R0          - distance to the planet
    Rsoft       - softening length
    inclination - planet inclination
    time        - time between snapshots
    number      - snapshot number
    rho_0       - density normalization
    zeta        - density profile slope
    rho_floor   - density floor
)

Contains the following functions:
MigrationTorque()
Sigma()
R_force()
N_force()
T_force()
inc_angle()
torque() -- function which calculates three components of the force acting on a planet as well as angular momentum flux through the disk.
angles() -- function which calculates inclination and precession angles of the disk
"""


# Python Modules
import numpy as np
import struct
import math
import numpy
import sys
import os
import time
from multiprocessing import Pool
from athena_read import athdf as athdf

class Writer():
  def __init__(self, location):
    self.location = location
    def write(self, string):
      with term.location(*self.location):
        print(string)

class SimulationSetup(object):
  def __init__(self, GM=1, GM1=1.0e-4, R0=1.0, Rsoft = 0.025, inclination=0.0, time=0.0, number = 1, rho_0=1.0, zeta=2.0, rho_floor=1.0e-9, hor = 0.1):
    self.GM = GM
    self.GM1 = GM1
    self.R0 = R0
    self.Rsoft = Rsoft
    self.inclination = inclination
    self.time = time
    self.number = number
    self.rho_0 = rho_0
    self.zeta = zeta
    self.rho_floor = rho_floor
    self.hor = hor


# Function for creating output names
def out_name(filename, string):
  parts = filename.split('.')
  p_length = len(parts)
  n_parts = parts
  n_parts[p_length-1] = string
  output = ".".join(n_parts)
  return output

##########################################################################

def torque(filename, setup):
  
  # 1. Read data from .athdf file

  array = athdf("athdf/" + filename)
  x1=array['x1f']
  x2=array['x2f']
  x3=array['x3f']
  data=array['dens']
  n1=len(x1)
  n2=len(x2)
  n3=len(x3)

  # 2. Transform to cell-centered coordinates

  y1=np.zeros((n1-1))
  y2=np.zeros((n2-1))
  y3=np.zeros((n3-1))
  d1=np.zeros((n1-1))
  d2=np.zeros((n2-1))
  d3=np.zeros((n3-1))

  for ir in range(0,n1-1):
    y1[ir] = 0.75*(x1[ir+1]+math.pow(x1[ir],3)/(math.pow(x1[ir],2) + x1[ir]*x1[ir+1] + math.pow(x1[ir+1],2)))
    d1[ir] = x1[ir+1] - x1[ir]
  for ith in range(0,n2-1):
    y2[ith] = (x2[ith]*np.cos(x2[ith]) - x2[ith+1]*np.cos(x2[ith+1]) - np.sin(x2[ith]) + np.sin(x2[ith+1]))/(np.cos(x2[ith])-np.cos(x2[ith+1]))
    d2[ith] = x2[ith+1] - x2[ith]
  for iphi in range(0,n3-1):
    y3[iphi] = 0.5*(x3[iphi]+x3[iphi+1])
    d3[iphi] = x3[iphi+1] - x3[iphi]

  # 3. Reshape

  a1 = y1.reshape(1, 1, n1-1)
  a2 = y2.reshape(1, n2-1, 1)
  a3 = y3.reshape(n3-1, 1, 1)
  b1 = d1.reshape(1, 1, n1-1)
  b2 = d2.reshape(1, n2-1, 1)
  b3 = d3.reshape(n3-1, 1, 1)

  # 4. Simulation setup

  GM0 = setup.GM
  GM1 = setup.GM1
  R0 = setup.R0
  pl_inc = setup.inclination
  Rsoft = setup.Rsoft
  time = setup.number*setup.time
  phase = math.pow((GM0+GM1)/R0/R0/R0,0.5)*time
  h = setup.hor
  zeta = setup.zeta
  rho0 = setup.rho_0

  # 5. Calculate surface density

  result = data*np.sin(a2)*a1*a1*b2*b3
  sigma = result.sum(axis = 0).sum(axis = 0)

  output = out_name(filename,"sigma.dat")
  out = open("sigma/"+output, 'w+')
  for ir in range(0,n1-1):
    out.write(str(y1[ir]))
    out.write("\t")
    out.write(str(sigma[ir]))
    out.write("\t")
    out.write(str(d1[ir]))
    out.write("\n")
  out.close()

  del result, sigma

  c1 = y1.reshape(1, 1, n1-1)
  c2 = y2.reshape(1, n2-1, 1)
  c3 = y3.reshape(n3-1, 1, 1)
  e1 = (c1*np.sin(c2)*np.cos(c3) - np.ones((n3-1,n2-1,n1-1))*R0*np.cos(pl_inc)*np.cos(phase))
  e2 = (c1*np.sin(c2)*np.sin(c3) - np.ones((n3-1,n2-1,n1-1))*R0*np.sin(phase))
  e3 = (c1*np.cos(c2) - np.ones((n3-1,n2-1,n1-1))*R0*np.sin(pl_inc)*np.cos(phase))
  dist2 = e1*e1 + e2*e2 + e3*e3
  data0 = rho0*np.power(a1*np.sin(a2), -zeta)*np.exp((1.0 - 1.0/np.sin(a2))/h/h/a1)

  # 6. Calculate T-force

  arrayT = (dist2*dist2 + 3.5*dist2*Rsoft*Rsoft+np.ones((n3-1,n2-1,n1-1))*4.375*Rsoft*Rsoft*Rsoft*Rsoft)/(dist2+np.ones((n3-1,n2-1,n1-1))*Rsoft*Rsoft)/(dist2+np.ones((n3-1,n2-1,n1-1))*Rsoft*Rsoft)/(dist2+np.ones((n3-1,n2-1,n1-1))*Rsoft*Rsoft)/np.sqrt(dist2+np.ones((n3-1,n2-1,n1-1))*Rsoft*Rsoft)*(-e1*np.cos(pl_inc)*np.sin(phase) + e2*np.cos(phase) - e3*np.sin(pl_inc)*np.sin(phase))
  resultT = arrayT*data*np.sin(a2)*a1*a1*b2*b3
  forceT = resultT.sum(axis = 0).sum(axis = 0)

  output = out_name(filename,"Tforce.dat")
  out = open("Tforce/"+output, 'w+')
  for ir in range(0,n1-1):
    out.write(str(y1[ir]))
    out.write("\t")
    out.write(str(forceT[ir]))
    out.write("\t")
    out.write(str(d1[ir]))
    out.write("\n")
  out.close()

  del forceT, resultT

  # 7. Calculate R-force

  arrayR = (dist2*dist2 + 3.5*dist2*Rsoft*Rsoft+np.ones((n3-1,n2-1,n1-1))*4.375*Rsoft*Rsoft*Rsoft*Rsoft)/(dist2+np.ones((n3-1,n2-1,n1-1))*Rsoft*Rsoft)/(dist2+np.ones((n3-1,n2-1,n1-1))*Rsoft*Rsoft)/(dist2+np.ones((n3-1,n2-1,n1-1))*Rsoft*Rsoft)/np.sqrt(dist2+np.ones((n3-1,n2-1,n1-1))*Rsoft*Rsoft)*(e1*np.cos(pl_inc)*np.cos(phase) + e2*np.sin(phase) + e3*np.sin(pl_inc)*np.cos(phase))
  resultR = arrayR*data*np.sin(a2)*a1*a1*b2*b3
  forceR = resultR.sum(axis = 0).sum(axis = 0)

  output = out_name(filename,"Rforce.dat")
  out = open("Rforce/"+output, 'w+')
  for ir in range(0,n1-1):
    out.write(str(y1[ir]))
    out.write("\t")
    out.write(str(forceR[ir]))
    out.write("\t")
    out.write(str(d1[ir]))
    out.write("\n")
  out.close()
  
  del forceR, resultR

  resultRpert = arrayR*(data - data0)*np.sin(a2)*a1*a1*b2*b3
  forceRpert = resultRpert.sum(axis = 0).sum(axis = 0)

  output = out_name(filename,"RforcePert.dat")
  out = open("RforcePert/"+output, 'w+')
  for ir in range(0,n1-1):
    out.write(str(y1[ir]))
    out.write("\t")
    out.write(str(forceRpert[ir]))
    out.write("\t")
    out.write(str(d1[ir]))
    out.write("\n")
  out.close()
  
  del forceRpert, resultRpert

  # 8. Calculate N-force

  arrayN = (dist2*dist2 + 3.5*dist2*Rsoft*Rsoft+np.ones((n3-1,n2-1,n1-1))*4.375*Rsoft*Rsoft*Rsoft*Rsoft)/(dist2+np.ones((n3-1,n2-1,n1-1))*Rsoft*Rsoft)/(dist2+np.ones((n3-1,n2-1,n1-1))*Rsoft*Rsoft)/(dist2+np.ones((n3-1,n2-1,n1-1))*Rsoft*Rsoft)/np.sqrt(dist2+np.ones((n3-1,n2-1,n1-1))*Rsoft*Rsoft)*(-e1*np.sin(pl_inc) + e2*0.0 + e3*np.cos(pl_inc))
  resultN = arrayN*(data)*np.sin(a2)*a1*a1*b2*b3
  forceN = resultN.sum(axis = 0).sum(axis = 0)

  output = out_name(filename,"Nforce.dat")
  out = open("Nforce/"+output, 'w+')
  for ir in range(0,n1-1):
    out.write(str(y1[ir]))
    out.write("\t")
    out.write(str(forceN[ir]))
    out.write("\t")
    out.write(str(d1[ir]))
    out.write("\n")
  out.close()
  
  del forceN, resultN

  resultNpert = arrayN*(data - data0)*np.sin(a2)*a1*a1*b2*b3
  forceNpert = resultNpert.sum(axis = 0).sum(axis = 0)

  output = out_name(filename,"NforcePert.dat")
  out = open("NforcePert/"+output, 'w+')
  for ir in range(0,n1-1):
    out.write(str(y1[ir]))
    out.write("\t")
    out.write(str(forceNpert[ir]))
    out.write("\t")
    out.write(str(d1[ir]))
    out.write("\n")
  out.close()
  
  del forceNpert, resultNpert
"""
  # 9. Calculate angular momentum flux

  mom1=array['mom1']
  mom2=array['mom2']
  mom3=array['mom3']

  Mdot   = mom1 * np.sin(a2)*a1*a1*b2*b3
  Fj     = (data - data0)*(data - data0) * np.sin(a2)*a1*a1*b2*b3
  stress = a1 * mom1 * (mom3/data - np.sqrt(GM0/a1/np.sin(a2)) * (1.0 - zeta * h*h*a1/np.sin(a2)/np.sin(a2))) * np.sin(a2)*a1*a1*b2*b3

  MdotF = Mdot.sum(axis = 0).sum(axis = 0)
  FjF = Fj.sum(axis = 0).sum(axis = 0)
  stressF = stress.sum(axis = 0).sum(axis = 0)


  output = out_name(filename,"AM.dat")
  out = open("".join(["angular_momentum/",output]), 'w+')
  for ir in range(0,n1-1):
    out.write(str(y1[ir]))
    out.write("\t")
    out.write(str(MdotF[ir]))
    out.write("\t")
    out.write(str(FjF[ir]))
    out.write("\t")
    out.write(str(stressF[ir]))
    out.write("\t")
    out.write(str(d1[ir]))
    out.write("\n")
  out.close()
"""
##########################################################################

def angles(filename):
  # Read the file
  array = athdf("athdf/"+filename)
  x1=array['x1f']
  x2=array['x2f']
  x3=array['x3f']
  data=array['dens']
  v1=array['mom1']
  v2=array['mom2']
  v3=array['mom3']
  n1=len(x1)
  n2=len(x2)
  n3=len(x3)
  output1 = out_name(filename,"inc_angle.dat")
  output2 = out_name(filename,"prec_angle.dat")
  y1=np.zeros((n1-1))
  y2=np.zeros((n2-1))
  y3=np.zeros((n3-1))
  d1=np.zeros((n1-1))
  d2=np.zeros((n2-1))
  d3=np.zeros((n3-1))
  # Thansform to cell-centered coordinates
  for ir in range(0,n1-1):
    y1[ir] = 0.75*(x1[ir+1]+math.pow(x1[ir],3)/(math.pow(x1[ir],2) + x1[ir]*x1[ir+1] + math.pow(x1[ir+1],2)))
    d1[ir] = x1[ir+1] - x1[ir]
  for ith in range(0,n2-1):
    y2[ith] = (x2[ith]*np.cos(x2[ith]) - x2[ith+1]*np.cos(x2[ith+1]) - np.sin(x2[ith]) + np.sin(x2[ith+1]))/(np.cos(x2[ith])-np.cos(x2[ith+1]))
    d2[ith] = x2[ith+1] - x2[ith]
  for iphi in range(0,n3-1):
    y3[iphi] = 0.5*(x3[iphi]+x3[iphi+1])
    d3[iphi] = x3[iphi+1] - x3[iphi]
    
  # Reshape
  a1 = y1.reshape(1, 1, n1-1)
  a2 = y2.reshape(1, n2-1, 1)
  a3 = y3.reshape(n3-1, 1, 1)
  b1 = d1.reshape(1, 1, n1-1)
  b2 = d2.reshape(1, n2-1, 1)
  b3 = d3.reshape(n3-1, 1, 1)
    
  # Calculate inclination angle
  Lz = np.sin(a2)*np.sin(a2)*v3
  Ly = np.cos(a3)*np.sin(a2)*v2 - np.cos(a2)*np.sin(a2)*np.sin(a3)*v3
  Lx = -np.sin(a3)*np.sin(a2)*v2 - np.cos(a2)*np.cos(a3)*np.sin(a2)*v3
  resLz = Lz.sum(axis = 0).sum(axis = 0)
  resLy = Ly.sum(axis = 0).sum(axis = 0)
  resLx = Lx.sum(axis = 0).sum(axis = 0)
  L = np.sqrt(resLx*resLx + resLy*resLy + resLz*resLz)
  inc = np.arccos(resLz/L)
  prec = np.arccos(-resLy/np.sqrt(resLy*resLy + resLx*resLx))
    
  # Write output
  out1 = open("inc_angle/"+output1, 'w+')
  for ir in range(0,n1-1):
    out1.write(str(y1[ir]))
    out1.write("\t")
    out1.write(str(inc[ir]))
    out1.write("\t")
    out1.write(str(d1[ir]))
    out1.write("\n")
  out1.close()
  out2 = open("prec_angle/"+output2, 'w+')
  for ir in range(0,n1-1):
    out2.write(str(y1[ir]))
    out2.write("\t")
    out2.write(str(prec[ir]))
    out2.write("\t")
    out2.write(str(d1[ir]))
    out2.write("\n")
  out2.close()

