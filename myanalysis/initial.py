"""
This file stores functions for different configurations for initial conditions
"""

#----------------------------
# Import
#---------------------------
import numpy as np
import random
import os
#import h5py

#----------------------------
# Functions
#---------------------------
# distribution:
#	1 -> random charges on the SURFACE of a sphere
#	2 -> random charges within a cube

#	test cases
#	-1 -> 4 particles to make a square, uniform charges
#	-2 ->uniform charges for a 2D circular ring
#
# this function sets the initial condition
#	

#
# set_initial_position(): sets the initial position 
#

def set_initial_position(distribution, pos, length):    
    
  if distribution == 1:    
    # random charge with random distribution on SURFACE of sphere 
    sphere(pos, length)
    
  elif distribution == 2:
    # random charge with random distribution within a cube
    cube(pos, length)
     
  return
  
#
# cube(): randomly assigns position and charges within a cube of side length "side" centered at the origin 
# pos must be a 2D array
#

def cube(pos, side):
  
  # distribute molecules randomly in cube
  # the cube is centered at zero
  N = pos.shape[0]
  
  #center the random numbers so they go from [-0.5, 0.5)
  for i in range(N):
      pos[i,0] = (random.random() - 0.5) * side
      pos[i,1] = (random.random() - 0.5) * side
      pos[i,2] = (random.random() - 0.5) * side
      
  return
  
#  
# sphere(): randomly assigns position on the SURFACE of a sphere of radius "radius" centered at the origin 
# pos must be a 2D array
#

def sphere(pos, radius):
  
  # distribute particles randomly in a sphere
  # sphere is centered at zero
  
  N = pos.shape[0]    
  
  for i in range(N):
      theta = random.random() * np.pi
      phi = random.random() * 2 * np.pi
      pos[i, 0] = radius * np.sin(theta) * np.cos(phi)
      pos[i, 1] = radius * np.sin(theta) * np.sin(phi)
      pos[i, 2] = radius * np.cos(theta)
      
  return
  
#
# square_test(): set the initial condtion (position & velocity) for the sourcs. The velocity is zero
# position (sources)
#

def square_test(pos_src, pos_targ, length):
  
  xmax = 0.5
  ymax = 0.5
  zmax = 1.0

  # set the position of sources
  pos_src[0,:] = [-xmax, -ymax, 0.0]
  pos_src[1,:] = [xmax, -ymax, 0.0]
  pos_src[2,:] = [-xmax, ymax, 0.0]
  pos_src[3,:] = [xmax, ymax, 0.0]
  
  # set the initial position for the targets
  pos_targ[0,:] = [0.0, 0.0, -zmax]
  pos_targ[1,:] = [0.0, 0.0, 0.0]
  pos_targ[2,:] = [0.0, 0.0, zmax]
    
  pos_src 	= pos_src * length
  pos_targ	= pos_targ * length  

  return

#
# circle_test(): set the initial condtion position for the sources.
#
  
def circle_test(pos_src, pos_targ, length):
  
  nsources = pos_src.shape[0]
  
  # radius
  R = length
  
  # max z coordinate
  zmax = 1.0
  
  for i in range(nsources):
    pos_src[i,0] = R * np.cos(2*np.pi/nsources * i)
    pos_src[i,1] = R * np.sin(2*np.pi/nsources * i)
    pos_src[i,2] = 0.0
    
  # set the initial position for the targets
  pos_targ[0,:] = [0.0, 0.0, -zmax]
  pos_targ[1,:] = [0.0, 0.0, 0.0]
  pos_targ[2,:] = [0.0, 0.0, zmax]
         
  return
  
#
# set_charge(): sets the charge to be random numbers in the interval [-max_charge, max_charge]
#

def set_charge(charge, max_charge):
  
  N = len(charge)
  
  for i in range(N):
    charge[i] = (random.random() - 0.5) # keep real charges
    #charge[i] = (random.random() - 0.5) + (random.random() - 0.5)*1j
  
  charge = charge * max_charge
  return

#
# set_dipole(): sets the dipole 
#

def set_dipole(dipstr, dipvec, max_strength):
  nsources = len(dipstr)
  
  for i in range(nsources):
    #dipstr[i] = (random.random() - 0.5) + (random.random() - 0.5)*1j
    dipstr[i] = (random.random() - 0.5)
    dipstr[i] = dipstr[i] * max_strength
    
    # dipole direction is normalized to 1
    x = random.random() - 0.5
    y = random.random() - 0.5
    z = random.random() - 0.5
    
    # magnitude
    mag = np.sqrt(x**2 + y**2 + z**2)
    
    dipvec[i,0] = x/mag
    dipvec[i,1] = y/mag
    dipvec[i,2] = z/mag
  
  return
  
""" this function creates the timesteps depending on the pfasst level
Given a number of stepsizes of size dt, and number of lobatto nodes in each step,
the appropriate lobatto nodes will be made in range [0, dt*nsteps]
"""
def create_time(dt, nsteps, lnodes):
  
  # the lobatto nodes from [0, 1] for 5 or 3 nodes
  if lnodes == 3:
    lobatto = ( 1 + np.array((-1.0, 0.0, 1.0)) )/2
    
  elif lnodes == 5:
    lobatto = (1 + np.array((-1.0, -np.sqrt(3./7), 0, np.sqrt(3./7), 1.0)) )/2
  
  # the number of nodes
  nnodes = 1 + (len(lobatto) - 1) * nsteps
  
  # the temporal nodes
  t = np.zeros((nnodes, 1))
  
  # the timesteps between nodes
  dt_lobatto = np.zeros( len(lobatto) -1)
  
  # calculate the timestep difference in the lobatto nodes
  for i in range( len(dt_lobatto) ):
      dt_lobatto[i] = dt*(lobatto[i+1] - lobatto[i])
      
  for i in range(1, nnodes):
    j = np.mod(i-1, len(dt_lobatto))
    t[i] = t[i-1] + dt_lobatto[j]
    
  return t    
  