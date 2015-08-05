"""
 If I have problems, run this code as python2.6
 This code plots the initial (source) position 
"""

#-----------------------------------------------------------
# import
#-----------------------------------------------------------
import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import my_utils as mu

#-----------------------------------------------------------
# Function declarations
#-----------------------------------------------------------

def plot3D_timestep(y, nsources, step):
  fig = plt.figure()
  ax = Axes3D(fig)
  ax = fig.add_subplot(111, projection='3d')
  ax.scatter(y[step, 0, 0:nsources, 0], y[step, 0, 0:nsources, 1], y[step, 0, 0:nsources, 2])
  
  plt.show()
  return
  
#-----------------------------------------------------------
# Run code
#-----------------------------------------------------------

# the input directory the directory for input data

input_dir = '/home/namdi/Documents/School/UNC/Parallel_Time/Data/test_fmm/t_0'

# read parameters
fname = input_dir + str("/info.txt")
nsources, ntargets, distribution, dt, nsteps = mu.load_parameters(fname)

N = nsources + ntargets

# create data array
y 	= np.zeros( (nsteps+1, 2, N, 3) )

# load initial data
y0, charge_src, charge_targ, dipstr_src, dipvec_src   = mu.load_initial_data(input_dir)
y[0,:] = y0[:]

# plot the data
plot3D_timestep(y, nsources, 0)

uu =1
