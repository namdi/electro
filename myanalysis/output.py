"""
 If I have problems, run this code as python2.6
"""

#-----------------------------------------------------------
# import
#-----------------------------------------------------------
import numpy as np
import matplotlib.pylab as plt
import my_utils as mu

#-----------------------------------------------------------
# Function declarations
#-----------------------------------------------------------

#
# testing
#

# This function plots the position of the source particles
def plot_timestep(y, nsources, step):
  plt.plot(y[step, 0, 0:nsources, 0], y[step, 0, 0:nsources, 1], '*')
  plt.show()
  return
  
#-----------------------------------------------------------
# Run code
#-----------------------------------------------------------

# the input directory the directory for input data

input_dir = '/home/namdi/Documents/School/UNC/Parallel_Time/Data/test/t_0'

# the pfasst information for output

output_dir = '/home/namdi/Documents/School/UNC/Parallel_Time/Data/test'

# read parameters
fname = input_dir + str("/info.txt")
nsources, ntargets, distribution, dt, nsteps = mu.load_parameters(fname)

N = nsources + ntargets

# finest level is the highest level
level = 1

iteration = 5

# create data array
y 	= np.zeros( (nsteps+1, 2, N, 3) )

# load initial data
y0, charge_src, charge_targ, dipstr_src, dipvec_src   = mu.load_initial_data(input_dir)
y[0,:] = y0[:]

# load pfasst data
mu.load_pfasst_data(output_dir, y, level, iteration)

#-----------------------------------------------------------
# testing
#-----------------------------------------------------------

# magnitude of velocities for each particle 
vmag = np.zeros( (nsteps+1, N) )

for i in range(vmag.shape[0]):
  for j in range(N):
    vmag[i,j] = np.linalg.norm( y[i,1,j] )

#    
# for debugging    
#

uu =1
