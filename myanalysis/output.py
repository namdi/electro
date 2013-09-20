"""
 If I have problems, run this code as python2.6
"""

#-----------------------------
# import
#---------------------------
import numpy as np
import matplotlib.pylab as plt
import my_utils as mu

# testing

def plot_timestep(y, nsources, step):
  plt.plot(y[step, 0, 0:nsources, 0], y[step, 0, 0:nsources, 1], '*')
  plt.show()
  return
#-----------------------------
# Run code
#-----------------------------

# the input directory the directory for input data

input_dir = '/home/namdi/Documents/School/UNC/Parallel_Time/Data/test/t_0'

# the pfasst information for output

output_dir = '/home/namdi/Documents/School/UNC/Parallel_Time/Data/test'
#output_dir2 = '/home/namdi/Documents/School/UNC/Parallel_Time/Data/test/fmm'

# read parameters
fname = input_dir + str("/info.txt")
nsources, ntargets, distribution, dt, nsteps = mu.load_parameters(fname)

N = nsources + ntargets

# finest level is the highest level
level = 1

iteration = 4

# create data array
y 	= np.zeros( (nsteps+1, 2, N, 3) )

# load initial data
y0, charge_src, charge_targ, dipstr_src, dipvec_src   = mu.load_initial_data(input_dir)
y[0,:] = y0[:]

# load pfasst data
mu.load_pfasst_data(output_dir, y, level, iteration)

# testing
velx = np.abs(y[:, 1, 0, 0])
posx = np.abs(y[:, 0, 0, 0])

#print "velx: " 
#print velx

#print "posx: "
#print posx

plt.plot(velx)
plt.show()
uu =1
