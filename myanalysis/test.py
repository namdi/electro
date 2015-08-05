"""
Plot error-vs-time for Laplace equation equation to demonstrate that PFASST is
converging properly.

This assumes that you have already have made a reference solution and a non-reference solution

In KS example, make sure you run:
$ ./main.exe probin.ks.ref.nml
$ mpiexec -n 5 ./main.exe probin.ks.pf.nml

Make sure that the LIBPFASST directory is in your PYTHONPATH:

export PYTHONPATH=/path/to/libpfasst

position error: idx=0
velocity error: idx=1

for debugging /home/namdi/Documents/School/UNC/Parallel_Time/Code/fpfasst/libpfasst/pf/my_convergence.py
"""

#--------------------------------------
# Import
#--------------------------------------
import pylab
import pf
import my_utils as mu
import numpy as np
import pf.my_convergence as mc
import numpy.linalg as la
import matplotlib.pylab as plt

#---------------------------------------
# Run
#---------------------------------------
# the input directory the directory for input data

data_dir = '/home/namdi/Documents/School/UNC/Parallel_Time/Data'

input_dir = data_dir + "/Results/2014_02/19/sfmm_s256_a3/t_0"

# the pfasst information for output

output_dir = data_dir + "/Results/2014_02/19/sfmm_s256_a3"

# read parameters
fname = input_dir + str("/info.txt")
nsources, ntargets, distribution, dt, nsteps = mu.load_parameters(fname)

N = nsources + ntargets

# finest level is the highest level
level = 1
iteration = 1
niters = 4

# create data array
y 	= np.zeros( (nsteps+1, 2, N, 3) )

# load initial data
y0, charge_src, charge_targ, dipstr_src, dipvec_src   = mu.load_initial_data(input_dir)
y[0,:] = y0[:]

# load pfasst data
mu.load_pfasst_data(output_dir, y, level, iteration)

# load pfasst data
idx = 1 # for velocity (niters, nsteps, N, 3)
data = mc.get_data(output_dir, nsources, idx)

data2 = np.zeros((niters, nsteps+1, 2, N, 3))
mu.load_pfasst_data_all(output_dir, data2, level)

err = abs ( data  - data2[:,1:,1,:nsources] )

magerr = np.zeros((niters, nsteps, nsources))
maxmag = np.zeros((niters, nsteps))
p = np.zeros(maxmag.shape)

""""
for k in range(niters):
  for i in range(nsteps):  
    for j in range(nsources):
      magerr[k,i,j] = la.norm(err[k,i,j,:])
    maxmag[k,i] = max(magerr[k,i,:])
    p[k,i] = np.log(maxmag[k,i])/np.log(10)
"""
mom = np.zeros((niters, nsteps, 3))
magmom = np.zeros((niters, nsteps))

vel = data2[:,:,1,:nsources] #(niters, nsteps, nsources, 3)

for i in range(niters):
  for j in range(nsteps):  
    mom[i,j] = sum(vel[i,j])
    magmom[i,j] = la.norm(mom[i,j])        
    p[i,j] = np.log(magmom[i,j])/np.log(10)

print p[:, nsteps-1]  
plt.show()
