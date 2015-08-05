"""
This comapres the PFASST results from different runs

This code plots the CHANGE in total mommentum for each timestep per titeration. 

position error: idx=0
velocity error: idx=1

"""

#--------------------------------------
# Import
#--------------------------------------
import pylab
import pf
import my_utils as mu
import pf.my_convergence as mc
import numpy as np
import numpy.linalg as la
import matplotlib.pylab as plt
  
#---------------------------------------
# Function
#---------------------------------------
def calcL2norm(data):
  # assuming data.shape = (niters, nsteps, 3)
  niters	= data.shape[0]
  nsteps	= data.shape[1]
    
  L2		= np.zeros( (niters, nsteps) )
    
  for i in range(niters):
    for j in range(nsteps):
      L2[i,j] = la.norm(data[i,j,:])    
      
  return L2
  
def calc_mommentum(vel):
  # calculate the mommentum for 1 timstep
  # vel has the shape (nsteps, N, 3)  
  
  nsteps = vel.shape[0]
  N = vel.shape[1]
  
  mommentum = np.zeros((nsteps, 3))
  
  for i in range(nsteps):
    
    mommentum[i,:] = sum(vel[i,:])
  
  return mommentum
  
#---------------------------------------
# Run
#---------------------------------------

data_dir = "/home/namdi/Documents/School/UNC/Parallel_Time/Data"

fpath = data_dir + "/Results/2014_05/10" 
fpath = fpath + "/pmrfmm_s64_a3_f5"

# load parameters
fname = fpath + "/t_0/info.txt"
nsources, ntargets, distribution = mu.load_parameters2(fname)

# the data has shape (niters, nsteps, N, 3)
# the data does not contain the initial condition
idx = 1
vel = mc.get_data(fpath, nsources, idx)

# load initial data y0 (2, N, 3)
fpath = fpath + "/t_0"
y0, charge_src, charge_targ, dipstr_src, dipvec_src   = mu.load_initial_data(fpath)

# store the number of steps (not including the initial condition) and iterations
niters = vel.shape[0]
nsteps = vel.shape[1]

# the mommentum of the initial condition
mom0 = sum(y0[1,:])

mom = np.zeros((niters, nsteps,3))
mag = np.zeros((niters, nsteps))
p = np.zeros(mag.shape)

# calculate the change in mommentum
for i in range(niters):
  for j in range(nsteps):
    mom[i,j] = sum(vel[i,j]) - mom0
    mag[i,j] = la.norm(mom[i,j])
    p[i,j] = np.log(mag[i,j])/np.log(10)

# plotting
label=[]
plt.xlabel("timstep/ processor")
plt.ylabel("log of change of mommentum")
for i in range(niters):
  #print p[i,nsteps-1]
  plt.plot(p[i,:])  
  label.append( "iter: " + str(i+1) )

plt.legend(label, 2)  
plt.show()  


uu = 1