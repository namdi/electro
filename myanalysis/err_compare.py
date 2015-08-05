"""
This comapres the PFASST results from different runs

This code plots the error per iteration when 
comparing both data sources

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
# Run
#---------------------------------------

def calcL2norm(data):
  # assuming data.shape = (niters, nsteps, N, 3)  
  niters	= data.shape[0]
  nsteps	= data.shape[1]
  N		= data.shape[2]
  L2 		= np.zeros( (niters, nsteps) )
  
  temp = np.zeros((N))
  
  for itr in range(niters):    
    for step in range(nsteps):
      temp[:] = 0      
      for i in range(N):
	temp[i] = la.norm( data[itr, step, i, :], 2 )
      L2[itr, step] = temp.max()
      
  return L2
  
data_dir = "/home/namdi/Documents/School/UNC/Parallel_Time/Data"

fpath1 = data_dir + "/Results/2014_05/09" # for pfmm_s64_a3_f5
fpath2 = data_dir + "/Results/2014_10/13" 
fpath3 = data_dir + "/Results/2014_10/14"

#fpath3 = data_dir + "/kd"
info_dir = data_dir + "/Results/2014_02/11/pdirect_s64"

app1 = fpath3 + "/gpfmm_s128_a23_f5"
app2 = fpath3 + "/gpfmm_s128_a33_f5"

# load parameters
fname = info_dir + "/t_0/info.txt"
nsources, ntargets, distribution = mu.load_parameters2(fname)
idx = 1 # pos/vel index

plot_aerr = True

# the data has shape (niters, nsteps, N, 3)
# the data does not contain the initial condition
data1 = mc.get_data(app1, nsources, idx)
data2 = mc.get_data(app2, nsources, idx)

# store the number of steps (not including the initial condition) and iterations
niters1 = data1.shape[0]
nsteps1 = data1.shape[1]

niters2 = data2.shape[0]
nsteps2 = data2.shape[1]

nsteps = nsteps2
niters = niters2

accurate = data1
#accurate = data1[:,1:nsteps1:2] # when accurate has double the steps as approx
approx = data2

# absolute error
aerr = abs(accurate - approx)

# relative error
rerr = abs(aerr / accurate)

# the x-axis as steps
x = np.arange(1, nsteps+1)

for i in range(0,2):
  
  # alternate between  absolute and relative error
  if i == 0:
    err = aerr
  else:
    err = rerr
    
  # norm is in the shape (niters, nsteps)  
  norm = calcL2norm(err)
  logNorm = np.log(norm)/np.log(10) # shape(niters, nsteps)
  p = logNorm * -1
  
  # set the figure
  fig = plt.figure(i)
  
  # add the legend  
  label = []
  plt.xlabel("timestep/ processor")
  plt.ylabel("log10 of error")

  # plot
  for i in range(niters):  
    plt.plot(x, logNorm[i,:])
    label.append("iter: " + str(i+1))
  
  plt.legend(label, loc=2)

# show plots  
plt.show()
uu = 1