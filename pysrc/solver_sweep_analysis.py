"""
  This file does analysis between the data and the least squares polynomial approximations to the data.
  These are done to see how smooth the data is notabbly the far-field potential to the targets and the near-field
  potential.
"""
#--------------------------------------------
# import
#-------------------------------------------
import numpy as np
import os
import matplotlib.pylab as plt
import commands

#import module_path
import my_utils
import data_dir
import load_utils as lu
import numpy.linalg as la
#-------------------------------------------

"""----------------------------------------------------------"""
# run
#-------------------------------------------

# this only works if one ran $ source my_activate
x		= commands.getstatusoutput("echo $MY_DATA_DIR")
data_path	= x[1]
fpath		= data_path + "/test.d"

# read in the parameters
fname	= "/t_0/info.txt"
nsources, ntargets, distribution, dt, nsteps = my_utils.load_parameters(fpath + fname)

# fmm accuracy
acc = np.arange(-2,6)
nfmm = len(acc)

# all of the fields (FMM + direct)
field_all = np.zeros((nfmm+1, nsources, 3))

# absolute error and relative error
aerr = np.zeros((nfmm, nsources, 3))
rerr = np.zeros(aerr.shape)

# the norm, maxnorm, and the degree of precision with the maxnorm
norm = np.zeros((nfmm, nsources))
maxnorm = np.zeros(nfmm)
p = np.zeros(maxnorm.shape)

# load the runtime
fname = "/runtime.txt"
runtime = np.loadtxt(fpath + fname)

# load the field data
dtype = 'f8'
for i in range(nfmm + 1):
  
  if i < nfmm:
    fname 	= "/field_source." + str(acc[i]) + ".txt"    
  else:
    fname 	= "/field_source.d.txt"
          
  temp	= np.loadtxt( fpath + fname).view(dtype)
  field_all[i,:] = temp[:,0:3]
  
# calculate the error  
for i in range(nfmm):
  # absolute and relative error
  aerr[i,:] = abs(field_all[i,:] - field_all[nfmm,:])
  rerr[i,:] = abs(aerr[i,:] / field_all[nfmm,:])
  
# specify type of error  
err = rerr

# L2 norm for each particle
for i in range(nfmm):
  
  for j in range(nsources):  
    norm[i,j] = la.norm(err[i,j,:])
    
  maxnorm[i] = max(norm[i,:])
  
#digits of precision
p = np.log(maxnorm)/np.log(10)  * -1.0

plt.figure(0)
plt.plot(runtime, "r*")
plt.xlabel("solver type")
plt.ylabel("runtime [s]")

plt.figure(1)
plt.plot(acc, p, '.')
plt.xlabel("FMM accuracy")
plt.ylabel("digits of precision")

plt.show()