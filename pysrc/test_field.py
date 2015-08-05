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
fpath		= data_path + "/test"
fpath1 		= data_path + "/test.-1"
fpath2		= data_path + "/test.d"

fname	= "/t_0/info.txt"

# read in the parameters
nsources, ntargets, distribution, dt, nsteps = my_utils.load_parameters(fpath + fname)

field1 = np.zeros((nsources, 3))
field2 = np.zeros(field1.shape)

norm = np.zeros((nsources, 1))

# get the times
t = np.linspace(start=0.0, stop=dt*nsteps, num=nsteps+1, endpoint = True)

dtype = 'f8'

temp		= np.loadtxt( fpath1 + "/field_source.txt").view(dtype)
field1[:,:] 	= temp[:,0:3]

temp		= np.loadtxt( fpath2 + "/field_source.txt").view(dtype)
field2[:,:] 	= temp[:,0:3]

# absolute error
aerr = abs(field1 - field2)

# relative error
rerr = abs(aerr/field2)

err = rerr

# L2 norm for each particle
for i in range(nsources):
  norm[i] = la.norm(err[i,:])

#digits of precision
p = np.log(norm)/np.log(10)  * -1.0

max_val = max(norm)
min_p = min(p)

print "max err: " + str(max_val)
print "min p: " + str(min_p)