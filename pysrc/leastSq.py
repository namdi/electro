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
#-------------------------------------------

#-------------------------------------------
# Function declarations
#-------------------------------------------

def find_approx(coeff, t):
  
  """ This function finds the least squares approximation at the data points, t, by evaluating
  the least Sq coefficients, coeff, at the data points"""
  
  # the number of particles and number of timesteps
  npart = coeff.shape[1]
  nsteps = t.shape[0] - 1
  
  approx = np.zeros( (npart, nsteps+1), dtype=np.dtype('f8'))
  
  # evaluate the least sq coefficients at the data points
  for i in range(npart):
    approx[i,:] = np.polyval(coeff[:,i], t)
  
  return approx
  
def potential_analysis(potential, t, degree):
  """ calculate:
  pot_lsq:	the least squares approximation to the potential: 	
  maxerr:	the maxerr over all particles between the least squares approx 
		and the numerical calculation over for each time
  log_maxerr:	The amounts of precision (base 10) in the error
  
  """
  
  temp = np.transpose(potential.real)
  
  # find the coefficients of the least squares approximation
  pot_lsq_coeff	= np.polyfit(t, temp, degree)
  
  # make sure temp points to nothing
  temp = None
  
  # evaluate the least squares approx at the data points
  pot_lsq 	= find_approx(pot_lsq_coeff, t)
  
  # the error of the approximation for every particle at every data point
  
  # the absolute error between the actual value and the approximation, 
  aerr = abs( potential.real - pot_lsq)
  rerr = abs(aerr/ potential.real)
  
  # set the error
  err = aerr
  
  # find the maximum error
  maxerr = max_error(err)
  
  # find the precision 
  log_maxerr = np.log(maxerr) / np.log(10) * -1
  
  return pot_lsq, maxerr, log_maxerr

def max_error(err):
  
  """ Find the maximum error over all the particles at every time"""
  
  nsteps 	= err.shape[1] - 1
  maxerr	= np.zeros(nsteps+1)
  
  for i in range(nsteps+1):
    maxerr[i] = max(err[:,i])

  return maxerr
  

def smoothness(potential_far, potential_lcl, t, degree):
  """ This function calculates a degree "deg" least sq approx to the far field .
  Then it tries to find the equivalent least sq approx with about the same accuracy 
  to the potential in the local field"""
  
  nsteps = len(t) -1
  
  # the largest degree of the least sq approximation of the local field
  deg_max 	= 20
  
  # the current degree of the least sq approx
  deg = degree
  
  stop = False
  
  # get the data for the least sq approximation for the far and near field
  pot_lsq_far, maxerr_far, log_maxerr_far, prec_err_far = smooth_data(potential_far, t, deg)
  pot_lsq_lcl, maxerr_lcl, log_maxerr_lcl, prec_err_lcl = smooth_data(potential_lcl, t, deg)
  
  if (prec_err_lcl >= prec_err_far):
    stop = True
  
  print "deg: " + str(deg)
  print "prec_err_far: " + str(prec_err_far)
  print "prec_err_lcl: " + str(prec_err_lcl)
  
  # check accuarcy of the lcl approximation with increasing degrees of precision
  
  while stop is False and deg <= deg_max:
    deg = deg + 2
    
    pot_lsq_lcl, maxerr_lcl, log_maxerr_lcl, prec_err_lcl = smooth_data(potential_lcl, t, deg)
    
    print "deg: " + str(deg)
    print "prec_err_lcl: " + str(prec_err_lcl)
    
    # we achieved the same accuracy
    if prec_err_lcl == prec_err_far:
      stop = True

    # the approx may be more accuracte than the far field approx. Try least sq approx for degree = deg-1
    elif prec_err_lcl > prec_err_far:
      stop = True
      
      pot_lsq_lcl2, maxerr_lcl2, log_maxerr_lcl2, prec_err_lcl2 = smooth_data(potential_lcl, t, deg-1)
      
      if prec_err_lcl2 >= prec_err_far:
	prec_err_lcl 	= prec_err_lcl2
	pot_lsq_lcl	= pot_lsq_lcl2
	maxerr_lcl	= maxerr_lcl2
	log_maxerr_lcl	= log_maxerr_lcl2
	deg = deg-1
	
  return pot_lsq_far, maxerr_far, log_maxerr_far, prec_err_far, pot_lsq_lcl, maxerr_lcl, log_maxerr_lcl, prec_err_lcl, deg

def smooth_data(potential, t, deg):
  
  """" This function calculates the degree "deg" least sq approximation to potential due to info at t
  
  The function returns:
    pot:	the degree "deg" least sq approx
    maxerr:	the max error (per iteration)
    log_maxerr:	the degree of precision of the err(per iteration) 
    prec_err:	the mean of the degree of precision error (rounded)
  """
  
  # store the number of time steps
  nsteps = len(t) - 1
  
  # calculate the degree "deg" least sq approx, the max error (per iteration), the degree of precision of the err(per iteration) 
  pot, maxerr, log_maxerr 	= potential_analysis(potential, t, deg)
  
  # calculate the mean of the degree of precision of the error
  prec_err = np.sum(log_maxerr) /(nsteps+1)
  prec_err = np.round(prec_err)
  
  return pot, maxerr, log_maxerr, prec_err
  
  
""" plotting """

def plot_log_maxerr(log_maxerr, degree, title, fig_id,):
  
  plt.figure(fig_id)
  
  plt.xlabel('timestep')
  plt.ylabel('digits of precision')
  title = title + " deg = " + str(degree)
  plt.title(title)
  
  plt.plot(log_maxerr, '-o')
  
  return

def plot_log_maxerr_combo(log_maxerr1, log_maxerr2, label1, label2, title, fig_id,):
   # set up the plotter
  fig = plt.figure(fig_id)
  ax1 = fig.add_subplot(111)
  
  # plot the data
  ax1.plot(log_maxerr1, '-o')
  ax1.plot(log_maxerr2, '-o')

  # create a label for the legend
  label = []
  
  label.append(label1)
  label.append(label2)
  
  plt.xlabel('timestep')
  plt.ylabel('digits of precision')
  plt.legend(label)
  plt.title(title)
  
  return

"""  
def plot_log_maxerr_combo(log_maxerr1, log_maxerr2, label1, label2, degree, fig_id,):
   # set up the plotter
  fig = plt.figure(fig_id)
  ax1 = fig.add_subplot(111)
  
  # plot the data
  ax1.plot(log_maxerr1, '-o')
  ax1.plot(log_maxerr2, '-o')

  # create a label for the legend
  label = []
  
  label.append(label1)
  label.append(label2)
  
  plt.xlabel('timestep')
  plt.ylabel('digits of precision')
  plt.legend(label)

  
  title = "error of leastSq approx. deg = " + str(degree)
  plt.title(title)
  
  return
"""
def plot_potential_combo(pot1, pot2, title, label1, label2, degree, fig_id):
  # set up the plotter
  fig = plt.figure(fig_id)
  ax1 = fig.add_subplot(111)
  
  # plot the data
  ax1.plot(pot1, '-')
  ax1.plot(pot2, '-')

  # create a label for the legend
  label = []
  
  label.append(label1)
  label.append(label2)
  
  plt.xlabel('timestep')
  plt.ylabel('potential')
  plt.legend(label)
  plt.title(title)
  
  return

def plot_log_maxerr_dist(potential_far, potential_lcl, t, deg_far, deg_lcl, fig_id):
  """ The name stands for plot the log of the max error of a distribution of least square polynomials with varying degrees.
      This function plots the digits of precision in the error between the least squares approx and FMM calculation of the far and near-fields.
      The function prints out the error from the least squares approx of varying degrees for comparision. """
  
  # set up plotter
  fig = plt.figure(fig_id)
  ax1 = fig.add_subplot(111)
  
  plt.title("error in least squares approx")
  plt.xlabel("timestep")
  plt.ylabel("digits of precision")
      
  # create a label for the legend
  label = []
  
  # first plot the error due to the far field
  for deg in deg_far:
    pot_far, maxerr_far, log_maxerr_far =  potential_analysis(potential_far, t, deg)
    ax1.plot(log_maxerr_far, '-o')
    label.append("far, deg = " + str(deg))
    
  """ plot the error due to the near field"""
  for deg in deg_lcl:
    pot_lcl, maxerr_lcl, log_maxerr_lcl =  potential_analysis(potential_lcl, t, deg)
    ax1.plot(log_maxerr_lcl, '--x')
    label.append("near, deg = " + str(deg))
    
  # add the legend
  plt.legend(label)
  
  return
  
  
""" saving """
def save_pot_target_approx_deg(fpath, potential_far, potential_lcl, degree):
  
  nsteps	 = potential_far.shape[1] - 1
  npart 	 = potential_far.shape[0]
  
  #calculate the approximation, the max error of the approx, and the degree of precision (base 10) of the max error
  pot_far, maxerr_far, log_maxerr_far =  potential_analysis(potential_far, t, degree)
  pot_lcl, maxerr_lcl, log_maxerr_lcl =  potential_analysis(potential_lcl, t, degree)
  
  for i in range(nsteps+1):
    fpath_temp = fpath + "t_" + str(i) + "/"
    
    fname = "pot_target_far_approx_deg" + str(degree) + ".txt"
    np.savetxt( fpath_temp + fname, pot_far[:,i] )
    
    fname = "pot_target_lcl_approx_deg" + str(degree) + ".txt"
    np.savetxt( fpath_temp + fname, pot_lcl[:,i] )
    
  return pot_far, pot_lcl, maxerr_far, maxerr_lcl, log_maxerr_far, log_maxerr_lcl



"""----------------------------------------------------------"""
# run
#-------------------------------------------

# this only works if one ran $ source my_activate
x		= commands.getstatusoutput("echo $MY_DATA_DIR")
data_path	= x[1] 
#fpath 		= data_path + "/kd/test_fe"
fpath 		= data_path + "/Results/2014_10/14/least_sq_dt"
fname	= "/t_0/info.txt"

# read in the parameters
nsources, ntargets, distribution, dt, nsteps = my_utils.load_parameters(fpath + fname)

# get the times
t = np.linspace(start=0.0, stop=dt*nsteps, num=nsteps+1, endpoint = True)

save = False

doTarget = False
# this is the actual data
if doTarget:
  pot_far	= lu.load_cplex_potential(fpath, "/pot_far_target.txt", ntargets, nsteps)
  pot_near	= lu.load_cplex_potential(fpath, "/pot_near_target.txt", ntargets, nsteps)
else:
  pot_far	= lu.load_cplex_potential(fpath, "/pot_far_source.txt", nsources, nsteps)
  pot_near	= lu.load_cplex_potential(fpath, "/pot_near_source.txt", nsources, nsteps)


# set the degree of the approximating polynomial
deg = 3

# find out how many least sq coefficients are needed to have the far field and local field have the same accuarcy
pot_lsq_far, maxerr_far, log_maxerr_far, prec_err_far, pot_lsq_near, maxerr_near, log_maxerr_near, prec_err_near, deg_near = smoothness(pot_far, pot_near, t, deg)

#title = "error of leastSq approx"
#plot_log_maxerr_combo(log_maxerr_far, log_maxerr_near, "far-field"+str(deg), "near-field"+str(deg_near),title, 3)


# calculate the approximation, the max error of the approx, and the degree of precision (base 10) of the max error
pot_lsq_far, maxerr_far, log_maxerr_far 	=  potential_analysis(pot_far, t, deg)
pot_lsq_near, maxerr_near, log_maxerr_near 	=  potential_analysis(pot_near, t, deg)

title = "error of leastSq approx. deg = " + str(deg)
#plot_log_maxerr_combo(log_maxerr_far, log_maxerr_near, "far-field", "near-field", title, 4)


# plot the error over a distribution of degrees
deg_far = np.array( (0,1, 2), dtype=np.dtype('i4') )
deg_near = np.array( (2,3,4), dtype=np.dtype('i4') )

plot_log_maxerr_dist(pot_far, pot_near, t, deg_far, deg_near, 5)

"""
# Testing
idx = 0
title = "target far-field potential"
label1 = "fmm solution"
label2 = "L2 approx deg = " + str(deg)
plot_potential_combo(potential_far[idx,:].real, pot_far[idx,:], title, label1, label2, deg, 10)

title = "target local-field potential"
plot_potential_combo(potential_lcl[idx,:].real, pot_lcl[idx,:], title, label1, label2, deg, 11)
"""
plt.show()
