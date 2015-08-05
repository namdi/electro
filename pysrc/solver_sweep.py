#=========================================== 
# Explaination
#===========================================

"""
This code runs sweeep for the FMM field strength over multiple spatial accuracies.  It also runs the 
direct solver. The results are stored in the data filed $MY_DATA_DIR/test/field_source.1.txt (if 1 is the precison)

Run the code as:
$ python solver_sweep.py -t 4 -i $MY_DATA_DIR/test -o $MY_DATA_DIR/test

"""

#---------------------------
# Import
#---------------------------
import numpy as np
import time
import os

import my_utils as mu
import solver_utils as su
#---------------------------
  
#========================================
# Run
#========================================

show_time = True
#------------------------------------------
# Set up parameters
#------------------------------------------
  
# read the command line. The only thing of relevance is nsources, nthreads, accuracy, and data_path. 
# Everything else is overwritten
nsources, ntargets, nthreads, distribution, accuracy, nsteps, indir, outdir, solver, niters = mu.parse_cmd_line()

# the exact path to the initial data
input_dir =  indir + "/t_0"
  
# the info-file name
fname = input_dir + "/info.txt"
  
# store the actuall parameters from the info file
nsources, ntargets, distribution, dt, nsteps = mu.load_parameters(fname)

print nsources
# set the accuracy
acc = np.arange(-2,6)

# how many FMM trials to run
nfmm = len(acc)

# how many trials to run (FMM + direct)
runtime = np.zeros(nfmm + 1)

# output directory name
if outdir is None:
  output = None
else:
  output = outdir
  
# Load the INITIAL CONDITION
q0, charge_src, charge_targ, dipstr_src, dipvec_src = mu.load_initial_data(input_dir)

pos_src		= q0[0, 0:nsources,:]
pos_targ 	= q0[0, nsources:, :]

vel_src 	= q0[1, 0:nsources, :]
vel_targ	= q0[1, nsources:, :]

pot_src		= np.zeros((nsources,1), np.dtype('c16'))
pot_targ	= np.zeros((ntargets,1), np.dtype('c16'))

field_src	= np.zeros((nsources,3), np.dtype('c16'))
field_targ	= np.zeros((ntargets,3), np.dtype('c16'))

pot_far_src	= np.zeros(pot_src.shape, pot_src.dtype)
pot_near_src	= np.zeros(pot_src.shape, pot_src.dtype)

pot_far_targ	= np.zeros(pot_targ.shape, pot_targ.dtype)
pot_near_targ	= np.zeros(pot_targ.shape, pot_targ.dtype)

field_far_src	= np.zeros(field_src.shape, field_src.dtype)
field_near_src	= np.zeros(field_src.shape, field_src.dtype)

field_far_targ		= np.zeros(field_targ.shape, field_targ.dtype)
field_near_targ		= np.zeros(field_targ.shape, field_targ.dtype)

# flags if to indicate if we are calculating source and/or target particles
dosource	= True
dotarget	= False

# flags to indicate whether we are caclualting the potential and/or the field
dopot 		= False
dofield		= True

# flag to calculate monopole and dipole
docharge	= True
dodipole 	= False

# flags for multirate, fmm, or direct
dofmm_near, dofmm_multi, dofmm, dodirect = mu.get_solver(solver)  


field_all = np.zeros((nfmm +1, nsources, 3), field_src.dtype)

#-------------------------------------------------------
# Run the solvers for 1 timestep
#-------------------------------------------------------

# run the FMM solvers

for i in range(nfmm + 1):    
  
  # start timing
  t1 = time.time()

  # do the direct solver
    
  if i == nfmm:
  
    su.direct(pos_src, pos_targ, charge_src, dipstr_src, dipvec_src, pot_src, field_src, pot_targ, field_targ, 
	dosource, dotarget, dopot, dofield, docharge, dodipole, nthreads)
    
    print "direct solver ..."
    
  # do some version of the FMM(full, multirate, near)
  else:
    print "FMM accuracy: " + str(acc[i]) + str(" ...")
    su.fmm(pos_src, pos_targ, charge_src, dipstr_src, dipvec_src, acc[i],
	pot_src, field_src, pot_targ, field_targ,
	pot_far_src, 	pot_near_src,  field_far_src,  field_near_src,
	pot_far_targ,  pot_near_targ, field_far_targ, field_near_targ,
	dosource, dotarget, dofmm, dofmm_multi, dofmm_near,
	dopot, dofield, docharge, dodipole, nthreads)  

  t2 = time.time()

  # timing
  runtime[i] = t2 - t1
  
  # store the field
  field_all[i,:] = field_src[:]   
  
  field_src[:] = 0
  field_far_src[:] = 0
  field_near_src[:] = 0
  
  pot_src[:] = 0
  pot_far_src[:] = 0
  pot_near_src[:] = 0
  
  field_targ[:] = 0
  field_far_targ[:] = 0
  field_near_targ[:] = 0
  
  pot_targ[:] = 0
  pot_far_targ[:] = 0
  pot_near_targ[:] = 0
  
# save the data
if output is not None:
  
  # save the FMM E-fields
  for i in range(nfmm + 1):
    if i < nfmm:
      fname = "/field_source." + str(acc[i]) + ".txt"      
    else:
      fname = "/field_source.d.txt"
      
    np.savetxt( outdir + fname, np.column_stack([field_all[i,:].real, field_all[i,:].imag]) )
    
  # save the runtime
  fname = "/runtime.txt" 
  np.savetxt( outdir + fname, runtime )  
  
uu = 1