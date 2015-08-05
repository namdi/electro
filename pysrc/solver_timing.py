#===========================================
# Explaination
#===========================================

"""
The purpose of this code is to compare the solvers (multirate FMM, regular FMM, direct). The code currently compares the Electric field caclualtions.
The code runs niters steps for testing.
The code uses ntrials to know how many test should be run so that we can get an average
If the save flag is True, the output is saved to the directory (ex: test.d)

Run the code as the following:
$ python solver_test.py -i test -o test.d -t 4 -k 10 -s $FMM -a 3
This code runs solver_test with
input: {$MY_DATA_DIR}/test,
output: test/test.d (This doesn't make sense if k is > 1)
iterations: 10 to average over
solver: full FMM
accuracy: 3
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
  
# read the command line. The only thing of relevance is nthreads and data_path. 
# Everything else is overwritten
nsources, ntargets, nthreads, distribution, accuracy, nsteps, indir, outdir, solver, niters = mu.parse_cmd_line()

#nsources = 4
#ntargets = 3
#nthreads = 1
#distribution = -1
#accuracy = 3
#nsteps = 1
#indir = "/home/namdi/Documents/School/UNC/Parallel_Time/Data/test.d"
#outdir = None
#solver = 2
#niters = 1

# the exact path to the initial data
input_dir =  indir + "/t_0"
  
# the info-file name
fname = input_dir + "/info.txt"
  
# store the actuall parameters from the info file
nsources, ntargets, distribution, dt, nsteps = mu.load_parameters(fname)

tol = mu.set_tol(accuracy)

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

# how many trials to run
ntrials = niters
runtime = np.zeros(ntrials)

print "accuracy: " + str(accuracy)

#-------------------------------------------------------
# Run the solvers for 1 timestep
#-------------------------------------------------------


for i in range(ntrials):
  # start timing
  t1 = time.time()

  # do the direct solver
  if dodirect:
  
    su.direct(pos_src, pos_targ, charge_src, dipstr_src, dipvec_src, pot_src, field_src, pot_targ, field_targ, 
	dosource, dotarget, dopot, dofield, docharge, dodipole, nthreads)
	
  # do some version of the FMM(full, multirate, near)
  else:
  
    su.fmm(pos_src, pos_targ, charge_src, dipstr_src, dipvec_src, accuracy,
	pot_src, field_src, pot_targ, field_targ,
	pot_far_src, 	pot_near_src,  field_far_src,  field_near_src,
	pot_far_targ,  pot_near_targ, field_far_targ, field_near_targ,
	dosource, dotarget, dofmm, dofmm_multi, dofmm_near,
	dopot, dofield, docharge, dodipole, nthreads)  

  t2 = time.time()

  # timing
  runtime[i] = t2 - t1
  
avg = np.mean(runtime)  

print
print "suggested FMM-field tol:  %.2e " % (tol * 10**3)
print "ntrials: " + str(ntrials)
print "average runtime[s]: %.2e" % avg

#-------------------------------------------------------
# save
#-------------------------------------------------------
#dofmm_multi = dofmm_near
if outdir is not None:   
    if not os.path.exists(output):
      os.makedirs(output)
      
    mu.save_pot_field_new( output, dofmm, dodirect, dofmm_multi, dofmm_near, dosource, dotarget, dopot, dofield,
		    pot_src, field_src, pot_targ, field_targ,
		    pot_far_src, pot_near_src, 
		    field_far_src, field_near_src,
		    pot_far_targ, field_far_targ,
		    pot_near_targ, field_near_targ)
		    
"""==========================================="""
# End
#===========================================
uu = 1