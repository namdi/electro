""" This is the driver for the 3d laplace potential
This program's purpose is to show the multiple spacial scales for the far-field and the near-field.
This is done by using the fmm serially in time. However, the tree structure is FIXED on the 
tree structure based on the initial condition.

The idea is that after the run, we check the least-squares errors of the near/ far field of the target particles.

$ python data_script.py -n 1000 -m 4 -d 2 -o test 

"""


#============================================
# Import
#============================================
import numpy as np
import os
import time

import my_utils as mu
import st_fmm_test_fe as fe # test using forward Euler
#import st_fmm_test_verlet

"""============================================"""
# Run
#============================================

# get parameters from command line
nsources, ntargets, nthreads, distribution, accuracy, nsteps, indir, outdir, solver, niters = mu.parse_cmd_line()

# flag to ouput data at every iteration
save = True

# set the precision flag 
iprec = accuracy
tol = mu.set_tol(iprec)

# flags to do ANY calculations on the fmm/ direct calculation. These must be 1 for the program to wok.  Because the position calculations
# depend on the full FMM and the direct calc, respectively
iffmm    = 0
ifdirect = 0

# the iffmm_multi flag
iffmm_multi = 1

# set the source type flags and output flags
# tells wheather to calculate the source potential/field for the direct/fmm
# need the potential for least squares analysis
# need the field to move the particles
ifpot = 1
iffld = 1

ifcharge=1
ifdipole=0	# turn the dipole strength off

ifpottarg=1	# add a parameter that checks to see if that if info.txt says targets are 0. These are set to zero.
iffldtarg=1

# My flags to tell tell whether to calculate the local and farfield interaction
iffldsrc_far	= 1
iffldsrc_near	= 1

iffldtarg_far	= 1
iffldtarg_near	= 1

ifpotsrc_far	= 0	# default 1
ifpotsrc_near	= 0	# default 1

ifpottarg_far	= 1
ifpottarg_near	= 1

#------------------------------------
# LOAD INITIAL CONDITION
#------------------------------------

# initial condition 
fpath	= indir + "/t_0"
fname 	= fpath + "/info.txt"
  
# load parameters
nsources, ntargets, distribution, dt, nsteps = mu.load_parameters(fname)
    
# load data from data files
if ifcharge == 1:
  charge_src = np.loadtxt( fpath + "/charge.txt").view('c16')
else:
  charge_src = np.zeros(nsources, dtype=np.dtype('c16') )
    
if ifdipole == 1:
  dipstr_src = np.loadtxt( fpath + "/dipstr.txt").view('c16')  
else:
  dipstr_src = np.zeros(nsources, dtype=np.dtype('c16') )
dipvec_src		= np.loadtxt( fpath + "/dipvec.txt").view('f8')

# load position/ velocity for source and targets
pos_src		= np.loadtxt( fpath + "/pos_source.txt").view('f8')
vel_src 	= np.loadtxt( fpath + "/vel_source.txt").view('f8')

pos_targ	= np.loadtxt( fpath + "/pos_target.txt").view('f8')
vel_targ	= np.loadtxt( fpath + "/vel_target.txt").view('f8')

# load the masses
mass_src	= np.loadtxt( fpath + "/mass_source.txt").view('f8')
mass_targ	= np.ones(ntargets, dtype=np.dtype('f8'))

# for the calculations
pot_src		= np.zeros( nsources, dtype=np.dtype('c16') )
field_src	= np.zeros( (nsources, 3), dtype=np.dtype('c16') )

pot_targ	= np.zeros( (ntargets), dtype=np.dtype('c16') ) 
field_targ	= np.zeros( (ntargets, 3), dtype=np.dtype('c16') )


# for the direct calculation, copy info into the direct position/velocity
d_pos_src	= np.array(pos_src, dtype=np.dtype('f8') )
d_vel_src 	= np.array(vel_src, dtype=np.dtype('f8') )

d_pos_targ	= np.array(pos_targ, dtype=np.dtype('f8'))
d_vel_targ	= np.array(vel_targ, dtype=np.dtype('f8'))

#------------------------------------------
# create potential and field arrays
#------------------------------------------

d_pot_src	= np.zeros( pot_src.shape, dtype=np.dtype('c16') )
d_field_src	= np.zeros( field_src.shape, dtype=np.dtype('c16') )

d_pot_targ	= np.zeros( pot_targ.shape, dtype=np.dtype('c16') )
d_field_targ	= np.zeros( field_targ.shape, dtype=np.dtype('c16') )

# for the FMM
pot_far_src	= np.zeros(pot_src.shape, pot_src.dtype)
pot_near_src	= np.zeros(pot_src.shape, pot_src.dtype)

field_far_src 	= np.zeros(field_src.shape, field_src.dtype)
field_near_src 	= np.zeros(field_src.shape, field_src.dtype)

pot_far_targ	= np.zeros( pot_targ.shape, pot_targ.dtype )
pot_near_targ	= np.zeros( pot_targ.shape, pot_targ.dtype )

field_far_targ		= np.zeros( field_targ.shape, field_targ.dtype )
field_near_targ		= np.zeros( field_targ.shape, field_targ.dtype )

#------------------------------------------
# set up timestepping info    
#------------------------------------------
t0 = 0.0
tend = t0 + nsteps * dt
times = np.linspace(t0, tend, num=nsteps+1, endpoint=True)

# time to do the fmm and direct calculations for each iteration
  
tfmm = np.zeros(times.shape)	# the +1 is because the time loop goes through (nsteps + 1) iterations
tdirect = np.zeros(times.shape)

# initialize timing call (in seconds)
tick = time.time()

# save the autorun parameters (#-1 means N/A)
iterations="N/A"
nnodes="N/A"
levels="N/A"
refinement="N/A"
mpi_nodes=1
#mu.save_params(input_dir, iterations, nnodes, levels, refinement, mpi_nodes, nthreads, iprec)

#------------------------------------------
# run the FMM
#------------------------------------------
ier, nlev, nboxes, max_box_fine_lvl, lused7 = fe.run( iprec, ifpot, iffld, ifcharge, ifdipole, ifpottarg, iffldtarg, 
						pos_src, nsources, charge_src, vel_src, dipstr_src, dipvec_src, pot_src, field_src, mass_src, 
						pos_targ, ntargets, vel_targ, pot_targ, field_targ, mass_targ, nthreads, 
						times, dt, nsteps, d_pos_src, d_vel_src, d_pos_targ, d_vel_targ, d_pot_src, d_field_src,
						ifdirect, iffmm, iffmm_multi, d_pot_targ, d_field_targ, tfmm, tdirect,
						ifpotsrc_far,  ifpotsrc_near,
						ifpottarg_far, ifpottarg_near,
						iffldsrc_far,  iffldsrc_near,
						iffldtarg_far, iffldtarg_near,
						pot_far_src, pot_near_src, pot_far_targ, pot_near_targ,
						field_far_src, field_near_src,
						field_far_targ, field_near_targ, save, outdir )												 									
	
toc = time.time()
