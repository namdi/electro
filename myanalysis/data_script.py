"""
This program sets the initial condition and writes it to an input directory called indir/t_0.

The program uses only the following paramters:
nsoures, ntargets, distribution, accuracy, nsteps, indir, outdir.

distribution has the following layount:
1: random particles in a cube with random charges
2: random particles on the SURFACE of a sphere, with random charges

Tet cases:
-1: 4 particles on 2D square, with uniform charge
-2: particles uniformly distributed on the CIRCUMFERENCE of a circle, uniform charges
"""

#--------------------
# Import
# -------------------
import numpy as np
import os
#import h5py

import my_utils as mu
import initial

#--------------------------
# Run
#-------------------------

# read from the command line/ use default paramters when necessary
nsources, ntargets, nthreads, distribution, accuracy, nsteps, indir, outdir, solver, niters, nlevs, fnodes = mu.parse_cmd_line()

if outdir is None:
  raise NotImplementedError, 'need to specify an output directory'

pos_src 	= np.zeros( (nsources,3), dtype=np.dtype('f8') )
vel_src		= np.zeros( pos_src.shape, dtype=np.dtype('f8') )

pos_targ 	= np.zeros( (ntargets, 3), dtype=np.dtype('f8') )
vel_targ	= np.zeros( pos_targ.shape, dtype=np.dtype('f8') )

charge_src	= np.ones( nsources, dtype=np.dtype('c16') )
dipstr_src  	= np.ones( nsources, dtype=np.dtype('c16') )
dipvec_src	= np.zeros( pos_src.shape, dtype=np.dtype('f8') )
mass_src	= np.ones( nsources, dtype=np.dtype('f8') )

dt = 0.001
length = 1.0

max_charge = 1.0
max_strength = 1.0

if distribution > 0:
  
  # set the initial position for the sources
  initial.set_initial_position(distribution, pos_src, length)
  
  # set the dipoles for the sources 
  initial.set_dipole(dipstr_src, dipvec_src, max_strength)

  # set the charge for the sources
  initial.set_charge(charge_src, max_charge)
  
  # set the initial position for the targets
  initial.set_initial_position(distribution, pos_targ, length)

# test case: the 4 particles in a square test case, uniform charge
elif distribution == -1:
  
  nsources = 4
  ntargets = 3
  
  pos_src 	= np.zeros( (nsources,3), dtype=np.dtype('f8') )
  vel_src	= np.zeros( pos_src.shape, dtype=np.dtype('f8') )

  pos_targ 	= np.zeros( (ntargets, 3), dtype=np.dtype('f8') )
  vel_targ	= np.zeros( pos_targ.shape, dtype=np.dtype('f8') )

  charge_src	= np.ones( nsources, dtype=np.dtype('c16') )
  dipstr_src  	= np.ones( nsources, dtype=np.dtype('c16') )  
  dipvec_src	= np.zeros( pos_src.shape, dtype=np.dtype('f8') )
  mass_src	= np.ones( nsources, dtype=np.dtype('f8') )
  
  # set the square test case:  
  initial.square_test(pos_src, pos_targ, length)
  
  # set charge to be uniform 
  charge_src = charge_src * max_charge
  
  # set dipole to have random normalized directions but zero strength
  max_strength = 0.0
  initial.set_dipole(dipstr_src, dipvec_src, max_strength)
  
# test case: uniformly distributited particles on a 2D sphere with uniform charge
elif distribution == -2:
    
  ntargets = 3
  
  pos_src 	= np.zeros( (nsources,3), dtype=np.dtype('f8') )
  vel_src	= np.zeros( pos_src.shape, dtype=np.dtype('f8') )

  pos_targ 	= np.zeros( (ntargets, 3), dtype=np.dtype('f8') )
  vel_targ	= np.zeros( pos_targ.shape, dtype=np.dtype('f8') )

  charge_src	= np.ones( nsources, dtype=np.dtype('c16') )
  dipstr_src  	= np.ones( nsources, dtype=np.dtype('c16') )  
  dipvec_src	= np.zeros( pos_src.shape, dtype=np.dtype('f8') )
  mass_src	= np.ones( nsources, dtype=np.dtype('f8') )  
  
  # set the circle test case:  
  initial.circle_test(pos_src, pos_targ, length)
  
  # set charge to be uniform 
  charge_src = charge_src * max_charge
  
  # set dipole to have random normalized directions but zero strength
  max_strength = 0.0
  initial.set_dipole(dipstr_src, dipvec_src, max_strength)
  
#------------------------------------------------------
# save data to .txt files
#--------------------------------------------------------

# save some data to files
fpath = outdir + "/t_0"

if not os.path.exists(fpath):
  os.makedirs(fpath)

# save the initial data
np.savetxt( fpath + "/pos_source.txt", pos_src)
np.savetxt( fpath + "/vel_source.txt", vel_src)

np.savetxt( fpath + "/pos_target.txt", pos_targ)
np.savetxt( fpath + "/vel_target.txt", vel_targ)

np.savetxt( fpath + "/d_pos_source.txt", pos_src)
np.savetxt( fpath + "/d_vel_source.txt", vel_src)

np.savetxt( fpath + "/d_pos_target.txt", pos_targ)
np.savetxt( fpath + "/d_vel_target.txt", vel_targ)

np.savetxt( fpath + "/dipstr.txt", np.column_stack([dipstr_src.real, dipstr_src.imag]) )
np.savetxt( fpath + "/dipvec.txt", dipvec_src)
np.savetxt( fpath + "/charge.txt", np.column_stack([charge_src.real, charge_src.imag]) )
np.savetxt( fpath + "/mass_source.txt", mass_src)

# write basic info to info file
f = open(fpath + '/info.txt', 'w')
f.write("nsources: " + str(nsources) + "\n")
f.write("ntargets: " + str(ntargets) + "\n")
f.write("distribution: " + str(distribution) + "\n")
f.write("dt: " + str(dt) + "\n")
f.write("nsteps: " + str(nsteps) + "\n")
f.close()