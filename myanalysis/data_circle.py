"""" THis file creates a test initial condition to use of 4 targets and 4 sources.
The mass is one. The initial velocity is zero.  The sources are on the unit square on the plane z=0.
The targets are a diamond between the source points.
"""

import numpy as np
import os

def run(nsources, ntargets, nsteps, distribution, outdir):

  sources 	= np.zeros( (nsources,3), dtype=np.dtype('f8') )
  velocity	= np.zeros( (nsources, 3), dtype=np.dtype('f8') )

  targets 	= np.zeros( (ntargets, 3), dtype=np.dtype('f8') )
  veltarg	= np.zeros( targets.shape, dtype=np.dtype('f8') )

  dipstr  	= np.zeros( nsources, dtype=np.dtype('c16') )
  charge	= np.ones( nsources, dtype=np.dtype('c16') )
  dipvec	= np.zeros( (nsources, 3), dtype=np.dtype('f8') )
  mass		= np.ones( nsources, dtype=np.dtype('f8') )

  dt = 0.01
  R = 1.0	# radius

  max_charge = 1.0
  max_strength = 1.0

  # set the initial condtion (position & velocity) for the sourcs. The velocity is zero
  # position (sources)
  zmax = 1.0

  for i in range(nsources):
    sources[i,0] = R*np.cos(2*np.pi/nsources * i)
    sources[i,1] = R*np.sin(2*np.pi/nsources * i)

  # set the dipoles for the sources. It is zero

  # set the initial position for the targets
  targets[0,:] = [0.0, 0.0, -zmax]
  targets[1,:] = [0.0, 0.0, 0.0]
  targets[2,:] = [0.0, 0.0, zmax]

  # set the charge for the sources. It is one.

  # save some data to files
  fpath = outdir + "/t_0"

  #--------------------------------------------------------
  # save initial data 
  #--------------------------------------------------------

  q = np.zeros( (2, nsources + ntargets, 3), dtype=np.dtype('f8'))

  q[0, 0:nsources, :] 	= sources[:,:]
  q[1, 0:nsources, :] 	= velocity[:,:]
  q[0, nsources:, :]	= targets[:,:]
  q[1, nsources:, :]	= veltarg[:,:]

  # not using h5 files
  #h5 = h5py.File(data_path + '/electro.h5', 'w')
  #h5.create_dataset('q0', data=q)
  #h5.close()

  #------------------------------------------------------
  # save data to .txt files
  #--------------------------------------------------------
  
  if not os.path.exists(fpath):
    os.makedirs(fpath)

  # save the initial data
  np.savetxt( fpath + "/pos_source.txt", sources)
  np.savetxt( fpath + "/vel_source.txt", velocity)

  np.savetxt( fpath + "/pos_target.txt", targets)
  np.savetxt( fpath + "/vel_target.txt", veltarg)

  np.savetxt( fpath + "/d_pos_source.txt", sources)
  np.savetxt( fpath + "/d_vel_source.txt", velocity)

  np.savetxt( fpath + "/d_pos_target.txt", targets)
  np.savetxt( fpath + "/d_vel_target.txt", veltarg)

  np.savetxt( fpath + "/dipstr.txt", np.column_stack([dipstr.real, dipstr.imag]) )
  np.savetxt( fpath + "/dipvec.txt", dipvec)
  np.savetxt( fpath + "/charge.txt", np.column_stack([charge.real, charge.imag]) )
  np.savetxt( fpath + "/mass_source.txt", mass)

  # write basic info to info file
  f = open(fpath + '/info.txt', 'w')
  f.write("nsources: " + str(nsources) + "\n")
  f.write("ntargets: " + str(ntargets) + "\n")
  f.write("distribution: " + str(distribution) + "\n")
  f.write("dt: " + str(dt) + "\n")
  f.write("nsteps: " + str(nsteps) + "\n")
  f.close()
