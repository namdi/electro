""" 
This program's purpose is to show the multiple spacial scales for the far-field and the near-field.
This is done by using the fmm serially in time using FORWARD EULER. However, the tree structure 
is FIXED on the tree structure based on the initial condition. 
"""

#============================================
# Import
#============================================
import numpy as np
import os
import time

import my_utils as mu
import fortfmm
import solver_utils as su

#============================================
# function declarations
#============================================
 
def run(iprec, ifpot, iffld, ifcharge, ifdipole, ifpottarg, iffldtarg, 
	pos_src, nsources, charge_src, vel_src, dipstr_src, dipvec_src, pot_src, field_src, mass_src, 
	pos_targ, ntargets, vel_targ, pot_targ, field_targ, mass_targ, nthreads, 
	times, dt, nsteps, d_pos_src, d_vel_src, d_pos_targ, d_vel_targ, d_pot_src, d_field_src,
	ifdirect, iffmm, iffm_multi, d_pot_targ, d_field_targ, tfmm, tdirect,
	ifpotsrc_far, ifpotsrc_near,
	ifpottarg_far, ifpottarg_near,
	iffldsrc_far, iffldsrc_near,
	iffldtarg_far, iffldtarg_near,
	pot_far_src, pot_near_src,
	pot_far_targ, pot_near_targ,
	field_far_src, field_near_src,
	field_far_targ,field_near_targ,
	save, outdir  ):
	  
  
  # dealing with the tree structure
  ladder	= np.zeros( (200,2), dtype=np.dtype('i4')) #200 is the maximum amount of levels
  box 		= np.zeros(20, dtype=np.dtype('i4'))
  center_temp 	= np.zeros(3, dtype=np.dtype('f8'))
  corners_temp	= np.zeros( (8,3), dtype=np.dtype('f8'))

  # arrays for the accelartion for the fmm
  accel_src		= np.zeros((nsources, 3), dtype=np.dtype('f8'))
  accel_targ		= np.zeros((ntargets, 3), dtype=np.dtype('f8'))
  accel_old_src		= np.zeros(accel_src.shape, accel_src.dtype)
  accel_old_targ	= np.zeros(accel_targ.shape, accel_targ.dtype)
  
  # arrays for the acceleration for the direct solver
  d_accel_src		= np.zeros(accel_src.shape, accel_src.dtype)
  d_accel_old_src	= np.zeros(accel_src.shape, accel_src.dtype)
  d_accel_targ		= np.zeros(accel_targ.shape, accel_targ.dtype)
  d_accel_old_targ 	= np.zeros(accel_targ.shape, accel_targ.dtype)
  
  # this reshaping must be done so that python can use charge to do matrix-vector multiplication later
  charge_src_reshape = np.reshape(charge_src, (nsources, 1))
  
  # set some flags
  if ifpot == 1 or iffld == 1:
    bsource = True
  else:
    bsource = False

  if ifpottarg == 1 or iffldtarg == 1:
    btarget = True
  else:
    btarget = False

  # these are the id's of the sources in tree-order (aka box-order)
  isource_tree	= np.zeros(nsources, dtype=np.dtype('i4'))
  itarget_tree	= np.zeros(ntargets, dtype=np.dtype('i4'))
  
  # set number of source/ box and set fmm tolerance
  nbox = su.set_sources_per_box(iprec, nsources, ntargets)  # number of sources per box  
  epsfmm = su.set_fmm_tolerance(iprec)  
  
  ifnear_only = 0
  
  #--------------------------------------------------
  # create the tree structure
  #--------------------------------------------------
    
  # Could be a problem with when the tree structure is made and how Verlocity Verlet works!
  t1 = time.time()
    
  # get the length of the tree structure object for wlists (ntot)
  ier, ntot = fortfmm.pre_create_tree(iprec, pos_src, nsources, pos_targ, ntargets, nbox, epsfmm)
  
  # create wlists (it's a workspace)
  # wlists_original 	= np.zeros(ntot, dtype=np.dtype('f8'))
  wlists = np.zeros(ntot, dtype=np.dtype('f8'))
  
  # create the tree structure
  ier, nlev, nboxes, lused7_original, iisource, iwlists, lwlists, iitarget, size = fortfmm.create_tree(iprec, 
							  pos_src, nsources, pos_targ, ntargets, center_temp, ladder, 
							  wlists, ntot, nbox, epsfmm)
  """
    ladder - an integer array dimensioned [nlev, 2], describing the
      numbers of boxes on various levels of sybdivision, so that
      the first box on level (i) has sequence number laddr[i,1],
      and there are laddr[i,2] boxes on level i.
      ladder starts counting from level 0
         
    nlev - the maximum level number on which any boxes have 
      been created. the maximum number possible is 200. 
      it is recommended that the array laddr above be 
      dimensioned at least (200,2), in case the user underestimates
      the number of levels required.
         
    size: size of the box on level 0
  """
  # wlists = np.array(wlists_original, dtype=np.dtype('f8'))
  
  t2 = time.time()  
  #print "In time it took to create the tree structure: " + str(t2-t1) + " [s]"
  
  # Store information about the tree structure
  boxes 	= np.zeros( (nboxes, 20), dtype=np.dtype('i4'))
  centers	= np.zeros( (nboxes, 3), dtype=np.dtype('f8'))
  corners	= np.zeros( (nboxes, 8, 3), dtype=np.dtype('f8'))
  
  # moved the following out of the loop
  # these are the id's of the sources in tree-order (aka box-order)
  #isource_tree	= np.zeros(nsources, dtype=np.dtype('i4'))
  #itarget_tree	= np.zeros(ntargets, dtype=np.dtype('i4'))
  
  mu.set_boxes(boxes, corners, centers, isource_tree, itarget_tree, wlists, iisource, iitarget, iwlists)
  
  # save the boxes data
  # set iter->0 and nsteps->0 (ie setpath(0,0)) to avoid making a new directory
  
  # fpath = os.path.dirname( os.getcwd() ) # this should be /fmmlib3d
  fpath = outdir + "/t_0"
  #fpath = data_path + "/t_" + str(i)
  mu.save_boxes(fpath, boxes, centers, corners, isource_tree, itarget_tree, ladder)
    
  #--------------------------------------------------------------------------------------------------
  # finished tree structure
  #---------------------------------------------------------------------------------------------------
  # Entering time loop
  
  # =========================================================================================
  for i in range(nsteps+1): # loop through nsteps + 1 so I can get the velocity for i==nsteps
    
    # set to the original wlists
    lused7	= lused7_original            
    
    mu.reset_pot_field(iffmm, ifdirect, ifpot, iffld, ifpottarg, iffldtarg,
		  ifpotsrc_far, ifpotsrc_near,
		  iffldsrc_far, iffldsrc_near,
		  ifpottarg_far, ifpottarg_near,
		  iffldtarg_far, iffldtarg_near,
		  pot_src, field_src, pot_targ, field_targ, 
		  d_pot_src, d_field_src, d_pot_targ, d_field_targ,
		  pot_far_src, field_far_src,
		  pot_near_src, field_near_src,
		  pot_far_targ, field_far_targ,
		  pot_near_targ, field_near_targ)
   
    
    # do the FMM
    
    if (iffmm == 1 or iffm_multi == 1):
      
      t1 = time.time()
      
      ier, lused7 = fortfmm.fmm(iprec, nsources, pos_src, ifcharge, charge_src, ifdipole, dipstr_src, dipvec_src,
	  ifpot, pot_src, iffld,field_src, ntargets, pos_targ, ifpottarg, pot_targ, iffldtarg, field_targ,
	  ladder, nlev, nboxes, nbox, 
	  epsfmm, lused7, nthreads, iisource, iwlists, lwlists, iitarget, size, wlists, 
	  ifpotsrc_far, ifpotsrc_near,
	  iffldsrc_far, iffldsrc_near,
	  ifpottarg_far, ifpottarg_near,
	  iffldtarg_far, iffldtarg_near,
	  ifnear_only,
	  pot_far_src,  field_far_src,
	  pot_near_src, field_near_src,
	  pot_far_targ, field_far_targ,
	  pot_near_targ, field_near_targ)           
      
      t2 = time.time()
      tfmm[i] = t2 -t1
      
      # sum the far and near field contributions
      pot_src 	= pot_far_src + pot_near_src
      pot_targ	= pot_far_targ + pot_near_targ
      
      field_src	= field_far_src + field_near_src
      field_targ = field_far_targ + field_near_targ                 

    # do the direct calculation
    if (ifdirect == 1):
      #print "\nDoing the direct calculation. nthreads = " + str(nthreads)
      
      t1 = time.time()
      
      fortfmm.direct(d_pos_src, ifcharge,charge_src, ifdipole, dipstr_src, dipvec_src, ifpot, d_pot_src, iffld, 
		      d_field_src, d_pos_targ, ifpottarg, d_pot_targ, iffldtarg, d_field_targ, nthreads)
      
      t2 = time.time()
      tdirect[i] = t2 - t1
           
    # calculate the new acceleration
    accel_src[:,:] = field_src.real[:,:] * charge_src.real[:]
    accel_targ[:,:] = field_targ.real[:,:] * 1.0
      
    # update velocity using FORWARD EULER
    mu.forwardEuler(vel_src, accel_src, dt)
    mu.forwardEuler(vel_targ, accel_targ, dt)
      
    # update position using FORWARD EULER
    mu.forwardEuler(pos_src, vel_src, dt)
    mu.forwardEuler(pos_targ, vel_targ, dt)
    
    if save is True:
      
      # set the path
      newpath, oldpath = mu.set_path(outdir, i, nsteps)
            
      mu.save_pot_field( oldpath, iffmm, ifdirect, ifpot, iffld, ifpottarg, iffldtarg, 
		    ifpotsrc_far, ifpotsrc_near,
		    iffldsrc_far, iffldsrc_near,
		    ifpottarg_far, ifpottarg_near,
		    iffldtarg_far, iffldtarg_near,
		    pot_src, pot_targ, field_src, field_targ,
		    d_pot_src, d_field_src, d_pot_targ, d_field_targ,
		    pot_far_src, pot_near_src, 
		    field_far_src, field_near_src,
		    pot_far_targ, field_far_targ, 
		    pot_near_targ, field_near_targ,
		    accel_src, accel_targ, d_accel_src, d_accel_targ)
     
      # save the velocity (v_n). Do not save velocity for i == 0. Already known.
      if i > 0:
	mu.save_vel(oldpath, vel_src, vel_targ, d_vel_src, d_vel_targ)
	mu.save_pos(oldpath, pos_src, pos_targ, d_pos_src, d_pos_targ)
    
  # =========================================================================================
  
  """------------------- Outside the time loop. Need to update the velocities ------------------""" 
  
  # testing
  #printStuff(iisource, iitarget, iwlists, lwlists, lused7, nbox, nboxes, ladder, size, epsfmm, nlev)
  
  if save is True:   
    # save potentials and fields
    mu.save_pot_field( oldpath, iffmm, ifdirect, ifpot, iffld, ifpottarg, iffldtarg, 
		    ifpotsrc_far, ifpotsrc_near,
		    iffldsrc_far, iffldsrc_near,
		    ifpottarg_far, ifpottarg_near,
		    iffldtarg_far, iffldtarg_near,
		    pot_src, pot_targ, field_src, field_targ,
		    d_pot_src, d_field_src, d_pot_targ, d_field_targ,
		    pot_far_src, pot_near_src, 
		    field_far_src, field_near_src,
		    pot_far_targ, field_far_targ, 
		    pot_near_targ, field_near_targ,
		    accel_src, accel_targ, d_accel_src, d_accel_targ)
    
  return ier, nlev, nboxes, nbox, lused7