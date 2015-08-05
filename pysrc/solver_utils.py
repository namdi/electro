#===========================================
# Explaination
#===========================================

"""
The purpose of this code is to hold functions for the the solvers (multirate FMM, regular FMM, direct).

fmm() runs the full FMM and the multirate FMM
direct() runs the direct solver
"""

#---------------------------
# Import
#---------------------------
import numpy as np
import time # testing

import fortfmm	
import my_utils as mu
#---------------------------

#===========================================
# Function Definitions
#===========================================
def set_sources_per_box(accuracy, nsources, ntargets):
  
  nbox = 0
  
  if (accuracy == -2): nbox = 40
  if (accuracy == -1): nbox = 50
  if (accuracy == 0): nbox = 80
  if (accuracy == 1): nbox = 160
  if (accuracy == 2): nbox = 400
  if (accuracy == 3): nbox = 800
  if (accuracy == 4): nbox = 1200
  if (accuracy == 5): nbox = 1400
  if (accuracy == 6): nbox = nsources + ntargets
    
  return nbox
  
def set_fmm_tolerance(accuracy): 
  
  epsfmm = 0.0
  
  if (accuracy == -2): epsfmm = 0.5
  if (accuracy == -1): epsfmm = 0.5 * 10**-1
  if (accuracy == 0): epsfmm = 0.5 * 10**-2
  if (accuracy == 1): epsfmm = 0.5 * 10**-3
  if (accuracy == 2): epsfmm = 0.5 * 10**-6
  if (accuracy == 3): epsfmm = 0.5 * 10**-9
  if (accuracy == 4): epsfmm = 0.5 * 10**-12
  if (accuracy == 5): epsfmm = 0.5 * 10**-15
  if (accuracy == 6): epsfmm = 0.0
    
  return epsfmm  
  
  
# the direct solver
def direct(pos_src, pos_targ, charge_src, dipstr_src, dipvec_src, pot_src, field_src, pot_targ, field_targ, 
	dosource, dotarget, dopot, dofield, docharge=True, dodipole=False, nthreads=1):
	
  iffld 	= 0
  ifpot 	= 0
  iffldtarg 	= 0
  ifpottarg 	= 0
  
  # set flags	   
  if dosource:
    if dofield:
      iffld = 1
      
    if dopot:
      ifpot = 1
      
  if dotarget:
    if dofield:
      iffldtarg = 1
      
    if dopot:
      ifpottarg = 1
  
  # do the direct calculation
  fortfmm.direct(pos_src, int(docharge), charge_src, int(dodipole), dipstr_src, dipvec_src, ifpot, pot_src, iffld, field_src, 
	   pos_targ, ifpottarg, pot_targ, iffldtarg, field_targ, nthreads)

  return

#  Greengard's FMM

def fmm(pos_src, pos_targ, charge_src, dipstr_src, dipvec_src, iprec,
	pot_src, field_src, pot_targ, field_targ,
	pot_far_src, 	pot_near_src,  	field_far_src,	field_near_src,
	pot_far_targ, 	pot_near_targ, 	field_far_targ,	field_near_targ,
	dosource, dotarget, dofmm, dofmm_multi, dofmm_near,
	dopot, dofield, docharge=True, dodipole=False, nthreads=1):

   # set flags for the potential
  (ifpot, ifpottarg, ifpotsrc_far, ifpotsrc_near, ifpottarg_far, ifpottarg_near, ifnear_only) = mu.set_flags_pot(dopot, dofmm, dofmm_multi, dofmm_near, False, dosource, dotarget)
  
  # set flags for the electric-field
  (iffld, iffldtarg, iffldsrc_far, iffldsrc_near, iffldtarg_far, iffldtarg_near, ifnear_only) = mu.set_flags(dofmm, dofmm_multi, dofmm_near, False, dosource, dotarget)  
  
  # set flags
  ifcharge = int(docharge)
  ifdipole = int(dodipole)
  
  nsources = pos_src.shape[0]
  ntargets = pos_targ.shape[0]  
    
  # set number of source/ box and set fmm tolerance
  nbox = set_sources_per_box(iprec, nsources, ntargets)  # number of sources per box  
  epsfmm = set_fmm_tolerance(iprec)  
  
  # testing to see the far field calculation only
  #iffldsrc_near = 0
  
  # testing to see the near field calculation only
  #iffldsrc_far = 0
  
  # call the FMM  
  ier, nlev, nboxes, nbox, lused7 = fmm_run(iprec, ifpot, iffld, ifcharge, ifdipole, ifpottarg, iffldtarg, 
	pos_src, nsources, charge_src, dipstr_src, dipvec_src, pot_src, field_src, 
	pos_targ, ntargets, pot_targ, field_targ, nthreads, nbox, epsfmm,
	ifpotsrc_far, ifpotsrc_near,
	iffldsrc_far, iffldsrc_near,
	ifpottarg_far, ifpottarg_near,
	iffldtarg_far, iffldtarg_near, ifnear_only,
	pot_far_src,  field_far_src,  pot_near_src,  field_near_src,
	pot_far_targ, field_far_targ, pot_near_targ, field_near_targ)
	
  return ier, nlev, nboxes, nbox, lused7
  
# fmm_run() actually runs the FMM

def fmm_run(iprec, ifpot, iffld, ifcharge, ifdipole, ifpottarg, iffldtarg, 
	pos_src, nsources, charge_src, dipstr_src, dipvec_src, pot_src, field_src, 
	pos_targ, ntargets, pot_targ, field_targ, nthreads, nbox, epsfmm,
	ifpotsrc_far,  ifpotsrc_near,
	iffldsrc_far,  iffldsrc_near,
	ifpottarg_far, ifpottarg_near,
	iffldtarg_far, iffldtarg_near, ifnear_only,
	pot_far_src,  field_far_src,  pot_near_src,  field_near_src,
	pot_far_targ, field_far_targ, pot_near_targ, field_near_targ):
	  
  
  # dealing with the tree structure
  ladder	= np.zeros( (200,2), dtype=np.dtype('i4')) #200 is the maximum amount of levels
  box 		= np.zeros(20, dtype=np.dtype('i4'))
  center_temp 	= np.zeros(3, dtype=np.dtype('f8'))
  corners_temp	= np.zeros( (8,3), dtype=np.dtype('f8'))
  
  # this reshaping must be done so that python can use charge to do matrix-vector multiplication later
  charge_src_reshape = np.reshape(charge_src, (nsources, 1))
   
  # these are the id's of the sources in tree-order (aka box-order)
  isource_tree	= np.zeros(nsources, dtype=np.dtype('i4'))
  itarget_tree	= np.zeros(ntargets, dtype=np.dtype('i4'))
   
  #--------------------------------------------------
  # create the tree structure
  #--------------------------------------------------
  
    
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
  
  # Store information about the tree structure
  boxes 	= np.zeros( (nboxes, 20), dtype=np.dtype('i4'))
  centers	= np.zeros( (nboxes, 3), dtype=np.dtype('f8'))
  corners	= np.zeros( (nboxes, 8, 3), dtype=np.dtype('f8'))
  
  mu.set_boxes(boxes, corners, centers, isource_tree, itarget_tree, wlists, iisource, iitarget, iwlists)
  
  # save the boxes data
  # set iter->0 and nsteps->0 (ie setpath(0,0)) to avoid making a new directory
  
  # fpath = os.path.dirname( os.getcwd() ) # this should be /fmmlib3d
  # fpath = fpath + "/mydata/t_0/"
    
  #--------------------------------------------------------------------------------------------------
  # finished tree structure
  #---------------------------------------------------------------------------------------------------
    
  # set to the original wlists
  lused7	= lused7_original
    
  ier, lused7 = fortfmm.fmm(iprec, nsources, pos_src, ifcharge, charge_src, ifdipole, dipstr_src, dipvec_src,ifpot, pot_src, iffld, 
	  field_src, ntargets, pos_targ, ifpottarg, pot_targ, iffldtarg, field_targ, ladder, nlev, nboxes, nbox, 
	  epsfmm, lused7, nthreads, iisource, iwlists, lwlists, iitarget, size, wlists, 
	  ifpotsrc_far, ifpotsrc_near,
	  iffldsrc_far, iffldsrc_near,
	  ifpottarg_far, ifpottarg_near,
	  iffldtarg_far, iffldtarg_near, ifnear_only,
	  pot_far_src,  field_far_src,
	  pot_near_src, field_near_src,
	  pot_far_targ, field_far_targ,
	  pot_near_targ, field_near_targ)   
  
  return ier, nlev, nboxes, nbox, lused7
  
  