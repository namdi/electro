""" Python wrapper for some fortran FMM routines.
This module acts as a wrapper for the fortran FMM routines from leslie Greengard's FMM.

"""
#============================================
# Import
#============================================
from ctypes import *
import numpy as np
from numpy.ctypeslib import ndpointer

#import data_dir
import os

#============================================
# load libraries
#============================================

# creates a python library/module called fort to access function in _my_fortfmmlib.so
# when _my_fortfmmlib.so is made with a make file, one does not need the initial _ in _my_fortfmmlib.so

obj_path	= os.environ["PY_OBJ_DIR"]
lib_name 	= obj_path + "/_my_fortfmmlib.so"

fort = CDLL(lib_name)

#Notice: That I do not need to add the "./" in CDLL if I add the entire path.
#The "./" adds the entire path to the shared object file name
#fort = CDLL("./_my_fortfmmlib.so")

#============================================
# function prototypes
#============================================

float64 = POINTER(c_double)
integer = POINTER(c_int)
# pass complex arrays as float64 (in memory a complex array is just the real part of the array followed
# by the complex part. It's an array twice as long.

fort.pre_create_tree_wrap.restype  = None
fort.pre_create_tree_wrap.argtypes = [c_int, ndpointer(dtype=np.dtype('f8')), c_int, ndpointer(dtype=np.dtype('f8')), 
				      c_int, ndpointer(dtype=np.dtype('i4')), c_int, c_int, c_double] 

fort.create_tree_wrap.restype  = None
fort.create_tree_wrap.argtypes = [ c_int, ndpointer(dtype=np.dtype('f8')), c_int, ndpointer(dtype=np.dtype('f8')), 
				c_int, ndpointer(dtype=np.dtype('f8')), ndpointer(dtype=np.dtype('i4')),
				ndpointer(dtype=np.dtype('i4')), c_int, ndpointer(dtype=np.dtype('f8')),
				c_int, ndpointer(dtype=np.dtype('f8')), c_int, c_int, c_double ]
				
fort.fmm_wrap.restype  = None
fort.fmm_wrap.argtypes = [ c_int, c_int, ndpointer(dtype=np.dtype('f8')), c_int, ndpointer(dtype=np.dtype('c16')),
			  c_int, ndpointer(dtype=np.dtype('c16')), ndpointer(dtype=np.dtype('f8')), c_int, 
			  ndpointer(dtype=np.dtype('c16')), c_int, ndpointer(dtype=np.dtype('c16')), #12
			  c_int, ndpointer(dtype=np.dtype('f8')), c_int, ndpointer(dtype=np.dtype('c16')),
			  c_int, ndpointer(dtype=np.dtype('c16')), ndpointer(dtype=np.dtype('i4')),
			  ndpointer(dtype=np.dtype('i4')), c_int, ndpointer(dtype=np.dtype('f8')), #22
			  c_int, ndpointer(dtype=np.dtype('f8')), c_int, #25
			  c_int, c_int, 
			  c_int, c_int,
			  c_int, c_int,
			  c_int, c_int, # 33
			  ndpointer(dtype=np.dtype('c16')), ndpointer(dtype=np.dtype('c16')), 
			  ndpointer(dtype=np.dtype('c16')), ndpointer(dtype=np.dtype('c16')),
			  ndpointer(dtype=np.dtype('c16')), ndpointer(dtype=np.dtype('c16')),
			  ndpointer(dtype=np.dtype('c16')), ndpointer(dtype=np.dtype('c16')), 	
			  c_int ] # 42

fort.direct_wrap.restype  = None
fort.direct_wrap.argtypes = [ ndpointer(dtype=np.dtype('i4')), c_int, ndpointer(dtype=np.dtype('f8')), 
			      ndpointer(dtype=np.dtype('c16')), 
			      ndpointer(dtype=np.dtype('c16')), ndpointer(dtype=np.dtype('f8')),
			      ndpointer(dtype=np.dtype('c16')), ndpointer(dtype=np.dtype('c16')),
			      ndpointer(dtype=np.dtype('f8')), ndpointer(dtype=np.dtype('c16')),
			      ndpointer(dtype=np.dtype('c16'))]
			      
fort.get_box_info_wrap.restype  = None
fort.get_box_info_wrap.argtypes = [  c_int, ndpointer(dtype=np.dtype('i4')), ndpointer(dtype=np.dtype('f8')), 
				    ndpointer(dtype=np.dtype('f8')), ndpointer(dtype=np.dtype('f8')), 
				    ndpointer(dtype=np.dtype('i4')), c_int, c_int ]
				    
fort.get_winfo_wrap.restype  = None
fort.get_winfo_wrap.argtypes = [ndpointer(dtype=np.dtype('i4')), ndpointer(dtype=np.dtype('f8')), c_int, c_int]

#=========================================================
# function declarations			  
#=========================================================

def pre_create_tree(iprec, pos_src, nsources, pos_targ, ntarget, nbox, epsfmm):
   
  idata_len = 2
  
  # idata[0] = ier	idata[1] = ntot
  
  idata = np.zeros(idata_len, dtype=np.dtype('i4'))
 
  # get value ntot from the tree
  fort.pre_create_tree_wrap(iprec, pos_src, nsources, pos_targ, ntarget, idata, idata_len, nbox, epsfmm)
			
  # unwrap the integer data
  ier		= idata[0]
  ntot 		= idata[1]
  
  return ier, ntot
			  
def create_tree(iprec, pos_src, nsources, pos_targ, ntarget, center, ladder, wlists, ntot, nbox, epsfmm):
  
  # idata is a hack to abstract integer information from the fortran functions
  # idata[0] = ier	idata[1] = nlev	idata[2] = nboxes	idata[3] = lused7
  
  idata_len = 9
  ddata_len = 2
  
  # idata[0] = ier	idata[1] = nlev		idata[2] = nboxes	idata[3] = nbox	
  ## idata[4] = lused7	idata[5] = iisourced	idata[6] = iwlists	idata[7] = lwlists
  # idata[8] = iitarget
  
  # ddata[0] = epsfmm	ddata[1] = size
  
  idata = np.zeros(idata_len, dtype=np.dtype('i4'))
  ddata = np.zeros(ddata_len, dtype=np.dtype('f8'))
  
  # create the tree
  fort.create_tree_wrap(iprec, pos_src, nsources, pos_targ, ntarget, center, ladder, idata, idata_len, ddata, 
			  ddata_len, wlists, ntot, nbox, epsfmm)
			
  # unwrap the integer data
  ier		= idata[0]
  nlev 		= idata[1]
  nboxes 	= idata[2]
  #nbox	 	= idata[3]
  lused7 	= idata[4]
  iisource 	= idata[5]
  iwlists 	= idata[6]
  lwlists 	= idata[7]
  iitarget 	= idata[8]
  
  #epsfmm = ddata[0]
  size   = ddata[1]
    
  return ier, nlev, nboxes, lused7, iisource, iwlists, lwlists, iitarget, size # create nbox, espfmm

"""
  This is the FMM. It assumes the tree has alread y been created.
"""

def fmm(iprec, nsources, pos_src, ifcharge, charge, ifdipole, dipstr, dipvec,ifpot, pot_src, iffld, 
	field_src, ntarget, pos_targ, ifpottarg, pot_targ, iffldtarg, field_targ, ladder, nlev, nboxes, nbox, 
	epsfmm, lused7, nthreads, iisource, iwlists, lwlists, iitarget, size, wlists,
	ifpotsrc_far, ifpotsrc_near,
	iffldsrc_far, iffldsrc_near,
	ifpottarg_far,ifpottarg_near,
	iffldtarg_far,iffldtarg_near, ifnear_only,
	pot_far_src,  field_far_src,  pot_near_src,  field_near_src,
	pot_far_targ, field_far_targ, pot_near_targ, field_near_targ):
  
  #print "ifnear_only (fortfmm.py): " + str(ifnear_only)
  idata_len = 10
  ddata_len = 2
  
  idata = np.zeros(idata_len, dtype=np.dtype('i4'))
  ddata = np.zeros(ddata_len, dtype=np.dtype('f8'))

  # idata[1] = ier	idata[1] = nlev	idata[2] = nboxes	idata[3] = nbox
  # idata[4] = lused7	idata[5] = nthread	idata[6] = iisource	idata[7] = iwlists
  # idata[8] = lwlists	idata[9] = iitarget

  # ddata[0] = epsfmm	ddata[1] = size
  
  idata[1] = nlev
  idata[2] = nboxes
  idata[3] = nbox
  idata[4] = lused7
  idata[5] = nthreads
  idata[6] = iisource
  idata[7] = iwlists
  idata[8] = lwlists
  idata[9] = iitarget
  
  ddata[0] = epsfmm
  ddata[1] = size
  
  ntot = len(wlists)
    
  #print "In fortfmm.py, before fort.fmm_wrap: field_near_src"
  
  fort.fmm_wrap(iprec, nsources, pos_src, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, pot_src, #10
		iffld, field_src, ntarget, pos_targ, ifpottarg, pot_targ, iffldtarg, field_targ,
		ladder, idata, idata_len, ddata, ddata_len, wlists, ntot, #25
		ifpotsrc_far,	ifpotsrc_near,
		iffldsrc_far,	iffldsrc_near,
		ifpottarg_far, 	ifpottarg_near,
		iffldtarg_far, 	iffldtarg_near, #33
		pot_far_src,  field_far_src,
		pot_near_src, field_near_src,
		pot_far_targ, field_far_targ,
		pot_near_targ, field_near_targ, ifnear_only ) # 42
		
  #print "In fortfmm.py, after fort.fmm_wrap: field_src"
  #print field_near_src
  
  # unwrap the integer output data
  ier = idata[0]
  lused7 = idata[4]
  
  return ier, lused7

"""
  This is the direct method
"""

def direct(pos_src, ifcharge, charge_src, ifdipole, dipstr, dipvec, ifpot, pot_src, iffld, field_src, 
	   pos_targ, ifpottarg, pot_targ, iffldtarg, field_targ, nthreads):
  
  idata_len = 9
  idata = np.zeros( idata_len, dtype=np.dtype('i4') )
  
  nsources = pos_src.shape[0]
  ntargets = pos_targ.shape[0]
  
  idata[0] = ifcharge
  idata[1] = ifdipole
  idata[2] = ifpot
  idata[3] = iffld
  idata[4] = ifpottarg
  idata[5] = iffldtarg
  idata[6] = nthreads
  idata[7] = nsources
  idata[8] = ntargets
    
  fort.direct_wrap(idata, idata_len, pos_src, charge_src, dipstr, dipvec, pot_src, field_src, 
		   pos_targ, pot_targ, field_targ)
		   
  #fort.direct_wrap(nsources, pos_src, ifcharge, charge_src, ifdipole, dipstr, dipvec,
			  #ifpot, pot_src, iffld, field_src, ntargets, pos_targ, ifpottarg, pot_targ,
			  #iffldtarg, fieldtarg, nthreads)
  
  return

"""
  This function takes out the information from boxes. From the fortran code.
"""

def get_box_info(ibox, box, center, corners, w, iwlists):
  """ This function takes in a box index and the workspace w and returns information
    about the box ibox. The following information is returned: 
    
    Input: ibox, w
    output: ier, box, center, conrners

        this entry returns to the user the characteristics of
        user-specified box  ibox.  

                     input parameters:

    ibox - the box number for which the information is desired
    w - storage area as created by the entry d3tstrcr (see above)

                     output parameters:

    ier - the error return code.
      ier=0 means successful execution
      ier=2 means that ibox is either greater than the number of boxes
           in the structure or less than 1.
    box - an integer array dimensioned box(20). its elements describe 
        the box number ibox, as follows:

      1. level - the level of subdivision on which this box 
             was constructed; 
      2, 3, 4  - the coordinates of this box among  all
             boxes on this level 
	      ! namdi: (x,y,z) coordinates? How is this an integer. 
	    Maybe 1 corresponds to 1 box over in a direction from an origin box?
	    !The units are by boxes themselves
      5 - the daddy of this box, identified by it address
             in array "boxes"
      6,7,8,9,10,11,12,13 - the  list of children of this box 
             (eight of them, and the child is identified by its address
             in the array "boxes"; if a box has only one child, only the
             first of the four child entries is non-zero, etc.)
      14 - the location in the array iz of the particles 
             living in this box
      15 - the number of particles living in this box
      16 - the location in the array iztarg of the targets
             living in this box
      17 - the number of targets living in this box
      18 - source box type: 0 - empty, 1 - leaf node, 2 - sub-divided
      19 - target box type: 0 - empty, 1 - leaf node, 2 - sub-divided
      20 - reserved for future use

    center - the center of the box number ibox 
    corners - the corners of the box number ibox 

    . . . return to the user all information about the box ibox
  """
  idata_len = 1
  ier = 0
  
  idata = np.zeros(idata_len, dtype=np.dtype('i4'))
  idata[0] = ier
  
  #print "In get_box_info() of fortfmm.py. About to enter get_box_info_wrap()" 
  fort.get_box_info_wrap(ibox, box, center, corners, w, idata, idata_len, iwlists)
  
  ier = idata[0]
  return ier

"""
  This function gets in information from wlists from fortran.
"""

def get_winfo(iarray, wlists, idx, n):
  """ this function gets the selected information from wlists (which is a linked list in fortran)
   Input: 
   idx: the array index for the selected info in wlists
   n:	the length of the array to be returned
   iarray: the data to be returned
  """
  fort.get_winfo_wrap(iarray, wlists, idx, n)
  
  return