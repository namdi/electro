"""
File: my_utils.py
"""

#-----------------------------------------
# Import 
#----------------------------------------

from optparse import OptionParser
import numpy as np
import time
import os

import fortfmm
#import data_dir
"""
The following parses the command line for the arguments:
  nsources, ntargets, nthreads, distribution, accuracy
"""

#============================================
# Function Declarations
#============================================

#----------------------------------------
# Command Line functions
#---------------------------------------

# parse the command line
def parse_cmd_line():
  
  parser = OptionParser()
  
  parser.add_option("-n", dest = "nsources", help="number of sources", metavar="INT")
  
  parser.add_option("-m", dest = "ntargets", help="number of targets", metavar="INT")
  
  parser.add_option("-t", dest = "nthreads", help="number of threads", metavar="INT")
  
  parser.add_option("-d", dest = "distribution", help="the initial distribution", metavar="INT")
  
  parser.add_option("-a", dest = "accuracy", help="The accuracy of the FMM solver.",
		    metavar = "INT")

  parser.add_option("-b", dest="nsteps", help="The number of timesteps", metavar="INT")

  parser.add_option("-o", dest="outdir", help="The directory in which the data is saved")
  
  parser.add_option("-i", dest="indir", help="The directory in which the initial data is made")
  
  parser.add_option("-s", dest="solver", help="The solver to solve the spatial problem")
  
  parser.add_option("-k", dest="niters", help="the number of PFASST iterations")
  
  (options, args) = parser.parse_args()
  
  # if an argument was not passed in assign default values
  if options.nsources is None:
    nsources = 10000
  else:
    nsources = int(options.nsources)  
  
  if options.ntargets is None:
    ntargets = 4
  else:
    ntargets = int(options.ntargets)
  
  if options.nthreads is None:
    nthreads = 1
  else:
    nthreads = int(options.nthreads)  
    
  if options.distribution is None:
    distribution = 2
  else:
    distribution = int(options.distribution)
  
  # read the precision of the FMM for calculating the potential
  if options.accuracy is None:
    accuracy = 3
  else:
   accuracy = int(options.accuracy)  
  
  tol = set_tol(accuracy)
  
  if  tol == -1 :
    raise NotImplementedError, 'need to enter a valid value for iprec: -2, -1, 0, 1, 2, 3, 4, 5'
    
  if options.nsteps is None:
    nsteps = 1
  else:
    nsteps = int(options.nsteps)
  
  if options.outdir is None:
    outdir = None
  else:
    #outdir = data_dir.data_path + "/" + str(options.outdir)
    outdir = str(options.outdir)
    
  if options.indir is None:
    indir = None
  else:
    #indir = data_dir.data_path + "/" + str(options.indir)
    indir = str(options.indir)
  if options.solver is None:
    solver = 2
  else:
    solver = int(options.solver)
    
    if solver_check(solver) is False:
      raise NotImplementedError, "need to enter a value of a valid solver"
    
  if options.niters is None:
    niters = 3
  else:
    niters = int(options.niters)
    
  return nsources, ntargets, nthreads, distribution, accuracy, nsteps, indir, outdir, solver, niters

#===============================================
#  autorun.py functions
#===============================================

# this function checks to make sure the solver has the right values
def solver_check(solver):
  
  # valid values for solver
  num_solvers = 4
  vals = np.linspace(0, num_solvers, num_solvers+1)
  
  check = (vals == solver)
  
  ok = False
  
  # we have a match of a valid solver
  if sum(check) == 1:
    ok = True
  
  return ok
    

def set_flags(dofmm, dofmm_multi, dofmm_near, dodirect, dosource, dotarget):
  
  iffld 	= 0
  iffldtarg 	= 0
  
  iffldsrc_far		= 0
  iffldsrc_near		= 0
  
  iffldtarg_far		= 0  
  iffldtarg_near	= 0
  
  ifnear_only		= 0
  
  if dofmm_near is True:
    iffld, iffldtarg, iffldsrc_far, iffldsrc_near, iffldtarg_far, iffldtarg_near = set_fmm_near_flags(dosource, dotarget)
    ifnear_only = 1
    
  if dofmm is True:
    iffld, iffldtarg, iffldsrc_far, iffldsrc_near, iffldtarg_far, iffldtarg_near = set_fmm_flags(dosource, dotarget)
    
  elif dofmm_multi is True:
    iffld, iffldtarg, iffldsrc_far, iffldsrc_near, iffldtarg_far, iffldtarg_near = set_fmm_multi_flags(dosource, dotarget)
    
  elif dodirect is True:
    iffld, iffldtarg, iffldsrc_far, iffldsrc_near, iffldtarg_far, iffldtarg_near = set_direct_flags(dosource, dotarget)
    
  return iffld, iffldtarg, iffldsrc_far, iffldsrc_near, iffldtarg_far, iffldtarg_near, ifnear_only
  

def set_fmm_multi_flags(dosource, dotarget):
  
  # turn off flags for regular fmm
  iffld		= 0
  iffldtarg	= 0
  
  iffldsrc_far		= int(dosource)  
  iffldsrc_near		= int(dosource)
    
  iffldtarg_far		= int(dotarget)  
  iffldtarg_near	= int(dotarget)
  
  return iffld, iffldtarg, iffldsrc_far, iffldsrc_near, iffldtarg_far, iffldtarg_near
  
def set_fmm_near_flags(dosource, dotarget):
  
  # turn off flags for regular fmm
  iffld		= 0
  iffldtarg	= 0
  
  iffldsrc_far		= 0  
  iffldsrc_near		= int(dosource)
    
  iffldtarg_far		= 0
  iffldtarg_near	= int(dotarget)
  
  return iffld, iffldtarg, iffldsrc_far, iffldsrc_near, iffldtarg_far, iffldtarg_near
    
def set_fmm_flags(dosource, dotarget):
  # flags for direct calculations
  iffld 	= int(dosource)
  iffldtarg 	= int(dotarget)
    
  iffldsrc_far 		= 0  
  iffldsrc_near		= 0
    
  iffldtarg_far		= 0
  iffldtarg_near	= 0
  
  return iffld, iffldtarg, iffldsrc_far, iffldsrc_near, iffldtarg_far, iffldtarg_near

def set_direct_flags(dosource, dotarget):
  return set_fmm_flags(dosource, dotarget)
  
def set_flags_pot(dopot, dofmm, dofmm_multi, dofmm_near, dodirect, dosource, dotarget):
  
  # falgs for regular FMM
  ifpot 	= 0
  ifpottarg 	= 0
  
  # these are the multriate FMM flags for source potentials
  ifpotsrc_far 		= 0  
  ifpotsrc_near		= 0
   
  # these are the multirate FMM flags for target potentials
  ifpottarg_far		= 0
  ifpottarg_near	= 0
  
  ifnear_only		= 0
  
  if dopot is True:
    if (dosource is True) and ( (dofmm is True) or (dodirect is True) ):
      ifpot = 1
    if (dotarget is True) and ( (dofmm is True) or (dodirect is True) ):
      ifpottarg = 1

    if dosource and dofmm_multi:
      ifpotsrc_far = 1      
      ifpotsrc_near = 1

    if dotarget and dofmm_multi:
      ifpottarg_far = 1      
      ifpottarg_near = 1
      
    if dosource and dofmm_near:
      ifpotsrc_far = 0      
      ifpotsrc_near = 1

    if dotarget and dofmm_near:
      ifpottarg_far = 0      
      ifpottarg_near = 1
      
    ifnear_only = int(dofmm_near)
    
  return ifpot, ifpottarg, ifpotsrc_far, ifpotsrc_near, ifpottarg_far, ifpottarg_near, ifnear_only
  
"""" This program assigns the flags that determine which solver will be used for the spatial calculations"""
def get_solver(solver):
  fmm_near	= False
  fmm_multi	= False
  fmm 		= False
  direct 	= False
  
  if solver == 1:
    direct = True
  if solver == 2:
    fmm = True
  elif solver == 3 or solver==4:
    fmm_multi = True


  return fmm_near, fmm_multi, fmm, direct

#----------------------------------------
# Loading Functions
#---------------------------------------

""" I will neeed to change this to make it work assuming that there are no targets and to have different
file names for data"""

# load the intiial data for autorun.py
def load_initial_data(fpath):
  pos_src	= np.loadtxt(fpath + "/pos_source.txt").view('f8')
  vel_src	= np.loadtxt(fpath + "/vel_source.txt").view(pos_src.dtype)
  
  pos_targ	= np.loadtxt(fpath + "/pos_target.txt").view(pos_src.dtype)
  vel_targ	= np.loadtxt(fpath + "/vel_target.txt").view(pos_src.dtype)
  
  nsources = pos_src.shape[0]
  ntargets = pos_targ.shape[0]
  
  # the initial data ( 2 -> position/velocity, nparticles-> [nsources + ntargets], 3 spatial dimensions )
  q0	= np.zeros( (2, nsources + ntargets, 3), pos_src.dtype)
  
  q0[0, :nsources, :] 	= pos_src[:,:]
  q0[0, nsources:, :]	= pos_targ[:,:]
  q0[1, :nsources, :]	= vel_src[:,:]
  q0[1, nsources:, :]	= vel_targ[:,:]
  
  # load the source charges and the target charges
  charge_src 	= np.loadtxt(fpath + "/charge.txt").view('c16')
  charge_targ 	= np.ones ( (ntargets,1), dtype=np.dtype('c16') )
  
  dipstr_src	= np.loadtxt(fpath + "/dipstr.txt").view('c16')
  dipvec_src	= np.loadtxt(fpath + "/dipvec.txt").view('f8')
  
  return q0, charge_src, charge_targ, dipstr_src, dipvec_src  
  
def load_parameters(fname):
  f = open(fname, "r")
  
  idx = len("nsources: ")
  line = f.readline()
  nsources = int(line[idx:-1])
  
  idx = len("ntargets: ")
  line = f.readline()
  ntargets = int(line[idx:-1])
  
  idx = len("distribution: ")
  line = f.readline()
  distribution = int(line[idx:-1])
  
  idx = len("dt: ")
  line = f.readline()
  dt = float(line[idx:-1])
  
  idx = len("nsteps: ")
  line = f.readline()
  nsteps = int(line[idx:-1])
  
  f.close()
  
  return  nsources, ntargets, distribution, dt, nsteps

def load_pfasst_parameters(fname):
  f 	= open(fname, "r")
  
  idx 		= len("iterations: ")
  line		= f.readline()
  iterations 	= int(line[idx:-1])
  
  idx 		= len("nnodes: ")
  line		= f.readline()
  nnodes 	= int(line[idx:-1])
  
  idx 		= len("levels: ")
  line		= f.readline()
  levels 	= int(line[idx:-1])
  
  idx 		= len("refine: ")
  line		= f.readline()
  refine 	= str(line[idx:-1])
  
  idx 		= len("mpi size: ")
  line		= f.readline()
  mpi_size 	= int(line[idx:-1])
  
  idx 		= len("omp size: ")
  line		= f.readline()
  omp_size 	= int(line[idx:-1])
  
  # skip 2 lines
  f.readline()
  f.readline()
  
  idx 		= len("iprec0 ")
  line		= f.readline()
  iprec0 	= int(line[idx:-1])
  
  idx 		= len("iprec1 ")
  line		= f.readline()
  iprec1 	= int(line[idx:-1])
  
  f.close()
  
  return iterations, nnodes, levels, refine, mpi_size, omp_size, iprec0, iprec1
  
#----------------------------------------
# Integration Functions
#---------------------------------------
  
def forwardEuler(pos, vel, dt):
  pos[:,:] = pos[:,:] + vel[:,:] * dt
  return
  
"""The following integration methods use the Velocity Verlet Formula """
def updatePos_verlet(pos, vel, accel, dt):
  pos[:,:] = pos[:,:] + vel[:,:] * dt + 0.5 * accel[:,:] * dt**2
  return
  
def updateVel_verlet(vel, accel, accel_old, dt):
  vel[:,:] = vel[:,:] + 0.5 * dt * (accel[:,:] + accel_old[:,:])
  return

def updatePos(bsource, btarget, pos_src, vel_src, accel_src, pos_targ, vel_targ, accel_targ, dt):  
  if bsource is True:
    updatePos_verlet(pos_src, vel_src, accel_src, dt)
  
  if btarget is True:
    updatePos_verlet(pos_targ, vel_targ, accel_targ, dt)
  return

def updateVel(bsource, btarget, vel_src, accel_src, accel_old_src, vel_targ, accel_targ,
		    accel_old_targ, dt):
  
  if bsource is True:
    updateVel_verlet(vel_src, accel_src, accel_old_src, dt)
    
  if btarget is True:
    updateVel_verlet(vel_targ, accel_targ, accel_old_targ, dt)
  return

def updateAccel(bsource, btarget, accel_src, accel_old_src, field_src, accel_targ, accel_old_targ, field_targ, charge_src):
  # calculate the acclearation a = qE/m
  # assume that mass=1.0
  # assume that the charge of the targets is 1.0
  # charge must be in a column vector as in charge.shape is (n,1)
      
  if bsource is True:
    # store the old acceleartion      
    accel_old_src[:,:] = accel_src[:,:]
	
    # calculate the new acceleration
    accel_src[:,:] = field_src.real[:,:] * charge_src.real[:]
	
  if btarget is True:
    # store the old accleration
    accel_old_targ[:,:] = accel_targ[:,:]
	
    # calculate the new acceleration
    accel_targ[:,:] = field_targ.real[:,:] * 1.0
    
  return

#----------------------------------------
# FMM tree functions
#---------------------------------------
  
def set_boxes(boxes, corners, centers, isource_tree, itarget_tree, wlists, iisource, iitarget, iwlists):
  """ 
  This function reads in the data to about boxes from wlist.
  It stores in information about all the boxes.
  It fills in arrays that show where the target particle and the source particles are (tree-structure formation)
  """
  nboxes 	= boxes.shape[0]
  nsources	= len(isource_tree)
  ntargets	= len(itarget_tree)
  
  box 		= np.zeros(20, dtype=np.dtype('i4'))
  center_temp	= np.zeros(3, dtype=np.dtype('f8'))
  corners_temp	= np.zeros( (8,3), dtype=np.dtype('f8'))
  
  # fill iz and iztarg. These are the arrays that show where the target particles and the source
  # particles are in the tree structure
  fortfmm.get_winfo(isource_tree, wlists, iisource, nsources)
  fortfmm.get_winfo(itarget_tree, wlists, iitarget, ntargets)
  
  # store the information for all boxes into boxes
  # ibox starts with index 1 in fortran
  for i in range(nboxes):
    ibox = i+1
    ier = fortfmm.get_box_info(ibox, box, center_temp, corners_temp, wlists, iwlists)
    boxes[i,:] 		= box[:]
    centers[i,:] 	= center_temp[:]
    corners[i,:,:]	= corners_temp[:]
    
  return

#----------------------------------------
# FMM Potential/ Field Functions
#---------------------------------------
  
def reset_pot_field(iffmm, ifdirect, ifpot, iffld, ifpottarg, iffldtarg,
		  ifpotsrc_far, ifpotsrc_near,
		  iffldsrc_far, iffldsrc_near,
		  ifpottarg_far, ifpottarg_near,
		  iffldtarg_far, iffldtarg_near,
		  pot_src, field_src, pot_targ, field_targ, 
		  d_pot_src, d_field_src, d_pot_targ, d_field_targ,
		  pot_far_src, field_far_src,
		  pot_near_src, field_near_src,
		  pot_far_targ, field_far_targ,
		  pot_near_targ, field_near_targ):
		    
  if iffmm == 1:
    
    if ifpot == 1:
      pot_src[:] = 0.0
      
    if iffld == 1:
      field_src[:,:] = 0.0
      
    if ifpottarg == 1:
      pot_targ[:] = 0.0
      
    if iffldtarg == 1:
      field_targ[:,:] = 0.0
      
  if ifdirect == 1:
    
    if ifpot == 1:
      d_pot_src[:] = 0.0
      
    if iffld == 1:
      d_field_src[:,:] = 0.0
      
    if ifpottarg == 1:
      d_pot_targ[:] = 0.0
      
    if iffldtarg == 1:
      d_field_targ[:,:] = 0.0
      
  # FMM multirate source
  if ifpotsrc_far == 1:
    pot_far_src[:] = 0.0
      
  if ifpotsrc_near == 1:
    pot_near_src[:] = 0.0
    
  if iffldsrc_far == 1:
    field_far_src[:,:] = 0.0
          
  if iffldsrc_near == 1:
    field_near_src[:,:] = 0.0
    
  # FMM multirate targets
  if ifpottarg_far == 1:
    pot_far_targ[:] = 0.0
      
  if ifpottarg_near == 1:
    pot_near_targ[:] = 0.0
    
  if iffldtarg_far == 1:
    field_far_targ[:,:] = 0.0        
    
  if iffldtarg_near == 1:
    field_near_targ[:,:] = 0.0

  return

#----------------------------------------
# Save Functions
#---------------------------------------

""" This sets the FMM tolerance of the potential for Greengard's FMM.
  Note that usually the E-Field tolerance is 3 orders less than the potential's
"""
def set_tol(iprec):
  
  """
  -2 => tolerance =.5d0
  -1 => tolerance =.5d-1
  0 => tolerance =.5d-2
  1 => tolerance =.5d-3
  2 => tolerance =.5d-6
  3 => tolerance =.5d-9
  4 => tolerance =.5d-12
  5 => tolerance =.5d-15
  """
  
  tol = -1
  
  if iprec == -2:
    tol = 0.5
    
  elif iprec == -1:
    tol = 0.5 * 10**-1
    
  elif iprec == 0:
    tol = 0.5 * 10**-2
    
  elif iprec == 1:
    tol = 0.5 * 10**-3
    
  elif iprec == 2:
    tol = 0.5 * 10**-6
    
  elif iprec == 3:
    tol = 0.5 * 10**-9
    
  elif iprec == 4:
    tol = 0.5 * 10**-12
    
  elif iprec == 5:
    tol = 0.5 * 10**-15
    
  return tol

# this function saves data. given a directory and a file name. 

def save_data(fpath, fname, data):
  
  if not os.path.exists(fpath):
    os.makedirs(fpath)
  
  fname = fpath + fname
  
  np.savetxt(fname, data ) 
  
def save_params(fpath,iterations, nnodes, levels, refinement, mpi_size, nthreads, iprec0, iprec1):
 
  """ The purpose of this code is to save the autorun parameters into a file"""
  
  tol0 = set_tol(iprec0)
  tol1 = set_tol(iprec1)
  
  # write basic info to info file
  f = open(fpath + '/pfasst_params.txt', 'w')
  f.write("iterations: " + str(iterations) + "\n")
  f.write("nnodes: " + str(nnodes) + "\n")
  f.write("levels: " + str(levels) + "\n")
  f.write("refine: " + str(refinement) + "\n")
  f.write("mpi size: " + str(mpi_size) + "\n")
  f.write("omp size: " + str(nthreads) + "\n")  
  f.write("err tolerance of the fine fmm-potential calculation: " + str(tol0) + "\n")
  f.write("err tolerance of the coarse fmm-potential calculation: " + str(tol1) + "\n")
  f.write("iprec0 " + str(iprec0) + "\n")
  f.write("iprec1 " + str(iprec1) + "\n")
  f.close()
  return
  
def save_boxes(fpath, boxes, centers, corners, isource_tree, itarget_tree, ladder):
  # save information about the boxes
  
  nboxes = corners.shape[0]
  temp = np.reshape( corners, (nboxes*8, 3) )
  
  np.savetxt( fpath + "/boxes.txt", boxes )
  np.savetxt( fpath + "/box_centers.txt", centers )
  
  # the corners array is reshaped from 3D of (n,8,3) to 2D of (n*8, 3)
  np.savetxt( fpath + "/box_corners.txt", temp)
  
  # save the particles sorted in box-order
  np.savetxt( fpath + "/isource_tree.txt", isource_tree )
  np.savetxt( fpath + "/itarget_tree.txt", itarget_tree )
  
  # save the ladder
  np.savetxt( fpath + "/ladder.txt", ladder )
  return


def save_pot_field2(fpath, pot_far_src, pot_near_src, field_far_src, field_near_src):
		      
  np.savetxt( fpath + "/pot_near_source.txt", np.column_stack([pot_near_src.real, pot_near_src.imag]) )  
  np.savetxt( fpath + "/pot_far_source.txt", np.column_stack([pot_far_src.real, pot_far_src.imag]) )
  
  np.savetxt( fpath + "/field_near_source.txt", np.column_stack([field_near_src.real, field_near_src.imag]) )  
  np.savetxt( fpath + "/field_far_source.txt", np.column_stack([field_far_src.real, field_far_src.imag]) )
  
  return

def save_pot_field( fpath, iffmm, ifdirect, ifpot, iffld, ifpottarg, iffldtarg, 
		    ifpotsrc_far, ifpotsrc_near,
		    iffldsrc_far, iffldsrc_near,
		    ifpottarg_far, ifpottarg_near,
		    iffldtarg_far, iffldtarg_near,
		    pot_src, pot_targ, field_src, field_targ,
		    d_pot_src, d_field_src, d_pot_targ, d_field_targ,
		    pot_far_src,  pot_near_src, 
		    field_far_src,  field_near_src,
		    pot_far_targ, field_far_targ, pot_near_targ, field_near_targ,
		    accel_src, accel_targ, d_accel_src, d_accel_targ ):

  if iffmm == 1:
    
    if ifpot == 1:
      np.savetxt( fpath + "/pot_source.txt", np.column_stack([pot_src.real, pot_src.imag]) )
    if iffld == 1:
      np.savetxt( fpath + "/field_source.txt", np.column_stack([field_src.real, field_src.imag]) )
    
    if ifpottarg == 1:
      np.savetxt( fpath + "/pot_target.txt", np.column_stack([pot_targ.real, pot_targ.imag]) )
    if iffldtarg == 1:
      np.savetxt( fpath + "/field_target.txt", np.column_stack([field_targ.real, field_targ.imag]) )

  if ifdirect == 1:
    
    if ifpot == 1:
      np.savetxt( fpath + "/d_pot_source.txt", np.column_stack([d_pot_src.real, d_pot_src.imag]) ) 
    if iffld == 1:
      np.savetxt( fpath + "/d_field_source.txt", np.column_stack([d_field_src.real, d_field_src.imag]) )

    if ifpottarg == 1:
      np.savetxt( fpath + "/d_pot_target.txt", np.column_stack([d_pot_targ.real, d_pot_targ.imag]) )
    if iffldtarg == 1:
      np.savetxt( fpath + "/d_field_target.txt", np.column_stack([d_field_targ.real, d_field_targ.imag]) )

  # FMM multirate source
  if ifpotsrc_near ==1:
    np.savetxt( fpath + "/pot_near_source.txt", np.column_stack([pot_near_src.real, pot_near_src.imag]) )  

  if ifpotsrc_far == 1:
    np.savetxt( fpath + "/pot_far_source.txt", np.column_stack([pot_far_src.real, pot_far_src.imag]) )

  if iffldsrc_near == 1:
    np.savetxt( fpath + "/field_near_source.txt", np.column_stack([field_near_src.real, field_near_src.imag]) )
      
  if iffldsrc_far == 1:
    np.savetxt( fpath + "/field_far_source.txt", np.column_stack([field_far_src.real, field_far_src.imag]) )
    
  # FMM multirate target
  if ifpottarg_far == 1:
    np.savetxt( fpath + "/pot_far_target.txt", np.column_stack([pot_far_targ.real, pot_far_targ.imag]) )
        
  if ifpottarg_near == 1:
    np.savetxt( fpath + "/pot_near_target.txt", np.column_stack([pot_near_targ.real, pot_near_targ.imag]) )

  if iffldtarg_far == 1:
      np.savetxt( fpath + "/field_far_target.txt", np.column_stack([field_far_targ.real, field_far_targ.imag]) )
          
  if iffldtarg_near == 1:
      np.savetxt( fpath + "/field_near_target.txt", np.column_stack([field_near_targ.real, field_near_targ.imag]) )  
      
  np.savetxt(fpath + "/accel_source.txt", accel_src)
  np.savetxt(fpath + "/accel_target.txt", accel_targ)
  np.savetxt(fpath + "/d_accel_source.txt", d_accel_src)
  np.savetxt(fpath + "/d_accel_target.txt", d_accel_targ)
  
  return
 
def save_pos_new(fpath, pos_src, pos_targ, dofmm, dodirect, dofmm_multi, dosource, dotarget):
  
  if (dofmm or dofmm_multi):
    if dosource:
      np.savetxt( fpath + "/pos_source.txt", pos_src)
      
    if dotarget:
      np.savetxt( fpath + "/pos_target.txt", pos_targ)

  elif dodirect:
    if dosource:
      np.savetxt( fpath + "/d_pos_source.txt", pos_src)
      
    if dotarget:
      np.savetxt( fpath + "/d_pos_target.txt", pos_targ)
      
  return
  
def save_vel_new(fpath, vel_src, vel_targ, dofmm, dodirect, dofmm_multi, dosource, dotarget):
  
  if (dofmm or dofmm_multi):
    if dosource:
      np.savetxt( fpath + "/vel_source.txt", vel_src)
      
    if dotarget:
      np.savetxt( fpath + "/vel_target.txt", vel_targ)

  elif dodirect:
    if dosource:
      np.savetxt( fpath + "/d_vel_source.txt", vel_src)
      
    if dotarget:
      np.savetxt( fpath + "/d_vel_target.txt", vel_targ)  
      
  return
  
def save_pot_field_new( fpath, dofmm, dodirect, dofmm_multi, dofmm_near, dosource, dotarget, dopot, dofield,
		    pot_src, field_src, pot_targ, field_targ,
		    pot_far_src, pot_near_src, 
		    field_far_src, field_near_src,
		    pot_far_targ, field_far_targ,
		    pot_near_targ, field_near_targ):

  # ----------------------------------------------
  # Full FMM
  # ----------------------------------------------
  if (dofmm or dofmm_near) and dosource:
    if dopot:
      np.savetxt( fpath + "/pot_source.txt", np.column_stack([pot_src.real, pot_src.imag]) )
      
    if dofield:
      np.savetxt( fpath + "/field_source.txt", np.column_stack([field_src.real, field_src.imag]) )
    
  if (dofmm or dofmm_near) and dotarget:
    
      if dopot:
	np.savetxt( fpath + "/pot_target.txt", np.column_stack([pot_targ.real, pot_targ.imag]) )
	
      if dofield:
	np.savetxt( fpath + "/field_target.txt", np.column_stack([field_targ.real, field_targ.imag]) )

  # ----------------------------------------------
  # Direct Solve
  # ----------------------------------------------
  if dodirect and dosource:
    
    if dopot:
      np.savetxt( fpath + "/d_pot_source.txt", np.column_stack([pot_src.real, pot_src.imag]) ) 
      
    if dofield:
      np.savetxt( fpath + "/d_field_source.txt", np.column_stack([field_src.real, field_src.imag]) )

  if dodirect and dotarget: 
  
    if dopot:
      np.savetxt( fpath + "/d_pot_target.txt", np.column_stack([pot_targ.real, pot_targ.imag]) )
      
    if dofield:
      np.savetxt( fpath + "/d_field_target.txt", np.column_stack([field_targ.real, field_targ.imag]) )

  # ----------------------------------------------
  # Multirate FMM
  # ----------------------------------------------
  
  if (dofmm_multi or dofmm_near) and dosource:
    
    if dopot:
      np.savetxt( fpath + "/pot_near_source.txt", np.column_stack([pot_near_src.real, pot_near_src.imag]) )      
      np.savetxt( fpath + "/pot_far_source.txt", np.column_stack([pot_far_src.real, pot_far_src.imag]) )

    if dofield:
      np.savetxt( fpath + "/field_near_source.txt", np.column_stack([field_near_src.real, field_near_src.imag]) )      
      np.savetxt( fpath + "/field_far_source.txt", np.column_stack([field_far_src.real, field_far_src.imag]) )
    
  if (dofmm_multi or dofmm_near) and dotarget:
    
    if dopot:
      np.savetxt( fpath + "/pot_far_target.txt", np.column_stack([pot_far_targ.real, pot_far_targ.imag]) )      
      np.savetxt( fpath + "/pot_near_target.txt", np.column_stack([pot_near_targ.real, pot_near_targ.imag]) )

    if dofield:
      np.savetxt( fpath + "/field_far_target.txt", np.column_stack([field_far_targ.real, field_far_targ.imag]) )      
      np.savetxt( fpath + "/field_near_target.txt", np.column_stack([field_near_targ.real, field_near_targ.imag]) )  
      
  return
  
  
def save_pos(fpath, source, target, d_source, d_target):
  np.savetxt( fpath + "/pos_source.txt", source )
  np.savetxt( fpath + "/pos_target.txt", target )
      
  np.savetxt( fpath + "/d_pos_source.txt", d_source )
  np.savetxt( fpath + "/d_pos_target.txt", d_target )
  return
  
def save_vel(fpath, vel, veltarget, d_vel, d_veltarget):
  np.savetxt(fpath + "/vel_source.txt", vel)
  np.savetxt(fpath + "/vel_target.txt", veltarget)
	
  np.savetxt(fpath + "/d_vel_source.txt", d_vel)
  np.savetxt(fpath + "/d_vel_target.txt", d_veltarget)
  return

def set_path(fpath, iteration, nsteps):
  #fpath = os.path.dirname( os.getcwd() ) # this should be /fmmlib3d
  #fpath = fpath + "/mydata/"
  
  newpath = fpath + "/t_" + str(iteration+1)
  oldpath = fpath + "/t_" + str(iteration)

  if not os.path.exists(newpath) and iteration < nsteps:
	os.mkdir(newpath)
	
  return newpath, oldpath

#----------------------------------------
# Printing for Tests
#---------------------------------------
def print_flags(dosource, dotarget, dofmm, dofmm_multi, dodirect, dopot):
  
  # set flags for the potential
  (ifpot, ifpottarg, ifpotsrc_far, ifpotsrc_near, ifpottarg_far, ifpottarg_near) = set_flags_pot(dopot, dofmm, dofmm_multi, dodirect, dosource, dotarget)
  
  # set flags for the electric-field
  (iffld, iffldtarg, iffldsrc_far,  iffldsrc_near, iffldtarg_far, iffldtarg_near) = set_flags(dofmm, dofmm_multi, dodirect, dosource, dotarget)  
  
  if dofmm_multi:
    msg = "multirate FMM"
  
  if dofmm:
    msg = "full FMM"

  if dodirect:
    msg = "direct"
    
  msg = msg + "\n"
  
  print "\n-----------------------------------------"
  print "solver: " + msg
  print_func(ifpot, iffld, "ifpot", "iffld")
  print_func(ifpotsrc_far, iffldsrc_far, "ifpotsrc_far", "iffldsrc_far")  
  print_func(ifpotsrc_near, iffldsrc_near, "ifpotsrc_near", "iffldsrc_near")
  print
  print_func(ifpottarg, iffldtarg, "ifpottarg", "iffldtarg")
  print_func(ifpottarg_far, iffldtarg_far, "ifpottarg_far", "iffldtarg_far")  
  print_func(ifpottarg_near, iffldtarg_near, "ifpottarg_near", "iffldtarg_near")
  print "-----------------------------------------\n"
  
  return

def print_func(num1, num2, fname1, fname2):
  
  print "(" + fname1 + ", " + fname2 + "): " + "(" + str(num1) + ", " + str(num2) + ")"
  return
  
def printStuff(iisource, iitarget, iwlists, lwlists, lused7, nbox, nboxes, ladder, size, epsfmm, nlev):
  print "iisource: " + str(iisource)
  print "iitarget: " + str(iitarget)
  print "iwlists: " + str(iwlists)
  print "lwlists: " + str(lwlists)
  print "lused7: " + str(lused7)
  print "nbox: " + str(nbox)
  print "nboxes: " + str(nboxes)
  print "size: " + str(size)
  print "epsfmm: " + str(epsfmm)
  print "nlev: " + str(nlev)
  #print "ladder"
  #print ladder
  return
  
def print_box_info(ibox, nboxes, nbox, nlev, ladder, box, center, corners,  w, ier):
    
  print "ier: " + str(ier)
  print "ibox: " + str(ibox)
  print "nlev: " + str(nlev)
  print "nboxes: " + str(nboxes)
  
  # the last entry in ladder should be empty
  # boxes start at level 2 to get well-separated boxes (starting from 0)
  print "ladder: "
  print str(ladder[0:nlev+3,:])
  
  print "\nRECALL the index addresses may be in fortran ordering\n"
  
  print "level of subdivision on which this box was constructed:"
  print  str(box[0])
  """ 
  The coordinates: (1,1,1) indicate that box is the first on its levels. All coordinate locations
  are referenced to this box for the other boxes on the same level.
  """
  print "coordinates of this box among all boxes on this level"
  print "( " + str(box[1]) + ", " + str(box[2]) + ", " + str(box[3]) + " )"
  
  print "father of this box (address in array 'boxes'):"
  print str(box[4])
  print "children of this box (address in array 'boxes'): " 
  print str(box[5:12])
  
  print "location of the array iz of sources living in this box: "
  print  str(box[13])
  print "number of particles (sources?) living in this box: "
  print str(box[14])
  
  print "the location in the array iztarg of the targets living in this box: "
  print str(box[15])
  print "the number of targets living in this box: "
  print str(box[16])
  
  print "source box type: 0 -> empty, 1 -> leaf node, 1 -> sub-divided: "
  print str(box[17])
  print "target box type: 0 -> empty, 1 -> leaf node, 1 -> sub-divided: "
  print str(box[18])
  print "reserved for future use: "
  print str(box[19])
  
  print ""
  print "center of this box:"
  print center
  print "corners of this box:"
  print corners
  print ""
  
  return   
  
"""============================================"""
# FMM Potential/ Field Functions
#============================================
  