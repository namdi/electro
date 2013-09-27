"""
File: my_utils.py
"""

#-----------------------------------------
# Import 
#----------------------------------------
from optparse import OptionParser

import numpy as np

"""
The following parses the command line for the arguments:
  nsources, ntargets, nthreads, distribution, accuracy
"""

#============================================
# Function Declarations
#============================================    


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
  
  parser.add_option("-l", dest="nlevs", help="the number of PFASST levels")
  
  parser.add_option("-f", dest="fnodes", help="the number of SDC nodes on the finest level")
  
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
    distribution = 0
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
    outdir = str(options.outdir)

  if options.indir is None:
    indir = None
  else:
    indir = str(options.indir)
    
  if options.solver is None:
    solver = 3
  else:
    solver = int(options.solver)
    
    if solver_check(solver) is False:
      raise NotImplementedError, "need to enter a value of a valid solver"
    
  if options.niters is None:
    niters = 3
  else:
    niters = int(options.niters)
    
  if options.fnodes is None:
    fnodes = 5
  else:
    fnodes = int(options.fnodes)
  
  if options.nlevs is None:
    nlevs = 2
  else:
    nlevs = int(options.nlevs)
    
  return nsources, ntargets, nthreads, distribution, accuracy, nsteps, indir, outdir, solver, niters, nlevs, fnodes
  
# this function checks to make sure the solver has the right values
def solver_check(solver):
  
  # valid values for solver
  num_solvers = 4
  vals = np.linspace(1, num_solvers, num_solvers)
  
  check = (vals == solver)
  
  ok = False
  
  # we have a match of a valid solver
  if sum(check) == 1:
    ok = True
  
  return ok  
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
  
  
def load_pfasst_data(fpath, q, level, iteration):
  
  # q.shape: [nsteps+1, 2, N, 3] 
  nsteps = q.shape[0] -1
  N = q.shape[2]
   
  for i in range(nsteps):    
    fname 	= fpath + "/s%05di%03dl%02d.npy" % (i, iteration, level)
    x 		= np.load(fname) # data has shape [3, N, 2]
    q[i+1, 0] 	= np.transpose( x[:,:,0] )
    q[i+1, 1]	= np.transpose( x[:,:,1] )
  return
    

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
  