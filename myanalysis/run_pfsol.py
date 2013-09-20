"""
This code is for error checking. This script makes the pfasst solution. This function is to be called from run_converge.py.
This code sets the initial condition and then runs the simulation. Since this is for testing, the test-cases usually are 
a 4 particle square, or a uniform circle.  I can make the condition in different directories and still compare results 
because these ICs do not have randomness.
"""

#------------------------------
# Import 
#------------------------------
import subprocess 

import data_circle
import data_sq

#------------------------------
# parse command line
#------------------------------
def run(nMPI, nthreads, nsources, ntargets, nsteps, niters, nlevs, fnodes, solver, accuracy, distribution, indir, outdir):
   
  #
  # Set the Initial Condition
  #
  
  # set the 4 particle square initial condition
  data_sq.run(nsteps, indir)

  # set the 2D ring initial condition
  #data_circle.run(nsources, ntargets, nsteps, distribution, indir)

  #
  # Call the fortran pfasst solver from the command line
  #
  
  # directory in which to load data from
  indir = indir + "/t_0"

  # get the directory of the executable
  exec_dir = "/home/namdi/Documents/School/UNC/Parallel_Time/Code/fpfasst/electro"
  fexec = exec_dir + "/main.exe"

  # create the command to be run from the command line
  command = "mpirun -n " + str(nMPI) + " " + str(fexec) + " " + str(nthreads) + " " + " " + str(solver) + " " + " " + str(accuracy)
  command = command + " " + str(niters) + " " + " " + str(nlevs) + " " + str(fnodes) + " " + str(indir) + " " + str(outdir)

  # call the pfasst executable
  p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
  stdout, stderr = p.communicate()
  
  return