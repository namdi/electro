"""
This code is for error checking. This script makes a reference solution.
These are some basic parameters for the reference solution:
$ python run_refsol.py -t 1 -b 20 -k 20 -s $DIRECT -a 3 -i $MY_DATA_DIR/ref -o $MY_DATA_DIR/ref

Make sure that SDC noodes are set to 9 in main.exe
Make sure reference solution is set to 1 level.
dt = 0.01
"""
#------------------------------
# Import 
#------------------------------
import subprocess 

import data_circle
import data_sq
import my_utils

#------------------------------
# parse command line
#------------------------------
nsources, ntargets, nthreads, distribution, accuracy, nsteps, indir, outdir, solver, niters, nlevs, fnodes = my_utils.parse_cmd_line()

#------------------------------
# Set the Initial Condition
#------------------------------
# set the 4 particle square initial condition
data_sq.run(nsteps, indir)

# set the 2D ring initial condition
#data_circle.run(nsources, ntargets, nsteps, distribution, indir)

#------------------------------
# Create Initial condition
#------------------------------

# numper of MPI processors
nMPI = 1

# directory in which to load data from
indir = indir + "/t_0"

# get the directory of the executable
exec_dir = "/home/namdi/Documents/School/UNC/Parallel_Time/Code/fpfasst/electro"
fexec = exec_dir + "/main.exe"

# create the command to be run from the command line
command = "mpirun -n " + str(nMPI) + " " + str(fexec) + " " + str(nthreads) + " " + " " + str(solver) + " " + " " + str(accuracy)
command = command + " " + str(niters) + " " + str(nlevs) + " " + str(fnodes) + " " + str(indir) + " " + str(outdir)

# call the command
p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

stdout, stderr = p.communicate()
