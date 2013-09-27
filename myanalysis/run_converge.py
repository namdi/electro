#
# This code runs the pfasst solution to an initia condition with the various parameters.
# This code make the initial condition and then runs the calculations calling the fortran pfasst solver in libfpasst.
#
# Make sure that the LIBPFASST directory is in your PYTHONPATH:
#
# export PYTHONPATH=/path/to/libpfasst
#
# Solvers:
# 1	# the Multirate FMM used in PFASST
# 2	# the FMM calculating far and local field separately 
# 3	# the full FMM
# 4	# the direct solver
#--------------------------------------------------------------------------------------

import run_pfsol	# run pfasst solution
#import run_refsol	# run reference solution

#nsources, ntargets, nthreads, distribution, accuracy, nsteps, indir, outdir, solver, niters, nlevs, fnodes = my_utils.parse_cmd_line()

FMM_MR		= 1 	# Multirate FMM
FMM_NEAR_FAR	= 2	
FMM		= 3
DIRECT 		= 4

nMPI		= 1
nthreads 	= 1

nsources 	= 4
ntargets 	= 3
nsteps		= 20
niters		= 4
nlevs		= 2
fnodes		= 5

solver 		= DIRECT
accuracy	= 3
distribution 	= -1 

indir 		= "/home/namdi/Documents/School/UNC/Parallel_Time/Data/fmm" 
outdir		="/home/namdi/Documents/School/UNC/Parallel_Time/Data/fmm"

# run the pfasst solution
run_pfsol.run(nMPI, nthreads, nsources, ntargets, nsteps, niters, nlevs, fnodes, solver, accuracy, distribution, indir, outdir)