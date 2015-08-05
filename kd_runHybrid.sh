#!/bin/csh                                                                                                                                            

#!/bin/csh                                                                                                                                            

# The order of command line paramters                                                                                                                 

# 1: The problem parameters

# activate parameters                                                                                                                                 
source /nas02/home/b/r/brandonn/fpfasst/electro/kd_activate

# run the program in serial
# bsub -n 1 -o out.%J -e err.%J -a openmp -q debug -R "span[ptile=1]" mpirun csh kd_runHybrid.sh $MY_DATA_DIR/test/params.nml
#
# where in params.nml
#
# nthreads_omp = 1
#
# run the program in parallel (4 MPI, 3 OpenMP per MPI)
# bsub -n 12 -o out.%J -e err.%J -a openmp -q debug -R "span[ptile=3]" mpirun csh kd_runHybrid.sh $MY_DATA_DIR/test/params.nml
#
# where in params.nml
#
# nthreads_omp = 3

# store the command line variables                                                                                                                    
set params = $1

# run                                                                                                                                                 
./main.exe $params

