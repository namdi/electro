# This script creates the following libraries that will be used to wrap fortran into python:
#	a library for the FMM to be used in fortfmm.py
#	a library for the direct solver to be used in direct.py

# This program creates the library that fortfmm.py will use to wrap fortran functions
# into a python module.
# fortfmm.f90 fortfmm_utils.f90 ../../fmmlib3d/src/lfmm3dpart.f  ../../fmmlib3d/src/lfmm3drouts.f

# DOES THIS APPLY TO THE PFASST CODE OR JUST TO GREENGARD'S FMM.
# Anyway, the PFASST compile_all.sh does make all of the object files make by running the make.

# Note. In order to have the object files from the original fortran files (ie, the ones that I did not create)
# one needs to run the following before hand:
# $ make cfmm3dpart-openmp (for example)
# This creates the object files like lfmm3dpart.o. The object files will be found in fmmlib3d/examples

# have variables for various file paths
INITIAL_DIR=$(pwd)				# store the current directory 
TEMP_DIR=$(dirname $INITIAL_DIR)		# go up a directory

#SRC_DIR=$TEMP_DIR/src				# store the path to the src folder
#MYSRC_DIR=$TEMP_DIR/mysrc			# store the path to the examples folder
#OBJ_DIR=$TEMP_DIR/myobj				# store#store the path to the src directory the path to the obj folder

OBJ_DIR=$TEMP_DIR/build				# store#store the path to the src directory the path to the obj folder
PY_OBJ_DIR=$TEMP_DIR/pyobj

#----------------------------------------
# Create Library to access fortran FMM functions
#----------------------------------------

# create the object files (Is the next line redundant)
#gfortran -g -c -fPIC -fopenmp -O3 ${MYSRC_DIR}/fortfmm.f90 	${MYSRC_DIR}/fortfmm_utils.f90 	${SRC_DIR}/my_func.f ${SRC_DIR}/driver.f

gfortran -g -shared -fopenmp -o	${PY_OBJ_DIR}/_my_fortfmmlib.so	${OBJ_DIR}/fortfmm.o 		${OBJ_DIR}/fortfmm_utils.o	\
				${OBJ_DIR}/d3tstrcr.o		${OBJ_DIR}/dlaran.o 		${OBJ_DIR}/l3dterms.o 		\
				${OBJ_DIR}/l3dtrans.o 		${OBJ_DIR}/laprouts3d.o		${OBJ_DIR}/lfmm3dpart.o 	\
				${OBJ_DIR}/lfmm3drouts.o 	${OBJ_DIR}/yrecursion.o		${OBJ_DIR}/d3mtreeplot.o	\
				${OBJ_DIR}/rotviarecur3.o	${OBJ_DIR}/prini.o		${OBJ_DIR}/prinm.o		\
				${OBJ_DIR}/my_func.o		
echo ""
echo "Finished creating _my_fortfmmlib.so"

#----------------------------------------
# Creating a library to the direct solvers
#----------------------------------------
#gfortran -g -c -fPIC -fopenmp -O3 ${MYSRC_DIR}/direct.f90

gfortran -g -shared -fopenmp -o ${PY_OBJ_DIR}/_directlib.so ${OBJ_DIR}/direct.o

echo ""
echo "finished creating _directlib.so"
echo ""

#-----------------------------------------
# creating test driver
#----------------------------------------
#export LD_LIBRARY_PATH=$OBJ_DIR
#gfortran -g -fopenmp ${OBJ_DIR}/driver.o ${OBJ_DIR}/_my_fortfmmlib.so -o my_test 
