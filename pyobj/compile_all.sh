# This script compiles all the source files with -fPIC


#store the path to the src directory
INITIAL_DIR=$(pwd)				# store the current directory 
TEMP_DIR=$(dirname $INITIAL_DIR)		# go up a directory

SRC_DIR=$TEMP_DIR/src				# store the path to the src folder (which holds data from Greengard's FMM)
MYSRC_DIR=$TEMP_DIR/mysrc			# store the path to the examples folder


# compile everything
# the -O3 flag is for optimization
gfortran -g -c -fPIC -fopenmp -O3 \
		  ${MYSRC_DIR}/dlaran.f 	${MYSRC_DIR}/hfmm3dpart_driver.f 	${MYSRC_DIR}/hfmm3dtria_driver.f 	\
		  ${MYSRC_DIR}/hkrand.f 	${MYSRC_DIR}/lfmm3dpart_driver.f 	${MYSRC_DIR}/lfmm3dtria_driver.f 	\
		  ${SRC_DIR}/cdjseval3d.f	${SRC_DIR}/helmrouts3d.f 		${SRC_DIR}/laprouts3d.f 		\
		  ${SRC_DIR}/prinm.f 		${SRC_DIR}/trilib.f 			${SRC_DIR}/d3mtreeplot.f  		\
		  ${SRC_DIR}/hfmm3dpart.f   	${SRC_DIR}/legeexps.f			${SRC_DIR}/projections.f		\
		  ${SRC_DIR}/triquadflatlib.f 	${SRC_DIR}/d3tstrcr.f 								\
		  ${SRC_DIR}/hfmm3drouts.f  	${SRC_DIR}/lfmm3dpart.f 		${SRC_DIR}/l3dtrans.f 			\
		  ${SRC_DIR}/rotproj.f  	${SRC_DIR}/yrecursion.f 							\
		  ${SRC_DIR}/dfft.f  		${SRC_DIR}/hfmm3dtria.f   		${SRC_DIR}/lfmm3drouts.f  		\
		  ${SRC_DIR}/rotviarecur3.f  	${SRC_DIR}/h3dterms.f  			${SRC_DIR}/l3dterms.f  			\
		  ${SRC_DIR}/lfmm3dtria.f  	${SRC_DIR}/triagauc.f 			${SRC_DIR}/h3dtrans.f  			\
		  ${SRC_DIR}/triahquad.f   	${SRC_DIR}/h3dtrirouts.f		${SRC_DIR}/l3dtrirouts.f  		\
		  ${SRC_DIR}/prini.f   		${SRC_DIR}/triasymq.f								

gfortran -g -c -fPIC -fopenmp -O3 \
		  ${MYSRC_DIR}/fortfmm.f90	${MYSRC_DIR}/fortfmm_utils.f90		${SRC_DIR}/my_func.f			\
		  ${MYSRC_DIR}/direct.f90	${SRC_DIR}/driver.f
echo ""
echo "Finished compiling EVERYTHING"
echo ""
