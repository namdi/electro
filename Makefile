#
# Makefile for electro example (bassed of off mpi-advection example).  This Makefile builds LIBPFASST
# inplace (instead of using the static library form of LIBPFASST) as
# well as the application specific code in the src directory.
#

LIBPFASST=/home/namdi/Documents/School/UNC/Parallel_Time/Code/fpfasst/libpfasst
PFDEV=/home/namdi/Documents/School/UNC/Parallel_Time/Code/fpfasst 

EXE = main.exe

include $(LIBPFASST)/Makefile.defaults

MYSRC_DIR=$(LIBPFASST)/myfmm
SRC_DIR=$(LIBPFASST)/fmm
LIB_DIR=$(LIBPFASST)/lib

F90SRC = $(wildcard src/*.f90 $(LIBPFASST)/src/*.f90 fmm/*.f90 myfmm/*.f90)  verlet.f90
FSRC   = $(wildcard fmm/*.f)
CSRC   = $(wildcard src/*.c $(LIBPFASST)/src/*.c)
OBJ    = $(addprefix build/, $(notdir $(F90SRC:.f90=.o) $(FSRC:.f=.o) $(CSRC:.c=.o)))

FFLAGS  += -g -fopenmp -O3 # do not optimize for debugging it changed the order of things
LDFLAGS += -g -fopenmp

VPATHS = $(LIBPFASST)/src src fmm myfmm $(PFDEV)

all: $(EXE)

include $(LIBPFASST)/Makefile.rules


build/%.o: %.f
	@mkdir -p build
	$(FC) $(FFLAGS) -c $< $(OUTPUT_OPTION)

vpath %.f $(VPATHS)

#
# dependencies
#

#build/fortfmm_utils.o:	lib/_my_fortfmmlib.so
build/my_type.o: build/pf_dtype.o build/pf_ndarray.o
build/tree.o: build/my_type.o build/fortfmm_utils.o
build/electro.o: build/tree.o build/my_type.o

build/feval.o: build/pfasst.o build/my_utils.o build/my_type.o build/electro.o build/tree.o build/print_mod.o build/fortfmm_utils.o # build/encap.o
build/hooks.o: build/pfasst.o build/probin.o
build/main.o: build/feval.o build/hooks.o build/transfer.o build/my_utils.o build/numpy.o build/verlet.o # build/encap.o

