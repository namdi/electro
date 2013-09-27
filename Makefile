#
# Makefile for electro example (bassed of off mpi-advection example).  This Makefile builds the 
# electro (Laplace/ Coulomb potential)solver to be used with pfasst (via libpfasst)
#

LIBPFASST=../libpfasst

EXE = main.exe

include $(LIBPFASST)/Makefile.defaults

F90SRC  = $(wildcard src/*.f90 $(LIBPFASST)/src/*.f90 fmm/*.f90 myfmm/*.f90)
FSRC    = $(wildcard fmm/*.f)
CSRC   += $(wildcard src/*.c)
OBJ     = $(addprefix build/, $(notdir $(F90SRC:.f90=.o) $(FSRC:.f=.o) $(CSRC:.c=.o)))

FFLAGS  += -g -fopenmp -O3 # do not optimize for debugging it changed the order of things
LDFLAGS += -g -fopenmp

VPATHS = $(LIBPFASST)/src src fmm myfmm

all: $(EXE)

include $(LIBPFASST)/Makefile.rules


build/%.o: %.f
	@mkdir -p build
	$(FC) $(FFLAGS) -c $< $(OUTPUT_OPTION)

vpath %.f $(VPATHS)

#
# dependencies
#

build/my_type.o: build/pf_dtype.o build/pf_ndarray.o
build/tree.o: build/my_type.o build/fortfmm_utils.o
build/electro.o: build/tree.o build/my_type.o

build/feval.o: build/pfasst.o build/my_utils.o build/my_type.o build/electro.o build/tree.o build/print_mod.o build/fortfmm_utils.o # build/encap.o
build/hooks.o: build/pfasst.o build/probin.o
build/main.o: build/feval.o build/hooks.o build/transfer.o build/my_utils.o build/numpy.o #build/verlet.o # build/encap.o

