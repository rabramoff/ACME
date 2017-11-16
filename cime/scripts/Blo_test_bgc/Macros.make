#
# COMPILER=gnu
# OS=Darwin
# MACH=perth
#
# Makefile Macros 
CPPDEFS+= -DFORTRANUNDERSCORE -DNO_R16 -DCPRGNU  -DSYSDARWIN   -DFORTRANUNDERSCORE -DNO_R16

SLIBS+=$(shell $(NETCDF_PATH)/bin/nf-config --flibs) -framework Accelerate

CFLAGS:= -mcmodel=medium 

CXX_LINKER:=FORTRAN

FC_AUTO_R8:= -fdefault-real-8 

FFLAGS:= -O -fconvert=big-endian -ffree-line-length-none -ffixed-line-length-none 

FFLAGS_NOOPT:= -O0 

FIXEDFLAGS:=  -ffixed-form 

FREEFLAGS:= -ffree-form 

HAS_F2008_CONTIGUOUS:=FALSE

MPICC:= /opt/local/bin/mpicc-mpich-gcc6  

MPICXX:= /opt/local/bin/mpicxx-mpich-gcc6 

MPIFC:= /opt/local/bin/mpif90-mpich-gcc6 

NETCDF_PATH:=/opt/local

PFUNIT_PATH:=$ENV{HOME}/local/pfunit/pfunit-sf.git.ae92605e8e

SCC:= /opt/local/bin/gcc-mp-6 

SCXX:= /opt/local/bin/gcc-mp-6 

SFC:= /opt/local/bin/gfortran-mp-6 

SUPPORTS_CXX:=TRUE

ifeq ($(DEBUG), FALSE) 
   FFLAGS +=  -O 
   CFLAGS +=  -O 
endif

ifeq ($(DEBUG), TRUE) 
   FFLAGS +=  -g -Wall 
   CFLAGS +=  -g -Wall -Og -fbacktrace -fcheck=bounds -ffpe-trap=invalid,zero,overflow
endif

ifeq ($(compile_threaded), true) 
   FFLAGS +=  -fopenmp 
   LDFLAGS +=  -fopenmp 
   CFLAGS +=  -fopenmp 
endif

ifeq ($(MODEL), cism) 
   CMAKE_OPTS +=  -D CISM_GNU=ON 
endif

ifeq ($(MODEL), clm) 
  ifeq ($(CLM_PFLOTRAN_COLMODE), TRUE) 
    ifeq ($(CLM_PFLOTRAN_COUPLED), TRUE) 
       CPPDEFS +=  -DCOLUMN_MODE 
    endif

  endif

  ifeq ($(CLM_PFLOTRAN_COUPLED), TRUE) 
     FFLAGS +=  -I$(CLM_PFLOTRAN_SOURCE_DIR) 
     CPPDEFS +=  -DCLM_PFLOTRAN 
  endif

endif

ifeq ($(MODEL), csm_share) 
   CFLAGS +=  -std=c99 
endif

ifeq ($(MODEL), driver) 
   LDFLAGS +=  -all_load 
  ifeq ($(CLM_PFLOTRAN_COUPLED), TRUE) 
     LDFLAGS +=  -L$(CLM_PFLOTRAN_SOURCE_DIR) -lpflotran $(PETSC_LIB) 
  endif

endif

ifeq ($(MODEL), pop) 
   CPPDEFS +=  -D_USE_FLOW_CONTROL 
endif

