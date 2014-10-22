# Change paths and associated libraries to compiled case directory  

CASEBLDDIR = /gscratch1/elroesl/acme_scratch/f1850c5_ne30/bld
CASEROOT = /home/elroesl/acme/ACME-interpic/ACME/scripts/f1850c5_ne30/
include $(CASEROOT)/Macros

NETCDFROOT = /projects/ccsm/netcdf4.3.2-intel12
NETCDF_PATH = /projects/ccsm/netcdf4.3.2-intel12
PNETCDF_PATH = /projects/ccsm/pnetcdf1.2.0-intel2

LIB_NETCDF = $(NETCDF_PATH)/lib
INC_NETCDF = $(NETCDF_PATH)/include
LIB_PNETCDF = $(PNETCDF_PATH)/lib

LIBDIR = $(CASEBLDDIR)/intel/openmpi/nodebug/nothreads

LIB_CSM_SHARE = $(LIBDIR)/MCT/noesmf/a1l1r1i1o1g1w1/csm_share
LIB_PIO = $(LIBDIR)/lib
LIB_MCT = $(LIBDIR)/lib
LIB_MPEU = $(LIBDIR)/lib


LIBS = -L$(LIB_CSM_SHARE) -lcsm_share -L$(LIB_PIO) -lpio -L$(LIB_MCT) -lmct -L$(LIB_MPEU) -lmpeu -L$(LIB_NETCDF) -lnetcdf -L$(LIB_PNETCDF) -lpnetcdf

MODDIR=$(LIBDIR)/include
## FC=ifort 
mkatmsrffile: mkatmsrffile.F90 
	$(FC) $(FFLAGS) -I$(MODDIR) -I$(INC_NETCDF) $< -o $@  $(LIBS) $(LDFLAGS) 
