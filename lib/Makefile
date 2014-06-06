# Makefile for Modules

include ../../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../Modules $(MOD_FLAG).

MODULES = \
asy.o \
bz_asy.o \
elastic_const.o \
piezoelectric_tensor.o \
gnuplot.o \
mp_asyn.o \
phdos_module.o \
ph_freq_module.o

all : thermo_lib.a

thermo_lib.a: 	$(MODULES)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *.F90 *.mod *.L 

include make.depend
