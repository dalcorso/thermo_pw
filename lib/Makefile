# Makefile for Modules

include ../../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../Modules $(MOD_FLAG)../../FFTXlib $(MOD_FLAG).

MODULES = \
asy.o \
debye_module.o \
bz_2d_form.o \
bz_asy.o \
elastic_const.o \
isoentropic.o \
latgen_2d.o \
gnuplot.o \
gnuplot_color.o \
lattices.o \
linear_solvers.o \
mp_asyn.o \
phdos_module.o \
ph_freq_module.o \
piezoelectric_tensor.o \
point_group.o \
quadratic_surfaces.o \
quartic_surfaces.o \
rotate.o \
space_groups.o \
strain.o \
xrdp.o 

all : thermo_lib.a

thermo_lib.a: 	$(MODULES)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *.F90 *.mod *.L 

include make.depend
