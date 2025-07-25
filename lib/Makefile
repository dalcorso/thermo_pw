# Makefile for thermo_pw/lib

include ../../make.inc

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG).

MODULES = \
asy.o \
atomic_pos.o \
bz_2d_form.o \
bz_asy.o \
cubic_surfaces.o \
color.o \
debye_module.o \
diagonalize.o \
elastic_const.o \
eldos_module.o \
eos.o \
gnuplot.o \
gnuplot_color.o \
isoentropic.o \
latgen_2d.o \
lattices.o \
linear_solvers.o \
linear_surfaces.o \
magnetic_point_group.o \
model_free_energy.o \
mp_asyn.o \
numerical_diff.o \
nye.o \
optical.o \
ph_freq_module.o \
phdos_module.o \
piezoelectric_tensor.o \
point_group.o \
polyfit.o \
polynomial.o \
quadratic_surfaces.o \
quartic_surfaces.o \
rotate.o \
space_groups.o \
strain.o \
thermodynamics.o \
vector.o \
voigt.o \
xrdp.o 

all : thermo_lib.a

thermo_lib.a: $(MODULES)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L 

veryclean: clean
	- /bin/rm -f ./make.depend

include make.depend
