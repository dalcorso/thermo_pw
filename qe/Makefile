# Makefile for Modules

include ../../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../FFTXlib \
          $(MOD_FLAG)../../PW/src $(MOD_FLAG)../../PHonon/PH \
          $(MOD_FLAG)../lib  $(MOD_FLAG).

MODULES = \
adddvepsi_us.o \
check_bands.o \
control_qe.o \
do_phonon.o \
electrons.o \
ev_xml.o \
phqscf.o \
prepare_q.o \
run_nscf.o \
read_dyn_from_file.o \
solve_linter.o 

all : thermo_qe.a

thermo_qe.a: 	$(MODULES)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *.F90 *.mod *.L 

include make.depend
