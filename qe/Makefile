# Makefile for Modules

include ../../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../FFTXlib \
          $(MOD_FLAG)../../PW/src $(MOD_FLAG)../../PHonon/PH \
          $(MOD_FLAG)../lib  $(MOD_FLAG).

MODULES = \
adddvepsi_us.o \
addusddenseq.o \
apply_ac.o \
check_bands.o \
compute_intq.o \
control_qe.o \
ccg_psi.o \
do_phonon.o \
dveqpsi_us.o \
dveqpsi_us_only.o \
electrons.o \
ev_xml.o \
initialize_ph.o \
orthogonalize_omega.o \
paw_add_symmetry.o \
phescf.o \
phqscf.o \
pheqscf.o \
phq_readin.o \
phq_summary.o \
polarizc.o \
prepare_q.o \
psyme_fpol.o \
psym_dmageq.o \
psymeq.o \
psym_dmage_fpol.o \
run_nscf.o \
read_dyn_from_file.o \
scal_prod.o \
set_intq_nc.o \
solve_e_fpolc.o \
solve_eq.o \
solve_linter.o  \
syme_fpol.o     \
symeq.o \
sym_dmage_fpol.o \
sym_dmageq.o \
transform_intq_nc.o \
transform_intq_so.o \
usnldiagq.o

all : thermo_qe.a

thermo_qe.a: 	$(MODULES)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *.F90 *.mod *.L 

include make.depend
