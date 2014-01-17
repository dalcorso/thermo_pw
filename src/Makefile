# Makefile for thermo_pw 

include ../../make.sys

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../PW/src $(MOD_FLAG)../../PHonon/PH \
          $(MOD_FLAG)../lib $(MOD_FLAG).
IFLAGS=

THERMOOBJS = \
bands_sub.o \
bcast_thermo_input.o \
collect_everything.o \
compute_cp.o \
compute_thermo.o  \
deallocate_asyn.o \
deallocate_thermo.o \
deallocate_thermo.o \
do_pwscf.o \
ev_sub.o \
ev_xml.o \
initialize_thermo_work.o \
initialize_thermo_master.o \
matdyn_sub.o \
plot_anhar.o \
plotband_sub.o \
plot_thermo.o \
q2r_sub.o \
run_thermo_asyncronously.o \
set_k_points.o \
set_paths_disp.o \
set_thermo_work_todo.o \
simple_plot.o \
solve_linter.o \
sym_band_sub.o \
thermo_mod.o \
thermo_readin.o \
thermo_pw.o \
write_anhar.o \
write_ev_input.o 
             

PWOBJS = ../../PHonon/PH/libph.a ../../PP/src/libpp.a ../../PW/src/libpw.a 
QEMODS = ../../Modules/libqemod.a

LIBOBJS	= ../lib/thermo_lib.a ../../PHonon/PH/libph.a ../../flib/ptools.a ../../flib/flib.a ../../clib/clib.a ../../iotk/src/libiotk.a

TLDEPS= bindir mods libs pw

all : libthermo_pw.a thermo_pw.x 

libs-ph : libph.a libphaux.a

thermo_pw.x : thermo_pw.o libthermo_pw.a $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ thermo_pw.o ../../PHonon/PH/rigid.o \
                 ../../PHonon/PH/write_eigenvectors.o libthermo_pw.a \
		 $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../THERMO_PW/src/ph.x . )

libthermo_pw.a : $(THERMOOBJS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *.F90 *.d *.mod *.i *.L
	- /bin/rm -rf ../../bin/thermo_pw.x

include make.depend
# DO NOT DELETE
