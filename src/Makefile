# Makefile for thermo_pw 

include ../../make.sys

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../PW/src $(MOD_FLAG)../../PHonon/PH \
          $(MOD_FLAG)../lib $(MOD_FLAG).
IFLAGS=

THERMOOBJS = \
addusdens1d_so.o \
allocate_thermo.o \
bands_sub.o \
bcast_thermo_input.o \
check_file_exists.o \
clean_dfft.o \
collect_everything.o \
compute_cp.o \
compute_gruneisen.o \
compute_omega_geo.o \
deallocate_asyn.o \
deallocate_thermo.o \
deallocate_thermo.o \
do_pwscf.o \
do_berry.o \
ev_sub.o \
ev_xml.o \
find_symmetry.o \
generate_k_along_lines_kz.o \
initialize_thermo_work.o \
initialize_thermo_master.o \
local_dos1d_so.o \
matdyn_sub.o \
mur.o \
plan_avg_sub.o \
plot_bz.o \
plot_anhar.o \
plotband_sub.o \
plot_2d_bz.o \
plot_e_ke.o \
plot_e_nk.o \
plot_thermo.o \
plot_mur.o \
polifit.o \
print_polarization.o \
proj_band_structure.o \
q2r_sub.o \
read_state_densities.o \
run_thermo_asyncronously.o \
set_2d_bz_path.o \
set_bz_path.o \
set_elastic_cons_work.o \
set_fft_mesh.o \
set_k_points.o \
set_paths_disp.o \
set_piezo_tensor_work.o \
set_temperature.o \
set_thermo_work_todo.o \
simple_nonanal.o \
simple_plot.o \
solve_linter.o \
sym_band_sub.o \
thermo_mod.o \
thermo_readin.o \
thermo_pw.o \
thermo_summary.o \
write_anhar.o \
write_e_nk.o \
write_e_ke.o  \
write_ev_input.o \
write_gruneisen_band.o \
write_read_energy.o \
write_thermo.o 

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
	- ( cd ../../bin ; ln -fs ../thermo_pw/src/thermo_pw.x . )

libthermo_pw.a : $(THERMOOBJS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *.F90 *.d *.mod *.i *.L
	- /bin/rm -rf ../../bin/thermo_pw.x 

include make.depend
# DO NOT DELETE
