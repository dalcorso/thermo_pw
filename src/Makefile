# Makefile for thermo_pw 

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../LAXlib  \
          $(MOD_FLAG)../../PW/src $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG)../../PHonon/PH $(MOD_FLAG)../lib $(MOD_FLAG)../qe $(MOD_FLAG).

THERMOOBJS = \
adjust_tau.o \
addusdens1d_so.o \
allocate_thermo.o \
bands_sub.o \
bcast_thermo_input.o \
check_existence.o \
check_file_exists.o \
clean_dfft.o \
collect_everything.o \
compute_cp.o \
compute_density.o \
compute_gruneisen.o \
compute_gruneisen_anis.o \
compute_omega_geo.o \
deallocate_asyn.o \
deallocate_thermo.o \
deallocate_thermo.o \
do_pwscf.o \
do_berry.o \
dos_sub.o \
elastic_constants_t.o \
ev_sub.o \
float_to_char.o \
find_symmetry.o \
generate_k_along_lines_kz.o \
initialize_thermo_work.o \
interpolate_cv.o \
local_dos1d_so.o \
matdyn_sub.o \
plan_avg_sub.o \
plot_bz.o \
plot_anhar.o \
plot_anhar_anis.o \
plotband_sub.o \
plot_2d_bz.o \
plot_dispersion.o \
plot_dos.o \
plot_e_ke.o \
plot_e_nk.o \
plot_formf.o \
plot_gruneisen_band_anis.o \
plot_elastic_t.o \
plot_epsilon_omega.o \
plot_multi_energy.o \
plot_mur.o \
plot_phdos.o \
plot_thermo.o \
plot_xrdp.o \
print_polarization.o \
proj_band_structure.o \
q2r_sub.o \
quadratic_fit.o \
readmodes.o \
read_state_densities.o \
run_thermo_asyncronously.o \
set_2d_bz_path.o \
set_bz_path.o \
set_dos_kpoints.o \
set_elastic_cons_work.o \
set_files_names.o \
set_fft_mesh.o \
set_k_points.o \
set_paths_disp.o \
set_piezo_tensor_work.o \
set_temperature.o \
set_thermo_work_todo.o \
simple_nonanal.o \
simple_plot.o \
sym_band_sub.o \
thermo_mod.o \
thermo_readin.o \
thermo_pw.o \
thermo_summary.o \
check_el_cons.o \
write_anhar.o \
write_anhar_anis.o \
write_e_nk.o \
write_e_ke.o  \
write_el_thermo.o \
write_ev_input.o \
write_elastic_t.o \
write_gruneisen_band.o \
write_gruneisen_band_anis.o \
write_gnuplot_energy.o \
write_minimum_energy_data.o \
write_mur.o \
write_phdos.o \
write_ph_freq.o \
write_ph_dispersions.o \
write_quadratic.o \
write_thermo.o 

LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PHonon/PH/libph.a ../../PP/src/libpp.a ../../PW/src/libpw.a 
QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a

LIBOBJS	= ../lib/thermo_lib.a ../qe/thermo_qe.a ../lib/thermo_lib.a ../../PHonon/PH/libph.a ../../clib/clib.a ../../iotk/src/libiotk.a

all : libthermo_pw.a thermo_pw.x 

libs-ph : libph.a libphaux.a

thermo_pw.x : thermo_pw.o libthermo_pw.a $(LRMODS) $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ thermo_pw.o ../../PHonon/PH/rigid.o \
                 ../../PHonon/PH/write_eigenvectors.o libthermo_pw.a \
         $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBS) ../../clib/clib.a
	- ( cd ../../bin ; ln -fs ../thermo_pw/src/thermo_pw.x . )

libthermo_pw.a : $(THERMOOBJS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *.F90 *.d *.mod *.i *.L
	- /bin/rm -rf ../../bin/thermo_pw.x 

include make.depend
# DO NOT DELETE
