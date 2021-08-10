# Makefile for thermo_pw 

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src  \
          $(MOD_FLAG)../../dft-d3  \
          $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG)../../PHonon/PH \
          $(MOD_FLAG)../lib \
          $(MOD_FLAG)../qe \
          $(MOD_FLAG).

#objects of libthermo_pw.a
THERMOOBJS = \
addusdens1d_so.o \
adjust_tau.o \
allocate_thermo.o \
bands_sub.o \
bcast_thermo_input.o \
b_factor_reduction.o \
check_all_geometries_done.o \
check_dynmat_all_geo_on_file.o \
check_el_cons.o \
check_el_cons_qha.o \
check_existence.o \
check_file_exists.o \
check_initial_status_all_geo.o \
check_thermo_all_geo.o \
clean_bfgs_history.o \
clean_dfft.o \
clean_ngeo.o \
clean_path.o \
collect_info.o \
collect_the_info.o \
compute_cp.o \
compute_density.o \
compute_el_comp_t.o \
compute_eref_band.o \
compute_omega_geo.o \
deallocate_asyn.o \
deallocate_thermo.o \
deallocate_q2r.o \
decorated_write.o \
do_pwscf.o \
do_berry.o \
do_phonon.o \
dos_sub.o \
dynmatrix_tpw.o \
elastic_constants_t.o \
ev_sub.o \
evaluate_celldm.o \
evaluate_vm.o \
fast_phq_readin.o \
find_band_sym_proj.o \
find_mode_sym_proj.o \
find_representations_mode_q.o \
find_symmetry.o \
fit_frequencies.o \
float_to_char.o \
gen_qpoints.o \
generalized_phdos.o \
generate_k_along_lines_kz.o \
initialize_elastic_cons.o \
initialize_thermo_work.o \
interpolate_bfact.o \
interpolate_el_cons.o \
interpolate_free_ener_strain.o \
interpolate_freq.o \
interpolate_thermo.o \
io_bands.o \
local_dos1d_so.o \
manage_anhar.o \
manage_anhar_routines.o \
manage_bands.o \
manage_collection.o \
manage_el_free_energy.o \
manage_elastic_cons.o \
manage_elastic_cons_qha.o \
manage_elastic_cons_qha_2.o \
manage_energy_minimum.o \
manage_piezo_tensor.o \
manage_ph.o \
manage_ph_all_geometries.o \
manage_ph_postproc.o \
manage_ph_run.o \
manage_ph_symmetry.o \
manage_surface_states.o \
manage_syn_1.o \
manage_syn_2.o \
manage_sync_pw.o \
manage_xrdp.o \
matdyn_sub.o \
plan_avg_sub.o \
plot_anhar.o \
plot_anhar_anis.o \
plot_bz.o \
plot_2d_bz.o \
plot_dispersion.o \
plot_dos.o \
plot_e_ke.o \
plot_e_nk.o \
plot_ef.o \
plot_elastic_t.o \
plot_epsilon_omega.o \
plot_formf.o \
plot_geo_p.o \
plot_gruneisen_band_anis.o \
plot_macro_el_t.o \
plot_multi_energy.o \
plot_mur.o \
plot_optical_const.o \
plot_phdos.o \
plot_thermo.o \
plot_xrdp.o \
plotband_sub.o \
prepare_do_phonon.o \
print_polarization.o \
proj_band_structure.o \
q2r_sub.o \
quadratic_fit.o \
read_elastic_qha.o \
read_minimal_info.o \
read_state_densities.o \
readmodes.o \
run_iwork.o \
run_thermo_asyncronously.o \
save_ph_variables.o \
set_2d_bz_path.o \
set_bz_path.o \
set_dos_kpoints.o \
set_elastic_constants_t.o \
set_fft_mesh.o \
set_files_names.o \
set_k_points.o \
set_paths_disp.o \
set_piezo_tensor_work.o \
set_temperature.o \
set_thermo_work_todo.o \
set_work_for_elastic_const.o \
simple_nonanal.o \
simple_plot.o \
surface.o \
sym_band_sub.o \
thermo_mod.o \
thermo_readin.o \
thermo_setup.o \
thermo_startup.o \
thermo_summary.o \
write_anhar.o \
write_anhar_anis.o \
write_e_ke.o  \
write_e_nk.o \
write_e_omega.o \
write_el_thermo.o \
write_elastic_t.o \
write_elastic_t_qha.o \
write_ev_input.o \
write_gnuplot_energy.o \
write_gruneisen_band.o \
write_gruneisen_band_anis.o \
write_minimum_energy_data.o \
write_mur.o \
write_ph_freq.o \
write_ph_dispersions.o \
write_phdos.o \
write_quadratic.o \
write_thermo.o 

#libraries of QE
PWOBJS = ../../PHonon/PH/libph.a \
         ../../PP/src/libpp.a \
         ../../LR_Modules/liblrmod.a \
         ../../PW/src/libpw.a 

#library modules of QE
QEMODS = ../../Modules/libqemod.a \
         ../../KS_Solvers/libks_solvers.a \
         ../../upflib/libupf.a    \
         ../../XClib/xc_lib.a     \
         ../../FFTXlib/libqefft.a \
         ../../LAXlib/libqela.a   \
         ../../UtilXlib/libutil.a \
         ../../dft-d3/libdftd3qe.a

#libraries of thermo_pw
LIBTHERMO = ../qe/thermo_qe.a \
            ../lib/thermo_lib.a 

MODULES = $(LIBTHERMO) $(PWOBJS) $(QEMODS) $(LIBOBJS)

all : thermo_pw.x 

thermo_pw.x : thermo_pw.o libthermo_pw.a $(MODULES)
	$(LD) $(LDFLAGS) -o $@ ./thermo_pw.o ./libthermo_pw.a \
              ../../PHonon/PH/rigid.o ../../PHonon/PH/write_eigenvectors.o \
        $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/src/thermo_pw.x . )

libthermo_pw.a : $(THERMOOBJS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -rf ../../bin/thermo_pw.x 

include make.depend
# DO NOT DELETE
