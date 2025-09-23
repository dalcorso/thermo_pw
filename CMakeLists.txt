set(src_tpw_lib
    lib/asy.f90
    lib/atomic_pos.f90
    lib/bz_2d_form.f90
    lib/bz_asy.f90
    lib/color.f90
    lib/cubic_surfaces.f90
    lib/debye_module.f90
    lib/diagonalize.f90
    lib/elastic_const.f90
    lib/eldos_module.f90
    lib/eos.f90
    lib/gnuplot_color.f90
    lib/gnuplot.f90
    lib/isoentropic.f90
    lib/latgen_2d.f90
    lib/lattices.f90
    lib/linear_solvers.f90
    lib/linear_surfaces.f90
    lib/magnetic_point_group.f90
    lib/model_free_energy.f90
    lib/mp_asyn.f90
    lib/numerical_diff.f90
    lib/nye.f90
    lib/optical.f90
    lib/phdos_module.f90
    lib/ph_freq_module.f90
    lib/piezoelectric_tensor.f90
    lib/point_group.f90
    lib/polyfit.f90
    lib/polynomial.f90
    lib/quadratic_surfaces.f90
    lib/quartic_surfaces.f90
    lib/rotate.f90
    lib/space_groups.f90
    lib/strain.f90
    lib/thermodynamics.f90
    lib/vector.f90
    lib/voigt.f90
    lib/xrdp.f90)
qe_enable_cuda_fortran("${src_tpw_lib}")

qe_add_library(qe_tpw_lib ${src_tpw_lib})
target_link_libraries(qe_tpw_lib
    PRIVATE
        qe_pw
        qe_lr_modules
        qe_modules
        qe_fftx
        qe_upflib
        qe_xclib
        qe_openmp_fortran
        qe_mpi_fortran)

set(src_tpw_fftpack5.1
    fftpack5.1/cfft1b_dev.f90
    fftpack5.1/cfft1f_dev.f90
    fftpack5.1/cfft1i_dev.f90
    fftpack5.1/c1f2kb_dev.f90
    fftpack5.1/c1f2kf_dev.f90
    fftpack5.1/c1f3kb_dev.f90
    fftpack5.1/c1f3kf_dev.f90
    fftpack5.1/c1f4kb_dev.f90
    fftpack5.1/c1f4kf_dev.f90
    fftpack5.1/c1f5kf_dev.f90
    fftpack5.1/c1f5kb_dev.f90
    fftpack5.1/c1fgkb_dev.f90
    fftpack5.1/c1fgkf_dev.f90
    fftpack5.1/c1fm1b_dev.f90
    fftpack5.1/c1fm1f_dev.f90
    fftpack5.1/r8_factor_dev.f90
    fftpack5.1/r8_mcfti1_dev.f90
    fftpack5.1/r8_tables_dev.f90
    fftpack5.1/xerfft_dev.f90)
qe_enable_cuda_fortran("${src_tpw_fftpack5.1}")

qe_add_library(qe_tpw_fftpack5.1 ${src_tpw_fftpack5.1})
target_link_libraries(qe_tpw_fftpack5.1)
target_include_directories(qe_tpw_fftpack5.1 PRIVATE fftpack5.1)

set(src_tpw_lapack
    lapack/dcabs1_dev.f90
    lapack/dcopy_dev.f90
    lapack/disnan_dev.f90
    lapack/dladiv1_dev.f90
    lapack/dladiv2_dev.f90
    lapack/dladiv_dev.f90
    lapack/dlae2_dev.f90
    lapack/dlaebz_dev.f90
    lapack/dlaev2_dev.f90
    lapack/dlagtf_dev.f90
    lapack/dlagts_dev.f90
    lapack/dlaisnan_dev.f90
    lapack/dlamc3_dev.f90
    lapack/dlamch_dev.f90
    lapack/dlanst_dev.f90
    lapack/dlapy2_dev.f90
    lapack/dlapy3_dev.f90
    lapack/dlarnv_dev.f90
    lapack/dlartg_dev.f90
    lapack/dlaruv_dev.f90
    lapack/dlascl_dev.f90
    lapack/dlasrt_dev.f90
    lapack/dlassq_dev.f90
    lapack/dnrm2_dev.f90
    lapack/dscal_dev.f90
    lapack/dstebz_dev.f90
    lapack/dsterf_dev.f90
    lapack/dznrm2_dev.f90
    lapack/idamax_dev.f90
    lapack/ieeeck_dev.f90
    lapack/ilaenv_dev.f90
    lapack/ilazlc_dev.f90
    lapack/ilazlr_dev.f90
    lapack/iparmq_dev.f90
    lapack/la_constants_dev.f90
    lapack/la_xisnan_dev.f90
    lapack/lsame_dev.f90
    lapack/xerbla_dev.f90
    lapack/zaxpy_dev.f90
    lapack/zcopy_dev.f90
    lapack/zdotc_dev.f90
    lapack/zdscal_dev.f90
    lapack/zgemm_dev.f90
    lapack/zgemv_dev.f90
    lapack/zgerc_dev.f90
    lapack/zheev_dev.f90
    lapack/zheevx_dev.f90
    lapack/zhegs2_dev.f90
    lapack/zhegst_dev.f90
    lapack/zhegv_dev.f90
    lapack/zhegvx_dev.f90
    lapack/zhemm_dev.f90
    lapack/zhemv_dev.f90
    lapack/zher2_dev.f90
    lapack/zher2k_dev.f90
    lapack/zherk_dev.f90
    lapack/zhetd2_dev.f90
    lapack/zhetrd_dev.f90
    lapack/zlacgv_dev.f90
    lapack/zlacpy_dev.f90
    lapack/zladiv_dev.f90
    lapack/zlanhe_dev.f90
    lapack/zlarfb_dev.f90
    lapack/zlarf_dev.f90
    lapack/zlarfg_dev.f90
    lapack/zlarft_dev.f90
    lapack/zlartg_dev.f90
    lapack/zlascl_dev.f90
    lapack/zlaset_dev.f90
    lapack/zlasr_dev.f90
    lapack/zlassq_dev.f90
    lapack/zlatrd_dev.f90
    lapack/zpotrf2_x_dev.f90
    lapack/zpotrf_dev.f90
    lapack/zscal_dev.f90
    lapack/zstein_dev.f90
    lapack/zsteqr_dev.f90
    lapack/zswap_dev.f90
    lapack/ztrmm_dev.f90
    lapack/ztrmv_dev.f90
    lapack/ztrsm_dev.f90
    lapack/ztrsv_dev.f90
    lapack/zung2l_dev.f90
    lapack/zung2r_dev.f90
    lapack/zungql_dev.f90
    lapack/zungqr_dev.f90
    lapack/zungtr_dev.f90
    lapack/zunm2l_dev.f90
    lapack/zunm2r_dev.f90
    lapack/zunmql_dev.f90
    lapack/zunmqr_dev.f90
    lapack/zunmtr_dev.f90)

qe_enable_cuda_fortran("${src_tpw_lapack}")

qe_add_library(qe_tpw_lapack ${src_tpw_lapack})

if (QE_ENABLE_CUDA) 
   target_compile_options(qe_tpw_lapack PUBLIC "-cuda=charstring")
endif()
target_include_directories(qe_tpw_lapack PRIVATE lapack)

set(src_tpw_qe
    qe/add_dkmds_tpw.f90
    qe/add_dvscf_rhs.f90
    qe/adddvscf_tran.f90
    qe/addnlcc_zstar_eu_us_tpw.f90
    qe/add_vuspsi_tpw.f90
    qe/add_zstar_ue_tpw.f90
    qe/add_zstar_us.f90
    qe/allocate_cg.f90
    qe/allocate_lanczos.f90
    qe/allocate_phq_tpw.f90
    qe/apply_ac.f90
    qe/apply_trev_r.f90
    qe/calc_chi.f90
    qe/c_bands_tpw.f90
    qe/c_bands_many_k.f90
    qe/ccg_psi.f90
    qe/cg_psi_many_k.f90
    qe/ch_psi_all_many_k.f90
    qe/cgsolve_all_many_k.f90
    qe/check_bands.f90
    qe/check_initial_geometry.f90
    qe/check_initial_status_tpw.f90
    qe/ch_psi_all_dev.f90
    qe/cegterg_vk.f90
    qe/cegterg_vk_dev.f90
    qe/cgsolve_all_dev.f90
    qe/clean_all_pw.f90
    qe/clean_pw_tpw.f90
    qe/collect_grid_files.f90
    qe/compute_deff_dev.f90
    qe/compute_dvloc_tpw.f90
    qe/compute_drhous_nc_tpw.f90
    qe/compute_drhous_tpw.f90
    qe/compute_int3_coeff.f90
    qe/compute_intq.f90
    qe/control_qe.f90
    qe/deallocate_ph_tpw.f90
    qe/dielec_tpw.f90
    qe/diago_dev.f90
    qe/divide_class.f90
    qe/divide_class_so.f90
    qe/do_cg_e.f90
    qe/do_cg_ph.f90
    qe/do_lanczos.f90
    qe/do_lanczos_psh.f90
    qe/drhodvus.f90
    qe/drho_tpw.f90
    qe/dveqpsi_us.f90
    qe/dveqpsi_us_only.f90
    qe/dvkb3_tpw.f90
    qe/dvqpsi_us_only_tpw.f90
    qe/dvqpsi_us_many_k.f90
    qe/dvqpsi_dev.f90
    qe/dv_of_drho_clf.f90
    qe/dv_of_drho_tran.f90
    qe/dvqpsi_us_tpw.f90
    qe/dyndiag.f90
    qe/ef_shift.f90
    qe/electrons_tpw.f90
    qe/ev_xml.f90
    qe/find_band_sym.f90
    qe/find_info_group.f90
    qe/find_inverse_s.f90
    qe/find_irrep_tpw.f90
    qe/find_mode_sym.f90
    qe/float_to_char.f90
    qe/g_psi_tpw.f90
    qe/g_psi_dev.f90
    qe/h_pcg_step.f90
    qe/h_psi_dev.f90
    qe/h_psi_tpw.f90
    qe/h_psi_gpu_tpw.f90
    qe/hdiag_dev.f90
    qe/incdrhoscf_dev.f90
    qe/initialize_ph_tpw.f90
    qe/init_representations_tpw.f90
    qe/init_run_tpw.f90
    qe/init_us_2_kernel.f90
    qe/irrek.f90
    qe/ke_g2kin_dev.f90
    qe/kpoint_grid.f90
    qe/kpoint_grid_serial_tpw.f90
    qe/lanczos_write_restart.f90
    qe/lr_dot_tpw.f90
    qe/lr_restart.f90
    qe/lr_sm1_psi_tpw.f90
    qe/lr_sym_mod_tpw.f90
    qe/lr_us.f90
    qe/manage_mixing.f90
    qe/many_k.f90
    qe/many_k_ph.f90
    qe/mix_pot.f90
    qe/nh_lanczos_step.f90
    qe/non_scf_tpw.f90
    qe/openfilq_tpw.f90
    qe/orthogonalize_omega.f90
    qe/orthogonalize_tpw.f90
    qe/orthog_dev.f90
    qe/paw_add_onecenter.f90
    qe/paw_add_symmetry.f90
    qe/pheqscf.f90
    qe/phescf_tpw.f90
    qe/phq_init_tpw.f90
    qe/phq_readin_tpw.f90
    qe/phqscf_tpw.f90
    qe/phq_setup_tpw.f90
    qe/phq_summary_tpw.f90
    qe/polarizc.f90
    qe/prepare_q_tpw.f90
    qe/prepare_sym_analysis.f90
    qe/print_gpu_memory.f90
    qe/print_symmetries.f90
    qe/psh_lanczos_step.f90
    qe/q2qstar_ph.f90
    qe/q_points_tpw.f90
    qe/read_dyn_from_file.f90
    qe/read_lanczos_chain.f90
    qe/realspace_grid_init.f90
    qe/rotate_all_psi.f90
    qe/rotate_all_psi_r.f90
    qe/rotate_and_add_dyn.f90
    qe/rotate_charges.f90
    qe/run_nscf_tpw.f90
    qe/s_psi_ch.f90
    qe/s_psi_tpw.f90
    qe/scal_prod.f90
    qe/select_c_bands.f90
    qe/set_defaults_pw.f90
    qe/set_fxc_tran.f90
    qe/set_int3q.f90
    qe/set_intq_nc.f90
    qe/set_irr_sym_tpw.f90
    qe/set_irr_tpw.f90
    qe/set_kplusq.f90
    qe/set_kplusq_nc.f90
    qe/set_kup_and_kdw.f90
    qe/setup_tpw.f90
    qe/solve_e_fpolc.f90
    qe/solve_eq.f90
    qe/solve_eq_tran.f90
    qe/solve_e_tpw.f90
    qe/solve_linear_system.f90
    qe/solve_linear_system_many_k.f90
    qe/solve_linter_tpw.f90
    qe/solve_linter_many_k.f90
    qe/s_psi_gpu_tpw.f90
    qe/s_psi_ch_gpu.f90
    qe/star_q.f90
    qe/summarize_alpha.f90
    qe/sym_and_write_magnetic_charges.f90
    qe/sym_def.f90
    qe/symmetrize_drho.f90
    qe/transform_intq_nc.f90
    qe/transform_intq_so.f90
    qe/tra_write_zstar.f90
    qe/vloc_psi_tpw.f90
    qe/vloc_psi_gpu_tpw.f90
    qe/write_group_info_ph.f90
    qe/write_rec_tpw.f90
    qe/ylm2_dev.f90
    qe/ylmr2_dev.f90
    qe/zstar_eu_loc.f90
    qe/zstar_eu_tpw.f90
    qe/zstar_eu_us_tpw.f90)
qe_enable_cuda_fortran("${src_tpw_qe}")

qe_add_library(qe_tpw_qe ${src_tpw_qe})
target_include_directories(qe_tpw_qe PRIVATE qe lapack fftpack5.1)
target_link_libraries(qe_tpw_qe
    PRIVATE
        qe_tpw_fftpack5.1
        qe_tpw_lapack
        qe_tpw_lib
        qe_pw
        qe_phonon_ph
        qe_lr_modules
        qe_modules
        qe_fftx
        qe_openmp_fortran
        qe_mpi_fortran
        qe_dftd3
        qe_devxlib
        qe_xclib
        qe_mbd
        qe_lax
    PUBLIC
        qe_openacc_fortran)

set(src_tpw
    src/addusdens1d_so.f90
    src/adjust_tau.f90
    src/allocate_thermo.f90
    src/bands_sub.f90
    src/bcast_thermo_input.f90
    src/b_factor_reduction.f90
    src/check_all_geometries_done.f90
    src/check_dynmat_all_geo_on_file.f90
    src/check_el_cons.f90
    src/check_el_cons_qha.f90
    src/check_existence.f90
    src/check_file_exists.f90
    src/check_initial_status_all_geo.f90
    src/check_thermo_all_geo.f90
    src/clean_bfgs_history.f90
    src/clean_dfft.f90
    src/clean_ngeo.f90
    src/clean_path.f90
    src/collect_info.f90
    src/collect_the_info.f90
    src/compute_cp.f90
    src/compute_density.f90
    src/compute_el_comp_t.f90
    src/compute_eref_band.f90
    src/compute_omega_geo.f90
    src/deallocate_asyn.f90
    src/deallocate_q2r.f90
    src/deallocate_thermo.f90
    src/decorated_write.f90
    src/do_berry.f90
    src/do_phonon_tpw.f90
    src/do_pwscf.f90
    src/dos_ph.f90
    src/dos_sub.f90
    src/dynmatrix_tpw.f90
    src/empirical_free_energy.f90
    src/evaluate_celldm.f90
    src/evaluate_vm.f90
    src/ev_mod.f90
    src/ev_sub.f90
    src/fast_phq_readin.f90
    src/find_band_sym_proj.f90
    src/find_mode_sym_proj.f90
    src/find_representations_mode_q.f90
    src/find_t_debye.f90
    src/find_symmetry.f90
    src/fit_frequencies.f90
    src/generalized_phdos.f90
    src/generate_k_along_lines_kz.f90
    src/gen_qpoints.f90
    src/geometry_file.f90
    src/initialize_elastic_cons.f90
    src/initialize_thermo_work.f90
    src/interpolate_bfact.f90
    src/interpolate_el_cons.f90
    src/interpolate_free_ener_strain.f90
    src/interpolate_freq.f90
    src/interpolate_thermo.f90
    src/io_bands.f90
    src/local_dos1d_so.f90
    src/manage_anhar.f90
    src/manage_anhar_routines.f90
    src/manage_bands.f90
    src/manage_collection.f90
    src/manage_elastic_cons.f90
    src/manage_elastic_cons_qha_2.f90
    src/manage_elastic_cons_qha.f90
    src/manage_el_free_energy.f90
    src/manage_energy_minimum.f90
    src/manage_hot_electrons.f90
    src/manage_ph_all_geometries.f90
    src/manage_ph.f90
    src/manage_ph_postproc.f90
    src/manage_ph_run.f90
    src/manage_ph_symmetry.f90
    src/manage_piezo_tensor.f90
    src/manage_plot_anhar.f90
    src/manage_surface_states.f90
    src/manage_syn_1.f90
    src/manage_syn_2.f90
    src/manage_sync_pw.f90
    src/manage_xrdp.f90
    src/matdyn_sub.f90
    src/plan_avg_sub.f90
    src/plot_2d_bz.f90
    src/plot_anhar_anis.f90
    src/plot_anhar.f90
    src/plotband_sub.f90
    src/plot_bz.f90
    src/plot_dispersion.f90
    src/plot_dos.f90
    src/plot_anhar_dw.f90
    src/plot_ef.f90
    src/plot_e_ke.f90
    src/plot_elastic_t.f90
    src/plot_e_nk.f90
    src/plot_epsilon_omega.f90
    src/plot_formf.f90
    src/plot_geo_p.f90
    src/plot_gruneisen_band_anis.f90
    src/plot_macro_el_t.f90
    src/plot_multi_energy.f90
    src/plot_mur.f90
    src/plot_optical_const.f90
    src/plot_phdos.f90
    src/plot_thermo.f90
    src/plot_xrdp.f90
    src/prepare_do_phonon.f90
    src/print_thermo_clocks.f90
    src/print_polarization.f90
    src/proj_band_structure.f90
    src/q2r_sub.f90
    src/quadratic_fit.f90
    src/read_minimal_info.f90
    src/readmodes.f90
    src/read_state_densities.f90
    src/redefine_energies.f90
    src/run_iwork.f90
    src/run_thermo_asyncronously.f90
    src/save_ph_variables.f90
    src/set_2d_bz_path.f90
    src/set_bz_path.f90
    src/set_dos_kpoints.f90
    src/set_elastic_constants_t.f90
    src/set_fft_mesh.f90
    src/set_files_names.f90
    src/set_geometry_el_cons.f90
    src/set_k_points.f90
    src/set_paths_disp.f90
    src/set_piezo_tensor_work.f90
    src/set_temperature.f90
    src/set_thermo_work_todo.f90
    src/set_work_for_elastic_const.f90
    src/simple_nonanal.f90
    src/simple_plot.f90
    src/summarize_anhar.f90
    src/surface.f90
    src/sym_band_sub.f90
    src/thermo_mod.f90
    src/thermo_pw.f90
    src/thermo_readin.f90
    src/thermo_setup.f90
    src/thermo_startup.f90
    src/thermo_summary.f90
    src/write_anhar_anis.f90
    src/write_anhar_aux.f90
    src/write_anhar.f90
    src/write_anhar_glob.f90
    src/write_e_ke.f90
    src/write_elastic_qha.f90
    src/write_elastic_t.f90
    src/write_elastic_t_qha.f90
    src/write_el_thermo.f90
    src/write_e_nk.f90
    src/write_e_omega.f90
    src/write_ev_input.f90
    src/write_free_energy.f90
    src/write_gnuplot_energy.f90
    src/write_gruneisen_band_anis.f90
    src/write_gruneisen_band.f90
    src/write_minimum_energy_data.f90
    src/write_mur.f90
    src/write_ph_dispersions.f90
    src/write_phdos.f90
    src/write_ph_freq.f90
    src/write_quadratic.f90
    src/write_thermo.f90)
qe_enable_cuda_fortran("${src_tpw}")

qe_add_library(qe_tpw ${src_tpw})
target_link_libraries(qe_tpw
    PRIVATE
        qe_tpw_qe
        qe_tpw_lib
        qe_pw
        qe_phonon_ph
        qe_modules
        qe_lr_modules
	qe_scalapack
        qe_fftx
        qe_upflib
        qe_openmp_fortran
        qe_mpi_fortran
        qe_devxlib
        qe_xclib
        qe_mbd
    PUBLIC
        qe_openacc_fortran)


###########################################################
# thermo_pw.x
###########################################################
set(src_tpw_x src/thermo_pw.f90)
qe_enable_cuda_fortran("${src_tpw_x}")
qe_add_executable(qe_tpw_exe ${src_tpw_x})
set_target_properties(qe_tpw_exe 
    PROPERTIES 
        OUTPUT_NAME thermo_pw.x)
target_link_libraries(qe_tpw_exe
    PRIVATE
        qe_tpw_lapack
        qe_tpw_fftpack5.1
        qe_tpw
        qe_tpw_qe
        qe_tpw_lib
        qe_lr_modules
        qe_phonon_phaux
        qe_phonon_ph
        qe_pw
        qe_pp
        qe_modules
        qe_mbd
    PUBLIC
        qe_openacc_fortran)


###########################################################
# space_groups.x
###########################################################
set(src_space_groups_x tools/space_groups.f90)
qe_enable_cuda_fortran("${src_space_groups_x}")
qe_add_executable(qe_tpw_space_groups_exe ${src_space_groups_x})
set_target_properties(qe_tpw_space_groups_exe 
    PROPERTIES 
        OUTPUT_NAME space_groups.x)
target_link_libraries(qe_tpw_space_groups_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# crystal_point_groups.x
###########################################################
set(src_crystal_point_group_x tools/crystal_point_group.f90)
qe_enable_cuda_fortran("${src_crystal_point_group_x}")
qe_add_executable(qe_tpw_crystal_point_group_exe ${src_crystal_point_group_x})
set_target_properties(qe_tpw_crystal_point_group_exe 
    PROPERTIES 
        OUTPUT_NAME crystal_point_group.x)
target_link_libraries(qe_tpw_crystal_point_group_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# mag_point_groups.x
###########################################################
set(src_mag_point_group_x tools/mag_point_group.f90)
qe_enable_cuda_fortran("${src_mag_point_group_x}")
qe_add_executable(qe_tpw_mag_point_group_exe ${src_mag_point_group_x})
set_target_properties(qe_tpw_mag_point_group_exe 
    PROPERTIES 
        OUTPUT_NAME mag_point_group.x)
target_link_libraries(qe_tpw_mag_point_group_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# units.x
###########################################################
set(src_units_x tools/units.f90)
qe_enable_cuda_fortran("${src_units_x}")
qe_add_executable(qe_tpw_units_exe ${src_units_x})
set_target_properties(qe_tpw_units_exe 
    PROPERTIES 
        OUTPUT_NAME units.x)
target_link_libraries(qe_tpw_units_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# hex_trig.x
###########################################################
set(src_hex_trig_x tools/hex_trig.f90)
qe_enable_cuda_fortran("${src_hex_trig_x}")
qe_add_executable(qe_tpw_hex_trig_exe ${src_hex_trig_x})
set_target_properties(qe_tpw_hex_trig_exe 
    PROPERTIES 
        OUTPUT_NAME hex_trig.x)
target_link_libraries(qe_tpw_hex_trig_exe
    PRIVATE
        qe_modules)

###########################################################
# bravais_lattices.x
###########################################################
set(src_bravais_lattices_x tools/bravais_lattices.f90)
qe_enable_cuda_fortran("${src_bravais_lattices_x}")
qe_add_executable(qe_tpw_bravais_lattices_exe ${src_bravais_lattices_x})
set_target_properties(qe_tpw_bravais_lattices_exe 
    PROPERTIES 
        OUTPUT_NAME bravais_lattices.x)
target_link_libraries(qe_tpw_bravais_lattices_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# kovalev.x
###########################################################
set(src_kovalev_x tools/kovalev.f90)
qe_enable_cuda_fortran("${src_kovalev_x}")
qe_add_executable(qe_tpw_kovalev_exe ${src_kovalev_x})
set_target_properties(qe_tpw_kovalev_exe 
    PROPERTIES 
        OUTPUT_NAME kovalev.x)
target_link_libraries(qe_tpw_kovalev_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# supercell.x
###########################################################

set(src_supercell_x tools/supercell.f90)
qe_enable_cuda_fortran("${src_supercell_x}")
qe_add_executable(qe_tpw_supercell_exe ${src_supercell_x})
set_target_properties(qe_tpw_supercell_exe 
    PROPERTIES 
        OUTPUT_NAME supercell.x)
target_link_libraries(qe_tpw_supercell_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# translate.x
###########################################################

set(src_translate_x tools/translate.f90)
qe_enable_cuda_fortran("${src_translate_x}")
qe_add_executable(qe_tpw_translate_exe ${src_translate_x})
set_target_properties(qe_tpw_translate_exe 
    PROPERTIES 
        OUTPUT_NAME translate.x)
target_link_libraries(qe_tpw_translate_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# elastic.x
###########################################################

set(src_elastic_x tools/elastic.f90)
qe_enable_cuda_fortran("${src_elastic_x}")
qe_add_executable(qe_tpw_elastic_exe ${src_elastic_x})
set_target_properties(qe_tpw_elastic_exe 
    PROPERTIES 
        OUTPUT_NAME elastic.x)
target_link_libraries(qe_tpw_elastic_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# optical.x
###########################################################

set(src_optical_x tools/optical.f90)
qe_enable_cuda_fortran("${src_optical_x}")
qe_add_executable(qe_tpw_optical_exe ${src_optical_x})
set_target_properties(qe_tpw_optical_exe 
    PROPERTIES 
        OUTPUT_NAME optical.x)
target_link_libraries(qe_tpw_optical_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# plot_surf_states.x
###########################################################

set(src_plot_sur_states_x tools/plot_sur_states.f90)
qe_enable_cuda_fortran("${src_plot_sur_states_x}")
qe_add_executable(qe_tpw_plot_sur_states_exe ${src_plot_sur_states_x})
set_target_properties(qe_tpw_plot_sur_states_exe 
    PROPERTIES 
        OUTPUT_NAME plot_sur_states.x)
target_link_libraries(qe_tpw_plot_sur_states_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# gener_2d_slab.x
###########################################################

set(src_gener_2d_slab_x tools/gener_2d_slab.f90)
qe_enable_cuda_fortran("${src_gener_2d_slab_x}")
qe_add_executable(qe_tpw_gener_2d_slab_exe ${src_gener_2d_slab_x})
set_target_properties(qe_tpw_gener_2d_slab_exe 
    PROPERTIES 
        OUTPUT_NAME gener_2d_slab.x)
target_link_libraries(qe_tpw_gener_2d_slab_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# gener_3d_slab.x
###########################################################

set(src_gener_3d_slab_x tools/gener_3d_slab.f90)
qe_enable_cuda_fortran("${src_gener_3d_slab_x}")
qe_add_executable(qe_tpw_gener_3d_slab_exe ${src_gener_3d_slab_x})
set_target_properties(qe_tpw_gener_3d_slab_exe 
    PROPERTIES 
        OUTPUT_NAME gener_3d_slab.x)
target_link_libraries(qe_tpw_gener_3d_slab_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# gener_nanowire.x
###########################################################

set(src_gener_nanowire_x tools/gener_nanowire.f90)
qe_enable_cuda_fortran("${src_gener_nanowire_x}")
qe_add_executable(qe_tpw_gener_nanowire_exe ${src_gener_nanowire_x})
set_target_properties(qe_tpw_gener_nanowire_exe 
    PROPERTIES 
        OUTPUT_NAME gener_nanowire.x)
target_link_libraries(qe_tpw_gener_nanowire_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# debye.x
###########################################################

set(src_debye_x tools/debye.f90)
qe_enable_cuda_fortran("${src_debye_x}")
qe_add_executable(qe_tpw_debye_exe ${src_debye_x})
set_target_properties(qe_tpw_debye_exe 
    PROPERTIES 
        OUTPUT_NAME debye.x)
target_link_libraries(qe_tpw_debye_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# density.x
###########################################################

set(src_density_x tools/density.f90)
qe_enable_cuda_fortran("${src_density_x}")
qe_add_executable(qe_tpw_density_exe ${src_density_x})
set_target_properties(qe_tpw_density_exe 
    PROPERTIES 
        OUTPUT_NAME density.x)
target_link_libraries(qe_tpw_density_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# rotate_tensor.x
###########################################################

set(src_rotate_tensors_x tools/rotate_tensors.f90)
qe_enable_cuda_fortran("${src_rotate_tensors_x}")
qe_add_executable(qe_tpw_rotate_tensors_exe ${src_rotate_tensors_x})
set_target_properties(qe_tpw_rotate_tensors_exe 
    PROPERTIES 
        OUTPUT_NAME rotate_tensors.x)
target_link_libraries(qe_tpw_rotate_tensors_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# test_colors.x
###########################################################

set(src_test_colors_x tools/test_colors.f90)
qe_enable_cuda_fortran("${src_test_colors_x}")
qe_add_executable(qe_tpw_test_colors_exe ${src_test_colors_x})
set_target_properties(qe_tpw_test_colors_exe 
    PROPERTIES 
        OUTPUT_NAME test_colors.x)
target_link_libraries(qe_tpw_test_colors_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# epsilon_tpw.x
###########################################################

set(src_epsilon_tpw_x tools/epsilon_tpw.f90)
qe_enable_cuda_fortran("${src_epsilon_tpw_x}")
qe_add_executable(qe_tpw_epsilon_tpw_exe ${src_epsilon_tpw_x})
set_target_properties(qe_tpw_epsilon_tpw_exe 
    PROPERTIES 
        OUTPUT_NAME epsilon_tpw.x)
target_link_libraries(qe_tpw_epsilon_tpw_exe
    PRIVATE
        qe_tpw_lib
        qe_pw
        qe_pp
        qe_lr_modules
        qe_modules)

###########################################################
# pdec.x
###########################################################

set(src_pdec_x tools/pdec.f90)
qe_enable_cuda_fortran("${src_pdec_x}")
qe_add_executable(qe_tpw_pdec_exe ${src_pdec_x})
set_target_properties(qe_tpw_pdec_exe 
    PROPERTIES 
        OUTPUT_NAME pdec.x)
target_link_libraries(qe_tpw_pdec_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# average_grun.x
###########################################################
set(src_average_grun_x tools/average_grun.f90)
qe_enable_cuda_fortran("${src_average_grun_x}")
qe_add_executable(qe_tpw_average_grun_exe ${src_average_grun_x})
set_target_properties(qe_tpw_average_grun_exe
    PROPERTIES
        OUTPUT_NAME average_grun.x)
target_link_libraries(qe_tpw_average_grun_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# change_dynmat_name.x
###########################################################
set(src_change_dynmat_name_x tools/change_dynmat_name.f90)
qe_enable_cuda_fortran("${src_change_dynmat_name_x}")
qe_add_executable(qe_tpw_change_dynmat_name_exe ${src_change_dynmat_name_x})
set_target_properties(qe_tpw_change_dynmat_name_exe
    PROPERTIES
        OUTPUT_NAME change_dynmat_name.x)
target_link_libraries(qe_tpw_change_dynmat_name_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# change_e_name.x
###########################################################
set(src_change_e_name_x tools/change_e_name.f90)
qe_enable_cuda_fortran("${src_change_e_name_x}")
qe_add_executable(qe_tpw_change_e_name_exe ${src_change_e_name_x})
set_target_properties(qe_tpw_change_e_name_exe
    PROPERTIES
        OUTPUT_NAME change_e_name.x)
target_link_libraries(qe_tpw_change_e_name_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# change_tel_name.x
###########################################################
set(src_change_tel_name_x tools/change_tel_name.f90)
qe_enable_cuda_fortran("${src_change_tel_name_x}")
qe_add_executable(qe_tpw_change_tel_name_exe ${src_change_tel_name_x})
set_target_properties(qe_tpw_change_tel_name_exe
    PROPERTIES
        OUTPUT_NAME change_tel_name.x)
target_link_libraries(qe_tpw_change_tel_name_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# emp_f.x
###########################################################
set(src_emp_f_x tools/emp_f.f90)
qe_enable_cuda_fortran("${src_emp_f_x}")
qe_add_executable(qe_tpw_emp_f_exe ${src_emp_f_x})
set_target_properties(qe_tpw_emp_f_exe
    PROPERTIES
        OUTPUT_NAME emp_f.x)
target_link_libraries(qe_tpw_emp_f_exe
    PRIVATE
        qe_tpw_lib
        qe_tpw_qe
        qe_modules)

###########################################################
# emp_g.x
###########################################################
set(src_emp_g_x tools/emp_g.f90)
qe_enable_cuda_fortran("${src_emp_g_x}")
qe_add_executable(qe_tpw_emp_g_exe ${src_emp_g_x})
set_target_properties(qe_tpw_emp_g_exe
    PROPERTIES
        OUTPUT_NAME emp_g.x)
target_link_libraries(qe_tpw_emp_g_exe
    PRIVATE
        qe_tpw_lib
        qe_tpw_qe
        qe_modules)

###########################################################
# test_eos.x
###########################################################
set(src_test_eos_x tools/test_eos.f90)
qe_enable_cuda_fortran("${src_test_eos_x}")
qe_add_executable(qe_tpw_test_eos_exe ${src_test_eos_x})
set_target_properties(qe_tpw_test_eos_exe
    PROPERTIES
        OUTPUT_NAME test_eos.x)
target_link_libraries(qe_tpw_test_eos_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# merge_interp.x
###########################################################
set(src_merge_interp_x tools/merge_interp.f90)
qe_enable_cuda_fortran("${src_merge_interp_x}")
qe_add_executable(qe_tpw_merge_interp_exe ${src_merge_interp_x})
set_target_properties(qe_tpw_merge_interp_exe
    PROPERTIES
        OUTPUT_NAME merge_interp.x)
target_link_libraries(qe_tpw_merge_interp_exe
    PRIVATE
        qe_tpw_lib
        qe_modules)

###########################################################
# wurtzite.x
###########################################################
set(src_wurtzite_x tools/wurtzite.f90)
qe_enable_cuda_fortran("${src_wurtzite_x}")
qe_add_executable(qe_tpw_wurtzite_exe ${src_wurtzite_x})
set_target_properties(qe_tpw_wurtzite_exe 
    PROPERTIES 
        OUTPUT_NAME wurtzite.x)
target_link_libraries(qe_tpw_wurtzite_exe
    PRIVATE
        qe_modules)

###########################################################
# analyze_dynmat.x
###########################################################
set(src_analyze_dynmat_x tools/analyze_dynmat.f90)
qe_enable_cuda_fortran("${src_analyze_dynmat_x}")
qe_add_executable(qe_tpw_analyze_dynmat_exe ${src_analyze_dynmat_x})
set_target_properties(qe_tpw_analyze_dynmat_exe 
    PROPERTIES 
        OUTPUT_NAME analyze_dynmat.x)
target_link_libraries(qe_tpw_analyze_dynmat_exe
    PRIVATE
        qe_modules)

###########################################################

qe_install_targets(
    qe_tpw_lapack
    qe_tpw_fftpack5.1
    qe_tpw_lib
    qe_tpw_qe
    qe_tpw
    qe_tpw_exe
    qe_tpw_analyze_dynmat_exe
    qe_tpw_average_grun_exe
    qe_tpw_bravais_lattices_exe
    qe_tpw_change_dynmat_name_exe
    qe_tpw_change_e_name_exe
    qe_tpw_change_tel_name_exe
    qe_tpw_crystal_point_group_exe
    qe_tpw_debye_exe
    qe_tpw_density_exe
    qe_tpw_elastic_exe
    qe_tpw_emp_f_exe
    qe_tpw_emp_g_exe
    qe_tpw_test_eos_exe
    qe_tpw_epsilon_tpw_exe
    qe_tpw_gener_2d_slab_exe
    qe_tpw_gener_3d_slab_exe
    qe_tpw_gener_nanowire_exe
    qe_tpw_hex_trig_exe
    qe_tpw_kovalev_exe
    qe_tpw_mag_point_group_exe
    qe_tpw_merge_interp_exe
    qe_tpw_optical_exe
    qe_tpw_pdec_exe
    qe_tpw_plot_sur_states_exe
    qe_tpw_rotate_tensors_exe
    qe_tpw_supercell_exe
    qe_tpw_space_groups_exe
    qe_tpw_test_colors_exe
    qe_tpw_translate_exe
    qe_tpw_units_exe
    qe_tpw_wurtzite_exe)
