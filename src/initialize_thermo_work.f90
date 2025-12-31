!
! Copyright (C) 2013-2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE initialize_thermo_work(nwork, part)
  !-----------------------------------------------------------------------
  !
  !  This routine is called by all images and initializes the global 
  !  variables that select the tasks to do. All images must have the same
  !  information.
  !  In addition to the variables that control the run it allocates, for
  !  the tasks that require it:
  !  ibrav_geo, celldm_geo, at_geo, omega_geo, energy_geo
  !  lpwscf, lstress, lberry, lphonon
  !
  USE kinds,          ONLY : DP
  USE thermo_mod,     ONLY : what, energy_geo, ef_geo, celldm_geo, ibrav_geo, &
                             omega_geo, tot_ngeo, no_ph, start_geometry,    &
                             last_geometry, iwho, tau_geo, at_geo, &
                             tot_ngeo_eos, no_ph_eos, uint_geo,    &
                             epsilon_infty_geo, zeu_geo, freq_geo, z_geo, &
                             epsilon_zero_geo, epsilon_zerom1_geo
  USE control_thermo, ONLY : lpwscf, lpwband, lphonon, lev_syn_1, lev_syn_2, &
                             lph, lef, lpwscf_syn_1, lbands_syn_1, lq2r,   &
                             ltherm, lconv_ke_test, lconv_nk_test,    &
                             lstress, lelastic_const, lpiezoelectric_tensor, &
                             lberry, lpolarization, lpart2_pw, do_scf_relax, &
                             ldos_syn_1, ltherm_dos, ltherm_freq, after_disp, &
                             lectqha, ltau_from_file, lpiezotqha
  USE control_pwrun,  ONLY : do_punch
  USE control_conv,   ONLY : nke, ke, deltake, nkeden, deltakeden, keden, &
                             nnk, nk_test, deltank, nsigma, sigma_test,   &  
                             deltasigma, ncutoffene
  USE equilibrium_conf, ONLY : celldm0, omega0
  USE control_atomic_pos, ONLY : max_nint_var
  USE initial_conf,   ONLY : ibrav_save, start_geometry_save,             &
                             last_geometry_save, omega_save
  USE piezoelectric_tensor, ONLY : allocate_piezo
  USE control_elastic_constants, ONLY : rot_mat, ngeom, use_free_energy,   &
                                   elalgen, work_base, start_geometry_qha, &
                                   last_geometry_qha, elastic_algorithm,   &
                                   el_con_omega_geo, el_con_geo,           &
                                   all_geometry_done_geo, found_dos_ec,    &
                                   found_ph_ec, min_y_t, min_yf_t, dyde,   &
                                   ngeo_strain, stype, tau_save_ec,        &
                                   el_con_celldm_geo, nint_var_ec,         &
                                   iconstr_internal_ec
  USE control_piezoelectric_tensor, ONLY : g_piezo_tensor_geo, &
                                   eg_piezo_tensor_geo, e_piezo_tensor_geo, &
                                   e_piezo_tensor_fi_geo,                   &
                                   d_piezo_tensor_geo, polar0_geo
  USE control_epsilon_infty,        ONLY : lepsilon_infty_geo, lzeu_geo
  USE temperature,   ONLY : ntemp
  USE control_eldos, ONLY : lel_free_energy
  USE gvecw,          ONLY : ecutwfc
  USE gvect,          ONLY : ecutrho
  USE control_quadratic_energy, ONLY : nvar
  USE ions_base,      ONLY : nat
  USE lattices,       ONLY : crystal_parameters
  USE start_k,        ONLY : nk1, nk2, nk3
  USE klist,          ONLY : degauss, lgauss, ltetra
  USE clib_wrappers,       ONLY : f_mkdir_safe
  USE io_global,      ONLY : meta_ionode
  !
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: nwork
  INTEGER, INTENT(IN)  :: part

  INTEGER  :: igeom, igeom_qha, ike, iden, icount, ink, isigma, iwork, &
              iwork_tot, start_geometry1, last_geometry1, istep, nstep, ios
  INTEGER  :: count_energies, iq, irr, start_omega, i, j
  REAL(DP) :: compute_omega_geo, dual, kev, kedenv
!
!  here initialize the work to do and sets to true the flags that activate
!  the different parts of thermo_pw. Sets also the number of tasks for each
!  what.
!
  nwork=0
  IF (part == 1) THEN
!
!   the restart directory is used in all cases
!
     iwho=0
     ios=0
     IF (meta_ionode) ios = f_mkdir_safe( 'restart' )
     ngeom=1
     tot_ngeo=1
     start_geometry=1
     last_geometry=1
     start_geometry1=1
     last_geometry1=1
     SELECT CASE (TRIM(what))
!
!   In these cases we do not do any asynchronous work in the first part
!
        CASE ( 'plot_bz', ' ') 
        CASE ( 'scf') 
           lpwscf_syn_1=.TRUE.
        CASE ('scf_bands') 
           lpwscf_syn_1=.TRUE.
           lbands_syn_1=.TRUE.
           IF (meta_ionode) ios = f_mkdir_safe( 'band_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
        CASE ('scf_2d_bands')
           lpwscf_syn_1=.TRUE.
           lbands_syn_1=.TRUE.
           IF (meta_ionode) ios = f_mkdir_safe( 'band_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
        CASE ('scf_dos') 
           lpwscf_syn_1=.TRUE.
           lbands_syn_1=.TRUE.
           ldos_syn_1=.TRUE.
           IF (meta_ionode) ios = f_mkdir_safe( 'band_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'therm_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
           CALL allocate_el_thermodynamics(1)
        CASE ('scf_ph') 
           lpwscf_syn_1=.NOT.after_disp
           lph=.TRUE.
           tot_ngeo=1
           tot_ngeo_eos=tot_ngeo
           ALLOCATE(no_ph(tot_ngeo))
           ALLOCATE(no_ph_eos(tot_ngeo))
           ALLOCATE(epsilon_infty_geo(3,3,tot_ngeo))
           ALLOCATE(zeu_geo(3,3,nat,tot_ngeo))
           ALLOCATE(lepsilon_infty_geo(tot_ngeo))
           ALLOCATE(lzeu_geo(tot_ngeo))
           ALLOCATE(omega_geo(tot_ngeo))
           omega_geo(1)=omega_save
           lepsilon_infty_geo=.FALSE.
           lzeu_geo=.FALSE.
           ALLOCATE(epsilon_zero_geo(3,3,tot_ngeo))
           ALLOCATE(epsilon_zerom1_geo(3,3,tot_ngeo))
           ALLOCATE(freq_geo(3*nat,tot_ngeo))
           ALLOCATE(z_geo(3*nat,3*nat,tot_ngeo))
           no_ph(1)=.FALSE.
           no_ph_eos(1)=.FALSE.
           IF (meta_ionode) ios = f_mkdir_safe( 'elastic_constants' )
           IF (meta_ionode) ios = f_mkdir_safe( 'dynamical_matrices' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
        CASE ('scf_disp')
           lpwscf_syn_1=.NOT.after_disp
           lph=.TRUE.
           tot_ngeo=1
           tot_ngeo_eos=tot_ngeo
           ALLOCATE(no_ph(tot_ngeo))
           ALLOCATE(no_ph_eos(tot_ngeo))
           ALLOCATE(epsilon_infty_geo(3,3,tot_ngeo))
           ALLOCATE(zeu_geo(3,3,nat,tot_ngeo))
           ALLOCATE(lepsilon_infty_geo(tot_ngeo))
           ALLOCATE(lzeu_geo(tot_ngeo))
           lepsilon_infty_geo=.FALSE.
           lzeu_geo=.FALSE.
           ALLOCATE(epsilon_zero_geo(3,3,tot_ngeo))
           ALLOCATE(epsilon_zerom1_geo(3,3,tot_ngeo))
           ALLOCATE(freq_geo(3*nat,tot_ngeo))
           ALLOCATE(z_geo(3*nat,3*nat,tot_ngeo))
           ALLOCATE(omega_geo(tot_ngeo))
           omega_geo(1)=omega_save
           no_ph(1)=.FALSE.
           no_ph_eos(1)=.FALSE.
           lq2r = .TRUE.
           ltherm = ltherm_dos .OR. ltherm_freq
           CALL allocate_thermodynamics()
           IF (meta_ionode) ios = f_mkdir_safe( 'elastic_constants' )
           IF (meta_ionode) ios = f_mkdir_safe( 'dynamical_matrices' )
           IF (meta_ionode) ios = f_mkdir_safe( 'phdisp_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'therm_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
!
!   In these cases we make asynchronous work in the first part
!
        CASE ( 'scf_ke') 
           nwork= count_energies(ecutwfc, ecutrho, deltake, deltakeden, nke,&
                                                                      nkeden)
           ncutoffene=nwork
           ALLOCATE(ke(nwork))
           ALLOCATE(keden(nwork))
           ALLOCATE(energy_geo(nwork))
           ALLOCATE(ef_geo(nwork))
           start_geometry1=MAX(1,start_geometry_save)
           last_geometry1=MIN(nwork, last_geometry_save)
           icount=0
           DO iden=1, nkeden
              kedenv = ecutrho + (iden-1) * deltakeden
              DO ike = 1, nke
                 kev = ecutwfc + (ike-1) * deltake
                 dual=kedenv/kev
                 IF ( dual > 3.9999_dp ) THEN
                    icount = icount + 1
                    ke(icount) = kev
                    keden(icount) = kedenv
                 ENDIF
              ENDDO
           ENDDO
           energy_geo=0.0_DP
           lconv_ke_test=.TRUE.
           do_punch=.FALSE.
           IF (meta_ionode) ios = f_mkdir_safe( 'energy_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
        CASE ( 'scf_nk' ) 
           nwork= nnk * nsigma
           ALLOCATE(nk_test(3,nwork))
           ALLOCATE(sigma_test(nwork))
           ALLOCATE(energy_geo(nwork))
           ALLOCATE(ef_geo(nwork))
           start_geometry1=MAX(1,start_geometry_save)
           last_geometry1=MIN(nwork, last_geometry_save)
           tot_ngeo=0
           icount=0
           DO isigma=1, nsigma
              DO ink = 1, nnk
                 icount = icount + 1
                 nk_test(1,icount) = nk1 + (ink - 1) * deltank(1)
                 nk_test(2,icount) = nk2 + (ink - 1) * deltank(2)
                 nk_test(3,icount) = nk3 + (ink - 1) * deltank(3)
                 sigma_test(icount) = degauss + (isigma - 1) * deltasigma
              ENDDO
           ENDDO
           energy_geo=0.0_DP
           lconv_nk_test=.TRUE.
           do_punch=.FALSE.
           IF (meta_ionode) ios = f_mkdir_safe( 'energy_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
!
!   in part 1 these cases do nothing
!
        CASE ('scf_elastic_constants') 
           lpart2_pw=.TRUE.
           tot_ngeo=1
           tot_ngeo_eos=tot_ngeo
           start_geometry_qha=1
           last_geometry_qha=1
           IF (meta_ionode) ios = f_mkdir_safe( 'therm_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'elastic_constants' )
        CASE ('scf_piezoelectric_tensor')
           lpart2_pw=.TRUE.
           tot_ngeo=1
           tot_ngeo_eos=tot_ngeo
           IF (meta_ionode) ios = f_mkdir_safe( 'elastic_constants' )
        CASE ('scf_polarization') 
           lpart2_pw=.TRUE.
           tot_ngeo=1
!
!   here all the cases that require the determination of the minimization
!   of the energy to find the equilibrium crystal parameters
!
        CASE ('mur_lc')
           lpart2_pw=lel_free_energy
           do_punch=lel_free_energy
           lev_syn_1=.TRUE.
           lpwscf_syn_1=do_scf_relax
           CALL initialize_mur(nwork)
           ALLOCATE(no_ph(nwork))
           ALLOCATE(no_ph_eos(nwork))
           no_ph=.TRUE.
           no_ph_eos=.TRUE.
           CALL summarize_geometries(nwork)
           tot_ngeo=nwork
           start_geometry=MAX(1,start_geometry_save)
           last_geometry=MIN(nwork, last_geometry_save)
           last_geometry1=nwork
           IF (meta_ionode) ios = f_mkdir_safe( 'energy_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
           IF (lel_free_energy) THEN
              IF (meta_ionode) ios = f_mkdir_safe( 'therm_files' )
              IF (meta_ionode) ios = f_mkdir_safe( 'anhar_files' )
              CALL allocate_el_thermodynamics(tot_ngeo)
           ENDIF
        CASE ('mur_lc_bands') 
           do_punch=.FALSE.
           lev_syn_1=.TRUE.
           lpwscf_syn_1=.TRUE.
           lbands_syn_1=.TRUE.
           CALL initialize_mur(nwork)
           last_geometry1=nwork
           IF (meta_ionode) ios = f_mkdir_safe( 'energy_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'band_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
        CASE ('mur_lc_dos') 
           do_punch = .FALSE.
           lev_syn_1=.TRUE.
           lpwscf_syn_1 = .TRUE.
           lbands_syn_1 = .TRUE.
           ldos_syn_1 = .TRUE.
           CALL initialize_mur(nwork)
           last_geometry1=nwork
           IF (meta_ionode) ios = f_mkdir_safe( 'energy_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'therm_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
           CALL allocate_el_thermodynamics(1)
        CASE ('mur_lc_ph') 
           do_punch=.FALSE.
           lev_syn_1=.TRUE.
           lpwscf_syn_1=.NOT.after_disp
           CALL initialize_mur(nwork)
           last_geometry1=nwork
           start_geometry_qha=1
           last_geometry_qha=1
           lph=.TRUE.
           tot_ngeo=1
           tot_ngeo_eos=tot_ngeo
           ALLOCATE(no_ph(tot_ngeo))
           ALLOCATE(no_ph_eos(tot_ngeo))
           no_ph(1)=.FALSE.
           no_ph_eos(1)=.FALSE.
           IF (meta_ionode) ios = f_mkdir_safe( 'energy_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'dynamical_matrices' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
        CASE ('mur_lc_disp')
           do_punch=.FALSE.
           lev_syn_1=.TRUE.
           lpwscf_syn_1=.NOT.after_disp
           CALL initialize_mur(nwork)
           last_geometry1=nwork
           start_geometry_qha=1
           last_geometry_qha=1
           lph=.TRUE.
           tot_ngeo=1
           tot_ngeo_eos=tot_ngeo
           ALLOCATE(no_ph(tot_ngeo))
           ALLOCATE(no_ph_eos(tot_ngeo))
           no_ph(1)=.FALSE.
           no_ph_eos(1)=.FALSE.
           lq2r = .TRUE.
           ltherm = ltherm_dos .OR. ltherm_freq
           CALL allocate_thermodynamics()
           IF (meta_ionode) ios = f_mkdir_safe( 'energy_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'dynamical_matrices' )
           IF (meta_ionode) ios = f_mkdir_safe( 'phdisp_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'therm_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
        CASE ('mur_lc_elastic_constants') 
           do_punch=.FALSE.
           lev_syn_1=.TRUE.
           lpwscf_syn_1=do_scf_relax
           CALL initialize_mur(nwork)
           last_geometry1=nwork
           start_geometry_qha=1
           last_geometry_qha=1
           lpart2_pw=.TRUE.
           tot_ngeo=nwork
           tot_ngeo_eos=tot_ngeo
           IF (meta_ionode) ios = f_mkdir_safe( 'therm_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'energy_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'elastic_constants' )
        CASE ('mur_lc_piezoelectric_tensor') 
           do_punch=.FALSE.
           lev_syn_1=.TRUE.
           lpwscf_syn_1=do_scf_relax
           CALL initialize_mur(nwork)
           last_geometry1=nwork
           start_geometry=MAX(1, start_geometry_save)
           last_geometry=MIN(nwork, last_geometry_save)
           lpart2_pw=.TRUE.
           tot_ngeo=1
           tot_ngeo_eos=tot_ngeo
           IF (meta_ionode) ios = f_mkdir_safe( 'energy_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'elastic_constants' )
        CASE ('mur_lc_polarization')
           do_punch=.FALSE.
           lev_syn_1=.TRUE.
           lpwscf_syn_1=do_scf_relax
           CALL initialize_mur(nwork)
           last_geometry1=nwork
           lpart2_pw=.TRUE.
           tot_ngeo=1
           IF (meta_ionode) ios = f_mkdir_safe( 'energy_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
!
!    Here all the cases that compute the free energy and minimize it
!
        CASE ('mur_lc_t')
           lev_syn_1=.TRUE.
           CALL initialize_mur(nwork)
           lph = .TRUE.
           lq2r = .TRUE.
           ltherm = .TRUE.
           lev_syn_2=.TRUE.
           tot_ngeo=nwork
           tot_ngeo_eos=tot_ngeo
           start_geometry=MAX(1,start_geometry_save)
           last_geometry=MIN(nwork, last_geometry_save)
           start_geometry1=start_geometry
           last_geometry1=last_geometry
           start_geometry_qha=1
           last_geometry_qha=1
           ALLOCATE(no_ph(tot_ngeo))
           ALLOCATE(no_ph_eos(tot_ngeo))
           ALLOCATE(found_dos_ec(tot_ngeo))
           ALLOCATE(found_ph_ec(tot_ngeo))
           ALLOCATE(all_geometry_done_geo(tot_ngeo))
           CALL initialize_no_ph(no_ph, no_ph_eos, tot_ngeo, ibrav_save)
           CALL summarize_geometries(nwork)
           nvar=crystal_parameters(ibrav_save)
           CALL allocate_thermodynamics()
           CALL allocate_anharmonic()
           IF (ANY(stype)) ALLOCATE(dyde(max_nint_var,21,tot_ngeo,ntemp))
           IF (meta_ionode) ios = f_mkdir_safe( 'energy_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'anhar_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'therm_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'dynamical_matrices' )
           IF (meta_ionode) ios = f_mkdir_safe( 'phdisp_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'elastic_constants' )
        CASE ('elastic_constants_geo')
           lectqha=.TRUE.
           lph=use_free_energy
           IF (.NOT.lph) lpart2_pw=lel_free_energy
           lq2r = use_free_energy
           ltherm = use_free_energy
           do_punch=(use_free_energy.OR.lel_free_energy)
           CALL initialize_mur_qha(ngeom)
           IF (start_geometry_qha<1) start_geometry_qha=1
           IF (last_geometry_qha>ngeom) last_geometry_qha=ngeom
           CALL initialize_elastic_cons(ngeom, nwork)
           !
           !  Check that start_geometry and last geometry are within the
           !  range of computed geometries otherwise disregard them
           !
           start_geometry=start_geometry_save
           last_geometry=last_geometry_save
           IF ((start_geometry < ((start_geometry_qha-1)*work_base+1)).OR. &
               (start_geometry > last_geometry_qha*work_base))             &
               start_geometry = (start_geometry_qha-1)*work_base+1
           IF ((last_geometry < start_geometry) .OR.                       &
               (last_geometry > last_geometry_qha*work_base))              &
               last_geometry = last_geometry_qha*work_base
           start_geometry=MAX((start_geometry_qha-1)*work_base+1, &
                                                       start_geometry)
           last_geometry=MIN(last_geometry_qha*work_base, last_geometry)
           start_geometry1=start_geometry
           last_geometry1=last_geometry
           tot_ngeo=nwork
           tot_ngeo_eos=tot_ngeo
           ALLOCATE(energy_geo(tot_ngeo))
           ALLOCATE(tau_save_ec(3,nat,tot_ngeo))
           ALLOCATE(all_geometry_done_geo(ngeom))
           IF (.NOT.lph) ALLOCATE(el_con_geo(6,6,ngeom))
           ALLOCATE(no_ph(tot_ngeo))
           ALLOCATE(no_ph_eos(tot_ngeo))
           IF (lel_free_energy) ALLOCATE(ef_geo(tot_ngeo))
           energy_geo=0.0_DP
           no_ph=.FALSE.
           no_ph_eos=.FALSE.
           IF (ltau_from_file) THEN
              DO iwork=1,nwork
                 CALL check_geometry_exist(iwork,1,iwho)
              ENDDO
           ENDIF
           IF (elastic_algorithm=='advanced'.OR.elastic_algorithm=='energy') &
                                                                          THEN
              ALLOCATE(omega_geo(nwork))
              DO iwork=1,nwork
                 omega_geo(iwork)=compute_omega_geo(ibrav_geo(iwork), &
                                  celldm_geo(:,iwork))
              ENDDO
              CALL summarize_geometries_ec(nwork)
           ENDIF
           IF (use_free_energy.OR.lel_free_energy) THEN
              CALL allocate_thermodynamics()
              CALL allocate_anharmonic()
              ALLOCATE(min_y_t(max_nint_var,ngeo_strain,21,ngeom,ntemp))
              ALLOCATE(min_yf_t(max_nint_var,ngeo_strain,21,ngeom,ntemp))
              ALLOCATE(dyde(max_nint_var,21,ngeom,ntemp))
           ELSE
              no_ph=.TRUE.
              no_ph_eos=.TRUE.
           ENDIF
           CALL allocate_debye()
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'elastic_constants' )
           IF (meta_ionode) ios = f_mkdir_safe( 'therm_files' )
           IF (use_free_energy) THEN
              IF (meta_ionode) ios = f_mkdir_safe( 'dynamical_matrices' )
              IF (meta_ionode) ios = f_mkdir_safe( 'phdisp_files' )
              IF (meta_ionode) ios = f_mkdir_safe( 'anhar_files' )
           ENDIF
        CASE ('piezoelectric_tensor_geo')
           do_punch=.TRUE.
           lpiezoelectric_tensor=.TRUE.
           lpiezotqha=use_free_energy
           lph=use_free_energy
           lq2r = use_free_energy
           ltherm = use_free_energy
           CALL initialize_mur_qha(ngeom)
           IF (start_geometry_qha<1) start_geometry_qha=1
           IF (last_geometry_qha>ngeom) last_geometry_qha=ngeom
           CALL set_piezo_tensor_work(ngeom, nwork)
           !
           !  Check that start_geometry and last geometry are within the
           !  range of computed geometries otherwise disregard them
           !
           start_geometry=start_geometry_save
           last_geometry=last_geometry_save
           IF ((start_geometry < ((start_geometry_qha-1)*work_base+1)).OR. &
               (start_geometry > last_geometry_qha*work_base))             &
               start_geometry = (start_geometry_qha-1)*work_base+1
           IF ((last_geometry < start_geometry) .OR.                       &
               (last_geometry > last_geometry_qha*work_base))              &
               last_geometry = last_geometry_qha*work_base
           start_geometry=MAX((start_geometry_qha-1)*work_base+1, &
                                                       start_geometry)
           last_geometry=MIN(last_geometry_qha*work_base, last_geometry)
           tot_ngeo=nwork
           tot_ngeo_eos=tot_ngeo
           ALLOCATE(all_geometry_done_geo(ngeom))
           ALLOCATE(no_ph(tot_ngeo))
           ALLOCATE(no_ph_eos(tot_ngeo))
           no_ph=.FALSE.
           no_ph_eos=.FALSE.
           ALLOCATE(polar0_geo(3,ngeom))
           ALLOCATE(g_piezo_tensor_geo(3,6,ngeom))
           ALLOCATE(eg_piezo_tensor_geo(3,6,ngeom))
           ALLOCATE(e_piezo_tensor_geo(3,6,ngeom))
           ALLOCATE(d_piezo_tensor_geo(3,6,ngeom))
           ALLOCATE(energy_geo(nwork))
           IF (use_free_energy) THEN
              CALL allocate_thermodynamics()
              CALL allocate_anharmonic()
              ALLOCATE(min_y_t(max_nint_var,ngeo_strain,21,ngeom,ntemp))
              ALLOCATE(min_yf_t(max_nint_var,ngeo_strain,21,ngeom,ntemp))
              ALLOCATE(dyde(max_nint_var,21,ngeom,ntemp))
           ELSE
              no_ph=.TRUE.
              no_ph_eos=.TRUE.
           ENDIF
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'elastic_constants' )
           IF (meta_ionode) ios = f_mkdir_safe( 'therm_files' )
           IF (use_free_energy) THEN
              IF (meta_ionode) ios = f_mkdir_safe( 'dynamical_matrices' )
              IF (meta_ionode) ios = f_mkdir_safe( 'phdisp_files' )
              IF (meta_ionode) ios = f_mkdir_safe( 'therm_files' )
              IF (meta_ionode) ios = f_mkdir_safe( 'anhar_files' )
           ENDIF
        CASE ('polarization_geo') 
           lpart2_pw=.TRUE.
           do_punch=.TRUE.
           tot_ngeo=1
           IF (meta_ionode) ios = f_mkdir_safe( 'elastic_constants' )
        CASE DEFAULT
           CALL errore('initialize_thermo_work','what not recognized',1)
     END SELECT
!
!  part 2
!
  ELSEIF (part == 2 ) THEN
!
!  In this part we do the phonon calculations
!
     SELECT CASE (TRIM(what))
        CASE ('scf_ph',        &
              'scf_disp',      &
              'mur_lc_ph',     &
              'mur_lc_disp')   
           CALL initialize_ph_work(nwork)
        CASE ('mur_lc_t')
           CALL initialize_ph_work(nwork)
        CASE ('elastic_constants_geo')
           IF (use_free_energy) THEN
              CALL initialize_ph_work(nwork)
           ELSEIF (lel_free_energy) THEN
              nwork=tot_ngeo
           ENDIF
        CASE ('piezoelectric_tensor_geo')
           IF (use_free_energy) THEN
              CALL initialize_ph_work(nwork)
           ENDIF
        CASE ('mur_lc')
             nwork=tot_ngeo
        CASE ('scf_elastic_constants', 'mur_lc_elastic_constants')

           IF (ALLOCATED(ibrav_geo))  DEALLOCATE(ibrav_geo)
           IF (ALLOCATED(celldm_geo)) DEALLOCATE(celldm_geo)
           IF (ALLOCATED(at_geo))     DEALLOCATE(at_geo)
           IF (ALLOCATED(tau_geo))    DEALLOCATE(tau_geo)
           IF (ALLOCATED(uint_geo))   DEALLOCATE(uint_geo)
           IF (ALLOCATED(energy_geo)) DEALLOCATE(energy_geo)
           IF (ALLOCATED(omega_geo))  DEALLOCATE(omega_geo)
!
!     initialise_elastic_cons allocates ibrav_geo and celldm_geo
!
           iwho=2
           CALL initialize_elastic_cons(ngeom,nwork)
           CALL set_unperturbed_geometry(ngeom)
           ALLOCATE(energy_geo(nwork))
           start_geometry=MAX(1, start_geometry_save)
           last_geometry=MIN(nwork, last_geometry_save)
           IF (ltau_from_file) THEN
              DO iwork=1,nwork
                  CALL check_geometry_exist(iwork,part,iwho)
              ENDDO
           ENDIF
           tot_ngeo=nwork
           energy_geo=0.0_DP
           IF (elalgen) THEN
              omega0= compute_omega_geo(ibrav_save, celldm0)
              el_con_omega_geo(1)=omega0
           ENDIF
           lelastic_const=.TRUE.
           do_punch=.FALSE.
           CALL allocate_debye()
        CASE ('scf_piezoelectric_tensor', 'mur_lc_piezoelectric_tensor')
           IF (ALLOCATED(energy_geo)) DEALLOCATE(energy_geo)
           IF (ALLOCATED(ibrav_geo)) DEALLOCATE(ibrav_geo)
           IF (ALLOCATED(celldm_geo)) DEALLOCATE(celldm_geo)
           IF (ALLOCATED(at_geo)) DEALLOCATE(at_geo)
           IF (ALLOCATED(tau_geo)) DEALLOCATE(tau_geo)
           IF (ALLOCATED(omega_geo)) DEALLOCATE(omega_geo)
           IF (ALLOCATED(uint_geo)) DEALLOCATE(uint_geo)
           iwho=2
           CALL set_unperturbed_geometry(ngeom)
           start_geometry_qha=1
           last_geometry_qha=1
           CALL set_piezo_tensor_work(1,nwork) 
           start_geometry=MAX(1, start_geometry_save)
           last_geometry=MIN(nwork, last_geometry_save)
           ALLOCATE(energy_geo(nwork))
           ALLOCATE(omega_geo(nwork))
           ALLOCATE(polar0_geo(3,ngeom))
           ALLOCATE(g_piezo_tensor_geo(3,6,ngeom))
           ALLOCATE(eg_piezo_tensor_geo(3,6,ngeom))
           ALLOCATE(e_piezo_tensor_geo(3,6,ngeom))
           ALLOCATE(d_piezo_tensor_geo(3,6,ngeom))
           energy_geo=0.0_DP
           lpiezoelectric_tensor=.TRUE.
           do_punch=.TRUE.
        CASE ('scf_polarization', 'mur_lc_polarization')
           IF (ALLOCATED(energy_geo)) DEALLOCATE(energy_geo)
           IF (ALLOCATED(omega_geo)) DEALLOCATE(omega_geo)
           nwork=1
           lpolarization=.TRUE.
           CALL allocate_piezo(nwork)
           ALLOCATE(energy_geo(nwork))
           ALLOCATE(omega_geo(nwork))
           do_punch=.TRUE.
        CASE ('polarization_geo')
           CALL initialize_mur(nwork)
           ALLOCATE(no_ph(nwork))
           ALLOCATE(no_ph_eos(nwork))
           no_ph=.TRUE.
           no_ph_eos=.TRUE.
           CALL summarize_geometries(nwork)
           start_geometry=MAX(1, start_geometry_save)
           last_geometry=MIN(nwork, last_geometry_save)
           CALL allocate_piezo(nwork)
           lpolarization=.TRUE.
     END SELECT
  ELSE
     CALL errore('initialize_thermo_work','unknown part',1)
  ENDIF

  IF ( nwork == 0 ) RETURN

  ALLOCATE( lpwscf(nwork) )
  ALLOCATE( lpwband(nwork) )
  ALLOCATE( lef(nwork) )
  ALLOCATE( lstress(nwork) )
  ALLOCATE( lberry(nwork) )
  ALLOCATE( lphonon(nwork) )
  lpwscf  = .FALSE.
  lpwband = .FALSE.
  lef = .FALSE.
  lstress = .FALSE.
  lberry  = .FALSE.
  lphonon = .FALSE.

  IF (part == 1) THEN
     SELECT CASE (TRIM(what))
!
!  The cases in which there is pwscf asynchronous work in the first part
!
        CASE ('scf_ke',                      &
              'scf_nk',                      &
              'mur_lc',                      &
              'mur_lc_bands',                &
              'mur_lc_dos',                  &
              'mur_lc_ph',                   &
              'mur_lc_disp',                 &
              'mur_lc_t',                    &
              'mur_lc_elastic_constants',    &
              'mur_lc_piezoelectric_tensor', &
              'mur_lc_polarization' )
           DO iwork=start_geometry1,last_geometry1
              lpwscf(iwork)=.TRUE.
              lef(iwork)=(lgauss.OR.ltetra)
           ENDDO
        CASE ('elastic_constants_geo')  
           DO igeom_qha=start_geometry_qha, last_geometry_qha
              DO iwork=1,work_base
                 iwork_tot= (igeom_qha-1)*work_base + iwork
                 IF (iwork_tot<start_geometry1.OR.iwork_tot>last_geometry1) CYCLE
                 lpwscf(iwork_tot)=.TRUE.
              ENDDO
           ENDDO
           lstress(1:nwork)=.NOT.elalgen
           IF (lel_free_energy) lef(1:nwork)=(lgauss.OR.ltetra)
        CASE ('piezoelectric_tensor_geo')
           DO iwork=start_geometry,last_geometry
              lpwscf(iwork)=.TRUE.
              lberry(iwork)=.TRUE.
           ENDDO
        CASE ('polarization_geo')
        CASE DEFAULT
          CALL errore('initialize_thermo_work','unknown what',1)
     END SELECT
  ENDIF 

  IF ( part == 2 ) THEN
     SELECT CASE (TRIM(what))
!
!   Here the cases in which there are phonon calculations in the second part
! 
        CASE ('scf_ph',           &
              'scf_disp',         &
              'mur_lc_ph',        &
              'mur_lc_disp',      &
              'mur_lc_t')         
             
           CALL initialize_flags_for_ph(nwork)

!
!  Here the case in which lel_free_energy is .true. In this case we must
!  do a el_dos calculation
!
        CASE( 'mur_lc' )
           DO iwork=start_geometry, last_geometry
              lpwband(iwork)=.TRUE.
           ENDDO
!
!  Here the case of elastic_constants_geo. In this case we must
!  do a el_dos calculation if use_free_energy is .FALSE. and
!  lel_free_energy is .TRUE. or a phonon calculation is 
!  use_free_energy is  .TRUE.
!
        CASE('elastic_constants_geo')
           IF (use_free_energy) THEN
              CALL initialize_flags_for_ph(nwork)
           ELSE
              DO igeom_qha=start_geometry_qha, last_geometry_qha
                 DO iwork=1,work_base
                    iwork_tot= (igeom_qha-1)*work_base + iwork
                    IF (iwork_tot.GE.start_geometry.AND.iwork_tot &
                              .LE.last_geometry) lpwband(iwork_tot)=.TRUE.
                 ENDDO
              ENDDO
           ENDIF
!
!  Here the cases in which there are pwscf calculation in the second part
!
        CASE ('scf_elastic_constants', 'mur_lc_elastic_constants')
           DO iwork=start_geometry,last_geometry
              lpwscf(iwork)=.TRUE.
              IF (.NOT.elalgen) lstress(iwork)=.TRUE.
           ENDDO
!
!  Here the cases in which there are pwscf and berry phase calculation 
!  in the second part
!
        CASE ('scf_piezoelectric_tensor', 'mur_lc_piezoelectric_tensor', &
              'scf_polarization', 'mur_lc_polarization','polarization_geo')
           DO iwork=start_geometry,last_geometry
              lpwscf(iwork)=.TRUE.
              lberry(iwork)=.TRUE.
           ENDDO
        CASE('piezoelectric_tensor_geo') 
           CALL initialize_flags_for_ph(nwork)
     END SELECT
  ENDIF

  RETURN
END SUBROUTINE initialize_thermo_work

!-----------------------------------------------------------------------
SUBROUTINE set_celldm_geo(celldm_geo, nwork)
!-----------------------------------------------------------------------
!
!   This routine sets the grid of values on celldm_geo.
!
USE kinds,         ONLY : DP
USE control_thermo, ONLY : lgeo_from_file
USE constants,     ONLY : pi
USE thermo_mod,    ONLY : step_ngeo, ngeo
USE initial_conf,  ONLY : celldm_save
USE geometry_file, ONLY : set_celldm_geo_from_file
USE control_atomic_pos, ONLY : ninternal, linternal_thermo

IMPLICIT NONE
INTEGER, INTENT(IN) :: nwork
REAL(DP), INTENT(INOUT) :: celldm_geo(6,nwork)

INTEGER  :: igeo1, igeo2, igeo3, igeo4, igeo5, igeo6
INTEGER  :: iwork, i, total_work, iw, pos
REAL(DP) :: angle1, angle2, angle3, delta(6), celldm_geo_eff(6,nwork)

IF (lgeo_from_file) THEN
   CALL set_celldm_geo_from_file(celldm_geo, ngeo(1))
   RETURN
ENDIF

delta=0.0_DP
DO i=1,6
   IF (MOD(ngeo(i),2)==0) delta(i)=step_ngeo(i)/2.0_DP
ENDDO

iwork=0
celldm_geo=0.0_DP
total_work=0
DO igeo6 = -ngeo(6)/2, (ngeo(6)-1)/2
   angle3 = ACOS(celldm_save(6))+(igeo6*step_ngeo(6))*pi/180.0_DP
   DO igeo5 = -ngeo(5)/2, (ngeo(5)-1)/2
      angle2 = ACOS(celldm_save(5))+(igeo5*step_ngeo(5))*pi/180.0_DP
      DO igeo4 = -ngeo(4)/2, (ngeo(4)-1)/2
         angle1 = ACOS(celldm_save(4))+(igeo4*step_ngeo(4))*pi/180.0_DP
         DO igeo3 = -ngeo(3)/2, (ngeo(3)-1)/2
            DO igeo2 = -ngeo(2)/2, (ngeo(2)-1)/2
               DO igeo1 = -ngeo(1)/2, (ngeo(1)-1)/2
                  total_work=total_work+1
                  IF (total_work>nwork/ninternal) &
                     CALL errore('set_celldm_geo','nwork too low', 1)
                  iwork=iwork+1
                  celldm_geo(1,iwork)=celldm_save(1)+igeo1*step_ngeo(1)
                  celldm_geo(2,iwork)=celldm_save(2)+igeo2*step_ngeo(2)
                  celldm_geo(3,iwork)=celldm_save(3)+igeo3*step_ngeo(3)
                  celldm_geo(1,iwork)=celldm_geo(1,iwork)+delta(1)
                  celldm_geo(2,iwork)=celldm_geo(2,iwork)+delta(2)
                  celldm_geo(3,iwork)=celldm_geo(3,iwork)+delta(3)
                  angle1=angle1+delta(4)*pi/180.0_DP
                  angle2=angle2+delta(5)*pi/180.0_DP
                  angle3=angle3+delta(6)*pi/180.0_DP
                  celldm_geo(4,iwork)=COS(angle1)
                  celldm_geo(5,iwork)=COS(angle2)
                  celldm_geo(6,iwork)=COS(angle3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

!
!  if linternal_thermo is true the number of works is multiplied by
!  the number of atomic configurations studied for each set of 
!  crystal parameters.
!
IF (linternal_thermo) THEN
   IF (iwork*ninternal/=nwork) &
      CALL errore('set_celldm_geo','problem with nwork',1)
   pos=0
   DO iw=1,iwork
      DO i=1,ninternal
         pos=pos+1
         celldm_geo_eff(:,pos)=celldm_geo(:,iw)
      ENDDO
   ENDDO
   celldm_geo(:,1:nwork)=celldm_geo_eff(:,1:nwork)
ENDIF


RETURN
END SUBROUTINE set_celldm_geo
!
!-----------------------------------------------------------------------
SUBROUTINE set_tau_geo(ibrav_geo, celldm_geo, tau_geo, at_geo, uint_geo, &
                                 nwork, nat, nint_var)
!-----------------------------------------------------------------------
!
!   This routine sets the grid of values on uint_geo.
!   Data are arranged so that the ninternal displacements for the
!   first set of external parameters is written first. Then there 
!   are the ninternal displacements for the second set of external
!   parameters, etc.. ninternal parameters first vary the first parameter,
!   then the second etc.
!
USE kinds,         ONLY : DP
USE control_atomic_pos, ONLY : ninternal, int_ngeo, int_step_ngeo, &
                       iconstr_internal, uint0, uint_eq, tau_eq
USE cell_base,     ONLY : celldm, at, ibrav
USE ions_base,     ONLY : tau

IMPLICIT NONE
INTEGER, INTENT(IN) :: nwork, nat, nint_var, ibrav_geo(nwork)
REAL(DP), INTENT(IN) :: celldm_geo(6,nwork), at_geo(3,3,nwork)
REAL(DP), INTENT(INOUT) :: tau_geo(3,nat,nwork)
REAL(DP), INTENT(INOUT) :: uint_geo(nint_var,nwork)

INTEGER  :: igeo(nint_var), id(nint_var)
INTEGER  :: iwork, i, total_work, pos, iinternal, nexternal, iexternal, ivar
REAL(DP) :: delta(nint_var)


ALLOCATE(uint_eq(nint_var))
ALLOCATE(tau_eq(3,nat))
!
!   First determine the value of uint0 of the geometry given as
!   input to thermo_pw
!
ALLOCATE(uint0(nint_var))

CALL internal_to_tau(ibrav, celldm, tau, at, uint0, nat, nint_var, &
                                           iconstr_internal, 2)
!
!  then generate a mesh of dimensions nint_var of values of u centered
!  in uint0
!
delta=0.0_DP
DO ivar=1,nint_var
   IF (MOD(int_ngeo(ivar),2)==0) delta(ivar)=int_step_ngeo(ivar)/2.0_DP
   id(ivar)= int_ngeo(ivar)/2 + 1 
ENDDO
!
!  compute the internal parameters for each work to do
!
uint_geo=0.0_DP
total_work=0
nexternal=nwork/ninternal
DO iexternal=1,nexternal
   DO iinternal = 1, ninternal
      CALL find_ipoint(iinternal, nint_var, int_ngeo, igeo)
      total_work=total_work+1
      DO ivar=1,nint_var
         uint_geo(ivar,total_work)=(igeo(ivar)-id(ivar))*int_step_ngeo(ivar) &
                                                 + delta(ivar) + uint0(ivar)
      ENDDO
   ENDDO
ENDDO
!
!  and determine for each value of u the atomic coodinates
!
DO iwork=1,nwork
   CALL internal_to_tau(ibrav_geo(iwork), celldm_geo(1,iwork), &
                        tau_geo(1,1,iwork), at_geo(1,1,iwork), &
                       uint_geo(1,iwork), nat, nint_var, iconstr_internal, 1)
ENDDO

RETURN
END SUBROUTINE set_tau_geo
!
!-----------------------------------------------------------------------
SUBROUTINE summarize_geometries(nwork)
!-----------------------------------------------------------------------
USE constants,     ONLY : bohr_radius_si
USE ions_base,     ONLY : nat
USE thermo_mod,    ONLY : ibrav_geo, omega_geo, celldm_geo, tau_geo, no_ph, &
                          uint_geo
USE initial_conf,  ONLY : celldm_save, ibrav_save, tau_save, atm_save, &
                          ityp_save
USE control_atomic_pos, ONLY : linternal_thermo, nint_var, uint0 
USE io_global,     ONLY : stdout

IMPLICIT NONE
INTEGER :: nwork
INTEGER :: igeo, phcount, i, ipol, na

phcount=0
DO igeo=1,nwork
   IF(.NOT.no_ph(igeo)) phcount=phcount+1
ENDDO

WRITE(stdout,'(/,5x,70("-"))')
WRITE(stdout,'(/,5x,"Number of geometries to compute:", i5)') nwork
WRITE(stdout,'(5x,"Phonons computed in:", 12x, i5, " geometries")') phcount
WRITE(stdout,'(/,5x,"Input ibrav and celldm:")')
WRITE(stdout,'(12x, i3, 6f10.5,/)') ibrav_save, celldm_save(:)
IF (linternal_thermo) THEN
   WRITE(stdout,'(/,5x,"Input atomic positions: uint=",2f15.8)') (uint0(i), &
                                                   i=1, nint_var)
   WRITE( stdout, '(/,3x,"Cartesian axes")')
   WRITE( stdout, '(/,5x,"site n.     atom                  positions (alat units)")')
   WRITE( stdout, '(6x,i4,8x,a6," tau(",i4,") = (",3f12.7,"  )")') &
        (na, atm_save(ityp_save(na)), na, (tau_save(ipol,na), ipol=1,3), &
                                                               na=1,nat)
   WRITE( stdout, *)
ENDIF

WRITE(stdout,'(/,5x,"List of the geometries to compute:")')

DO igeo=1,nwork
   WRITE(stdout,'(5x,i5,": ", i3,6f10.5,l2)') igeo, ibrav_geo(igeo), &
                                   celldm_geo(:,igeo), .NOT.no_ph(igeo)
   IF (linternal_thermo) THEN
      WRITE( stdout,'(5x,"Atomic positions: uint=",2f15.8)') &
                                            (uint_geo(i,igeo), i=1, nint_var)
      WRITE( stdout, '(6x,i4,8x,a6," tau(",i4,") = (",3f12.7,"  )")')       &
           (na, atm_save(ityp_save(na)), na, (tau_geo(ipol,na,igeo),        &
                                                   ipol=1,3),  na=1,nat)
      WRITE( stdout, *)
   ENDIF
ENDDO
WRITE(stdout,'(/,5x,"Volumes: ",10x,"(a.u.)^3",10x,"(A)^3")')
DO igeo=1,nwork
   WRITE(stdout,'(5x,i5,2f20.10)') igeo, omega_geo(igeo), &
                            omega_geo(igeo)*(bohr_radius_si)**3/1.D-30
ENDDO
WRITE(stdout,'(/,5x,70("-"))')

RETURN
END SUBROUTINE summarize_geometries
!
!-----------------------------------------------------------------------
SUBROUTINE summarize_geometries_ec(nwork)
!-----------------------------------------------------------------------
USE thermo_mod,    ONLY : ibrav_geo, omega_geo, celldm_geo, no_ph, uint_geo
USE initial_conf,  ONLY : celldm_save, ibrav_save, atm_save, ityp_save
USE control_elastic_constants, ONLY : ngeo_strain, nstep_ec, ngeom, tau_acc, &
                                      nint_var_ec, stypec, stype, nmove,     &
                                      ninternal_ec
USE io_global,     ONLY : stdout
USE ions_base,     ONLY : nat

IMPLICIT NONE
INTEGER :: nwork
INTEGER :: igeo, igeom, istep, iwork, phcount, ipol, na, ivar, imove, &
                                               iinternal
INTEGER :: compute_nwork_ph

phcount=compute_nwork_ph(no_ph, nwork)

WRITE(stdout,'(/,5x,70("-"))')
WRITE(stdout,'(/,5x,"Number of geometries to compute:", i5)') nwork
WRITE(stdout,'(5x,"Phonons computed in:", 12x, i5, " geometries")') phcount
WRITE(stdout,'(/,5x,"Input ibrav and celldm:")')
WRITE(stdout,'(12x, i3, 6f10.5,/)') ibrav_save, celldm_save(:)

WRITE(stdout,'(/,5x,"List of the geometries to compute:")')

IF (ANY(stypec)) THEN
   iwork=0
   DO igeom=1, ngeom
      DO istep=1,nstep_ec
         DO igeo=1,ngeo_strain
            IF (stypec(istep)) THEN
               DO iinternal=1,ninternal_ec(istep)
                  iwork=iwork+1
                  WRITE(stdout,'(5x,i5,": ", i3,6f10.5,l2)') iwork, &
                      ibrav_geo(iwork), celldm_geo(:,iwork), .NOT.no_ph(iwork)
                  WRITE( stdout,'(5x,"Atomic positions: uint=",5f15.8)') &
                         (uint_geo(ivar,iwork), ivar=1, nint_var_ec(istep))
                  WRITE( stdout, '(6x,i4,8x,a6," tau(",i4,") = &
                                                         &(",3f12.7,"  )")') &
                      (na, atm_save(ityp_save(na)), na, &
                                (tau_acc(ipol,na,iwork), ipol=1,3),  na=1,nat)
                  WRITE( stdout, *)
               ENDDO
            ELSE
               iwork=iwork+1
               WRITE(stdout,'(5x,i5,": ", i3,6f10.5,l2)') iwork, &
                      ibrav_geo(iwork), celldm_geo(:,iwork), .NOT.no_ph(iwork)
            ENDIF
         ENDDO
      ENDDO 
   ENDDO
ENDIF

IF (ANY(stype)) THEN
   iwork=0
   DO igeom=1, ngeom
      DO istep=1,nstep_ec
         IF (stype(istep)) THEN
            DO igeo=1,ngeo_strain
               DO imove=1, nmove
                  iwork=iwork+1
                  WRITE(stdout,'(5x,i5,": ", i3,6f10.5,l2)') iwork, &
                  ibrav_geo(iwork), celldm_geo(:,iwork), .NOT.no_ph(iwork)
                  WRITE( stdout,'(5x,"Atomic positions:")')          
                  WRITE( stdout, '(6x,i4,8x,a6," tau(",i4,") &
                     &= (",3f12.7,"  )")') (na, atm_save(ityp_save(na)), na, &
                        (tau_acc(ipol,na,iwork), ipol=1,3),  na=1,nat)
                  WRITE( stdout, *)
               ENDDO
            ENDDO
         ELSE
            DO igeo=1,ngeo_strain
               iwork=iwork+1
               WRITE(stdout,'(5x,i5,": ", i3,6f10.5,l2)') iwork, &
                  ibrav_geo(iwork), celldm_geo(:,iwork), .NOT.no_ph(iwork)
            ENDDO
         ENDIF
      ENDDO
   ENDDO
ENDIF
RETURN
END SUBROUTINE summarize_geometries_ec
!
!-----------------------------------------------------------------------
INTEGER FUNCTION compute_nwork()
!-----------------------------------------------------------------------
!
!  This function computes the number of tasks needed for energy minimization
!
USE thermo_mod,  ONLY : ngeo, central_geo
USE control_mur, ONLY : lmurn

IMPLICIT NONE

INTEGER :: i, auxgeo

auxgeo=ngeo(1)*ngeo(2)*ngeo(3)*ngeo(4)*ngeo(5)*ngeo(6)
IF (lmurn) auxgeo=ngeo(1)
central_geo=auxgeo/2
IF (MOD(auxgeo,2)==1) central_geo=central_geo+1

compute_nwork=auxgeo

RETURN
END FUNCTION compute_nwork
!
!-----------------------------------------------------
INTEGER FUNCTION compute_nwork_ph(no_ph,ndatatot)
!-----------------------------------------------------
!
!  This routines computes the number of geometries in which
!  phonon modes have to be calculated.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: ndatatot
LOGICAL, INTENT(IN) :: no_ph(ndatatot)

INTEGER :: idata, counter_ndata

counter_ndata=0
DO idata=1,ndatatot
   IF (.NOT. no_ph(idata)) counter_ndata=counter_ndata+1
ENDDO
compute_nwork_ph=counter_ndata

RETURN
END FUNCTION compute_nwork_ph

!-----------------------------------------------------------------------
SUBROUTINE initialize_mur(nwork)
!-----------------------------------------------------------------------
!
!   This routine sets the number of tasks needed to make an energy
!   minimization. It allocates the variables ibrav_geo, celldm_geo, 
!   energy_geo, ef_geo, and omega_geo. It allocates uint_geo and tau_geo
!   Moreover it sets the values of celldm_geo, in each task.
!   It allocates p0, omega_p, density_p, celldm_p, p2_p, p4_p 
!   for pressure dependent quantities.
!   It allocates p4_eq and p2_eq the polynomials that interpolate 
!   uint_geo and uint_p and tau_p the values of uint_geo as a function
!   of pressure.
!   If ntemp_plot>0 or npress_plot>0 it allocates space for the 
!   eos interpolation of the enthalpy at several pressures:
!   vmin_p, b0_p, b01_p, b02_p, and emin_p.
!
USE kinds,         ONLY : DP
USE thermo_mod,    ONLY : ibrav_geo, ngeo, celldm_geo, energy_geo, ef_geo, &
                          omega_geo, at_geo, tau_geo, tot_ngeo_eos, iwho,  &
                          uint_geo, epsilon_infty_geo, zeu_geo, z_geo,     &
                          freq_geo, epsilon_zero_geo, epsilon_zerom1_geo
USE control_atomic_pos, ONLY : linternal_thermo, nint_var, ninternal, &
                               int_ngeo, p2_eq, p4_eq, tau_p, uint_p
USE control_epsilon_infty, ONLY : lepsilon_infty_geo, lzeu_geo
USE ions_base,     ONLY : nat
USE control_vol,   ONLY : nvol, vmin_input, vmax_input, deltav
USE initial_conf,  ONLY : ibrav_save
USE temperature,   ONLY : ntemp_plot
USE geometry_file, ONLY : read_geometry_file
USE uniform_pressure, ONLY : omega_p, density_p, celldm_p, p2_p, p4_p
USE control_thermo, ONLY : lgeo_from_file, ltau_from_file
USE control_pressure, ONLY : npress, npress_plot
USE control_mur_p, ONLY : vmin_p, b0_p, b01_p, b02_p, emin_p
USE control_mur,   ONLY : p0

IMPLICIT NONE

INTEGER, INTENT(OUT) :: nwork

INTEGER              :: igeom, ivar, iwork
INTEGER              :: compute_nwork
REAL(DP)             :: compute_omega_geo

iwho=1
IF (lgeo_from_file) CALL read_geometry_file(ngeo)
nwork=compute_nwork()
tot_ngeo_eos=nwork
ninternal=1
IF (linternal_thermo) then
   DO ivar=1,nint_var
      ninternal=ninternal*int_ngeo(ivar)
   ENDDO
   nwork=nwork*ninternal
ENDIF
ALLOCATE(ibrav_geo(nwork))
ALLOCATE(celldm_geo(6,nwork))
ALLOCATE(at_geo(3,3,nwork))
ALLOCATE(energy_geo(nwork))
ALLOCATE(ef_geo(nwork))
ALLOCATE(tau_geo(3,nat,nwork))
ALLOCATE(omega_geo(nwork))
ALLOCATE(epsilon_infty_geo(3,3,nwork))
ALLOCATE(zeu_geo(3,3,nat,nwork))
ALLOCATE(lepsilon_infty_geo(nwork))
ALLOCATE(lzeu_geo(nwork))
ALLOCATE(epsilon_zero_geo(3,3,nwork))
ALLOCATE(epsilon_zerom1_geo(3,3,nwork))
ALLOCATE(freq_geo(3*nat,nwork))
ALLOCATE(z_geo(3*nat,3*nat,nwork))
ALLOCATE(uint_geo(nint_var,nwork))
lepsilon_infty_geo=.FALSE.
lzeu_geo=.FALSE.
CALL set_celldm_geo(celldm_geo, nwork)
IF (linternal_thermo) THEN
   CALL set_tau_geo(ibrav_geo, celldm_geo, tau_geo, at_geo, uint_geo, &
                    nwork, nat, nint_var)
ELSEIF (ltau_from_file) THEN
   DO iwork=1,nwork
      CALL check_geometry_exist(iwork,1,iwho)
   ENDDO
ENDIF

DO igeom = 1, nwork
   omega_geo(igeom)=compute_omega_geo(ibrav_save,celldm_geo(:,igeom))
ENDDO
energy_geo=0.0_DP
ibrav_geo=ibrav_save

IF (vmin_input == 0.0_DP) vmin_input=omega_geo(1) * 0.98_DP
IF (vmax_input == 0.0_DP) vmax_input=omega_geo(ngeo(1)*ninternal) * 1.02_DP
IF (nvol > 1) THEN
   deltav = (vmax_input - vmin_input)/(nvol-1)
ELSE
   IF (deltav > 0.0_DP) THEN
      nvol = NINT ( ( vmax_input - vmin_input ) / deltav + 1.51d0 )
   ELSE
      nvol = 51
      deltav = (vmax_input - vmin_input)/(nvol-1)
   ENDIF
ENDIF
ALLOCATE(p0(nvol))


ALLOCATE(omega_p(npress))
ALLOCATE(density_p(npress))
ALLOCATE(celldm_p(6,npress))
ALLOCATE(p2_p(npress))
ALLOCATE(p4_p(npress))

IF (linternal_thermo) THEN
   ALLOCATE(p4_eq(nint_var))
   ALLOCATE(p2_eq(nint_var))
   ALLOCATE(tau_p(3,nat,npress))
   ALLOCATE(uint_p(nint_var,npress))
ENDIF

IF (ntemp_plot>0.OR.npress_plot>0) THEN
   ALLOCATE(vmin_p(npress)) 
   ALLOCATE(b0_p(npress)) 
   ALLOCATE(b01_p(npress)) 
   ALLOCATE(b02_p(npress)) 
   ALLOCATE(emin_p(npress)) 
ENDIF

RETURN
END SUBROUTINE initialize_mur

!-----------------------------------------------------------------------
SUBROUTINE initialize_mur_qha(nwork)
!-----------------------------------------------------------------------
!
!   This routine sets the unperturbed geometries for the computation
!   of the elastic constants within the quasi-harmonic approximation.
!   It allocates the variables el_con_ibrav_geo, el_con_celldm_geo,
!   el_con_omega_geo, and el_con_tau_crys_geo.
!
USE kinds,         ONLY : DP
USE thermo_mod,    ONLY : ngeo, iwho
USE ions_base,     ONLY : nat
USE control_atomic_pos, ONLY : ninternal
USE control_elastic_constants, ONLY : el_con_ibrav_geo, el_con_celldm_geo, &
                                      el_con_tau_crys_geo, el_con_omega_geo, &
                                      el_con_tau_geo, el_con_at_geo
USE geometry_file, ONLY : read_geometry_file
USE control_thermo, ONLY : lgeo_from_file, ltau_el_cons_from_file
USE initial_conf,  ONLY : ibrav_save, tau_save_crys, atm_save, ityp_save
USE io_global, ONLY : stdout

IMPLICIT NONE

INTEGER, INTENT(OUT) :: nwork

INTEGER              :: igeom, iwork, ipol, na
INTEGER              :: compute_nwork
REAL(DP)             :: compute_omega_geo

iwho=2
ninternal=1
IF (lgeo_from_file) CALL read_geometry_file(ngeo)
nwork=compute_nwork()
ALLOCATE(el_con_ibrav_geo(nwork))
ALLOCATE(el_con_celldm_geo(6,nwork))
ALLOCATE(el_con_at_geo(3,3,nwork))
ALLOCATE(el_con_omega_geo(nwork))
ALLOCATE(el_con_tau_crys_geo(3,nat,nwork))
ALLOCATE(el_con_tau_geo(3,nat,nwork))
IF (ltau_el_cons_from_file) THEN
!
!  read from file geometry and atomic coordinates
!
   DO iwork=1,nwork
      CALL check_geometry_el_cons_exist(iwork,1)
   ENDDO
ELSE
   el_con_celldm_geo=0.0_DP
   CALL set_celldm_geo(el_con_celldm_geo(:,:), nwork)
   el_con_ibrav_geo(:)=ibrav_save
   DO iwork=1,nwork
      el_con_tau_crys_geo(:,:,iwork)=tau_save_crys(:,:)
!      el_con_omega_geo(iwork)=compute_omega_geo(ibrav_save, &
!                                    el_con_celldm_geo(:,iwork))
      CALL latgen(ibrav_save, el_con_celldm_geo(1,iwork), &
              el_con_at_geo(1,1,iwork), el_con_at_geo(1,2,iwork), &
              el_con_at_geo(1,3,iwork), el_con_omega_geo(iwork)) 
      el_con_at_geo(:,:,iwork)=el_con_at_geo(:,:,iwork)/   &
                                          el_con_celldm_geo(1,iwork)
   ENDDO
ENDIF

WRITE(stdout,'(/,5x,70("-"))')
WRITE(stdout,'(/,5x,"The celldm of the",i5," unperturbed geometries &
                                              &are:",/)') nwork
DO iwork=1,nwork
   WRITE(stdout,'(5x,i5,": ", 6f10.5)') iwork, &
                                   el_con_celldm_geo(1:6,iwork)
   IF (ltau_el_cons_from_file) THEN
      WRITE( stdout, '(6x,i4,8x,a6," tau(",i4,") = (",3f12.7,"  )")')      &
         (na, atm_save(ityp_save(na)), na, (el_con_tau_geo(ipol,na,iwork), &
                                              ipol=1,3), na=1,nat)
      WRITE( stdout, *)
   ENDIF
ENDDO
WRITE(stdout,'(/,5x,70("-"))')

RETURN
END SUBROUTINE initialize_mur_qha

!-----------------------------------------------------------------------
SUBROUTINE initialize_ph_work(nwork)
!-----------------------------------------------------------------------
!
!  This routine receives as input:
!  start_geometry
!  last_geometry
!  the first and last geometry computed in this run.
!  For each geometry it receives a structure collect_info_save(igeom).
!  This structure has the following information:
!  nqs    the number of q points to compute
!  irr_iq(iq) for each q point how many irreducible representations.  
!
!  This routine counts the number of tasks for one or many phonon 
!  For each phonon task it sets the geometry, the q point and the irr.
!  When use_ph_images=.TRUE. the q point and irr are not used 
!  and set to zero.
!
USE control_thermo, ONLY : all_geometries_together, geometry, iqw, irrw, &
                        comp_irr_iq_iw, comp_iq_iw, done_irr_iq_iw, done_iq_iw
USE thermo_mod,  ONLY : start_geometry, last_geometry
USE grid_irr_iq, ONLY : irr_iq
USE ions_base,   ONLY : nat
USE disp,        ONLY : nqs
USE freq_ph,     ONLY : nfs, fpol
USE images_omega,ONLY : omega_group
USE collect_info, ONLY : read_collect_info
USE initial_conf, ONLY : collect_info_save
USE control_ph,  ONLY : epsil, trans
USE control_qe,  ONLY : use_ph_images
USE mp_images,   ONLY : nimage

IMPLICIT NONE
INTEGER, INTENT(OUT) :: nwork

INTEGER :: iq, irr, igeom, iwork, nqs_max, ngeom, image, stge, lage, nqsx
INTEGER, ALLOCATABLE :: nqs_loc(:), irr_iq_loc(:,:)

nwork=0
IF (trans) THEN
!
!  find the number of geometries and the maximum number of q points
!  among geometries
!
   IF (all_geometries_together) THEN
      stge = start_geometry
      lage = last_geometry
      ngeom = lage - stge + 1
      nqsx = MAXVAL(collect_info_save(stge:lage)%nqs)
   ELSE
      stge = 1
      lage = 1
      ngeom = 1
      nqsx = nqs
   ENDIF
!
!   uniformize the case of many geometries and one geometry
!
   ALLOCATE(nqs_loc( stge : lage ))
   ALLOCATE(irr_iq_loc(nqsx, stge : lage))
   IF (all_geometries_together) THEN
      DO igeom = stge, lage
         nqs_loc(igeom)=collect_info_save(igeom)%nqs
         DO iq=1,nqs_loc(igeom)
            irr_iq_loc(iq,igeom)=collect_info_save(igeom)%irr_iq(iq)
         ENDDO
      ENDDO
   ELSE
      nqs_loc(1)=nqs
      DO iq=1,nqs
         irr_iq_loc(iq,1)=irr_iq(iq)
      ENDDO
   ENDIF
!
!  now count the number of works. 
!
   IF (use_ph_images) THEN
      nwork=nimage * ngeom
   ELSE
      nwork=0
      DO igeom = stge, lage
         DO iq = 1, nqs_loc(igeom)
            DO irr = 0, irr_iq_loc(iq,igeom)
               nwork = nwork + 1
            ENDDO
         ENDDO
      ENDDO
   ENDIF
!
!   then for each work set the geometry, iq, and irr.
!
   ALLOCATE(geometry(nwork))
   ALLOCATE(iqw(nwork))
   ALLOCATE(irrw(nwork))
   ALLOCATE(comp_irr_iq_iw(0:3*nat,nqsx,nwork))
   ALLOCATE(comp_iq_iw(nqsx,nwork))
   ALLOCATE(done_irr_iq_iw(0:3*nat,nqsx,nwork))
   ALLOCATE(done_iq_iw(nqsx,nwork))

   iqw=0
   irrw=0
   comp_irr_iq_iw=.FALSE.
   comp_iq_iw=.FALSE.
   done_irr_iq_iw=.FALSE.
   done_iq_iw=.FALSE.
   geometry=1
   IF (use_ph_images) THEN
      DO iwork=1, nwork
         geometry(iwork)= (iwork-1) / nimage + stge
         image=MOD(iwork-1,nimage)+1
         igeom=geometry(iwork)
         CALL read_collect_info(collect_info_save(igeom), nqs, &
              nat, image, comp_irr_iq_iw(0,1,iwork), comp_iq_iw(1,iwork), &
                          done_irr_iq_iw(0,1,iwork), done_iq_iw(1,iwork))
         iqw(iwork) = image
         irrw(iwork) = 0
      ENDDO
   ELSE
      nwork=0
      DO igeom = stge, lage
         DO iq=1, nqs_loc(igeom)
            DO irr=0, irr_iq_loc(iq,igeom)
               nwork = nwork + 1
               geometry(nwork) = igeom
               iqw(nwork) = iq
               irrw(nwork) = irr
               comp_irr_iq_iw(irr,iq,nwork)=.TRUE.
               comp_iq_iw(iq,nwork)=.TRUE.
               IF (collect_info_save(igeom)%done_irr_iq(irr,iq,1)==1) &
                                       done_irr_iq_iw(irr,iq,nwork)=.TRUE.
               IF (collect_info_save(igeom)%done_iq(iq,1)==1) &
                                       done_iq_iw(iq,nwork)=.TRUE.
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   DEALLOCATE(nqs_loc)
   DEALLOCATE(irr_iq_loc)
ELSEIF (fpol) THEN
   IF (nimage>1) THEN
      nwork=nfs/omega_group
      IF (nwork*omega_group /= nfs ) nwork=nwork+1
   ELSE
      nwork=1
   ENDIF
ELSEIF (epsil) THEN
   nwork=1
ELSE
   CALL errore('initialize_ph_work','Both trans and epsil are .FALSE.',1)
ENDIF

RETURN
END SUBROUTINE initialize_ph_work

!-----------------------------------------------------------------------
SUBROUTINE initialize_no_ph(no_ph, no_ph_eos, tot_ngeo, ibrav)
!-----------------------------------------------------------------------
!
!   This routine is used to skip phonon calculations in those geometries
!   that are not on the grid used to compute the vibrational energy.
!   When lgruneisen_gen is used, instead, phonons are computed only 
!   in a few geometries about a central geometry. According to 
!   ggrun_recipes these allow to compute a linear or a quadratic 
!   polynomial that fits the free energy.
!
USE thermo_mod, ONLY : ngeo, fact_ngeo, ngeo_ph
USE control_thermo, ONLY : lgruneisen_gen
USE control_gen_gruneisen, ONLY : icenter_grun, xngeo, ggrun_recipe, ind_rec3
USE lattices, ONLY : compress_int_vect, crystal_parameters
USE control_atomic_pos, ONLY : ninternal, linternal_thermo
IMPLICIT NONE

INTEGER, INTENT(IN)    :: tot_ngeo, ibrav
LOGICAL, INTENT(INOUT) :: no_ph(tot_ngeo), no_ph_eos(tot_ngeo)

LOGICAL :: todo(6)
INTEGER :: igeo1, igeo2, igeo3, igeo4, igeo5, igeo6, start, count_ngeo, &
           aux_deg(tot_ngeo), iwork, nvar, icount, iint_geom
INTEGER, ALLOCATABLE :: ipoint(:), inde(:)
LOGICAL :: select_ph_to_do
!
!   Here we set no_ph_eos taking a grid of
!   2*nvar+1 (ggrun_recipe=1)
!   3**nvar (ggrun_recipe=2) 
!   (nvar*(nvar+3))/2+1 (ggrun_recipe=3)
!   points about a geometry selected with icenter_grun. 
!   The derivatives of the free energy are then calculated numerically.
!   When ggrun_recipe=1 the free energy is modeled as a polynomial of 
!   first degree, when ggrun_recipe=2 or ggrun_recipe=3 the free energy
!   is modeled with a polynomial of second degree in the crystal parameters.
!
IF (lgruneisen_gen) THEN
   nvar=crystal_parameters(ibrav)
   ALLOCATE(ipoint(nvar))
   ALLOCATE(inde(nvar))
   CALL find_ipoint(icenter_grun, nvar, xngeo, ipoint )
   icount=0
   DO iwork=1,tot_ngeo
      CALL find_ipoint(iwork,nvar,xngeo,inde)
      no_ph_eos(iwork) = .NOT. select_ph_to_do(iwork,nvar,inde,ipoint, &
                                           ggrun_recipe)
      IF (ggrun_recipe==3.AND..NOT.no_ph_eos(iwork)) THEN
         icount=icount+1
         ind_rec3(:,icount)=inde(:)-ipoint(:)
      ENDIF
   ENDDO
   DEALLOCATE(inde)
   DEALLOCATE(ipoint)
   GOTO 100
ENDIF
!
!  test the compatibility of ngeo and fact_ngeo
!
IF (MOD(ngeo(1)-1, fact_ngeo(1))/=0 .OR. &
    MOD(ngeo(2)-1, fact_ngeo(2))/=0 .OR. &
    MOD(ngeo(3)-1, fact_ngeo(3))/=0 .OR. &
    MOD(ngeo(4)-1, fact_ngeo(4))/=0 .OR. &
    MOD(ngeo(5)-1, fact_ngeo(5))/=0 .OR. &
    MOD(ngeo(6)-1, fact_ngeo(6))/=0 ) CALL errore('initialize_no_ph', &
                                          'fact_ngeo incompatible with ngeo',1) 
no_ph_eos=.FALSE.
!
!  Here we set no_ph_eos based on ngeo_ph. Note that fact_ngeo should all
!  be 1 in this case
!
count_ngeo=0
DO igeo6 = 1, ngeo(6)
   todo(6)= (MOD(igeo6-1, fact_ngeo(6))==0)
   start = (ngeo(6)-ngeo_ph(6))/2
   todo(6)= todo(6).AND.(igeo6>start).AND.(igeo6<=start+ngeo_ph(6))
   DO igeo5 = 1, ngeo(5)
      todo(5)= (MOD(igeo5-1, fact_ngeo(5))==0)
      start = (ngeo(5)-ngeo_ph(5))/2
      todo(5)= todo(5).AND.(igeo5>start).AND.(igeo5<=start+ngeo_ph(5))
      DO igeo4 = 1, ngeo(4)
         todo(4)= (MOD(igeo4-1, fact_ngeo(4))==0)
         start = (ngeo(4)-ngeo_ph(4))/2
         todo(4)= todo(4).AND.(igeo4>start).AND.(igeo4<=start+ngeo_ph(4))
         DO igeo3 = 1, ngeo(3)
            todo(3)= (MOD(igeo3-1, fact_ngeo(3))==0)
            start = (ngeo(3)-ngeo_ph(3))/2
            todo(3)= todo(3).AND.(igeo3>start).AND.(igeo3<=start+ngeo_ph(3))
            DO igeo2 = 1, ngeo(2)
               todo(2)= (MOD(igeo2-1, fact_ngeo(2))==0)
               start = (ngeo(2)-ngeo_ph(2))/2
               todo(2)= todo(2).AND.(igeo2>start).AND.(igeo2<=start+ngeo_ph(2))
               DO igeo1 = 1, ngeo(1)
                  todo(1)= (MOD(igeo1-1, fact_ngeo(1))==0)
                  start = (ngeo(1)-ngeo_ph(1))/2
                  todo(1)= todo(1).AND.(igeo1>start).AND.&
                                       (igeo1<=start+ngeo_ph(1))
                  count_ngeo = count_ngeo + 1       
                  no_ph_eos(count_ngeo)= .NOT. (todo(1).AND.todo(2).AND.todo(3) &
                                     .AND.todo(4).AND.todo(5).AND.todo(6))
               END DO
            END DO
         END DO
      END DO
   END DO
END DO
100 CONTINUE
!
! If there are internal degrees of freedom we compute the phonons in 
! all geometries whose external geometries have to be calculated
! otherwise copy no_ph_eos in no_ph.
!
IF (linternal_thermo) THEN
   DO icount=1, tot_ngeo/ninternal
      DO iint_geom=1, ninternal
         no_ph((icount-1)*ninternal+iint_geom)=no_ph_eos(icount)
      ENDDO
   ENDDO
ELSE
   no_ph(1:tot_ngeo)=no_ph_eos(1:tot_ngeo)
ENDIF

RETURN
END SUBROUTINE initialize_no_ph
!
!-----------------------------------------------------------------------
INTEGER FUNCTION count_energies(ecutwfc, ecutrho, deltake, deltakeden, &
                                                             nke, nkeden)
!-----------------------------------------------------------------------
!
!   Given a set of cut-offs for the wavefunctions and the charge
!   densities, this functions counts how many possible calculations
!   can be done removing all those for which the cut-off for the
!   charge is smaller than 4 times the one for the wavefunctions
!
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER :: nke, nkeden
REAL(DP) :: ecutwfc, ecutrho, deltake, deltakeden

INTEGER :: icount, iden, ike
REAL(DP) :: keden, ke

icount=0
DO iden=1, nkeden
   keden=ecutrho + (iden-1) * deltakeden
   DO ike = 1, nke
      ke = ecutwfc + (ike-1) * deltake
      IF (keden/ke > 3.9999_DP) icount = icount + 1
   ENDDO
ENDDO
count_energies=icount

RETURN
END FUNCTION count_energies

!-----------------------------------------------------------------------
SUBROUTINE find_central_geo(ngeo, no_ph, central_geo)
!-----------------------------------------------------------------------
!
!
IMPLICIT NONE
INTEGER :: ngeo
LOGICAL :: no_ph(ngeo)
INTEGER :: central_geo

INTEGER :: igeo

central_geo=ngeo/2
IF (MOD(ngeo,2)==1) central_geo=central_geo+1
IF (no_ph(central_geo)) THEN
   DO igeo=1,ngeo/2
      central_geo=central_geo-igeo
      IF (.NOT. no_ph(central_geo)) EXIT
      central_geo=central_geo+2*igeo
      IF (.NOT. no_ph(central_geo)) EXIT
      central_geo=central_geo-igeo
   ENDDO
ENDIF

RETURN
END SUBROUTINE find_central_geo

!---------------------------------------------------------------------
SUBROUTINE initialize_flags_for_ph(nwork)
!---------------------------------------------------------------------
!
USE control_thermo, ONLY : lphonon, all_geometries_together, iqw,  &
                           irrw, comp_f_iw, geometry
USE initial_conf,   ONLY : collect_info_save
USE collect_info,   ONLY : something_to_do_all
USE freq_ph,        ONLY : fpol, nfs
USE images_omega,   ONLY : omega_group
USE mp_images,      ONLY : nimage
USE optical,        ONLY : start_freq, last_freq
USE control_ph,     ONLY : trans


IMPLICIT NONE
INTEGER, INTENT(IN) :: nwork
INTEGER :: i, j, iwork, igeom, iq, irr, start_omega

lphonon(1:nwork)=.TRUE.
IF (trans) THEN
   IF (all_geometries_together) THEN
      DO iwork=1,nwork
         igeom=geometry(iwork)
         iq=iqw(iwork)
         irr=irrw(iwork)
         lphonon(iwork)=something_to_do_all( &
                               collect_info_save(igeom),iwork,iq,irr)
      ENDDO
   ENDIF
ELSEIF (fpol) THEN
   ALLOCATE (comp_f_iw(1:nfs,nwork))
   IF (nimage>1) THEN
      comp_f_iw=.FALSE.
      DO iwork=1,nwork
         start_omega=(iwork-1)*omega_group
         DO i=1, omega_group
            j=MIN(start_omega+i, nfs)
            comp_f_iw(j,iwork)=.TRUE.
         ENDDO
         lphonon(iwork)=.FALSE.
         DO i=start_freq,last_freq
            lphonon(iwork)=lphonon(iwork).OR.comp_f_iw(i,iwork)
         ENDDO
      ENDDO
   ELSE
      comp_f_iw(:,1)=.TRUE.
   ENDIF
ENDIF
RETURN
END SUBROUTINE initialize_flags_for_ph
!
!---------------------------------------------------------------------
LOGICAL FUNCTION select_ph_to_do(iwork,nvar,inde,ipoint,ggrun_recipe)
!---------------------------------------------------------------------
!
!  This routine selects on which points to compute the phonon
!  dispersions in the case lgruneisen_gen=.TRUE.. It selects the 
!  point ipoint and all the 3^nvar
!  points about it (grun_recipe=2) or all the 2*nvar
!  points about it obtained by increasing or decreasing only 
!  one coordinate (grun_recipe=1) or the nvar * (nvar+3) /2 
!  points obtained by increasing or decreasing one coordinate
!  or increasing two of them (grun_recipe=3). 
!  Note that the user must provide an ipoint inside the mesh.  
!  If the ipoint is on the border this routine works but it 
!  will set to .TRUE. less points.
!  
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar, iwork, ggrun_recipe
INTEGER, INTENT(IN) :: ipoint(nvar)
INTEGER, INTENT(IN) :: inde(nvar)

INTEGER :: ivar, iaux
LOGICAL :: laux, laux1

iaux=0
DO ivar=1,nvar
   iaux=iaux+ABS(inde(ivar)-ipoint(ivar)) 
ENDDO
IF (ggrun_recipe==1) THEN
   laux=((iaux==0).OR.(iaux==1))
ELSEIF (ggrun_recipe==2) THEN
   laux=.TRUE.
   DO ivar=1,nvar
      laux=laux.AND.( ((inde(ivar)-ipoint(ivar))==0) .OR.  &
                      ((inde(ivar)-ipoint(ivar))==1) .OR.  &
                      ((inde(ivar)-ipoint(ivar))==-1) )
   ENDDO
ELSEIF (ggrun_recipe==3) THEN
   laux=((iaux==0).OR.(iaux==1))
   iaux=0
   DO ivar=1, nvar
      IF ((inde(ivar)-ipoint(ivar))==1) &
         iaux=iaux+(inde(ivar)-ipoint(ivar)) 
   ENDDO
   laux1=(iaux==2)
!
!  remove the points with indeces lower that ipoint 
!  or that in any direction are more distant than one point
!
   DO ivar=1, nvar
      IF (inde(ivar)< ipoint(ivar)) laux1=.FALSE.
      IF (ABS(inde(ivar)-ipoint(ivar))>1) laux1=.FALSE.
   ENDDO
   laux=laux.OR.laux1
ENDIF
select_ph_to_do=laux

RETURN
END FUNCTION select_ph_to_do
!
!---------------------------------------------------------------------
SUBROUTINE find_ipoint(iwork,nvar,nd,ipoint)
!---------------------------------------------------------------------
!
!  This routine receives the dimensions, nd(nvar), of a grid of nvar 
!  variables, one point iwork and sets the array ipoint with the coordinate
!  of iwork in the grid.
!  The grid is supposed to be ordered as
!  iwork
!  1     (1,1,...,1) (nvar variables)
!  2     (2,1,...,1)
!  ....
!  last  (nd(1),nd(2),...,nd(nvar)) 
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar, iwork
INTEGER, INTENT(IN) :: nd(nvar)
INTEGER, INTENT(OUT) :: ipoint(nvar)

INTEGER :: iaux, ivar
INTEGER :: inde(nvar), ncoef(nvar)

ncoef(1)=nd(1)
DO ivar=2,nvar
   ncoef(ivar)=ncoef(ivar-1)*nd(ivar)
ENDDO
IF (iwork<1.OR.iwork>ncoef(nvar)) CALL errore('find_ipoint','wrong iwork',1)

iaux=iwork
DO ivar=nvar,2,-1
   inde(ivar)= (iaux - 1) / ncoef(ivar-1) + 1
   iaux = MOD(iaux, ncoef(ivar-1)) 
   IF (iaux==0) iaux=iaux + ncoef(ivar-1)
ENDDO
inde(1)=iaux

ipoint(1:nvar)=inde(1:nvar)

RETURN
END SUBROUTINE find_ipoint
!
!---------------------------------------------------------------------
SUBROUTINE set_unperturbed_geometry(ngeom)
!---------------------------------------------------------------------
!
! This routine sets the unperturbed geometry for the elastic constants
! for the cases what='scf_elastic_constants' and 
! what='mur_lc_elastic_constants'. It allocates el_con_ibrav_geo and
! el_con_celldm_geo and copy on them the geometry read in input.
!
USE kinds, ONLY : DP
USE control_elastic_constants, ONLY : el_con_ibrav_geo, el_con_celldm_geo, &
                                      el_con_omega_geo, el_con_at_geo
USE equilibrium_conf, ONLY : celldm0, omega0, at0
USE initial_conf, ONLY : ibrav_save, at_save
IMPLICIT NONE
INTEGER, INTENT(OUT) :: ngeom
REAL(DP) :: compute_omega_geo

ngeom=1

ALLOCATE(el_con_ibrav_geo(1))
ALLOCATE(el_con_celldm_geo(6,1))
ALLOCATE(el_con_at_geo(3,3,1))
ALLOCATE(el_con_omega_geo(1))
el_con_ibrav_geo(1)=ibrav_save
el_con_celldm_geo(:,1)=celldm0(:)
el_con_at_geo(:,:,1)=at0(:,:)
IF (ibrav_save/=0) THEN
   omega0= compute_omega_geo(ibrav_save, celldm0)
ELSE
   CALL volume (1.0_dp, at0(1,1), at0(1,2), at0(1,3), omega0)
ENDIF
el_con_omega_geo(1)=omega0

RETURN
END SUBROUTINE set_unperturbed_geometry

!-----------------------------------------------------------------------
SUBROUTINE dtau_dinternal(celldm_geo, dtau_dint, nat, nint_var, &
                                                   iconstr_internal)
!-----------------------------------------------------------------------
!
!   This routine set the vector dtau_dint, the derivative of the atomic
!   positions with respect to the internal parameters
!
USE kinds,         ONLY : DP

IMPLICIT NONE
INTEGER, INTENT(IN) :: nat, nint_var
INTEGER, INTENT(IN) :: iconstr_internal
REAL(DP), INTENT(IN) :: celldm_geo(6)
REAL(DP), INTENT(INOUT) :: dtau_dint(3,nat,nint_var)

!
!  This routine provides the derivative of tau with respect to u
!
dtau_dint(:,:,:)=0.0_DP
IF ( iconstr_internal==1) THEN
!
!  In this constraint u is one dimensional and it is the parameter that
!  determines the atomic positions of the wurtzite structure. Atomic
!  coordinates on output are cartesian in units of celldm(1)
!
    dtau_dint(3,3,1)= celldm_geo(3)
    dtau_dint(3,4,1)= celldm_geo(3)
ELSEIF ( iconstr_internal==2) THEN
!
!  In this constraint the hcp structure is distorted with a strain
! (\epsilon,0,0,0,0,0). In this case the routine receives the displacement
! u along y (one dimensional) and gives as output the displacements of
! the coordinates of the two atoms of the hcp structure
!
    dtau_dint(2,1,1) = - 1.0_DP / 2.0_DP
    dtau_dint(2,2,1) =   1.0_DP / 2.0_DP
ELSEIF ( iconstr_internal==3) THEN
!
!  In this constraint the wurtzite structure is distorted with a strain
! (0,0,0,0,\epsilon,0). In this case the routine receives the displacement
! u along x (one dimensional) and gives as output the displacements of
! the cation moves along -x while the anion along +x
!
    dtau_dint(1,1,1) =   1.0_DP 
    dtau_dint(1,2,1) =   1.0_DP 
    dtau_dint(1,3,1) =  -1.0_DP 
    dtau_dint(1,4,1) =  -1.0_DP 
ENDIF

RETURN
END SUBROUTINE dtau_dinternal
