!
! Copyright (C) 2013-2019 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE set_thermo_work_todo(iwork, part, iq_point, irr_value, auxdyn_loc)
  !-----------------------------------------------------------------------
!
!  This routine receives from the asynchronous driver the work to do in
!  the integer iwork and sets the input variables of pwscf or of the
!  phonon according to iwork before performing the actual calculation.
!  On output it gives the q point and the irreducible representation, if
!  this correspond to a phonon calculation.
!
  USE kinds,            ONLY : DP
  USE thermo_mod,       ONLY : what, celldm_geo, start_geometry, &
                               last_geometry
  USE control_thermo,   ONLY : outdir_thermo, lstress, lphonon, &
                               all_geometries_together, geometry, irrw, &
                               iqw, comp_f_work
  USE control_elastic_constants, ONLY : frozen_ions
  USE control_conv,     ONLY : ke, keden, nk_test, sigma_test
  USE initial_conf,     ONLY : ibrav_save, tau_save_crys, collect_info_save
  USE equilibrium_conf, ONLY : at0, tau0
  USE images_omega,     ONLY : omega_group
  USE control_qe,       ONLY : use_ph_images
  USE collect_info,     ONLY : copy_collect_info
!
!  the library modules
!
  USE elastic_constants, ONLY : epsilon_geo
  USE strain_mod,        ONLY : apply_strain, print_strain
!
!  the pw variables that are set here or used to set the input
!
  USE input_parameters, ONLY : electron_maxstep, k_points, xk, wk, k1, k2, &
                               k3, nkstot, etot_conv_thr, forc_conv_thr
  USE control_flags,    ONLY : lbfgs, nstep, niter 
  USE cell_base,   ONLY : cell_base_init, at
  USE ions_base,   ONLY : tau, nat
  USE gvecw,       ONLY : ecutwfc
  USE gvect,       ONLY : ecutrho
  USE gvecs,       ONLY : dual
  USE force_mod,   ONLY : lforce, lstres
  USE relax,       ONLY : epse, epsf
  USE start_k,     ONLY : init_start_k
  USE klist,       ONLY : degauss
  USE freq_ph,     ONLY : fpol, nfs
  USE io_files,    ONLY : tmp_dir, wfc_dir, check_tempdir
!
!   the phonon variables set here or used to set the input
!
  USE grid_irr_iq, ONLY : irr_iq, comp_irr_iq, done_irr_iq
  USE disp,        ONLY : nqs, comp_iq, done_iq
  USE control_ph,       ONLY : epsil, trans, recover
  USE images_omega, ONLY : comp_f
  USE mp_images,    ONLY : nimage

  USE io_global,   ONLY : stdout
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iwork, part
  INTEGER, INTENT(OUT) :: iq_point, irr_value
  CHARACTER(LEN=256), INTENT(INOUT) :: auxdyn_loc

  INTEGER :: jwork, irr, iq, i, j, ia, nk1, nk2, nk3, ibrav, start_omega, &
             igeom, startge, lastge, target_nqs, target_irr_iq, image
  REAL(DP) :: rd_ht(3,3), zero, celldm(6)
  CHARACTER(LEN=10) :: cell_units
  CHARACTER(LEN=6) :: int_to_char
  CHARACTER(LEN=256) :: outdir
  LOGICAL :: exst, parallelfs, trd_ht
  !
  iq_point=0
  irr_value=0
  zero=0.0_DP

  IF (part == 1) THEN
     SELECT CASE (TRIM(what))
!
!  first the cases where there is nothing to do
!
        CASE ('scf',              &
              'scf_bands',        &
              'scf_dos',          &
              'scf_ph',           &
              'scf_disp') 
!
!  then the cases in which we set the kinetic energy and the k points
!
        CASE ('scf_ke')
           ecutwfc = ke(iwork)
           ecutrho = keden(iwork)
           dual = ecutrho / ecutwfc
           CALL set_fft_mesh()
           outdir=TRIM(outdir_thermo)//'/ke'//TRIM(int_to_char(iwork))//'/'
           tmp_dir = TRIM ( outdir )
           wfc_dir = tmp_dir
           CALL check_tempdir ( tmp_dir, exst, parallelfs )
        CASE ('scf_nk')
           degauss = sigma_test(iwork)
           nk1=nk_test(1,iwork)
           nk2=nk_test(2,iwork)
           nk3=nk_test(3,iwork)
           IF (TRIM(k_points) /='automatic') &
              CALL errore('set_thermo_work_todo', &
                               'kpoint test requires automatic k point',1)
!
!   for the shift we use the same parameters read in input
!
           CALL init_start_k ( nk1, nk2, nk3, k1, k2, k3, k_points, &
                               nkstot, xk, wk )
           CALL set_fft_mesh()
           outdir=TRIM(outdir_thermo)//'/nk'//TRIM(int_to_char(iwork))//'/'
           tmp_dir = TRIM ( outdir )
           wfc_dir = tmp_dir
           CALL check_tempdir ( tmp_dir, exst, parallelfs )
!
!   then all the cases that require many energies calculations at
!   different geometries. The geometry is set here in the pwscf variables
!
        CASE ('mur_lc',                      &
              'mur_lc_bands',                &
              'mur_lc_dos',                  &
              'mur_lc_ph',                   &
              'mur_lc_disp',                 &
              'mur_lc_t',                    &
              'mur_lc_elastic_constants',    &
              'mur_lc_piezoelectric_tensor', &
              'mur_lc_polarization')

           IF (frozen_ions) THEN
              lbfgs=.FALSE.
           ELSE
              lforce=.TRUE.
              lstres=.TRUE.
              lbfgs = .TRUE.
              nstep = 20
              epse = etot_conv_thr
              epsf = forc_conv_thr
           ENDIF
!
!   now set the celldm. 
!
           celldm(:)=celldm_geo(:,iwork)
           rd_ht=0.0_DP
           CALL cell_base_init ( ibrav_save, celldm, zero, zero, zero, zero, &
                                     zero, zero, .FALSE., rd_ht, ' ' )
           CALL set_fft_mesh()
!
! strain uniformly the coordinates to the new celldm
!
           tau=tau_save_crys
           CALL cryst_to_cart( nat, tau, at, 1 )

           outdir=TRIM(outdir_thermo)//'/g'//TRIM(int_to_char(iwork))//'/'
           tmp_dir = TRIM ( outdir )
           wfc_dir = tmp_dir
           CALL check_tempdir ( tmp_dir, exst, parallelfs )

           IF (.NOT.frozen_ions) CALL clean_bfgs_history()
!
!    the case of quasi-harmonic elastic constants
!
        CASE ('elastic_constants_t')
           WRITE(stdout,'(/,2x,76("-"))')
           niter = electron_maxstep
           IF (frozen_ions) THEN
              lstres=lstress(iwork)
              lbfgs=.FALSE.
           ELSE
              lforce=.TRUE.
              lstres=lstress(iwork)
              lbfgs = .TRUE.
              nstep = 20
              epse = etot_conv_thr
              epsf = forc_conv_thr
           ENDIF
           CALL set_work_for_elastic_const(iwork)
     
        CASE DEFAULT
           CALL errore('set_thermo_work_todo','unknown what',1)
     END SELECT
  ELSE IF (part==2) THEN
     SELECT CASE (TRIM(what))
!
!   here we set the representation and the q point to calculate in the
!   phonon calculation for the present geometry
!
        CASE ('scf_ph',         &
              'scf_disp',       &
              'mur_lc_ph',      &
              'mur_lc_disp',    &
              'mur_lc_t',       &
              'elastic_constants_t')
           IF (.NOT.lphonon(iwork)) RETURN
           igeom=1
           iq_point=1
           irr_value=0
           IF (all_geometries_together) THEN
               igeom=geometry(iwork)
               CALL initialize_ph_geometry(igeom, auxdyn_loc)
           ENDIF
           comp_irr_iq=.FALSE.
           comp_iq=.FALSE.
           comp_f=.FALSE.
           IF (trans) THEN
              IF (use_ph_images) THEN
                 image=MOD(iwork-1,nimage)+1
                 CALL copy_collect_info(collect_info_save(igeom), nqs, &
                            nat, nimage, image, comp_irr_iq, done_irr_iq, &
                            comp_iq, done_iq, irr_iq)
                 iq_point=image
              ELSE
                 comp_irr_iq(irrw(iwork),iqw(iwork))=.TRUE.
                 comp_iq(iqw(iwork))=.TRUE.
                 iq_point=iqw(iwork)
                 irr_value=irrw(iwork)
              ENDIF
           ELSEIF (fpol) THEN
              comp_iq(1)=.TRUE.
              comp_irr_iq(0,1)=.TRUE.
              comp_f(:)=comp_f_work(:,iwork)
           ELSEIF (epsil) THEN
              comp_iq(1)=.TRUE.
              comp_irr_iq(0,1)=.TRUE.
           ENDIF
!
!    Here the elastic constant calculation
!
        CASE ('scf_elastic_constants', 'mur_lc_elastic_constants')
           niter = electron_maxstep
           IF (frozen_ions) THEN
              lstres=lstress(iwork)
              lbfgs=.FALSE.
           ELSE
              lforce=.TRUE.
              lstres=lstress(iwork)
              lbfgs = .TRUE.
              nstep = 20
              epse = etot_conv_thr
              epsf = forc_conv_thr
           ENDIF
           CALL set_work_for_elastic_const(iwork)
        CASE ('scf_piezoelectric_tensor', 'mur_lc_piezoelectric_tensor')
           ibrav=0
           niter = electron_maxstep
           DO i=1, 3
              CALL apply_strain(at0(1,i), at(1,i), epsilon_geo(1,1,iwork))
           ENDDO
           DO ia=1,nat
              CALL apply_strain(tau0(1,ia), tau(1,ia), epsilon_geo(1,1,iwork))
           ENDDO
           WRITE(stdout,'(/,2x,76("-"))')
           CALL print_strain(epsilon_geo(:,:,iwork))
           IF (frozen_ions) THEN
              lstres=.TRUE.
              lbfgs=.FALSE.
           ELSE
              lforce=.TRUE.
              lstres=.TRUE.
              lbfgs = .TRUE.
              nstep = 20
              epse = etot_conv_thr
              epsf = forc_conv_thr
           ENDIF
           rd_ht = TRANSPOSE( at ) 
           trd_ht=.TRUE.
           cell_units='alat'
           CALL cell_base_init ( ibrav, celldm, zero, zero, zero, zero, &
                         zero, zero, trd_ht, rd_ht, cell_units )
           CALL set_fft_mesh()

           outdir=TRIM(outdir_thermo)//'/g'//TRIM(int_to_char(iwork))//'/'
           tmp_dir = TRIM ( outdir )
           wfc_dir = tmp_dir
           CALL check_tempdir ( tmp_dir, exst, parallelfs )

           IF (.NOT.frozen_ions) CALL clean_bfgs_history()

        CASE ('scf_polarization','mur_lc_polarization')
     END SELECT
  ELSE
     CALL errore('set_thermo_work_todo','unknown part',1)
  END IF
  !
  RETURN
  !
END SUBROUTINE set_thermo_work_todo
