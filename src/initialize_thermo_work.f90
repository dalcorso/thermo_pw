!
! Copyright (C) 2013 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE initialize_thermo_work(nwork, part)
  !-----------------------------------------------------------------------
  !
  !  This is called by all images and initializes the global variables
  !  necessary to select the job to do. All images must have the same
  !  information.
  !
  USE kinds,      ONLY : DP
  USE thermo_mod, ONLY : what, alat_geo, step_ngeo, energy_geo, ngeo, &
                         celldm_geo, omega_geo, ibrav_geo, tot_ngeo
  USE control_thermo, ONLY : lpwscf, lbands, lphonon, lev_syn_1, lev_syn_2, &
                             lph, lpwscf_syn_1, lbands_syn_1, ldos, lq2r,   &
                             lmatdyn, ltherm, lconv_ke_test, lconv_nk_test, &
                             lstress, lelastic_const, &
                             lpiezoelectric_tensor, lberry, lpolarization
  USE control_conv,   ONLY : nke, ke, deltake, nkeden, deltakeden, keden, &
                             nnk, nk_test, deltank, nsigma, sigma_test, &  
                             deltasigma
  USE control_pwrun,  ONLY : celldm_save
  USE piezoelectric_tensor, ONLY : polar_geo
  USE control_elastic_constants, ONLY : elastic_algorithm, omega0, ibrav_save
  USE control_mur,    ONLY : vmin, vmin_input, vmax_input
  USE wvfct,          ONLY : ecutwfc
  USE gvect,          ONLY : ecutrho
  USE input_parameters, ONLY : nk1, nk2, nk3
  USE klist,          ONLY : degauss
  USE cell_base,  ONLY : alat
  USE grid_irr_iq, ONLY : irr_iq
  USE disp, ONLY : nqs
  !
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: nwork
  INTEGER, INTENT(IN) :: part
  INTEGER :: igeom, iq, irr, ike
  INTEGER :: iden, icount, ink, isigma
  REAL(DP) :: alat_new, celldm(6)
  REAL(DP) :: compute_omega_geo
  !
  nwork=0
!
!  In this part we do: all calculations for murnaghan
!
  IF (part == 1) THEN
     SELECT CASE (TRIM(what))
        CASE ( 'scf', 'scf_bands', 'scf_2d_bands', 'scf_ph', 'scf_disp' )
           ALLOCATE(alat_geo(1))
           ALLOCATE(energy_geo(1))
           tot_ngeo=0
           lpwscf_syn_1=.TRUE.
           lev_syn_1=.FALSE.
           IF ( TRIM(what)=='scf_ph'.OR. TRIM(what)=='scf_disp' ) THEN
              lph=.TRUE.
              tot_ngeo=1
           ENDIF
           IF ( TRIM(what)=='scf_disp' ) THEN 
              lq2r = .TRUE.
              ldos = .TRUE.
              lmatdyn = .TRUE.
              ltherm = .TRUE.
              CALL allocate_thermodynamics()
           ENDIF
           IF (what=='scf_bands'.OR.what=='scf_2d_bands') lbands_syn_1=.TRUE.
        CASE ( 'scf_ke') 
           nwork= nke * nkeden
           ALLOCATE(ke(nwork))
           ALLOCATE(keden(nwork))
           ALLOCATE(energy_geo(nwork))
           tot_ngeo=0
           icount=0
           DO iden=1, nkeden
              DO ike = 1, nke
                 icount = icount + 1
                 ke(icount) = ecutwfc + (ike-1) * deltake
                 keden(icount) = ecutrho + (iden-1) * deltakeden
              ENDDO
           ENDDO
           energy_geo=0.0_DP
           lconv_ke_test=.TRUE.
         CASE ( 'scf_nk' ) 
           nwork= nnk * nsigma
           ALLOCATE(nk_test(nwork))
           ALLOCATE(sigma_test(nwork))
           ALLOCATE(energy_geo(nwork))
           tot_ngeo=0
           icount=0
           DO isigma=1, nsigma
              DO ink = 1, nnk
                 icount = icount + 1
                 nk_test(icount)=nk1 + (ink - 1) * deltank
                 sigma_test(icount) = degauss + (isigma - 1) * deltasigma
              ENDDO
           ENDDO
           energy_geo=0.0_DP
           lconv_nk_test=.TRUE.
        CASE ('mur_lc', 'mur_lc_bands', 'mur_lc_ph', 'mur_lc_disp', &
              'mur_lc_elastic_constants', 'mur_lc_piezoelectric_tensor', &
              'mur_lc_polarization')
           nwork=ngeo(1)
           ALLOCATE(alat_geo(nwork))
           ALLOCATE(energy_geo(nwork))
           ALLOCATE(omega_geo(nwork))
           tot_ngeo=1
           DO igeom = 1, nwork
              alat_geo(igeom)=celldm_save(1)+(igeom-(ngeo(1)+1.0_DP)/2.0_DP)&
                                 *step_ngeo(1)
              celldm=0.0_DP
              celldm(1)=alat_geo(igeom)
              omega_geo(igeom)=compute_omega_geo(ibrav_save,celldm)
           ENDDO
           energy_geo=0.0_DP
           lev_syn_1=.TRUE.
           IF (TRIM(what)/='mur_lc'.AND.&
               TRIM(what)/='mur_lc_elastic_constants'.AND.&
               TRIM(what)/='mur_lc_piezoelectric_tensor'.AND.&
               TRIM(what)/='mur_lc_polarization') lpwscf_syn_1=.TRUE.
           IF ( TRIM(what)=='mur_lc_ph' .OR. TRIM(what)=='mur_lc_disp') THEN
              lph=.TRUE.
              tot_ngeo=1
           ENDIF
           IF ( TRIM(what)=='mur_lc_disp' ) THEN
              lq2r = .TRUE.
              ldos = .TRUE.
              lmatdyn = .TRUE.
              ltherm = .TRUE.
              CALL allocate_thermodynamics()
           ENDIF
           IF (what=='mur_lc_bands') lbands_syn_1=.TRUE.
        CASE ('mur_lc_t')
           nwork=ngeo(1)
           ALLOCATE(alat_geo(nwork))
           ALLOCATE(energy_geo(nwork))
           ALLOCATE(omega_geo(nwork))
           tot_ngeo=nwork
           DO igeom = 1, ngeo(1)
              alat_geo(igeom) = celldm_save(1) + &
                                (igeom-(ngeo(1)+1.0_DP)/2.0_DP)*step_ngeo(1)
              celldm=0.0_DP
              celldm(1)=alat_geo(igeom)
              omega_geo(igeom)=compute_omega_geo(ibrav_save,celldm)
           ENDDO
           energy_geo=0.0_DP
           lph = .TRUE.
           lq2r = .TRUE.
           ldos = .TRUE.
           ltherm = .TRUE.
           lmatdyn = .TRUE.
           lev_syn_1=.TRUE.
           lev_syn_2=.TRUE.
           CALL allocate_thermodynamics()
           CALL allocate_anharmonic()
        CASE ('elastic_constants', 'piezoelectric_tensor', 'polarization') 
!
!   in part 1 this case does nothing
!
        CASE DEFAULT
           CALL errore('initialize_thermo_work','what not recognized',1)
     END SELECT
  ELSE IF (part == 2 ) THEN
!
!  In this part we do the phonon calculations
!
     SELECT CASE (TRIM(what))
        CASE ( 'scf_ph', 'scf_disp', 'mur_lc_ph', 'mur_lc_disp', 'mur_lc_t' )
           nwork=0
           DO iq=1,nqs
              DO irr=0, irr_iq(iq)
                 nwork=nwork+1
              ENDDO
           ENDDO
        CASE ('elastic_constants', 'mur_lc_elastic_constants')
           IF (ALLOCATED(alat_geo)) DEALLOCATE(alat_geo)
           IF (ALLOCATED(energy_geo)) DEALLOCATE(energy_geo)
           IF (ALLOCATED(omega_geo)) DEALLOCATE(omega_geo)
           IF (elastic_algorithm=='standard') THEN
              CALL set_elastic_cons_work(nwork)
              ALLOCATE(celldm_geo(6,nwork))
              ALLOCATE(ibrav_geo(nwork))
           ELSEIF (elastic_algorithm=='advanced'.OR. &
                                       elastic_algorithm=='energy') THEN
!
!      celldm_geo and ibrav_geo are allocated inside the routine
!
              CALL set_elastic_cons_work_adv(nwork)
           ELSE
              CALL errore('initialize_thermo_work',&
                                        'elastic_algorithm unknown',1)
           ENDIF
           ALLOCATE(alat_geo(nwork))
           ALLOCATE(energy_geo(nwork))
           ALLOCATE(omega_geo(nwork))
           energy_geo=0.0_DP
           IF (elastic_algorithm=='energy') THEN
              DO igeom = 1, nwork
                 omega_geo(igeom)=compute_omega_geo(ibrav_geo(igeom),&
                                                    celldm_geo(1,igeom))
              ENDDO
              omega0 = compute_omega_geo(ibrav_save, celldm_save)
           ENDIF
           lelastic_const=.TRUE.
        CASE ('piezoelectric_tensor', 'mur_lc_piezoelectric_tensor')
           IF (ALLOCATED(alat_geo)) DEALLOCATE(alat_geo)
           IF (ALLOCATED(energy_geo)) DEALLOCATE(energy_geo)
           IF (ALLOCATED(omega_geo)) DEALLOCATE(omega_geo)
           CALL set_piezo_tensor_work(nwork) 
           ALLOCATE(alat_geo(nwork))
           ALLOCATE(energy_geo(nwork))
           ALLOCATE(omega_geo(nwork))
           lpiezoelectric_tensor=.TRUE.
        CASE ('polarization', 'mur_lc_polarization')
           IF (ALLOCATED(alat_geo)) DEALLOCATE(alat_geo)
           IF (ALLOCATED(energy_geo)) DEALLOCATE(energy_geo)
           IF (ALLOCATED(omega_geo)) DEALLOCATE(omega_geo)
           nwork=1
           lpolarization=.TRUE.
           ALLOCATE(polar_geo(3,nwork))
           ALLOCATE(alat_geo(nwork))
           ALLOCATE(energy_geo(nwork))
           ALLOCATE(omega_geo(nwork))
           polar_geo=0.0_DP
     END SELECT
  ELSE
     CALL errore('initialize_thermo_work','unknown part',1)
  END IF

  IF ( nwork == 0 ) RETURN

  ALLOCATE( lpwscf(nwork) )
  ALLOCATE( lstress(nwork) )
  ALLOCATE( lbands(nwork) )
  ALLOCATE( lberry(nwork) )
  ALLOCATE( lphonon(nwork) )
  lpwscf  = .FALSE.
  lstress = .FALSE.
  lbands  = .FALSE.
  lberry  = .FALSE.
  lphonon = .FALSE.

  IF (part == 1) THEN
     SELECT CASE (TRIM(what))
        CASE ('scf', 'scf_bands', 'scf_2d_bands', 'scf_ph', 'scf_disp', &
            'elastic_constants', 'piezoelectric_tensor','polarization')
        CASE ('scf_ke', 'scf_nk', 'mur_lc', 'mur_lc_bands', 'mur_lc_ph', &
              'mur_lc_disp', 'mur_lc_t', 'mur_lc_elastic_constants', &
                          'mur_lc_piezoelectric_tensor', &
                          'mur_lc_polarization' )
           lpwscf(1:nwork)=.TRUE.
        CASE DEFAULT
          CALL errore('initialize_thermo_work','unknown what',1)
     END SELECT
  ELSE IF ( part == 2 ) THEN
     SELECT CASE (TRIM(what))
        CASE ('scf_ph', 'scf_disp', 'mur_lc_ph', 'mur_lc_disp','mur_lc_t')
           lphonon(1:nwork)=.TRUE.
        CASE ('elastic_constants', 'mur_lc_elastic_constants')
           lpwscf(1:nwork)=.TRUE.
           lstress(1:nwork)=.TRUE.
        CASE ('piezoelectric_tensor', 'mur_lc_piezoelectric_tensor',&
              'polarization', 'mur_lc_polarization')
           lpwscf(1:nwork)=.TRUE.
           lberry(1:nwork)=.TRUE.
     END SELECT
  END IF

  RETURN
END SUBROUTINE initialize_thermo_work
