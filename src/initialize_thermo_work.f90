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
  USE thermo_mod, ONLY : what, alat_geo, step_ngeo, energy_geo, ngeo, omega_geo
  USE control_thermo, ONLY : lpwscf, lbands, lphonon, lev_syn_1, lev_syn_2, &
                             lph, lpwscf_syn_1, lbands_syn_1, ldos, lq2r,   &
                             lmatdyn, ltherm, lconv_ke_test, lconv_nk_test, &
                             compute_lc, lstress, lelastic_const
  USE control_conv,   ONLY : nke, ke, deltake, nkeden, deltakeden, keden, &
                             nnk, nk_test, deltank, nsigma, sigma_test, &  
                             deltasigma
  USE control_mur,    ONLY : vmin
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
  REAL(DP) :: alat_new
  REAL(DP) :: compute_omega_geo
  !
  nwork=0
!
!  In this part we do: all calculations for murnaghan
!
  compute_lc=.TRUE.
  IF (part == 1) THEN
     SELECT CASE (TRIM(what))
        CASE ( 'scf', 'scf_bands', 'scf_ph', 'scf_disp' )
           IF (ALLOCATED(energy_geo)) compute_lc=.FALSE.
           IF (compute_lc) THEN
              ALLOCATE(alat_geo(ngeo))
              ALLOCATE(energy_geo(ngeo))
              lpwscf_syn_1=.TRUE.
              lev_syn_1=.FALSE.
              IF ( TRIM(what)=='scf_ph'.OR. TRIM(what)=='scf_disp' ) lph=.TRUE.
              IF ( TRIM(what)=='scf_disp' ) THEN 
                 lq2r = .TRUE.
                 ldos = .TRUE.
                 lmatdyn = .TRUE.
                 ltherm = .TRUE.
                 CALL allocate_thermodynamics()
              ENDIF
              IF (what=='scf_bands') lbands_syn_1=.TRUE.
           ENDIF
        CASE ( 'scf_ke') 
           nwork= nke * nkeden
           IF (ALLOCATED(energy_geo)) compute_lc=.FALSE.
           IF (compute_lc) THEN
              ALLOCATE(ke(nwork))
              ALLOCATE(keden(nwork))
              ALLOCATE(energy_geo(nwork))
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
           ENDIF
         CASE ( 'scf_nk' ) 
           nwork= nnk * nsigma
           IF (ALLOCATED(energy_geo)) compute_lc=.FALSE.
           IF (compute_lc) THEN
              ALLOCATE(nk_test(nwork))
              ALLOCATE(sigma_test(nwork))
              ALLOCATE(energy_geo(nwork))
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
           ENDIF
        CASE ('mur_lc', 'mur_lc_bands', 'mur_lc_ph', 'mur_lc_disp', &
              'mur_lc_elastic_constants')
           nwork=ngeo
           IF ( ALLOCATED(energy_geo) ) THEN
              compute_lc=.FALSE.
              alat_new=alat_geo(ngeo/2+1)*( vmin / omega_geo(ngeo/2+1) ) &
                                         ** (1.0_DP /3.0_DP)
              IF ( ABS(alat_geo(ngeo/2+1)-alat_new) > step_ngeo ) THEN
                 alat=alat_new
                 compute_lc=.TRUE.
              ENDIF
           ELSE
              ALLOCATE(alat_geo(ngeo))
              ALLOCATE(energy_geo(ngeo))
              ALLOCATE(omega_geo(ngeo))
           ENDIF
           IF (compute_lc) THEN
              DO igeom = 1, ngeo
                alat_geo(igeom)=alat+(igeom-(ngeo+1.0_DP)/2.0_DP)*step_ngeo
                omega_geo(igeom)=compute_omega_geo(alat_geo(igeom))
              ENDDO
              energy_geo=0.0_DP
              lev_syn_1=.TRUE.
              IF (TRIM(what)/='mur_lc'.AND.&
                  TRIM(what)/='mur_lc_elastic_constants') lpwscf_syn_1=.TRUE.
              IF ( TRIM(what)=='mur_lc_ph' .OR. TRIM(what)=='mur_lc_disp') lph=.TRUE.
              IF ( TRIM(what)=='mur_lc_disp' ) THEN
                 lq2r = .TRUE.
                 ldos = .TRUE.
                 lmatdyn = .TRUE.
                 ltherm = .TRUE.
                 CALL allocate_thermodynamics()
              ENDIF
              IF (what=='mur_lc_bands') lbands_syn_1=.TRUE.
           ENDIF
        CASE ('mur_lc_t')
           nwork=ngeo
           IF ( ALLOCATED(energy_geo) ) THEN
              compute_lc=.FALSE.
              alat_new=alat_geo(ngeo/2+1)*( vmin / omega_geo(ngeo/2+1) ) &
                                         ** (1.0_DP /3.0_DP)
              IF ( ABS(alat_geo(ngeo/2+1)-alat_new) > step_ngeo ) THEN
                 alat=alat_new
                 compute_lc=.TRUE.
              ENDIF
           ELSE
              ALLOCATE(alat_geo(ngeo))
              ALLOCATE(energy_geo(ngeo))
              ALLOCATE(omega_geo(ngeo))
           ENDIF
           IF (compute_lc) THEN
              DO igeom = 1, ngeo
                 alat_geo(igeom) = alat + (igeom-(ngeo+1.0_DP)/2.0_DP)*step_ngeo
                 omega_geo(igeom)=compute_omega_geo(alat_geo(igeom))
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
           ENDIF
        CASE ('elastic_constants') 
!
!   in part 1 this case does nothing
!
             compute_lc=.FALSE.
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
           CALL set_elastic_cons_work(nwork) 
           ALLOCATE(alat_geo(nwork))
           ALLOCATE(energy_geo(nwork))
           ALLOCATE(omega_geo(nwork))
           lelastic_const=.TRUE.
     END SELECT
  ELSE
     CALL errore('initialize_thermo_work','unknown part',1)
  END IF

  IF ( nwork == 0 .OR. .NOT. compute_lc ) RETURN

  ALLOCATE( lpwscf(nwork) )
  ALLOCATE( lstress(nwork) )
  ALLOCATE( lbands(nwork) )
  ALLOCATE( lphonon(nwork) )
  lpwscf  = .FALSE.
  lstress = .FALSE.
  lbands  = .FALSE.
  lphonon = .FALSE.

  IF (part == 1) THEN
     SELECT CASE (TRIM(what))
        CASE ('scf', 'scf_bands', 'scf_ph', 'scf_disp', 'elastic_constants')
        CASE ('mur_lc', 'mur_lc_bands', 'mur_lc_ph', 'mur_lc_disp', &
              'mur_lc_t', 'mur_lc_elastic_constants' )
           lpwscf(1:ngeo)=.TRUE.
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
     END SELECT
  END IF

  RETURN
END SUBROUTINE initialize_thermo_work
