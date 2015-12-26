!
! Copyright (C) 2013-2015 Andrea Dal Corso
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
  USE thermo_mod, ONLY : what, step_ngeo, energy_geo, ngeo, &
                         celldm_geo, omega_geo, ibrav_geo, tot_ngeo, no_ph
  USE control_thermo, ONLY : lpwscf, lbands, lphonon, lev_syn_1, lev_syn_2, &
                             lph, lpwscf_syn_1, lbands_syn_1, ldos, lq2r,   &
                             lmatdyn, ltherm, lconv_ke_test, lconv_nk_test, &
                             lstress, lelastic_const, lpiezoelectric_tensor,&
                             lberry, lpolarization, lpart2_pw, do_scf_relax
  USE control_conv,   ONLY : nke, ke, deltake, nkeden, deltakeden, keden, &
                             nnk, nk_test, deltank, nsigma, sigma_test, &  
                             deltasigma
  USE control_pwrun,  ONLY : celldm_save, ibrav_save, do_punch
  USE piezoelectric_tensor, ONLY : polar_geo
  USE control_elastic_constants, ONLY : elastic_algorithm, omega0
  USE control_mur,    ONLY : vmin, vmin_input, vmax_input
  USE wvfct,          ONLY : ecutwfc
  USE gvect,          ONLY : ecutrho
  USE control_quadratic_energy, ONLY : nvar, degree
  USE input_parameters, ONLY : nk1, nk2, nk3
  USE klist,          ONLY : degauss
  USE cell_base,  ONLY : alat
  !
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: nwork
  INTEGER, INTENT(IN) :: part
  INTEGER :: igeom, iq, irr, ike
  INTEGER :: iden, icount, ink, isigma
  REAL(DP) :: compute_omega_geo
  !
  nwork=0
  IF (part == 1) THEN
     SELECT CASE (TRIM(what))
!
!   In these cases we do not do any asyncronous work in the first part
!
        CASE ( 'scf') 
           ALLOCATE(energy_geo(1))
           lpwscf_syn_1=.TRUE.
        CASE ('scf_bands') 
           ALLOCATE(energy_geo(1))
           lpwscf_syn_1=.TRUE.
           lbands_syn_1=.TRUE.
        CASE ('scf_2d_bands')
           ALLOCATE(energy_geo(1))
           lpwscf_syn_1=.TRUE.
           lbands_syn_1=.TRUE.
        CASE ('scf_ph') 
           ALLOCATE(energy_geo(1))
           lpwscf_syn_1=.TRUE.
           lph=.TRUE.
           tot_ngeo=1
           ALLOCATE(no_ph(tot_ngeo))
           no_ph(1)=.FALSE.
        CASE ('scf_disp')
           ALLOCATE(energy_geo(1))
           lpwscf_syn_1=.TRUE.
           lph=.TRUE.
           tot_ngeo=1
           ALLOCATE(no_ph(tot_ngeo))
           no_ph(1)=.FALSE.
           lq2r = .TRUE.
           lmatdyn = .TRUE.
           ldos = .TRUE.
           ltherm = .TRUE.
           CALL allocate_thermodynamics()
!
!   In these cases we make asyncronous work in the first part
!
        CASE ( 'scf_ke') 
           nwork= nke * nkeden
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
           do_punch=.FALSE.
        CASE ( 'scf_nk' ) 
           nwork= nnk * nsigma
           ALLOCATE(nk_test(3,nwork))
           ALLOCATE(sigma_test(nwork))
           ALLOCATE(energy_geo(nwork))
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
!
!   in part 1 these cases do nothing
!
        CASE ('scf_elastic_constants') 
           lpart2_pw=.TRUE.
        CASE ('scf_piezoelectric_tensor')
           lpart2_pw=.TRUE.
        CASE ('scf_polarization') 
           lpart2_pw=.TRUE.
!
!   here all the cases that require the determination of the minimization
!   of the energy to find the equilibrium crystal parameters
!
        CASE ('mur_lc')
           CALL initialize_mur(nwork)
           do_punch=.FALSE.
        CASE ('mur_lc_bands') 
           CALL initialize_mur(nwork)
           do_punch=.FALSE.
           lpwscf_syn_1=.TRUE.
           lbands_syn_1=.TRUE.
        CASE ('mur_lc_ph') 
           CALL initialize_mur(nwork)
           do_punch=.FALSE.
           lpwscf_syn_1=.TRUE.
           lph=.TRUE.
           tot_ngeo=1
           ALLOCATE(no_ph(tot_ngeo))
           no_ph(1)=.FALSE.
        CASE ('mur_lc_disp')
           CALL initialize_mur(nwork)
           do_punch=.FALSE.
           lpwscf_syn_1=.TRUE.
           lph=.TRUE.
           tot_ngeo=1
           ALLOCATE(no_ph(tot_ngeo))
           no_ph(1)=.FALSE.
           lq2r = .TRUE.
           ldos = .TRUE.
           lmatdyn = .TRUE.
           ltherm = .TRUE.
           CALL allocate_thermodynamics()
        CASE ('mur_lc_elastic_constants') 
           CALL initialize_mur(nwork)
           lpwscf_syn_1=do_scf_relax
           lpart2_pw=.TRUE.
           do_punch=.FALSE.
        CASE ('mur_lc_piezoelectric_tensor') 
           CALL initialize_mur(nwork)
           lpwscf_syn_1=do_scf_relax
           lpart2_pw=.TRUE.
           do_punch=.FALSE.
        CASE ('mur_lc_polarization')
           CALL initialize_mur(nwork)
           lpwscf_syn_1=do_scf_relax
           lpart2_pw=.TRUE.
           do_punch=.FALSE.
!
!    Here all the cases that compute the free energy and minimize it
!

        CASE ('mur_lc_t')
           CALL initialize_mur(nwork)
           lph = .TRUE.
           tot_ngeo=nwork
           ALLOCATE(no_ph(tot_ngeo))
           CALL initialize_no_ph(no_ph, tot_ngeo)
           lq2r = .TRUE.
           ldos = .TRUE.
           ltherm = .TRUE.
           lmatdyn = .TRUE.
           lev_syn_2=.TRUE.
           CALL compute_degree(ibrav_save, degree, nvar)
           CALL allocate_thermodynamics()
           CALL allocate_anharmonic()
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
        CASE ('scf_ph')
           CALL initialize_ph_work(nwork)
        CASE ('scf_disp') 
           CALL initialize_ph_work(nwork)
        CASE ('mur_lc_ph') 
           CALL initialize_ph_work(nwork)
        CASE ('mur_lc_disp') 
           CALL initialize_ph_work(nwork)
        CASE ('mur_lc_t')
           CALL initialize_ph_work(nwork)
        CASE ('scf_elastic_constants', 'mur_lc_elastic_constants')
           IF (ALLOCATED(ibrav_geo))  DEALLOCATE(ibrav_geo)
           IF (ALLOCATED(celldm_geo)) DEALLOCATE(celldm_geo)
           IF (ALLOCATED(energy_geo)) DEALLOCATE(energy_geo)
           IF (ALLOCATED(omega_geo))  DEALLOCATE(omega_geo)
!
!     set_elastic_constant_work or set_elastic_constants_work_adv 
!     allocate ibrav_geo and celldm_geo
!
           IF (elastic_algorithm=='standard') THEN
              CALL set_elastic_cons_work(nwork)
           ELSEIF (elastic_algorithm=='advanced'.OR. &
                                       elastic_algorithm=='energy') THEN
              CALL set_elastic_cons_work_adv(nwork)
           ELSE
              CALL errore('initialize_thermo_work',&
                                        'elastic_algorithm unknown',1)
           ENDIF
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
           do_punch=.FALSE.
           CALL allocate_debye()
        CASE ('scf_piezoelectric_tensor', 'mur_lc_piezoelectric_tensor')
           IF (ALLOCATED(energy_geo)) DEALLOCATE(energy_geo)
           IF (ALLOCATED(omega_geo)) DEALLOCATE(omega_geo)
           CALL set_piezo_tensor_work(nwork) 
           ALLOCATE(energy_geo(nwork))
           ALLOCATE(omega_geo(nwork))
           lpiezoelectric_tensor=.TRUE.
           do_punch=.TRUE.
        CASE ('scf_polarization', 'mur_lc_polarization')
           IF (ALLOCATED(energy_geo)) DEALLOCATE(energy_geo)
           IF (ALLOCATED(omega_geo)) DEALLOCATE(omega_geo)
           nwork=1
           lpolarization=.TRUE.
           ALLOCATE(polar_geo(3,nwork))
           ALLOCATE(energy_geo(nwork))
           ALLOCATE(omega_geo(nwork))
           polar_geo=0.0_DP
           do_punch=.TRUE.
     END SELECT
  ELSE
     CALL errore('initialize_thermo_work','unknown part',1)
  END IF

  IF ( nwork == 0 ) RETURN

  ALLOCATE( lpwscf(nwork) )
  ALLOCATE( lbands(nwork) )
  ALLOCATE( lstress(nwork) )
  ALLOCATE( lberry(nwork) )
  ALLOCATE( lphonon(nwork) )
  lpwscf  = .FALSE.
  lbands  = .FALSE.
  lstress = .FALSE.
  lberry  = .FALSE.
  lphonon = .FALSE.

  IF (part == 1) THEN
     SELECT CASE (TRIM(what))
!
!  First the cases in which there is no asynchronous work in the first part
!
        CASE ( 'scf',                        &
               'scf_bands',                  &
               'scf_2d_bands',               &
               'scf_ph',                     & 
               'scf_disp',                   &
               'scf_elastic_constants',      &
               'scf_piezoelectric_tensor',   &
               'scf_polarization' )
!
!  Then the cases in which there is pwscf asynchronous work in the first part
!
        CASE ('scf_ke',                      &
              'scf_nk',                      &
              'mur_lc',                      &
              'mur_lc_bands',                &
              'mur_lc_ph',                   &
              'mur_lc_disp',                 &
              'mur_lc_t',                    &
              'mur_lc_elastic_constants',    &
              'mur_lc_piezoelectric_tensor', &
              'mur_lc_polarization' )
           lpwscf(1:nwork)=.TRUE.
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
           lphonon(1:nwork)=.TRUE.
!
!  Here the cases in which there are pwscf calculation in the second part
!
        CASE ('scf_elastic_constants', 'mur_lc_elastic_constants')
           lpwscf(1:nwork)=.TRUE.
           IF (elastic_algorithm=='standard'.OR.elastic_algorithm=='advanced')&
              lstress(1:nwork)=.TRUE.
!
!  Here the cases in which there are pwscf and berry phase calculation 
!  in the second part
!
        CASE ('scf_piezoelectric_tensor', 'mur_lc_piezoelectric_tensor', &
              'scf_polarization', 'mur_lc_polarization')
           lpwscf(1:nwork)=.TRUE.
           lberry(1:nwork)=.TRUE.
     END SELECT
  END IF

  RETURN
END SUBROUTINE initialize_thermo_work

SUBROUTINE clean_ngeo(ngeo,fact_ngeo,ibrav)

USE control_mur, ONLY : lmurn
IMPLICIT NONE
INTEGER, INTENT(IN) :: ibrav
INTEGER, INTENT(INOUT) :: ngeo(6), fact_ngeo(6)
INTEGER :: ngeo_aux(6), fact_ngeo_aux(6)

INTEGER :: i

ngeo_aux=1
ngeo_aux(1)=ngeo(1)
fact_ngeo_aux=1
fact_ngeo_aux(1)=fact_ngeo(1)

IF (.NOT.lmurn) THEN
   SELECT CASE (ibrav)
      CASE(1,2,3)
      CASE (4,6,7)
         IF (ngeo(3) /= 0) THEN
             ngeo_aux(3)=ngeo(3)
         ELSE
             ngeo_aux(3)=5
         ENDIF
         fact_ngeo_aux(3)=fact_ngeo(3)
      CASE (5)
         IF (ngeo(4) /= 0) THEN
            ngeo_aux(4)=ngeo(4)
         ELSE
            ngeo_aux(4)=ngeo(4)
         ENDIF
         fact_ngeo_aux(4)=fact_ngeo(4)
      CASE(8,9,-9,91,10,11)
         IF (ngeo(2) /= 0) THEN
            ngeo_aux(2)=ngeo(2)
         ELSE
            ngeo_aux(2)=5
         ENDIF
         IF (ngeo(3) /= 0) THEN
            ngeo_aux(3)=ngeo(3)
         ELSE
            ngeo_aux(3)=5
         ENDIF
         fact_ngeo_aux(2)=fact_ngeo(2)
         fact_ngeo_aux(3)=fact_ngeo(3)
      CASE(12,13)
         IF (ngeo(2) /= 0) THEN
            ngeo_aux(2)=ngeo(2)
         ELSE
            ngeo_aux(2)=5
         ENDIF
         IF (ngeo(3) /= 0) THEN
            ngeo_aux(3)=ngeo(3)
         ELSE
            ngeo_aux(3)=5
         ENDIF
         IF (ngeo(4) /= 0) THEN
            ngeo_aux(4)=ngeo(4)
         ELSE
            ngeo_aux(4)=5
         ENDIF
         fact_ngeo_aux(2)=fact_ngeo(2)
         fact_ngeo_aux(3)=fact_ngeo(3)
         fact_ngeo_aux(4)=fact_ngeo(4)
      CASE(-12,-13)   
         IF (ngeo(2) /= 0) THEN
            ngeo_aux(2)=ngeo(2)
         ELSE
            ngeo_aux(2)=5
         ENDIF
         IF (ngeo(3) /= 0) THEN
            ngeo_aux(3)=ngeo(3)
         ELSE
            ngeo_aux(3)=5
         ENDIF
         IF (ngeo(5) /= 0) THEN
            ngeo_aux(5)=ngeo(5)
         ELSE
            ngeo_aux(5)=5
         ENDIF
         fact_ngeo_aux(2)=fact_ngeo(2)
         fact_ngeo_aux(3)=fact_ngeo(3)
         fact_ngeo_aux(5)=fact_ngeo(5)
   CASE DEFAULT

!  If the Bravais lattice is unkown, 14 or 0 we let the user choose
!
         ngeo_aux=ngeo
         DO i=2,6
            IF (ngeo_aux(i)==0) ngeo_aux(i)=5
         ENDDO
         fact_ngeo_aux=fact_ngeo
   END SELECT
END IF

ngeo=ngeo_aux
fact_ngeo=fact_ngeo_aux

RETURN
END SUBROUTINE clean_ngeo

SUBROUTINE set_celldm_geo()
!
!   This routine sets the celldm to a grid of points.
!
USE kinds, ONLY : DP
USE thermo_mod, ONLY : what, step_ngeo, ngeo, celldm_geo, reduced_grid, &
                       start_geo, jump_geo
USE control_pwrun, ONLY : celldm_save
IMPLICIT NONE
INTEGER :: igeo1, igeo2, igeo3, igeo4, igeo5, igeo6
INTEGER :: iwork, igeo, jgeo, i, total_work
REAL(DP) :: angle1, angle2, angle3, angle

iwork=0
celldm_geo=0.0_DP
DO igeo6 = 1, ngeo(6)
   angle3 = ACOS(celldm_save(6)) +  &
            (igeo6-(ngeo(6)+1.0_DP)/2.0_DP)*step_ngeo(6)
   DO igeo5 = 1, ngeo(5)
      angle2 = ACOS(celldm_save(5)) +  &
             (igeo5-(ngeo(5)+1.0_DP)/2.0_DP)*step_ngeo(5)
      DO igeo4 = 1, ngeo(4)
         angle1 = ACOS(celldm_save(4)) +  &
              (igeo4-(ngeo(4)+1.0_DP)/2.0_DP)*step_ngeo(4)
         DO igeo3 = 1, ngeo(3)
            DO igeo2 = 1, ngeo(2)
               DO igeo1 = 1, ngeo(1)
                  total_work=total_work+1
                  IF (reduced_grid.AND.(total_work - start_geo) >= 0 .AND. &
                          MOD(total_work - start_geo, jump_geo)/=0) CYCLE
                  iwork=iwork+1
                  celldm_geo(1,iwork)=celldm_save(1)+&
                        (igeo1-(ngeo(1)+1.0_DP)/2.0_DP)*step_ngeo(1)
                  celldm_geo(2,iwork)=celldm_save(2)+&
                        (igeo2-(ngeo(2)+1.0_DP)/2.0_DP)*step_ngeo(2)
                  celldm_geo(3,iwork)=celldm_save(3)+&
                        (igeo3-(ngeo(3)+1.0_DP)/2.0_DP)*step_ngeo(3)
                  celldm_geo(4,iwork)= COS(angle1)
                  celldm_geo(5,iwork)= COS(angle2)
                  celldm_geo(6,iwork)= COS(angle3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

RETURN
END SUBROUTINE set_celldm_geo


INTEGER FUNCTION compute_nwork()

USE thermo_mod, ONLY : ngeo, reduced_grid, central_geo, start_geo, jump_geo
USE control_mur, ONLY : lmurn
IMPLICIT NONE
INTEGER :: iwork, auxgeo

auxgeo=ngeo(1)*ngeo(2)*ngeo(3)*ngeo(4)*ngeo(5)*ngeo(6)
IF (reduced_grid) auxgeo=(auxgeo+1)/jump_geo - (start_geo -1)/jump_geo
compute_nwork=auxgeo
IF (lmurn) compute_nwork=ngeo(1)
central_geo=compute_nwork/2
IF (MOD(compute_nwork,2)==1) central_geo=central_geo+1

RETURN
END FUNCTION compute_nwork

SUBROUTINE initialize_mur(nwork)
USE kinds, ONLY : DP
USE thermo_mod, ONLY : ngeo, ibrav_geo, celldm_geo, energy_geo, omega_geo, &
                       fact_ngeo
USE control_pwrun, ONLY : ibrav_save
USE control_thermo, ONLY : lev_syn_1

IMPLICIT NONE
INTEGER, INTENT(OUT) :: nwork
INTEGER :: igeom
INTEGER :: compute_nwork
REAL(DP) :: compute_omega_geo

   CALL clean_ngeo(ngeo,fact_ngeo,ibrav_save)
   nwork=compute_nwork()
   ALLOCATE(ibrav_geo(nwork))
   ALLOCATE(celldm_geo(6,nwork))
   ALLOCATE(energy_geo(nwork))
   ALLOCATE(omega_geo(nwork))
   CALL set_celldm_geo()
   DO igeom = 1, nwork
      omega_geo(igeom)=compute_omega_geo(ibrav_save,celldm_geo(:,igeom))
   ENDDO
   energy_geo=0.0_DP
   ibrav_geo=ibrav_save
   lev_syn_1=.TRUE.

RETURN
END SUBROUTINE initialize_mur

SUBROUTINE initialize_ph_work(nwork)
USE kinds, ONLY : DP
USE grid_irr_iq, ONLY : irr_iq
USE disp, ONLY : nqs

IMPLICIT NONE
INTEGER, INTENT(OUT) :: nwork

INTEGER :: iq, irr

   nwork=0
   DO iq=1,nqs
      DO irr=0, irr_iq(iq)
         nwork=nwork+1
      ENDDO
   ENDDO

RETURN
END SUBROUTINE initialize_ph_work

SUBROUTINE initialize_no_ph(no_ph, tot_ngeo)
USE thermo_mod, ONLY : ngeo, fact_ngeo, reduced_grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: tot_ngeo
LOGICAL, INTENT(INOUT) :: no_ph(tot_ngeo)
INTEGER :: count_ngeo

LOGICAL :: todo(6)
INTEGER :: igeo1, igeo2, igeo3, igeo4, igeo5, igeo6, i
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

IF (reduced_grid) THEN
   no_ph=.FALSE.
   RETURN
ENDIF

count_ngeo=0
DO igeo6 = 1, ngeo(6)
   todo(6)= (MOD(igeo6-1, fact_ngeo(6))==0)
   DO igeo5 = 1, ngeo(5)
      todo(5)= (MOD(igeo5-1, fact_ngeo(5))==0)
      DO igeo4 = 1, ngeo(4)
         todo(4)= (MOD(igeo4-1, fact_ngeo(4))==0)
         DO igeo3 = 1, ngeo(3)
            todo(3)= (MOD(igeo3-1, fact_ngeo(3))==0)
            DO igeo2 = 1, ngeo(2)
               todo(2)= (MOD(igeo2-1, fact_ngeo(2))==0)
               DO igeo1 = 1, ngeo(1)
                  todo(1)= (MOD(igeo1-1, fact_ngeo(1))==0)
                  count_ngeo = count_ngeo + 1       
                  no_ph(count_ngeo)= .NOT. (todo(1).AND.todo(2).AND.todo(3) &
                                     .AND.todo(4).AND.todo(5).AND.todo(6))
               END DO
            END DO
         END DO
      END DO
   END DO
END DO

RETURN
END SUBROUTINE initialize_no_ph

SUBROUTINE set_equilibrium_conf( celldm, tau, at, omega )
!
!  After the murnaghan equation this routine sets the input variables
!  as if they were the equilibrium crystal parameters. This
!  has to be changed soon, the input variables should not be changed
!  and this variables should have a different name.
!
USE kinds, ONLY : DP
USE ions_base, ONLY : nat
USE control_elastic_constants, ONLY : omega0, tau_save, at_save
USE control_mur, ONLY : celldm0
IMPLICIT NONE

REAL(DP), INTENT(IN) :: celldm(6), tau(3,nat), at(3,3), omega
REAL(DP) :: bg(3,3)

tau_save(:,:) = tau(:,:)
celldm0(:) = celldm(:)
at_save(:,:) = at(:,:)
omega0=omega
!
!  tau_save is given in crystal coordinates
!
CALL recips(at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
CALL cryst_to_cart( nat, tau_save, bg, -1 )

RETURN
END SUBROUTINE set_equilibrium_conf
