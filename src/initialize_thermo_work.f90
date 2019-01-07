!
! Copyright (C) 2013-2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE initialize_thermo_work(nwork, part, iaux)
  !-----------------------------------------------------------------------
  !
  !  This routine is called by all images and initializes the global 
  !  variables that select the tasks to do. All images must have the same
  !  information.
  !  In addition to the variables that control the run it allocates, for
  !  the tasks that require it:
  !  ibrav_geo, celldm_geo, omega_geo, energy_geo
  !  lpwscf, lstress, lberry, lphonon
  !
  USE kinds,          ONLY : DP
  USE thermo_mod,     ONLY : what, step_ngeo, energy_geo, ngeo, &
                             celldm_geo, omega_geo, ibrav_geo, tot_ngeo, no_ph,&
                             start_geometry, last_geometry
  USE control_thermo, ONLY : lpwscf, lphonon, lev_syn_1, lev_syn_2, &
                             lph, lpwscf_syn_1, lbands_syn_1, lq2r,   &
                             ltherm, lconv_ke_test, lconv_nk_test, &
                             lstress, lelastic_const, lpiezoelectric_tensor,&
                             lberry, lpolarization, lpart2_pw, do_scf_relax, &
                             ldos_syn_1, ltherm_dos, ltherm_freq, after_disp
  USE control_pwrun,  ONLY : do_punch
  USE control_conv,   ONLY : nke, ke, deltake, nkeden, deltakeden, keden, &
                             nnk, nk_test, deltank, nsigma, sigma_test, &  
                             deltasigma, ncutoffene
  USE equilibrium_conf, ONLY : celldm0, omega0
  USE initial_conf,   ONLY : celldm_save, ibrav_save
  USE piezoelectric_tensor, ONLY : polar_geo
  USE control_elastic_constants, ONLY : elastic_algorithm, rot_mat
  USE elastic_constants, ONLY : epsilon_voigt, sigma_geo, epsilon_geo
  USE gvecw,          ONLY : ecutwfc
  USE gvect,          ONLY : ecutrho
  USE control_quadratic_energy, ONLY : nvar, degree
  USE lattices,       ONLY : crystal_parameters
  USE quadratic_surfaces, ONLY : quadratic_var
  USE start_k,        ONLY : nk1, nk2, nk3
  USE klist,          ONLY : degauss
  USE wrappers,       ONLY : f_mkdir_safe
  USE io_global,      ONLY : meta_ionode
  !
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: nwork, iaux
  INTEGER, INTENT(IN) :: part

  INTEGER :: igeom, ike, iden, icount, ink, isigma, ios
  INTEGER :: count_energies
  REAL(DP) :: compute_omega_geo, dual, kev, kedenv
  !
!
!   the restart directory is used in all cases
!
  ios=0
  IF (meta_ionode) ios = f_mkdir_safe( 'restart' )
!
!  here initialize the work to do and sets to true the flags that activate
!  the different parts of thermo_pw. Sets also the number of tasks for each
!  what.
!
  nwork=0
  iaux=0
  IF (part == 1) THEN
     SELECT CASE (TRIM(what))
!
!   In these cases we do not do any asynchronous work in the first part
!
        CASE ( 'plot_bz', ' ') 
        CASE ( 'scf') 
           ALLOCATE(energy_geo(1))
           lpwscf_syn_1=.TRUE.
        CASE ('scf_bands') 
           ALLOCATE(energy_geo(1))
           lpwscf_syn_1=.TRUE.
           lbands_syn_1=.TRUE.
           IF (meta_ionode) ios = f_mkdir_safe( 'band_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
        CASE ('scf_2d_bands')
           ALLOCATE(energy_geo(1))
           lpwscf_syn_1=.TRUE.
           lbands_syn_1=.TRUE.
           IF (meta_ionode) ios = f_mkdir_safe( 'band_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
        CASE ('scf_dos') 
           ALLOCATE(energy_geo(1))
           lpwscf_syn_1=.TRUE.
           lbands_syn_1=.TRUE.
           ldos_syn_1=.TRUE.
           IF (meta_ionode) ios = f_mkdir_safe( 'band_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'therm_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
        CASE ('scf_ph') 
           ALLOCATE(energy_geo(1))
           lpwscf_syn_1=.NOT.after_disp
           lph=.TRUE.
           tot_ngeo=1
           ALLOCATE(no_ph(tot_ngeo))
           no_ph(1)=.FALSE.
           IF (meta_ionode) ios = f_mkdir_safe( 'dynamical_matrices' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
        CASE ('scf_disp')
           ALLOCATE(energy_geo(1))
           lpwscf_syn_1=.NOT.after_disp
           lph=.TRUE.
           tot_ngeo=1
           ALLOCATE(no_ph(tot_ngeo))
           no_ph(1)=.FALSE.
           lq2r = .TRUE.
           ltherm = ltherm_dos .OR. ltherm_freq
           CALL allocate_thermodynamics()
!
!   In these cases we make asynchronous work in the first part
!
           IF (meta_ionode) ios = f_mkdir_safe( 'dynamical_matrices' )
           IF (meta_ionode) ios = f_mkdir_safe( 'phdisp_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'therm_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
        CASE ( 'scf_ke') 
           nwork= count_energies(ecutwfc, ecutrho, deltake, deltakeden, nke,&
                                                                      nkeden)
           ncutoffene=nwork
           ALLOCATE(ke(nwork))
           ALLOCATE(keden(nwork))
           ALLOCATE(energy_geo(nwork))
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
           IF (meta_ionode) ios = f_mkdir_safe( 'therm_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'elastic_constants' )
        CASE ('scf_piezoelectric_tensor')
           lpart2_pw=.TRUE.
           tot_ngeo=1
        CASE ('scf_polarization') 
           lpart2_pw=.TRUE.
           tot_ngeo=1
!
!   here all the cases that require the determination of the minimization
!   of the energy to find the equilibrium crystal parameters
!
        CASE ('mur_lc')
           lev_syn_1=.TRUE.
           lpwscf_syn_1=do_scf_relax
           do_punch=.FALSE.
           CALL initialize_mur(nwork)
           IF (meta_ionode) ios = f_mkdir_safe( 'energy_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
        CASE ('mur_lc_bands') 
           do_punch=.FALSE.
           lev_syn_1=.TRUE.
           lpwscf_syn_1=.TRUE.
           lbands_syn_1=.TRUE.
           CALL initialize_mur(nwork)
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
           IF (meta_ionode) ios = f_mkdir_safe( 'energy_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'therm_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
        CASE ('mur_lc_ph') 
           do_punch=.FALSE.
           lpwscf_syn_1=.NOT.after_disp
           lev_syn_1=.TRUE.
           lph=.TRUE.
           CALL initialize_mur(nwork)
           tot_ngeo=1
           ALLOCATE(no_ph(tot_ngeo))
           no_ph(1)=.FALSE.
           IF (meta_ionode) ios = f_mkdir_safe( 'energy_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'dynamical_matrices' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
        CASE ('mur_lc_disp')
           do_punch=.FALSE.
           lev_syn_1=.TRUE.
           lpwscf_syn_1=.NOT.after_disp
           lph=.TRUE.
           CALL initialize_mur(nwork)
           tot_ngeo=1
           ALLOCATE(no_ph(tot_ngeo))
           no_ph(1)=.FALSE.
           lq2r = .TRUE.
           ltherm = ltherm_dos .OR. ltherm_freq
           CALL allocate_thermodynamics()
           IF (meta_ionode) ios = f_mkdir_safe( 'energy_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'dynamical_matrices' )
           IF (meta_ionode) ios = f_mkdir_safe( 'phdisp_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'therm_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
        CASE ('mur_lc_elastic_constants') 
           lev_syn_1=.TRUE.
           lpwscf_syn_1=do_scf_relax
           lpart2_pw=.TRUE.
           do_punch=.FALSE.
           CALL initialize_mur(nwork)
           tot_ngeo=1
           IF (meta_ionode) ios = f_mkdir_safe( 'therm_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'energy_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'elastic_constants' )
        CASE ('mur_lc_piezoelectric_tensor') 
           lev_syn_1=.TRUE.
           lpwscf_syn_1=do_scf_relax
           lpart2_pw=.TRUE.
           do_punch=.FALSE.
           CALL initialize_mur(nwork)
           tot_ngeo=1
           IF (meta_ionode) ios = f_mkdir_safe( 'energy_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
        CASE ('mur_lc_polarization')
           lev_syn_1=.TRUE.
           lpwscf_syn_1=do_scf_relax
           lpart2_pw=.TRUE.
           do_punch=.FALSE.
           CALL initialize_mur(nwork)
           tot_ngeo=1
           IF (meta_ionode) ios = f_mkdir_safe( 'energy_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
!
!    Here all the cases that compute the free energy and minimize it
!
        CASE ('mur_lc_t')
           lev_syn_1=.TRUE.
           lph = .TRUE.
           lq2r = .TRUE.
           ltherm = .TRUE.
           lev_syn_2=.TRUE.
           CALL initialize_mur(nwork)
           tot_ngeo=nwork
           ALLOCATE(no_ph(tot_ngeo))
           CALL initialize_no_ph(no_ph, tot_ngeo)
           degree=crystal_parameters(ibrav_save)
           nvar=quadratic_var(degree)
           CALL allocate_thermodynamics()
           CALL allocate_anharmonic()
           IF (meta_ionode) ios = f_mkdir_safe( 'energy_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'anhar_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'therm_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'dynamical_matrices' )
           IF (meta_ionode) ios = f_mkdir_safe( 'phdisp_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'elastic_constants' )
        CASE ('elastic_constants_t')
           lev_syn_1=.TRUE.
           lpart2_pw=.TRUE.
           CALL initialize_mur(nwork)
           tot_ngeo=nwork
           CALL init_elastic_constants_t()
           IF (meta_ionode) ios = f_mkdir_safe( 'therm_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'gnuplot_files' )
           IF (meta_ionode) ios = f_mkdir_safe( 'elastic_constants' )
           nwork=0
           lev_syn_1=.FALSE.
        CASE DEFAULT
           CALL errore('initialize_thermo_work','what not recognized',1)
     END SELECT
     IF (start_geometry < 1) start_geometry=1
     IF (last_geometry > tot_ngeo) last_geometry=tot_ngeo
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
        CASE ('scf_elastic_constants', 'mur_lc_elastic_constants',&
                                       'elastic_constants_t')

           IF (ALLOCATED(ibrav_geo))  DEALLOCATE(ibrav_geo)
           IF (ALLOCATED(celldm_geo)) DEALLOCATE(celldm_geo)
           IF (ALLOCATED(energy_geo)) DEALLOCATE(energy_geo)
           IF (ALLOCATED(omega_geo))  DEALLOCATE(omega_geo)
           IF (ALLOCATED(epsilon_voigt)) DEALLOCATE(epsilon_voigt)
           IF (ALLOCATED(epsilon_geo)) DEALLOCATE(epsilon_geo)
           IF (ALLOCATED(sigma_geo))  DEALLOCATE(sigma_geo)
           IF (ALLOCATED(rot_mat))    DEALLOCATE(rot_mat)
!
!     set_elastic_constant_work allocates ibrav_geo and celldm_geo
!
           CALL set_elastic_cons_work(nwork)

           ALLOCATE(energy_geo(nwork))
           ALLOCATE(omega_geo(nwork))
           energy_geo=0.0_DP
           IF (elastic_algorithm=='energy') THEN
              DO igeom = 1, nwork
                 omega_geo(igeom)=compute_omega_geo(ibrav_geo(igeom),&
                                                    celldm_geo(1,igeom))
              ENDDO
              omega0 = compute_omega_geo(ibrav_save, celldm0)
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
           iaux=10
           lpiezoelectric_tensor=.TRUE.
           do_punch=.TRUE.
        CASE ('scf_polarization', 'mur_lc_polarization')
           IF (ALLOCATED(energy_geo)) DEALLOCATE(energy_geo)
           IF (ALLOCATED(omega_geo)) DEALLOCATE(omega_geo)
           iaux=20 
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
  ENDIF

  IF ( nwork == 0 ) RETURN

  ALLOCATE( lpwscf(nwork) )
  ALLOCATE( lstress(nwork) )
  ALLOCATE( lberry(nwork) )
  ALLOCATE( lphonon(nwork) )
  lpwscf  = .FALSE.
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
               'scf_dos',                    &
               'scf_2d_bands',               &
               'scf_ph',                     & 
               'scf_disp',                   &
               'scf_elastic_constants',      &
               'elastic_constants_t',        &
               'scf_piezoelectric_tensor',   &
               'scf_polarization' )
!
!  Then the cases in which there is pwscf asynchronous work in the first part
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
        CASE ('scf_elastic_constants', 'mur_lc_elastic_constants', &
                                       'elastic_constants_t')
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
  ENDIF

  RETURN
END SUBROUTINE initialize_thermo_work

SUBROUTINE clean_ngeo(ngeo,fact_ngeo,ibrav)
!
!  This routine cleans the ngeo variable, setting 1 on all the values
!  that correspond to crystallographic parameters that are fixed in
!  a given ibrav. Moreover if ngeo is zero for some parameter needed
!  by ibrav, the routine sets a default value of ndefault, presently set
!  to 5
!
USE control_mur, ONLY : lmurn
IMPLICIT NONE

INTEGER, INTENT(IN)    :: ibrav
INTEGER, INTENT(INOUT) :: ngeo(6), fact_ngeo(6)

INTEGER :: i, ndefault, ngeo_aux(6), fact_ngeo_aux(6)

ngeo_aux=1
ngeo_aux(1)=ngeo(1)
fact_ngeo_aux=1
fact_ngeo_aux(1)=fact_ngeo(1)
ndefault=5

IF (.NOT.lmurn) THEN
   SELECT CASE (ibrav)
      CASE(1,2,3)
      CASE (4,6,7)
         IF (ngeo(3) /= 0) THEN
             ngeo_aux(3)=ngeo(3)
         ELSE
             ngeo_aux(3)=ndefault
         ENDIF
         fact_ngeo_aux(3)=fact_ngeo(3)
      CASE (5)
         IF (ngeo(4) /= 0) THEN
            ngeo_aux(4)=ngeo(4)
         ELSE
            ngeo_aux(4)=ndefault
         ENDIF
         fact_ngeo_aux(4)=fact_ngeo(4)
      CASE(8,9,-9,91,10,11)
         IF (ngeo(2) /= 0) THEN
            ngeo_aux(2)=ngeo(2)
         ELSE
            ngeo_aux(2)=ndefault
         ENDIF
         IF (ngeo(3) /= 0) THEN
            ngeo_aux(3)=ngeo(3)
         ELSE
            ngeo_aux(3)=ndefault
         ENDIF
         fact_ngeo_aux(2)=fact_ngeo(2)
         fact_ngeo_aux(3)=fact_ngeo(3)
      CASE(12,13)
         IF (ngeo(2) /= 0) THEN
            ngeo_aux(2)=ngeo(2)
         ELSE
            ngeo_aux(2)=ndefault
         ENDIF
         IF (ngeo(3) /= 0) THEN
            ngeo_aux(3)=ngeo(3)
         ELSE
            ngeo_aux(3)=ndefault
         ENDIF
         IF (ngeo(4) /= 0) THEN
            ngeo_aux(4)=ngeo(4)
         ELSE
            ngeo_aux(4)=ndefault
         ENDIF
         fact_ngeo_aux(2)=fact_ngeo(2)
         fact_ngeo_aux(3)=fact_ngeo(3)
         fact_ngeo_aux(4)=fact_ngeo(4)
      CASE(-12,-13)   
         IF (ngeo(2) /= 0) THEN
            ngeo_aux(2)=ngeo(2)
         ELSE
            ngeo_aux(2)=ndefault
         ENDIF
         IF (ngeo(3) /= 0) THEN
            ngeo_aux(3)=ngeo(3)
         ELSE
            ngeo_aux(3)=ndefault
         ENDIF
         IF (ngeo(5) /= 0) THEN
            ngeo_aux(5)=ngeo(5)
         ELSE
            ngeo_aux(5)=ndefault
         ENDIF
         fact_ngeo_aux(2)=fact_ngeo(2)
         fact_ngeo_aux(3)=fact_ngeo(3)
         fact_ngeo_aux(5)=fact_ngeo(5)
   CASE DEFAULT

!  If the Bravais lattice is unkown, 14 or 0 we let the user choose
!
         ngeo_aux=ngeo
         DO i=2,6
            IF (ngeo_aux(i)==0) ngeo_aux(i)=ndefault
         ENDDO
         fact_ngeo_aux=fact_ngeo
   END SELECT
ENDIF

ngeo=ngeo_aux
fact_ngeo=fact_ngeo_aux

RETURN
END SUBROUTINE clean_ngeo

SUBROUTINE set_celldm_geo()
!
!   This routine sets the grid of values on celldm_geo. 
!   It modifies only celldm_geo.
!
USE kinds,         ONLY : DP
USE constants,     ONLY : pi
USE thermo_mod,    ONLY : step_ngeo, ngeo, celldm_geo, reduced_grid
USE initial_conf,  ONLY : celldm_save

IMPLICIT NONE
INTEGER  :: igeo1, igeo2, igeo3, igeo4, igeo5, igeo6
INTEGER  :: iwork, total_work
REAL(DP) :: angle1, angle2, angle3

iwork=0
celldm_geo=0.0_DP
total_work=0
DO igeo6 = 1, ngeo(6)
   angle3 = ACOS(celldm_save(6)) +  &
            (igeo6-(ngeo(6)+1.0_DP)/2.0_DP)*step_ngeo(6)*pi/180.0_DP
   DO igeo5 = 1, ngeo(5)
      angle2 = ACOS(celldm_save(5)) +  &
             (igeo5-(ngeo(5)+1.0_DP)/2.0_DP)*step_ngeo(5)*pi/180.0_DP
      DO igeo4 = 1, ngeo(4)
         angle1 = ACOS(celldm_save(4)) +  &
              (igeo4-(ngeo(4)+1.0_DP)/2.0_DP)*step_ngeo(4)*pi/180.0_DP
         DO igeo3 = 1, ngeo(3)
            DO igeo2 = 1, ngeo(2)
               DO igeo1 = 1, ngeo(1)
                  total_work=total_work+1
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
!
!  This function computes the number of tasks needed for energy minimization
!
USE thermo_mod, ONLY : ngeo, reduced_grid, central_geo
USE control_mur, ONLY : lmurn

IMPLICIT NONE

INTEGER :: iwork, auxgeo

auxgeo=ngeo(1)*ngeo(2)*ngeo(3)*ngeo(4)*ngeo(5)*ngeo(6)
IF (lmurn) auxgeo=ngeo(1)

central_geo=auxgeo/2
IF (MOD(auxgeo,2)==1) central_geo=central_geo+1

compute_nwork=auxgeo

RETURN
END FUNCTION compute_nwork

SUBROUTINE initialize_mur(nwork)
!
!   this routine sets the number of tasks needed to make an energy
!   minimization. It allocates the variables ibrav_geo, celldm_geo, 
!   energy_geo, and omega_geo. Moreover it sets the values of celldm_geo,
!   in each task.
!
USE kinds,         ONLY : DP
USE thermo_mod,    ONLY : ibrav_geo, celldm_geo, energy_geo, omega_geo
USE initial_conf,  ONLY : ibrav_save

IMPLICIT NONE

INTEGER, INTENT(OUT) :: nwork

INTEGER              :: igeom
INTEGER              :: compute_nwork
REAL(DP)             :: compute_omega_geo

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

RETURN
END SUBROUTINE initialize_mur

SUBROUTINE initialize_ph_work(nwork)
!
!  This routine counts how many tasks there are in a phonon calculation.
!
USE control_thermo, ONLY : all_geometries_together
USE thermo_mod,  ONLY : tot_ngeo, start_geometry, last_geometry
USE ions_base,   ONLY : nat
USE grid_irr_iq, ONLY : irr_iq
USE disp,        ONLY : nqs
USE freq_ph,     ONLY : nfs, fpol
USE images_omega,ONLY : omega_group
USE initial_conf, ONLY : geometry, collect_info_save
USE control_ph,  ONLY : epsil, trans
USE control_qe,  ONLY : use_ph_images
USE mp_images,   ONLY : nimage
USE mp_asyn,     ONLY : with_asyn_images

IMPLICIT NONE
INTEGER, INTENT(OUT) :: nwork

INTEGER :: iq, irr, igeom, iwork, nqs_max, ngeometries

nwork=0
IF (trans) THEN
   IF (use_ph_images) THEN
      IF (all_geometries_together) THEN
         nwork=nimage*(last_geometry - start_geometry + 1)
         ALLOCATE(geometry(nwork))
         DO iwork=1, nwork
            geometry(iwork)= (iwork-1) / nimage + start_geometry
         ENDDO
      ELSE
         nwork=nimage
      ENDIF
   ELSE
      IF (all_geometries_together) THEN
         nqs_max=0
         DO igeom=start_geometry, last_geometry
            nqs_max=MAX(collect_info_save(igeom)%nqs, nqs_max)
         ENDDO
         ngeometries=last_geometry - start_geometry + 1
         ALLOCATE(geometry((3*nat+1)*nqs_max*ngeometries))
         nwork=0
         DO igeom=start_geometry, last_geometry
            DO iq=1, collect_info_save(igeom)%nqs
               DO irr=0, collect_info_save(igeom)%irr_iq(iq)
                  nwork=nwork+1
                  geometry(nwork)=igeom
               ENDDO
            ENDDO
         ENDDO
      ELSE
         DO iq=1,nqs
            DO irr=0, irr_iq(iq)
               nwork=nwork+1
            ENDDO
         ENDDO
      ENDIF
   ENDIF
ELSEIF (fpol) THEN
   IF (with_asyn_images) THEN
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

SUBROUTINE initialize_no_ph(no_ph, tot_ngeo)
USE thermo_mod, ONLY : ngeo, fact_ngeo, reduced_grid
!
!   This routine is used to skip phonon calculations in those geometries
!   that are not on the grid used to compute the vibrational energy.
!
IMPLICIT NONE

INTEGER, INTENT(IN)    :: tot_ngeo
LOGICAL, INTENT(INOUT) :: no_ph(tot_ngeo)

LOGICAL :: todo(6)
INTEGER :: igeo1, igeo2, igeo3, igeo4, igeo5, igeo6, count_ngeo
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
no_ph=.FALSE.
IF (reduced_grid) RETURN

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
!  This routine sets the equilibrium variables to those given
!  in input. The tau are in cartesian coordinates and this routine
!  computes the reciprocal lattice vectors and the atomic
!  coordinates in crystal basis.
!
USE kinds,            ONLY : DP
USE ions_base,        ONLY : nat
USE equilibrium_conf, ONLY : celldm0, omega0, tau0, tau0_crys, at0, bg0

IMPLICIT NONE

REAL(DP), INTENT(IN) :: celldm(6), tau(3,nat), at(3,3), omega

celldm0(:) = celldm(:)
omega0=omega
at0(:,:) = at(:,:)
!
! compute the reciprocal lattice vectors
!
CALL recips(at0(1,1),at0(1,2),at0(1,3),bg0(1,1),bg0(1,2),bg0(1,3))
!
! set the equlibrium tau and computes them in crystal coordinates
!
tau0(:,:) = tau(:,:)
tau0_crys=tau0
CALL cryst_to_cart( nat, tau0_crys, bg0, -1 )

RETURN
END SUBROUTINE set_equilibrium_conf

INTEGER FUNCTION count_energies(ecutwfc, ecutrho, deltake, deltakeden, &
                                                             nke, nkeden)
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

SUBROUTINE find_central_geo(ngeo, no_ph, central_geo)

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
