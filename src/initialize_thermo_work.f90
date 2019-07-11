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
                             start_geometry, last_geometry, reduced_grid, &
                             in_degree
  USE control_thermo, ONLY : lpwscf, lphonon, lev_syn_1, lev_syn_2, &
                             lph, lpwscf_syn_1, lbands_syn_1, lq2r,   &
                             ltherm, lconv_ke_test, lconv_nk_test, &
                             lstress, lelastic_const, lpiezoelectric_tensor,&
                             lberry, lpolarization, lpart2_pw, do_scf_relax, &
                             ldos_syn_1, ltherm_dos, ltherm_freq, after_disp, &
                             lectqha
  USE control_pwrun,  ONLY : do_punch
  USE control_conv,   ONLY : nke, ke, deltake, nkeden, deltakeden, keden, &
                             nnk, nk_test, deltank, nsigma, sigma_test, &  
                             deltasigma, ncutoffene
  USE equilibrium_conf, ONLY : celldm0, omega0
  USE initial_conf,   ONLY : celldm_save, ibrav_save
  USE piezoelectric_tensor, ONLY : polar_geo
  USE control_elastic_constants, ONLY : elastic_algorithm, rot_mat
  USE control_elastic_constants_qha, ONLY : ngeom, use_free_energy
  USE elastic_constants, ONLY : epsilon_voigt, sigma_geo, epsilon_geo
  USE gvecw,          ONLY : ecutwfc
  USE gvect,          ONLY : ecutrho
  USE control_quadratic_energy, ONLY : ncoeff, nvar
  USE lattices,       ONLY : crystal_parameters
  USE quadratic_surfaces, ONLY : quadratic_ncoeff
  USE start_k,        ONLY : nk1, nk2, nk3
  USE klist,          ONLY : degauss
  USE wrappers,       ONLY : f_mkdir_safe
  USE io_global,      ONLY : meta_ionode
  !
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: nwork, iaux
  INTEGER, INTENT(IN) :: part

  INTEGER :: igeom, ike, iden, icount, ink, isigma, iwork, ios
  INTEGER :: count_energies
  REAL(DP) :: compute_omega_geo, dual, kev, kedenv
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
           lev_syn_1=.NOT.reduced_grid
           lph = .TRUE.
           lq2r = .TRUE.
           ltherm = .TRUE.
           lev_syn_2=.TRUE.
           CALL initialize_mur(nwork)
           tot_ngeo=nwork
           ALLOCATE(no_ph(tot_ngeo))
           IF (reduced_grid) ALLOCATE(in_degree(tot_ngeo))
           CALL initialize_no_ph(no_ph, tot_ngeo, ibrav_save)
           CALL summarize_geometries(nwork)
           nvar=crystal_parameters(ibrav_save)
           ncoeff=quadratic_ncoeff(nvar)
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
           lectqha=.TRUE.
           lph=use_free_energy
           lq2r = use_free_energy
           ltherm = use_free_energy
           CALL initialize_mur_qha(ngeom)
           CALL set_elastic_cons_work(ngeom, nwork)
           tot_ngeo=nwork
           ALLOCATE(energy_geo(tot_ngeo))
           ALLOCATE(no_ph(tot_ngeo))
           energy_geo=0.0_DP
           no_ph=.FALSE.
           CALL summarize_geometries(nwork)
           IF (use_free_energy) THEN
              CALL allocate_thermodynamics()
              CALL allocate_anharmonic()
           ELSE
              no_ph=.TRUE.
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
        CASE ('elastic_constants_t')
           IF (use_free_energy) THEN
              CALL initialize_ph_work(nwork)
           ELSE
              nwork=0
           ENDIF
        CASE ('scf_elastic_constants', 'mur_lc_elastic_constants')

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
           CALL set_elastic_cons_work(1,nwork)

           ALLOCATE(energy_geo(nwork))
           ALLOCATE(omega_geo(nwork))
           energy_geo=0.0_DP
           IF (elastic_algorithm=='energy'.OR. &
                                   elastic_algorithm=='energy_std') THEN
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
           IF (reduced_grid) THEN
              DO iwork=1,nwork
                 IF (no_ph(iwork)) lpwscf(iwork)=.FALSE.
              ENDDO
           ENDIF
        CASE ('elastic_constants_t')  
           lpwscf(1:nwork)=.TRUE.
           lstress(1:nwork)=( elastic_algorithm=='standard' .OR. &
                              elastic_algorithm=='advanced')
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
              'mur_lc_t',         &
              'elastic_constants_t')

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
  ENDIF

  RETURN
END SUBROUTINE initialize_thermo_work

SUBROUTINE clean_ngeo(ngeo,fact_ngeo,ngeo_ph,ibrav)
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
INTEGER, INTENT(INOUT) :: ngeo(6), fact_ngeo(6), ngeo_ph(6)

INTEGER :: i, ndefault, ngeo_aux(6), fact_ngeo_aux(6)
LOGICAL :: clean_fact

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
!
!  if ngeo_ph has been set, check if it is compatible with ngeo and
!  clean fact_ngeo which is not used.
!
clean_fact=.FALSE.
DO i=1,6
   IF (ngeo_ph(i)<=0) ngeo_ph(i)=ngeo(i)
   IF (ngeo_ph(i)>ngeo(i)) ngeo_ph(i)=ngeo(i)
   IF (MOD(ngeo_ph(i),2) /= MOD(ngeo(i),2)) &
               CALL errore('clean_ngeo','ngeo_ph incompatible with ngeo',1)
   IF (ngeo_ph(i)/=ngeo(i)) clean_fact=.TRUE.
ENDDO
IF (clean_fact) fact_ngeo=1

RETURN
END SUBROUTINE clean_ngeo

SUBROUTINE set_celldm_geo(celldm_geo, nwork)
!
!   This routine sets the grid of values on celldm_geo.
!
USE kinds,         ONLY : DP
USE constants,     ONLY : pi
USE thermo_mod,    ONLY : step_ngeo, ngeo, reduced_grid
USE initial_conf,  ONLY : celldm_save

IMPLICIT NONE
INTEGER, INTENT(IN) :: nwork
REAL(DP), INTENT(INOUT) :: celldm_geo(6,nwork)

INTEGER  :: igeo1, igeo2, igeo3, igeo4, igeo5, igeo6
INTEGER  :: iwork, i, total_work
REAL(DP) :: angle1, angle2, angle3, delta(6)

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
                  IF (total_work>nwork) &
                     CALL errore('set_celldm_geo','nwork too low', 1)
                  iwork=iwork+1
                  celldm_geo(1,iwork)=celldm_save(1)+igeo1*step_ngeo(1)
                  celldm_geo(2,iwork)=celldm_save(2)+igeo2*step_ngeo(2)
                  celldm_geo(3,iwork)=celldm_save(3)+igeo3*step_ngeo(3)
                  IF (reduced_grid.AND.(ABS(igeo1)+ABS(igeo2)+ABS(igeo3)+  &
                            ABS(igeo4)+ABS(igeo5)+ABS(igeo6)==1)) THEN
                     IF (igeo1/=0) celldm_geo(1,iwork)=celldm_geo(1,iwork) &
                                                                  +delta(1)
                     IF (igeo2/=0) celldm_geo(2,iwork)=celldm_geo(2,iwork) &
                                                                  +delta(2)
                     IF (igeo3/=0) celldm_geo(3,iwork)=celldm_geo(3,iwork) &
                                                                  +delta(3)
                     IF (igeo4/=0) angle1=angle1+delta(4)*pi/180.0_DP
                     IF (igeo5/=0) angle2=angle2+delta(5)*pi/180.0_DP
                     IF (igeo6/=0) angle3=angle3+delta(6)*pi/180.0_DP
                  ELSEIF(.NOT.reduced_grid) THEN
                     celldm_geo(1,iwork)=celldm_geo(1,iwork)+delta(1)
                     celldm_geo(2,iwork)=celldm_geo(2,iwork)+delta(2)
                     celldm_geo(3,iwork)=celldm_geo(3,iwork)+delta(3)
                     angle1=angle1+delta(4)*pi/180.0_DP
                     angle2=angle2+delta(5)*pi/180.0_DP
                     angle3=angle3+delta(6)*pi/180.0_DP
                  ENDIF
                  celldm_geo(4,iwork)=COS(angle1)
                  celldm_geo(5,iwork)=COS(angle2)
                  celldm_geo(6,iwork)=COS(angle3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

RETURN
END SUBROUTINE set_celldm_geo

SUBROUTINE summarize_geometries(nwork)
USE thermo_mod,    ONLY : ibrav_geo, celldm_geo, no_ph
USE initial_conf,  ONLY : celldm_save, ibrav_save
USE io_global,     ONLY : stdout

IMPLICIT NONE
INTEGER :: nwork
INTEGER :: igeo, phcount

phcount=0
DO igeo=1,nwork
   IF(.NOT.no_ph(igeo)) phcount=phcount+1
ENDDO

WRITE(stdout,'(/,5x,70("-"))')
WRITE(stdout,'(/,5x,"Number of geometries to compute:", i5)') nwork
WRITE(stdout,'(5x,"Phonons computed in:", 12x, i5, " geometries")') phcount
WRITE(stdout,'(/,5x,"Input ibrav and celldm:")')
WRITE(stdout,'(12x, i3, 6f10.5,/)') ibrav_save, celldm_save(:)

DO igeo=1,nwork
   WRITE(stdout,'(5x,i5,": ", i3,6f10.5,l2)') igeo, ibrav_geo(igeo), &
                                   celldm_geo(:,igeo), .NOT.no_ph(igeo)
ENDDO
WRITE(stdout,'(/,5x,70("-"))')

RETURN
END SUBROUTINE summarize_geometries

INTEGER FUNCTION compute_nwork()
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
CALL set_celldm_geo(celldm_geo, nwork)
DO igeom = 1, nwork
   omega_geo(igeom)=compute_omega_geo(ibrav_save,celldm_geo(:,igeom))
ENDDO
energy_geo=0.0_DP
ibrav_geo=ibrav_save

RETURN
END SUBROUTINE initialize_mur

SUBROUTINE initialize_mur_qha(nwork)
!
!   this routine sets the unperturbed geometries for the computation
!   of the elastic constants within the quasiharmonic approximation.
!   It allocates the variables ibrav_save_qha and celldm0_qha. 
!
USE kinds,         ONLY : DP
USE ions_base,     ONLY : nat
USE control_elastic_constants_qha, ONLY : ibrav_save_qha, celldm0_qha, &
                                          tau_crys_qha, omega0_qha
USE initial_conf,  ONLY : ibrav_save, tau_save_crys
USE io_global, ONLY : stdout

IMPLICIT NONE

INTEGER, INTENT(OUT) :: nwork

INTEGER              :: igeom, iwork
INTEGER              :: compute_nwork
REAL(DP)             :: compute_omega_geo

nwork=compute_nwork()
ALLOCATE(ibrav_save_qha(nwork))
ALLOCATE(celldm0_qha(6,nwork))
ALLOCATE(omega0_qha(nwork))
ALLOCATE(tau_crys_qha(3,nat,nwork))
celldm0_qha=0.0_DP
CALL set_celldm_geo(celldm0_qha(:,:), nwork)
ibrav_save_qha(:)=ibrav_save
DO iwork=1,nwork
   tau_crys_qha(:,:,iwork)=tau_save_crys(:,:)
   omega0_qha(iwork)=compute_omega_geo(ibrav_save, celldm0_qha(:,iwork))
ENDDO

WRITE(stdout,'(/,5x,70("-"))')
WRITE(stdout,'(/,5x,"The celldm of the",i5," unperturbed geometries &
                                              &are:",/)') nwork
DO iwork=1,nwork
   WRITE(stdout,'(5x,i5,": ", 6f10.5)') iwork, &
                                   celldm0_qha(1:6,iwork)
ENDDO
WRITE(stdout,'(/,5x,70("-"))')

RETURN
END SUBROUTINE initialize_mur_qha

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

SUBROUTINE initialize_no_ph(no_ph, tot_ngeo, ibrav)
USE thermo_mod, ONLY : ngeo, fact_ngeo, ngeo_ph, reduced_grid, &
                       red_central_geo, in_degree
USE lattices, ONLY : compress_int_vect

!
!   This routine is used to skip phonon calculations in those geometries
!   that are not on the grid used to compute the vibrational energy.
!
IMPLICIT NONE

INTEGER, INTENT(IN)    :: tot_ngeo, ibrav
LOGICAL, INTENT(INOUT) :: no_ph(tot_ngeo)

LOGICAL :: todo(6), in_grid
INTEGER :: igeo1, igeo2, igeo3, igeo4, igeo5, igeo6, start, count_ngeo, &
           aux_deg(tot_ngeo)

IF (reduced_grid) THEN
   count_ngeo=0
   no_ph=.TRUE.
   aux_deg=0
   DO igeo6 = -ngeo(6)/2, (ngeo(6)-1)/2
      DO igeo5 = -ngeo(5)/2, (ngeo(5)-1)/2
         DO igeo4 = -ngeo(4)/2, (ngeo(4)-1)/2
            DO igeo3 = -ngeo(3)/2, (ngeo(3)-1)/2
               DO igeo2 = -ngeo(2)/2, (ngeo(2)-1)/2
                  DO igeo1 = -ngeo(1)/2, (ngeo(1)-1)/2
                     count_ngeo=count_ngeo+1
                     IF (in_grid(igeo1,igeo2,igeo3,igeo4,igeo5,igeo6)) THEN
                        no_ph(count_ngeo)=.FALSE. 
                        IF (igeo1/=0) aux_deg(count_ngeo)=1
                        IF (igeo2/=0) aux_deg(count_ngeo)=2
                        IF (igeo3/=0) aux_deg(count_ngeo)=3
                        IF (igeo4/=0) aux_deg(count_ngeo)=4
                        IF (igeo5/=0) aux_deg(count_ngeo)=5
                        IF (igeo6/=0) aux_deg(count_ngeo)=6
                     ELSEIF ((ABS(igeo1)+ABS(igeo2)+ABS(igeo3)+ABS(igeo4)+ &
                                     ABS(igeo5)+ABS(igeo6))==0) THEN
                        no_ph(count_ngeo)=.FALSE. 
                        red_central_geo=count_ngeo
                     ENDIF 
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   CALL compress_int_vect(aux_deg,in_degree,tot_ngeo,ibrav)
   RETURN
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
no_ph=.FALSE.
!
!  Here we set no_ph based on ngeo_ph. Note that fact_ngeo should all
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
!
USE thermo_mod, ONLY : reduced_grid, red_central_geo
!
IMPLICIT NONE
INTEGER :: ngeo
LOGICAL :: no_ph(ngeo)
INTEGER :: central_geo

INTEGER :: igeo

IF (reduced_grid) THEN
   central_geo=red_central_geo
   RETURN
ENDIF

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

FUNCTION in_grid(igeo1, igeo2, igeo3, igeo4, igeo5, igeo6)
!
!  This logical function returns .TRUE. if only one of the six igeo 
!  is different from zero. These are the geometries at which the phonon
!  dispersions are calculated when reduced_grid=.TRUE.
!
IMPLICIT NONE
LOGICAL :: in_grid
INTEGER :: igeo1, igeo2, igeo3, igeo4, igeo5, igeo6, nzero

nzero=0
IF (igeo1 /= 0) nzero=nzero+1
IF (igeo2 /= 0) nzero=nzero+1
IF (igeo3 /= 0) nzero=nzero+1
IF (igeo4 /= 0) nzero=nzero+1
IF (igeo5 /= 0) nzero=nzero+1
IF (igeo6 /= 0) nzero=nzero+1

in_grid=(nzero==1)

RETURN
END FUNCTION in_grid
