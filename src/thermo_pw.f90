!
! Copyright (C) 2013-2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
PROGRAM thermo_pw
  !-----------------------------------------------------------------------
  !
  ! ... This is a driver for the calculation of thermodynamic quantities,
  ! ... using the harmonic and/or quasiharmonic approximation and the
  ! ... plane waves pseudopotential method.
  ! ... It reads the input of pw.x and an input that specifies
  ! ... the calculations to do and the parameters for these calculations.
  ! ... It checks the scratch directories to see what has been already
  ! ... calculated. The info for the quantities that have been already
  ! ... calculated is read inside the code. The others tasks are scheduled
  ! ... and distributed to the image driver.
  ! ... If there are several available images the different tasks are
  ! ... carried out in parallel. This driver can carry out a scf 
  ! ... calculation, a non scf calculation to determine the band structure,
  ! ... or a linear response calculation at a given q and for a given
  ! ... representation. Finally several post processing tasks are carried
  ! ... out in parallel or by the root image. 
  ! ... The type of calculations currently implemented 
  ! ... are described in the user's guide and in the developer's guide.
  ! ...

  USE kinds,            ONLY : DP

  USE thermo_mod,       ONLY : what, ngeo
  USE control_thermo,   ONLY : lev_syn_1, lev_syn_2, lpwscf_syn_1,         &
                               lph, lconv_ke_test, lconv_nk_test,    &
                               lelastic_const, lectqha,            &
                               lpiezoelectric_tensor, lpolarization,       &
                               lpart2_pw

  USE control_elastic_constants, ONLY : ngeom
!
!  variables of pw or phonon used here
!
  USE check_stop,       ONLY : max_seconds, check_stop_init
!
!  parallelization control 
!
  USE mp_asyn,          ONLY : stop_signal_activated
  !
  IMPLICIT NONE
  !
  INTEGER  :: part, nwork, iaux
  REAL(DP) :: polar(3)
  CHARACTER (LEN=9)   :: code = 'THERMO_PW'
  CHARACTER (LEN=256) :: auxdyn=' '
  !
  !  Initialize mpi, the clocks and all the parallelism of QE
  !
  CALL thermo_startup(code)
  !
  ! ... and begin with the initialization part
  ! ... readin input
  !
  CALL thermo_readin()
  !
  ! ... setup common variables needed for many tasks
  !
  CALL thermo_setup()
  !
  ! ... make a summary of what will be computed
  !
  CALL thermo_summary()
  !
  ! ... initialize the timer that stops the calculation after max_seconds
  !
  CALL check_stop_init(max_seconds)
  !
  !... setup the work to do for the given task
  !
  part = 1
  !
  CALL initialize_thermo_work(nwork, part, iaux)
  !
  !  In this part the images work asynchronously. No communication is
  !  allowed except though the master-workers mechanism
  !
  CALL run_thermo_asynchronously(nwork, part, iaux, auxdyn)
  !
  !  In this part all images are synchronized and can communicate 
  !  their results thought the world_comm communicator.
!
!  In the kinetic energy test write the results
!
  IF (lconv_ke_test) THEN
     CALL write_e_ke()
     CALL plot_e_ke()
  ENDIF
!
! In the k-point test write the results
!
  IF (lconv_nk_test) THEN
     CALL write_e_nk()
     CALL plot_e_nk()
  ENDIF
!
!  In a Murnaghan equation calculation determine the lattice constant,
!  bulk modulus and its pressure derivative and write the results.
!  Otherwise interpolate the energy with a quadratic or quartic polynomial.
!
  IF (lev_syn_1) CALL manage_energy_minimum(nwork)
!
!  When computing the elastic constants as a function of temperature
!  here we have the energy for each strained geometry and can compute
!  the T=0 K elastic constants for all the reference geometries.
!
  IF (lectqha) CALL manage_elastic_cons(nwork,ngeom)

  CALL deallocate_asyn()
!
!  This part makes a self consistent pw.x calculation followed by
!  a nonselfconsisten one if required. Only one image does the calculation.
!  This part is used also for surface band structure and to identify 
!  surface states.
!
  IF (lpwscf_syn_1) CALL manage_sync_pw()
!
!   here we make another asynchronous calculation with many runs
!   of pw.x. It can be used for instance to compute the T=0 K elastic
!   elastic constants at the minimum of the murnaghan computed in the
!   first part. 
!
  IF (lpart2_pw) THEN
     !
     part=2
     !
     CALL initialize_thermo_work(nwork, part, iaux)
     !
     !  Asynchronous work starts again. No communication is
     !  allowed except though the master workers mechanism
     !
     CALL run_thermo_asynchronously(nwork, part, 1, auxdyn)
     !
     ! here we return synchronized and calculate the elastic constants 
     ! from energy or stress 
     !
     IF (lelastic_const) CALL manage_elastic_cons(nwork, 1)
     !
     IF (lpiezoelectric_tensor) CALL manage_piezo_tensor(nwork)
     !
     IF (lpolarization) CALL print_polarization(polar, .TRUE. )

     CALL deallocate_asyn()
  ENDIF

  IF (what(1:8) /= 'mur_lc_t') ngeo=1
!
!   This part makes one or several phonon calculations, using 
!   images and running asynchronously.
!   After each calculation one can plot the phonon dispersions,
!   compute the density of phonon states and the harmonic thermodynamic
!   quantities. A single phonon dispersion made here can also provide
!   the dielectric constant, the eels spectrum, or the Born effective
!   charges.
!
  IF (lph) THEN

     CALL manage_ph_run()

     IF (stop_signal_activated) GOTO 1000
!
!     Here the Helmholtz free energy at each geometry is available.
!     We can write on file the free energy as a function of the volume at
!     any temperature. For each temperature we can fit the free energy
!     (or the Gibbs energy if we have a finite pressure) with a 
!     Murnaghan equation or with a quadratic or quartic polynomial. 
!     We save the minimum volume or the crystal parameters. 
!     With the Murnaghan fit we save also the bulk modulus and 
!     its pressure derivative for each temperature.
!     One can also compute here the Gruneisen parameters.
!
     IF (lev_syn_2) CALL manage_anhar_routines()
!
!   When lectqha=.TRUE. here we have computed the Helmholtz free energy for
!   all the strained geometries (for a single unpertubed geometry or
!   for a mesh of unperturbed geometries) and compute here the
!   temperature dependent elastic constants
!
     IF (lectqha) CALL manage_elastic_cons_qha()
  ENDIF
  !
1000  CALL thermo_end(code) 
  !
  STOP
  !
END PROGRAM thermo_pw
