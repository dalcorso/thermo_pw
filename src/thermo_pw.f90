!
! Copyright (C) 2013-2020 Andrea Dal Corso
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
  ! ... using the harmonic and/or quasi-harmonic approximation and the
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

  USE control_thermo,   ONLY : lev_syn_2, lph, lpwscf_syn_1, lectqha, &
                               lpart2_pw
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
  INTEGER  :: part, nwork
  CHARACTER (LEN=9)   :: code = 'THERMO_PW'
  !
  ! ... Initialize mpi, the clocks, and all the parallelism of QE
  !
  CALL thermo_startup(code)
  !
  ! ... and begin with the initialization part
  ! ... reading input
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
  !... setup the work to do in part 1 for the given task
  !
  part = 1
  !
  CALL initialize_thermo_work(nwork, part)
  !
  !  In this part the images work asynchronously. No communication is
  !  allowed except through the master-slave mechanism
  !
  CALL run_thermo_asynchronously(nwork, part, 1)
  !
  !  In this part all images are re-synchronized again and can communicate 
  !  through the world_comm communicator.
  !  Compute here the physical quantities that can be obtained
  !  from the first set of calculations.
  !
  CALL manage_syn_1(nwork)
  !
  !  This part makes a self consistent pw.x calculation followed by
  !  a non self-consistent one if required. Only one image does the 
  !  calculation. This part is used also for surface band structure 
  !  and to identify surface states.
  !
  IF (lpwscf_syn_1) CALL manage_sync_pw()
  !
  IF (lpart2_pw) THEN
     !
     !   here we make another asynchronous set of calculations with many runs
     !   of pw.x. 
     !
     !... setup the work to do in part 2 for the given task
     !
     part=2
     !
     CALL initialize_thermo_work(nwork, part)
     !
     !  Asynchronous work starts again. No communication is
     !  allowed except through the master-slave mechanism
     !
     CALL run_thermo_asynchronously(nwork, part, 1)
     !
     !  In this part all images are re-synchronized again and can communicate 
     !  through the world_comm communicator.
     !  Compute here the physical quantities that can be obtained
     !  from the second set of calculations.
     !
     CALL manage_syn_2(nwork)
     !
  ELSEIF (lph) THEN
     !
     !   Here we make one or several phonon calculations, using 
     !   images and running asynchronously.
     !   After each calculation one can plot the phonon dispersions,
     !   the density of phonon states and the harmonic thermodynamic
     !   quantities. A single phonon dispersion made here can also provide
     !   the dielectric constant, the eels spectrum, or the Born effective
     !   charges. The management of the phonon asynchronous work is made
     !   inside the manager.
     !
     CALL manage_ph_run()
     ! 
     !   The stop_signal stops all the calculation after each image has
     !   terminated its task even if there are still tasks to do.
     !   It is triggered by max_seconds.
     !
     IF (stop_signal_activated) GOTO 1000
     !
     !   Here the Helmholtz free energy at each geometry is available.
     !   We can write on file the free energy as a function of the volume at
     !   any temperature. For each temperature we can fit the free energy
     !   (or the Gibbs energy if we have a finite pressure) with an
     !   equation of state or with a quadratic or quartic polynomial. 
     !   We save the minimum volume or the crystal parameters. 
     !   With the equation of state we save also the bulk modulus and 
     !   its pressure derivative for each temperature.
     !   One can also compute here the Gruneisen parameters.
     !
     IF (lev_syn_2) CALL manage_anhar_routines()
     !
     !   When lectqha=.TRUE. here we have computed the Helmholtz free energy 
     !   for all the strained geometries (for a single unperturbed geometry or
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
