!
! Copyright (C) 2013-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE do_pwscf ( exit_status, lscf_ ) 
  !----------------------------------------------------------------------------
  !
  !! Author: Paolo Giannozzi
  !! License: GNU 
  !! Summary: Run an instance of the Plane Wave Self-Consistent Field code
  !!
  !! Run an instance of the Plane Wave Self-Consistent Field code 
  !! MPI initialization and input data reading is performed in the 
  !! calling code - returns in exit_status the exit code for pw.x, 
  !! returned in the shell. Values are:
  !! * 0: completed successfully
  !! * 1: an error has occurred (value returned by the errore() routine)
  !! * 2-127: convergence error
  !!    * 2: scf convergence error
  !!    * 3: ion convergence error
  !! * 128-255: code exited due to specific trigger
  !!    * 255: exit due to user request, or signal trapped,
  !!          or time > max_seconds
  !!     (note: in the future, check_stop_now could also return a value
  !!     to specify the reason of exiting, and the value could be used
  !!     to return a different value for different reasons)
  !! @Note
  !! 10/01/17 Samuel Ponce: Add Ford documentation
  !! @endnote
  !!
  !
  USE kinds,                ONLY : DP
  USE mp,                   ONLY : mp_bcast, mp_sum
  USE io_global,            ONLY : stdout, ionode, ionode_id
  USE parameters,           ONLY : ntypx, npk
  USE upf_params,           ONLY : lmaxx
  USE initial_param,        ONLY : ethr0
  USE cell_base,            ONLY : fix_volume, fix_area
  USE control_flags,        ONLY : conv_elec, gamma_only, ethr, lscf, treinit_gvecs
  USE control_flags,        ONLY : conv_ions, istep, nstep, restart, lmd, &
                                   lbfgs, io_level, lensemb, lforce=>tprnfor, &
                                   tstress
  USE cellmd,               ONLY : lmovecell
  USE command_line_options, ONLY : command_line
  USE force_mod,            ONLY : sigma, force
  USE check_stop,           ONLY : check_stop_init, check_stop_now
  USE control_pwrun,        ONLY : do_punch
  USE basis,                ONLY : starting_pot, starting_wfc, startingconfig
  USE mp_images,            ONLY : intra_image_comm
  USE extrapolation,        ONLY : update_file, update_pot
  USE scf,                  ONLY : rho
  USE lsda_mod,             ONLY : nspin
  USE fft_base,             ONLY : dfftp
  USE qmmm,                 ONLY : qmmm_initialization, qmmm_shutdown, &
                                   qmmm_update_positions, qmmm_update_forces
  USE qexsd_module,         ONLY : qexsd_set_status, qexsd_reset_steps
  USE xc_lib,               ONLY : xclib_dft_is, stop_exx, exx_is_active
  USE beef,                 ONLY : beef_energies
  USE ldaU,                 ONLY : lda_plus_u
  USE add_dmft_occ,         ONLY : dmft
  USE extffield,            ONLY : init_extffield, close_extffield
  USE input_parameters,     ONLY : nextffield

  USE device_fbuff_m,             ONLY : dev_buf
#if defined (__ENVIRON)
  USE plugin_flags,      ONLY : use_environ
  USE environ_pw_module, ONLY : is_ms_gcs, init_ms_gcs
#endif
#if defined (__OSCDFT)
  USE plugin_flags,      ONLY : use_oscdft
  USE oscdft_base,       ONLY : oscdft_ctx
  USE oscdft_functions,  ONLY : oscdft_run_pwscf
#endif
  !
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: exit_status
  !! Gives the exit status at the end
  LOGICAL :: lscf_
  LOGICAL, EXTERNAL :: matches
  ! checks if first string is contained in the second
  !
  ! ... local variables
  !
  INTEGER :: idone 
  ! counter of electronic + ionic steps done in this run
  INTEGER :: ions_status
  ! ions_status =  3  not yet converged
  ! ions_status =  2  converged, restart with nonzero magnetization
  ! ions_status =  1  converged, final step with current cell needed
  ! ions_status =  0  converged, exiting
  INTEGER :: ierr
  !
  LOGICAL :: optimizer_failed = .FALSE.
  !
  IF (lscf_) THEN
     CALL start_clock('tpw_scf_pw')
  ELSE
     CALL start_clock('tpw_nscf_pw')
  ENDIF
  exit_status = 0
  ions_status = 3
  IF ( ionode ) WRITE( UNIT = stdout, FMT = 9010 ) ntypx, npk, lmaxx
  !
  IF ( gamma_only ) WRITE( UNIT = stdout, &
     & FMT = '(/,5X,"gamma-point specific algorithms are used")' )

  IF (lscf_) THEN
     starting_pot ='atomic'
     startingconfig='input'
     lscf=.TRUE.
  ELSE
     starting_pot ='file'
     startingconfig='file'
     starting_wfc = 'atomic+random'
     lscf=.FALSE.
     lbfgs=.FALSE.
     lforce=.FALSE.
     tstress=.FALSE.
!
!   in the nscf case we save the wavefunctions to allow restart
!
     IF (io_level<1) io_level=1
  ENDIF
  ethr=ethr0
  istep=0

  !
  ! call to void routine for user defined / plugin patches initializations
  !
#if defined(__LEGACY_PLUGINS)
  CALL plugin_initialization()
#endif
#if defined (__ENVIRON)
  IF (use_environ) THEN
     IF (is_ms_gcs()) CALL init_ms_gcs()
  END IF
#endif
  !
  CALL check_stop_init()
  !
  CALL setup_tpw ()
  !
  CALL qmmm_update_positions()
  !
  ! ... dry run: code will stop here if called with exit file present
  ! ... useful for a quick and automated way to check input data
  !
  IF ( nstep == 0 .OR. check_stop_now() ) THEN
     CALL pre_init()
     CALL data_structure( gamma_only )
     CALL summary()
     CALL memory_report()
     CALL qexsd_set_status( exit_status )
     CALL punch( 'config-init' )
     exit_status = 255
     RETURN
  ENDIF
  !
  CALL init_run()
  !
  !  read external force fields parameters
  !
  IF ( nextffield > 0 .AND. ionode) THEN
     !
     CALL init_extffield( 'PW', nextffield )
     !
  END IF
  !
  IF ( check_stop_now() ) THEN
     CALL qexsd_set_status( exit_status )
     CALL punch( 'config' )
     exit_status = 255
     RETURN
  ENDIF
  !
  main_loop: DO idone = 1, nstep
     !
     ! ... electronic self-consistency or band structure calculation
     !
#if defined (__OSCDFT)
     IF (use_oscdft) THEN
        CALL oscdft_run_pwscf(oscdft_ctx)
     ELSE
#endif
     IF ( .NOT. lscf) THEN
        CALL non_scf()
     ELSE
        CALL electrons_tpw()
     END IF
#if defined (__OSCDFT)
     END IF
#endif

     !
     ! ... code stopped by user or not converged
     !
     IF ( check_stop_now() .OR. .NOT. conv_elec ) THEN
        IF ( check_stop_now() ) THEN
            exit_status = 255
        ELSE
           IF (dmft) THEN
              exit_status =  131
           ELSE
              exit_status = 2
           ENDIF
        ENDIF
        CALL qexsd_set_status(exit_status)
        IF(exx_is_active()) then
          CALL punch( 'all' )
        ELSE
          CALL punch( 'config' )
        ENDIF
        RETURN
     ENDIF
     !
     ! ... file in CASINO format written here if required
     !
     IF ( lmd ) THEN
        CALL pw2casino( istep )
     ELSE
        CALL pw2casino( 0 )
     END IF
     !
     ! ... ionic section starts here
     !
     CALL start_clock( 'ions' ); !write(*,*)' start ions' ; FLUSH(6)
     conv_ions = .TRUE.
     !
     ! ... force calculation
     !
     IF ( lforce ) CALL forces()
     !
     ! ... stress calculation
     !
     IF ( tstress ) CALL stress( sigma )
     !
     IF ( lmd .OR. lbfgs ) THEN
        !
        ! ... add information on this ionic step to xml file
        !
        CALL add_qexsd_step( idone )
        !
        IF (fix_volume) CALL impose_deviatoric_stress( sigma )
        IF (fix_area)   CALL impose_deviatoric_stress_2d( sigma )
        !
        ! ... save data needed for potential and wavefunction extrapolation
        !
        CALL update_file()
        !
        ! ... ionic step (for molecular dynamics or optimization)
        !
        CALL move_ions ( idone, ions_status, optimizer_failed )
        conv_ions = ( ions_status == 0 ) .OR. &
                    ( ions_status == 1 .AND. treinit_gvecs )
        !
        IF ( xclib_dft_is('hybrid') )  CALL stop_exx()
        !
        ! ... save restart information for the new configuration
        !
        IF ( idone <= nstep .AND. .NOT. conv_ions ) THEN
            exit_status = 255
            CALL qexsd_set_status( exit_status )
            CALL punch( 'config-only' )
        END IF
        !
     END IF
     !
     CALL stop_clock( 'ions' ); !write(*,*)' stop ions' ; FLUSH(6)
     !
     ! ... send out forces to MM code in QM/MM run
     !
     CALL qmmm_update_forces( force, rho%of_r, nspin, dfftp )
     !
     ! ... exit condition (ionic convergence) is checked here
     !
     IF ( conv_ions .OR. optimizer_failed ) EXIT main_loop
     !
     ! ... receive new positions from MM code in QM/MM run
     !
     CALL qmmm_update_positions()
     !
     ! ... terms of the hamiltonian depending upon nuclear positions
     ! ... are reinitialized here
     !
     IF ( lmd .OR. lbfgs ) THEN
        !
        IF ( ions_status == 1 ) THEN
           !
           ! ... final scf calculation with G-vectors for final cell
           !
           lbfgs=.FALSE.; lmd=.FALSE.
           WRITE( UNIT = stdout, FMT=9020 ) 
           !
           CALL reset_gvectors( )
           !
           ! ... read atomic occupations for DFT+U(+V)
           !
           IF ( lda_plus_u ) CALL read_ns()
           !
        ELSE IF ( ions_status == 2 ) THEN
           !
           ! ... check whether nonzero magnetization is real
           !
           CALL reset_magn()
           !
        ELSE
           !
           IF ( treinit_gvecs ) THEN
              !
              ! ... prepare for next step with freshly computed G vectors
              !
              IF ( lmovecell) CALL scale_h()
              CALL reset_gvectors ( )
              !
           ELSE
              !
              ! ... update the wavefunctions, charge density, potential
              ! ... update_pot initializes structure factor array as well
              !
              CALL update_pot()
              !
              ! ... re-initialize atomic position-dependent quantities
              !
              CALL hinit1()
              !
           END IF
           !
        END IF
        !
        !
     ENDIF
     ! ... Reset convergence threshold of iterative diagonalization for
     ! ... the first scf iteration of each ionic step (after the first)
     !
     ethr = 1.0D-6
     !
     CALL dev_buf%reinit( ierr )
     IF ( ierr .ne. 0 ) CALL infomsg( 'run_pwscf', 'Cannot reset GPU buffers! Some buffers still locked.' )
     !
     !
  ENDDO main_loop
  !
  ! Set correct exit_status
  !
  IF ( .NOT. conv_ions .OR. optimizer_failed ) THEN
      exit_status =  3
  ELSE
      ! All good
      exit_status = 0
  END IF
  !
  ! ... save final data file
  !
  CALL qexsd_set_status( exit_status )
  IF ( lensemb ) CALL beef_energies( )
  IF ( io_level > -1 .AND. do_punch) CALL punch( 'all' )
  IF (.NOT. do_punch) CALL qexsd_reset_steps()
  !
  CALL qmmm_shutdown()
  !
  CALL close_files(.TRUE.)
  !
  CALL laxlib_end()
  IF (lscf_) THEN
     CALL stop_clock('tpw_scf_pw')
  ELSE
     CALL stop_clock('tpw_nscf_pw')
  ENDIF
  CALL clean_pw( .FALSE. )
  !
  RETURN
  !
9010 FORMAT( /,5X,'Current dimensions of program PWSCF are:', &
           & /,5X,'Max number of different atomic species (ntypx) = ',I2,&
           & /,5X,'Max number of k-points (npk) = ',I6,       &
           & /,5X,'Max angular momentum in pseudopotentials (lmaxx) = ',i2)
9020 FORMAT( /,5X,'Final scf calculation at the relaxed structure.', &
          &  /,5X,'The G-vectors are recalculated for the final unit cell', &
          &  /,5X,'Results may differ from those at the preceding step.' )

  !
END SUBROUTINE do_pwscf


