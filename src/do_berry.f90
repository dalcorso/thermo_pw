!
! Copyright (C) 2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE do_berry ( exit_status, polar, tot_b_phase, nppl ) 
  !----------------------------------------------------------------------------
  !
  ! ... Run an instance of the Plane Wave Self-Consistent Field code 
  ! ... MPI initialization and input data reading is performed in the 
  ! ... calling code - returns in exit_status the exit code for pw.x, 
  ! ... returned in the shell. Values are:
  ! ... * 0: completed successfully
  ! ... * 1: an error has occurred (value returned by the errore() routine)
  ! ... * 2-127: convergence error
  ! ...   * 2: scf convergence error
  ! ...   * 3: ion convergence error
  ! ... * 128-255: code exited due to specific trigger
  !       * 255: exit due to user request, or signal trapped,
  !              or time > max_seconds
  ! ...     (note: in the future, check_stop_now could also return a value
  ! ...     to specify the reason of exiting, and the value could be used
  ! ..      to return a different value for different reasons)
  ! ... Will be eventually merged with NEB
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout, ionode
  USE parameters,       ONLY : ntypx, npk
  USE upf_params,       ONLY : lmaxx
  USE cell_base,        ONLY : at, omega, alat
  USE starting_scf,     ONLY : starting_pot, startingconfig
  USE control_flags,    ONLY : gamma_only, lscf, lbands, ethr, &
                               istep, nstep, lbfgs
  USE initial_conf,     ONLY : nosym_save
  USE initial_param,    ONLY : ethr0
  USE polarization_vector,  ONLY : mod_tot
  USE symm_base,        ONLY : nosym
  USE bp,               ONLY : pdl_tot, nppstr, gdir, lberry
  !
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: exit_status
  INTEGER, INTENT(IN) :: nppl
  REAL(DP), INTENT(OUT) :: polar(3), tot_b_phase(3) 

  REAL(DP) :: polar_at(3), atmod
  INTEGER :: idir
  !
  exit_status = 0
  IF ( ionode ) WRITE( unit = stdout, FMT = 9010 ) ntypx, npk, lmaxx
  !
  IF ( gamma_only ) WRITE( UNIT = stdout, &
     & FMT = '(/,5X,"gamma-point specific algorithms are used")' )
  !
  starting_pot ='file'
  startingconfig='file'
  lscf=.FALSE.
  lberry=.TRUE.
  nppstr=nppl

  ethr=ethr0
  istep=0

  DO idir=1,3 
     !
     gdir=idir
     !
     nosym=nosym_save
     !
     CALL setup ()
     !
     CALL init_run()
     !
     !
     ! ... band structure calculation and berry phase calculation
     !
     CALL non_scf ()
     !
     !   this is the phase, electronic+ionic. 
     !
     !
     !   the routine bp_c_phase does not bring the total phase
     !   in the standard interval [-1,1) (mod 2) or [-1/2,1/2). We do
     !   it here otherwise the derivatives of polarization might explode.
     !
     IF (mod_tot==2) then
        pdl_tot=pdl_tot-2.d0*nint(pdl_tot/2.d0)
     ELSE
        pdl_tot=pdl_tot-nint(pdl_tot)
     ENDIF
     !
     tot_b_phase(idir) = pdl_tot
     !
     !   We divide here by the volume that might change for each strain. 
     !
     polar_at(idir) = pdl_tot / omega
     !
     CALL close_files(.TRUE.)
     !
     CALL clean_pw_tpw( .FALSE. )
     !
  END DO 
  nosym=nosym_save
!
!  Here we compute polarization of this structure in cartesian coordinates
!  and in units of e/bohr**2
!
  polar(:) = polar_at(1)*at(:,1) + polar_at(2)*at(:,2) + polar_at(3)*at(:,3)
  polar(:) = polar * alat

  CALL print_polarization(polar(:), .FALSE.)
  !
  RETURN
  !
9010 FORMAT( /,5X,'Current dimensions of program PWSCF are:', &
           & /,5X,'Max number of different atomic species (ntypx) = ',I2,&
           & /,5X,'Max number of k-points (npk) = ',I9,&
           & /,5X,'Max angular momentum in pseudopotentials (lmaxx) = ',i2)
  !
END SUBROUTINE do_berry
