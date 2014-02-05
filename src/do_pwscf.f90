!
! Copyright (C) 2013 Andrea Dal Corso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE do_pwscf (lscf_) 
  !----------------------------------------------------------------------------
  !
  ! ... Run an instance of the Plane Wave Self-Consistent Field code 
  ! ... MPI initialization and input data reading is performed in the 
  ! ... calling code. This is a simplified version of run_pwscf that
  ! ... does only electronic scf or nscf run
  !
  USE force_mod,        ONLY : lforce, lstres, sigma, force
  USE control_flags,    ONLY : lscf
  USE basis,            ONLY : starting_pot, startingconfig
  !
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: lscf_
  !
  !
  IF (lscf_) THEN
     starting_pot ='atomic'
     startingconfig='input'
     lscf=.TRUE.
  ELSE
     starting_pot ='file'
     startingconfig='file'
     lscf=.FALSE.
  ENDIF
  
  CALL setup ()
  !
  CALL init_run(.TRUE.)
  !
  IF ( .NOT. lscf) THEN
     CALL non_scf ()
  ELSE
     CALL electrons()
  END IF
  IF ( lforce .AND. lscf ) CALL forces()
  !
  ! ... stress calculation
  !
  IF ( lstres .AND. lscf ) CALL stress ( sigma )
  !
  ! ... save final data file
  !
  CALL punch('all')
  !
  CALL close_files(.TRUE.)
  !
  CALL clean_pw( .FALSE. )

  !
  RETURN
  !
END SUBROUTINE do_pwscf
