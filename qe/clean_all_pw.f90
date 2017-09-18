
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! TB
! included deallocation of forcefield of gate 'forcegate'
!
!----------------------------------------------------------------------
SUBROUTINE clean_all_pw( )
  !----------------------------------------------------------------------
  !    
  ! ... This routine deallocates dynamically allocated arrays that
  ! ... clean_pw would clean only when lflag=.TRUE.
  ! ... if lflag=.TRUE.  all arrays are deallocated (end of calculation)
  ! ... if lflag=.FALSE. ion-related variables and arrays allocated
  ! ... at the very beginning of the calculation (routines iosys, read_file,
  ! ... setup, read_pseudo) are not deallocated; all others arrays are.
  ! ... This is used when a new calculation has to be performed (e.g. in neb,
  ! ... phonon, vc-relax). Beware: the new calculation should not call any
  ! ... of the routines mentioned above
  !
  USE force_mod,            ONLY : force
  USE symm_base,            ONLY : irt
  USE uspp_param,           ONLY : upf
  USE extfield,             ONLY : forcefield, forcegate
  USE atom,                 ONLY : msh, rgrid
  USE radial_grids,         ONLY : deallocate_radial_grid
  USE pseudo_types,         ONLY : deallocate_pseudo_upf
  USE ions_base,            ONLY : deallocate_ions_base

  !
  USE london_module,        ONLY : dealloca_london
  USE xdm_module,           ONLY : cleanup_xdm
  USE constraints_module,   ONLY : deallocate_constraint
  !
  !
  IMPLICIT NONE

  INTEGER :: nt
  !
  ! ... arrays allocated at the very beginning of the calculation
  !
  IF( ALLOCATED( upf ) ) THEN
     DO nt = 1, SIZE( upf )
        CALL deallocate_pseudo_upf( upf( nt ) )
     END DO
     DEALLOCATE( upf )
  END IF
  IF (ALLOCATED(msh)) DEALLOCATE (msh)
  CALL deallocate_radial_grid(rgrid)
  !
  CALL deallocate_ions_base()
  !
  IF ( ALLOCATED( force ) )      DEALLOCATE( force )
  IF ( ALLOCATED( forcefield ) ) DEALLOCATE( forcefield )
  IF ( ALLOCATED( forcegate ) )  DEALLOCATE( forcegate )
  IF ( ALLOCATED( irt ) )        DEALLOCATE( irt )
  !
  CALL dealloca_london()
  CALL cleanup_xdm()
  CALL deallocate_constraint()
  !
  RETURN
  !
END SUBROUTINE clean_all_pw
