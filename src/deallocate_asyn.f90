!
! Copyright (C) 2013 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE deallocate_asyn()
  !-----------------------------------------------------------------------
  !
  !  This routine deallocates the variables that control the thermo calculation
  !
  USE kinds, ONLY : DP
  USE control_thermo, ONLY : lpwscf, lbands, lphonon, lstress, lberry
  IMPLICIT NONE
  !
  IF (ALLOCATED(lpwscf))     DEALLOCATE(lpwscf) 
  IF (ALLOCATED(lbands))     DEALLOCATE(lbands) 
  IF (ALLOCATED(lberry))     DEALLOCATE(lberry) 
  IF (ALLOCATED(lphonon))    DEALLOCATE(lphonon) 
  IF (ALLOCATED(lstress))    DEALLOCATE(lstress) 
  ! 
  RETURN
  !
END SUBROUTINE deallocate_asyn
