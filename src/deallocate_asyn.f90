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
  USE control_thermo, ONLY : lpwscf, lphonon, lstress, lberry, geometry, &
                             iqw, irrw, comp_f_work
  IMPLICIT NONE
  !
  IF (ALLOCATED(lpwscf))     DEALLOCATE(lpwscf) 
  IF (ALLOCATED(lberry))     DEALLOCATE(lberry) 
  IF (ALLOCATED(lphonon))    DEALLOCATE(lphonon) 
  IF (ALLOCATED(lstress))    DEALLOCATE(lstress) 
  IF (ALLOCATED(geometry))   DEALLOCATE(geometry) 
  IF (ALLOCATED(iqw))        DEALLOCATE(iqw) 
  IF (ALLOCATED(irrw))       DEALLOCATE(irrw) 
  IF (ALLOCATED(comp_f_work)) DEALLOCATE(comp_f_work) 
  ! 
  RETURN
  !
END SUBROUTINE deallocate_asyn
