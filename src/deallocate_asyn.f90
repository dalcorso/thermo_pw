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
                             iqw, irrw, comp_irr_iq_iw, comp_iq_iw,      &
                             comp_f_iw, lpwband, lef
  IMPLICIT NONE
  !
  IF (ALLOCATED(lpwscf))     DEALLOCATE(lpwscf) 
  IF (ALLOCATED(lpwband))    DEALLOCATE(lpwband) 
  IF (ALLOCATED(lef))        DEALLOCATE(lef) 
  IF (ALLOCATED(lberry))     DEALLOCATE(lberry) 
  IF (ALLOCATED(lphonon))    DEALLOCATE(lphonon) 
  IF (ALLOCATED(lstress))    DEALLOCATE(lstress) 
  IF (ALLOCATED(geometry))   DEALLOCATE(geometry) 
  IF (ALLOCATED(iqw))        DEALLOCATE(iqw) 
  IF (ALLOCATED(irrw))       DEALLOCATE(irrw) 
  IF (ALLOCATED(comp_irr_iq_iw)) DEALLOCATE(comp_irr_iq_iw) 
  IF (ALLOCATED(comp_iq_iw))  DEALLOCATE(comp_iq_iw) 
  IF (ALLOCATED(comp_f_iw)) DEALLOCATE(comp_f_iw) 
  ! 
  RETURN
  !
END SUBROUTINE deallocate_asyn
