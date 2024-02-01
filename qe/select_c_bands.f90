! Copyright (C) 2023 Dal Corso Andrea
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
SUBROUTINE select_c_bands(iter)
!------------------------------------------------------------------------
!
USE control_qe, ONLY : many_k
IMPLICIT NONE
INTEGER, INTENT(IN) :: iter

IF (many_k) THEN
!
!  Call the routine that deals with many k points together
!
   CALL c_bands_many_k(iter)
ELSE
!
!  Standard case. Call the QE routine.
!
   CALL c_bands(iter)
ENDIF

RETURN
END SUBROUTINE select_c_bands
