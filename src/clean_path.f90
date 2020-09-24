!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE clean_path(nqaux, xqaux, wqaux)
!---------------------------------------------------------------------
!
! This routine receives a sets of reciprocal lattice points that
! determines a path. It removes from the paths the lines that are too
! short by setting to zero the number of point of that line.
! It computes the total lenght of the path and removes those line
! that are shorter than 1/20 of the total path.
!
! The points are in cartesian coordinates
!
USE kinds, ONLY : DP
IMPLICIT NONE

INTEGER, INTENT(IN) :: nqaux
INTEGER, INTENT(INOUT) :: wqaux(nqaux)
REAL(DP), INTENT(IN) :: xqaux(3, nqaux)

REAL(DP) :: total_length, average_length, fact, leng(nqaux)
INTEGER :: iq

IF (nqaux<2) RETURN
fact=5.0_DP

total_length=0.0_DP
DO iq=1, nqaux-1
   IF (wqaux(iq) > 0) THEN
      leng(iq)= SQRT( (xqaux(1,iq+1) - xqaux(1,iq))**2 + &
                      (xqaux(2,iq+1) - xqaux(2,iq))**2 + &
                      (xqaux(3,iq+1) - xqaux(3,iq))**2 )
      total_length=total_length + leng(iq)
   ENDIF
ENDDO
average_length = total_length / (nqaux-1)

DO iq=1, nqaux-1
   IF ( leng(iq) * fact < average_length ) wqaux(iq)=0
ENDDO

RETURN
END SUBROUTINE clean_path
