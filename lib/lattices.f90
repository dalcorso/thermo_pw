!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE lattices
!
!  This module contains variables and routines to deal with Bravais
!  lattices. 
!
!  Presently it has the routine: 
!
!  compute_conventional  : this routine receives the direct lattice
!                          vectors of a centered cell and computes the
!                          vectors of the corresponding primitive cell.
!                          atc in output are in the same units of the at
!                          in input. They are in cartesian coordinates.
!                          This routine must be matched to the routine
!                          latgen because it depends on the choice of the
!                          direct lattice vectors.
!
!
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC compute_primitive

CONTAINS

  SUBROUTINE compute_primitive(at, atp, ibrav)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ibrav
  REAL(DP), INTENT(IN) :: at(3,3)
  REAL(DP), INTENT(INOUT) :: atp(3,3)

  SELECT CASE (ibrav)
     CASE(2)
         atp(:,1)= -at(:,1) + at(:,2) - at(:,3)
         atp(:,2)= -at(:,1) + at(:,2) + at(:,3)
         atp(:,3)=  at(:,1) + at(:,2) - at(:,3)
     CASE(3)
         atp(:,1)=  at(:,1) - at(:,2) 
         atp(:,2)=  at(:,2) - at(:,3)
         atp(:,3)=  at(:,1) + at(:,3)
     CASE(7)
         atp(:,1)=  at(:,1) - at(:,3) 
         atp(:,2)= -at(:,1) + at(:,2)
         atp(:,3)=  at(:,2) + at(:,3)
     CASE(9)
         atp(:,1)=  at(:,1) - at(:,2) 
         atp(:,2)=  at(:,1) + at(:,2)
         atp(:,3)=  at(:,3)
     CASE(10)
         atp(:,1)=  at(:,1) + at(:,2) - at(:,3)
         atp(:,2)= -at(:,1) + at(:,2) + at(:,3)
         atp(:,3)=  at(:,1) - at(:,2) + at(:,3)
     CASE(11)
         atp(:,1)=  at(:,1) - at(:,2) 
         atp(:,2)=  at(:,2) - at(:,3)
         atp(:,3)=  at(:,1) + at(:,3)
     CASE(13)
         atp(:,1)=  at(:,1) + at(:,2) 
         atp(:,2)= -at(:,1) + at(:,2)
         atp(:,3)=  at(:,3)
     CASE(-13)
         atp(:,1)=  at(:,1) + at(:,3) 
         atp(:,2)=  at(:,2)
         atp(:,3)= -at(:,1) + at(:,3)
     CASE DEFAULT
!
!  If this is not a known centered lattice we simply copy the at in atp
!
         atp=at
  END SELECT   
  RETURN
  END SUBROUTINE compute_primitive

  END MODULE lattices
