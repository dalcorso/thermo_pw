! Copyright (C) 2019 C. Malica
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE expand_el_cons(el_cons_t, laue, ibrav)

USE kinds,      ONLY : DP
USE temperature,      ONLY : ntemp, temp

IMPLICIT NONE
REAL(DP), INTENT(INOUT) :: el_cons_t(6,6,ntemp)
INTEGER,  INTENT(IN) :: ibrav, laue
INTEGER :: itemp, i, j

SELECT CASE (laue)

   CASE(29,32)
!
!  cubic T_h (m-3), O_h (m-3m)
!
      DO itemp=2,ntemp-1
         el_cons_t(2,2,itemp)=el_cons_t(1,1,itemp)
         el_cons_t(3,3,itemp)=el_cons_t(1,1,itemp)
         el_cons_t(1,3,itemp)=el_cons_t(1,2,itemp)
         el_cons_t(2,3,itemp)=el_cons_t(1,2,itemp)
         el_cons_t(5,5,itemp)=el_cons_t(4,4,itemp)
         el_cons_t(6,6,itemp)=el_cons_t(4,4,itemp) 
      END DO

   CASE(25) 
!
!  trigonal D_3d (-3m)
!
      DO itemp=2,ntemp-1
         el_cons_t(2,2,itemp)=el_cons_t(1,1,itemp)
         el_cons_t(2,3,itemp)=el_cons_t(1,3,itemp)
         el_cons_t(2,4,itemp)=-el_cons_t(1,4,itemp)
         el_cons_t(5,5,itemp)=el_cons_t(4,4,itemp)
         el_cons_t(5,6,itemp)=el_cons_t(1,4,itemp)
         el_cons_t(6,6,itemp)=(el_cons_t(1,1,itemp)-&
                         el_cons_t(1,2,itemp))/2.0_DP
      END DO

   CASE (27)
!
!  trigonal S_6 (-3)
!
      DO itemp=2,ntemp-1
         el_cons_t(1,5,itemp)=-el_cons_t(2,5,itemp)
         el_cons_t(2,2,itemp)=el_cons_t(1,1,itemp)
         el_cons_t(2,3,itemp)=el_cons_t(1,3,itemp)
         el_cons_t(2,4,itemp)=-el_cons_t(1,4,itemp)
         el_cons_t(5,5,itemp)=el_cons_t(4,4,itemp)
         el_cons_t(5,6,itemp)=el_cons_t(1,4,itemp)
         el_cons_t(6,6,itemp)=(el_cons_t(1,1,itemp)-&
                         el_cons_t(1,2,itemp))/2.0_DP
      END DO

   CASE (19,23)
!
!  hexagonal C_6h (6/m), D_6h (6/mmm)
!
      DO itemp=2,ntemp-1
         el_cons_t(2,2,itemp)=el_cons_t(1,1,itemp)
         el_cons_t(2,3,itemp)=el_cons_t(1,3,itemp)
         el_cons_t(5,5,itemp)=el_cons_t(4,4,itemp)
         el_cons_t(6,6,itemp)=(el_cons_t(1,1,itemp)-&
                         el_cons_t(1,2,itemp))/2.0_DP
      END DO

   CASE(22)
!
!  tetragonal D_4h (4/mmm)
!
      DO itemp=2,ntemp-1
         el_cons_t(2,2,itemp)=el_cons_t(1,1,itemp)
         el_cons_t(2,3,itemp)=el_cons_t(1,3,itemp)
         el_cons_t(5,5,itemp)=el_cons_t(4,4,itemp)
      END DO

   CASE(20)
!
!  orthorhombic D_2h (mmm)
!
      !There are no other elastic constants 
      !dependent from those read.

   CASE(18)
!
!  tetragonal C_4h (4/m)
!
      DO itemp=2,ntemp-1
         el_cons_t(2,2,itemp)=el_cons_t(1,1,itemp)
         el_cons_t(2,3,itemp)=el_cons_t(1,3,itemp)
         el_cons_t(2,6,itemp)=-el_cons_t(1,6,itemp)
         el_cons_t(5,5,itemp)=el_cons_t(4,4,itemp)
      END DO

   CASE(16)
!
!    monoclinic case, class C_2h (2/m) 
!
!    There are no other elastic constants 
!    dependent from those read in both
!    b-unique and c-unique cases.

   CASE(2)
!
!    triclinic case or generic 
!
!    There are no other elastic constants 
!    dependent from those read.

END SELECT 

DO i=1, 6
   DO j=i+1, 6
      el_cons_t(j,i,:)=el_cons_t(i,j,:)
   END DO
END DO

RETURN
END SUBROUTINE expand_el_cons
