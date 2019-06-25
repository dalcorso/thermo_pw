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


! Se metto il loop che simetrizza nel caso cubico farebbe
! 9 cicli inutili perché sulla matrice ci sono zeri 
! Pero gli si può dire se trova un elemento nullo di 
! fare CYCLE

SELECT CASE (laue)
   CASE(29,32)
      DO itemp=2,ntemp-1
         el_cons_t(2,2,itemp)=el_cons_t(1,1,itemp)
         el_cons_t(3,3,itemp)=el_cons_t(1,1,itemp)
         el_cons_t(1,3,itemp)=el_cons_t(1,2,itemp)
         !el_cons_t(2,1,itemp)=el_cons_t(1,2,itemp)
         el_cons_t(2,3,itemp)=el_cons_t(1,2,itemp)
         !el_cons_t(3,1,itemp)=el_cons_t(1,2,itemp)
         !el_cons_t(3,2,itemp)=el_cons_t(1,2,itemp)
         el_cons_t(5,5,itemp)=el_cons_t(4,4,itemp)
         el_cons_t(6,6,itemp)=el_cons_t(4,4,itemp) 
      END DO
END SELECT 

DO itemp=2,ntemp-1
   DO i=1, 6
      DO j=i+1, 6
         el_cons_t(j,i,itemp)=el_cons_t(i,j,itemp)
      END DO
   END DO
END DO

RETURN
END SUBROUTINE expand_el_cons
