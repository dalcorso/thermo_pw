!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

SUBROUTINE polifit(x,y,ndati,a,m1)
USE kinds, ONLY : DP
USE quadratic_surfaces, ONLY : linsolvx
IMPLICIT NONE
INTEGER, INTENT(IN)   ::   ndati, &  ! number of data points
                           m1        ! polynomial coefficients 

REAL(DP), INTENT(IN)  :: x(ndati), y(ndati)
REAL(DP), INTENT(OUT) :: a(m1)

REAL(DP) :: amat(m1,m1), bvec(m1), eigv(m1,m1)
INTEGER  ::  i, j, k       ! counters

DO k=1,m1
   DO j=1,m1
      amat(k,j)=0.d0
      DO i=1,ndati
         amat(k,j)=amat(k,j)+x(i)**(j-1)*x(i)**(k-1)
      ENDDO
   ENDDO
   bvec(k)=0.d0
   DO i=1,ndati
      bvec(k)=bvec(k)+y(i)*x(i)**(k-1)
   ENDDO
ENDDO

CALL linsolvx(amat,m1,bvec,a)

RETURN
END SUBROUTINE polifit
!
