!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

SUBROUTINE polifit(x,y,ndati,a,m1)
USE kinds, ONLY : DP
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
!-------------------------------------------------------------------
SUBROUTINE linsolvx(hc,n,vc,alpha)
!-------------------------------------------------------------------
!
!    This routine is a driver for the correspondent lapack routines
!    which solve a linear system of equations with real coefficients. On 
!    input the matrix is contained in hc, and the known part in vc, on 
!    output the solution is on alpha.
!
!
USE kinds, ONLY : DP
IMPLICIT NONE

INTEGER, INTENT(IN)  :: n        ! input: logical dimension of hc

REAL(DP), INTENT(IN)  ::  hc(n,n),  &  ! input: the matrix to solve
                          vc(n)        ! input: the known part of the system

REAL(DP), INTENT(OUT) ::  alpha(n)     ! output: the solution

INTEGER, ALLOCATABLE    :: iwork(:)
INTEGER :: info

ALLOCATE(iwork(n))

CALL dgetrf(n,n,hc,n,iwork,info)
CALL errore('linsolvx','error in factorization',abs(info))
alpha=vc
CALL dgetrs('N',n,1,hc,n,iwork,alpha,n,info)
CALL errore('linsolvx','error in solving',abs(info))
 
DEALLOCATE( iwork )

RETURN
END SUBROUTINE linsolvx
