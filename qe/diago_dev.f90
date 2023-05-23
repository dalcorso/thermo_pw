!
! Copyright (C) 2023 Andrea Dal Corso and Xuejun Gong
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__CUDA)
!----------------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE diago_dev(ntime,nsizex,nsizek,nvec,outk,a,b,  &
                      w,x,work,lwork,adiag,bdiag,rwork,lrwork,iwork,liwork, &
                      info,ifail,m)
!----------------------------------------------------------------------------
!
! This routine diagonalizes simultaneously ntime generalized eigenvalues
! problems 
! A x = w B x
! using the lapack device routines zhegv (if all eigenvalues are required)
! or zhegvx (if a few eigenvalues are required). It is supposed to run on 
! a GPU and each thread diagonalizes one matrix. nvec contains the number 
! of bands required (the same for all threads). This routine must be called 
! using the chevron syntax
! CALL diago_dev<<<ntime,1>>>(...)
! On input, the matrices A and B contain the operators for all problems.
! On output, W and X contain the eigenvalues and eigenvectors for all problems.
! The matrices A and B can have different sizes nsizek(ntime) and
! nsizex is the maximum size.
! rwork is a real work space of size lrwork=7 * nsizex. 
! iwork is an integer work space of size lrwork=5 * nsizex. 
! adiag and bdiag are real vectors sufficient to 
! contain the diagonal of A and B for all problems, while work should have
! size lwork for all problems. lwork can be computed 
! outside the routine by calling the routine find_lwork provided below
! with the maximum size.
! The calling (host) routine must allocate work, rwork, iwork, adiag and
! bdiag on the device for all the problems.
! 
! On output, info(ntime) contains the info of each diagonalization.
! On output, ifail(nsizex,ntime) contains the eigenvalues not converged.
! On output, m(ntime) contains the number of eigenvalues computed.
! These arrays must be allocated on the device by the calling routine.
!
USE cudafor
USE kinds, ONLY : DP
IMPLICIT NONE

#include<zhegv_interf.f90>
#include<zhegvx_interf.f90>
INTEGER, VALUE :: nsizex
INTEGER, VALUE :: ntime
INTEGER, VALUE :: lrwork
INTEGER, VALUE :: liwork
INTEGER, VALUE :: lwork
INTEGER, VALUE :: nvec
LOGICAL, DEVICE :: outk(ntime)
INTEGER, DEVICE :: nsizek(ntime)
INTEGER, DEVICE :: m(ntime)
INTEGER, DEVICE :: ifail(nsizex,ntime)
INTEGER, DEVICE :: info(ntime)
INTEGER, DEVICE :: iwork(liwork,ntime)
COMPLEX(DP), DEVICE :: a(nsizex,nsizex,ntime), b(nsizex,nsizex,ntime), &
                      x(nsizex,nsizex,ntime), work(lwork,ntime)
REAL(DP), DEVICE :: rwork(lrwork,ntime)
REAL(DP), DEVICE :: w(nsizex,ntime), adiag(nsizex,ntime), bdiag(nsizex,ntime)

CHARACTER :: jobz, uplo, vrange
REAL(DP) :: vl, vu, abstol
INTEGER :: itype, lda, ldb, ldx, n, l, i, j, nsize, il, iu
INTEGER :: tn

tn = (Blockidx%x-1)*BlockDim%x+Threadidx%x
IF (tn>ntime) RETURN
IF (outk(tn)) RETURN
nsize=nsizek(tn)

jobz = 'V'
uplo = 'U'
itype = 1
lda = nsizex
ldb = nsizex
ldx = nsizex
n=nvec
!
!  call the lapack device routine
!
IF (n==nsize) THEN
!
!  copy a on the x matrix which will be overwritten with the eigenvectors
!
   DO i=1,nsize
      DO j=i,nsize
         x(i,j,tn) = a(i,j,tn)
      ENDDO
   ENDDO
!
!  save the diagonal part of B 
!
   DO i = 1, nsize
      bdiag(i,tn) = DBLE( b(i,i,tn) )
   END DO
   CALL ZHEGV_XG( itype, jobz, uplo, n, x(1,1,tn), lda, b(1,1,tn), ldb, &
            w(1,tn), work(1,tn), lwork, rwork(1,tn), info(tn) )
ELSE
   vl=0.0_DP
   vu=0.0_DP
   vrange='I'
   abstol=0.D-10
   il=1
   iu=n
!
!  save the diagonal part of A and of B
!
   DO i = 1, nsize
      adiag(i,tn) = DBLE( a(i,i,tn) )
      bdiag(i,tn) = DBLE( b(i,i,tn) )
   END DO
   CALL ZHEGVX_XG( itype, jobz, vrange, uplo, nsize, a(1,1,tn), lda, &
            b(1,1,tn), ldb, vl, vu, il, iu, abstol, m(tn), w(1,tn), &
            x(1,1,tn), ldx, work(1,tn), lwork, rwork(1,tn), iwork(1,tn), &
            ifail(1,tn), info(tn) )
!
! ... restore input a matrix from saved diagonal and lower triangle
!
   DO i = 1, nsize
      a(i,i,tn) = CMPLX( adiag(i,tn), 0.0_DP ,kind=DP)
      DO j = i + 1, nsize
         a(i,j,tn) = CONJG( a(j,i,tn) )
      END DO
      DO j = nsize + 1, lda
         a(j,i,tn) = ( 0.0_DP, 0.0_DP )
      END DO
   END DO
ENDIF
!
! ... restore input b matrix from saved diagonal and lower triangle
!
DO i = 1, nsize
   b(i,i,tn) = CMPLX( bdiag(i,tn), 0.0_DP ,kind=DP)
   DO j = i + 1, nsize
      b(i,j,tn) = CONJG( b(j,i,tn) )
   END DO
   DO j = nsize + 1, ldb
      b(j,i,tn) = ( 0.0_DP, 0.0_DP )
   END DO
END DO
!
END SUBROUTINE diago_dev
#endif
!--------------------------------------------------------------------
SUBROUTINE find_lwork(n, lwork)
!--------------------------------------------------------------------

IMPLICIT NONE
INTEGER :: n, lwork
INTEGER :: ILAENV
EXTERNAL :: ILAENV

INTEGER :: nb

nb = ILAENV( 1, 'ZHETRD', 'U', n, -1, -1, -1 )
!
IF ( nb < 1 .OR. nb >= n) THEN
   !
   lwork = 2*n
   !
ELSE
   !
   lwork = ( nb + 1 )*n
   !
END IF

RETURN
END SUBROUTINE find_lwork
