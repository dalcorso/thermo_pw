INTERFACE zherk_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZHERK_XG(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
#else
      SUBROUTINE ZHERK_XG(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
#endif
!
!  -- Reference BLAS level3 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDC,N
      CHARACTER TRANS,UPLO
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),C(LDC,*)
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: UPLO,TRANS,N,K,ALPHA,LDA,BETA,LDC
      ATTRIBUTES(DEVICE) ::  A,C
#endif
!
      END SUBROUTINE ZHERK_XG
END INTERFACE zherk_interf

