INTERFACE zhemm_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZHEMM_XG(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
#else
      SUBROUTINE ZHEMM_XG(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
#endif
!
!  -- Reference BLAS level3 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      COMPLEX*16 ALPHA,BETA
      INTEGER LDA,LDB,LDC,M,N
      CHARACTER SIDE,UPLO
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
!     ..
!

#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: SIDE,UPLO,M,N,LDA,LDB,LDC
      ATTRIBUTES(DEVICE) ::  A, B, C, ALPHA,BETA
#endif
!
      END SUBROUTINE ZHEMM_XG
END INTERFACE zhemm_interf

