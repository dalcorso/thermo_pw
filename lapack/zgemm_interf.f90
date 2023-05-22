INTERFACE zgemm_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZGEMM_XG(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
#else
      SUBROUTINE ZGEMM_XG(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
#endif
!
!  -- Reference BLAS level3 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      COMPLEX*16 ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
!     ..
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: TRANSA,TRANSB,M,N,K,LDA,LDB,LDC
      ATTRIBUTES(DEVICE) :: A, B, C,ALPHA,BETA
#endif
!
      END SUBROUTINE ZGEMM_XG
END INTERFACE zgemm_interf

