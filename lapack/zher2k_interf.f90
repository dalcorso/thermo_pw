INTERFACE zher2k_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZHER2K_XG(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
#else
      SUBROUTINE ZHER2K_XG(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
#endif
!
!  -- Reference BLAS level3 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      COMPLEX*16 ALPHA
      DOUBLE PRECISION BETA
      INTEGER K,LDA,LDB,LDC,N
      CHARACTER TRANS,UPLO
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
!     ..
!

#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: UPLO,TRANS,N,K,LDA,LDB,BETA,LDC
      ATTRIBUTES(DEVICE) ::  A, B, C,ALPHA
#endif

      END SUBROUTINE ZHER2K_XG
END INTERFACE zher2k_interf

