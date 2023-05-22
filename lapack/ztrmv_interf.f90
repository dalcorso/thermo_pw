INTERFACE ztrmv_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZTRMV_XG(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
#else
      SUBROUTINE ZTRMV_XG(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
#endif
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),X(*)
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  UPLO,TRANS,DIAG,N,LDA,INCX
      ATTRIBUTES(DEVICE) ::  A, X
#endif
!
      END  SUBROUTINE ZTRMV_XG
END INTERFACE ztrmv_interf
