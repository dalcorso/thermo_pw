INTERFACE ztrsm_interf

#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZTRSM_XG(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
#else
      SUBROUTINE ZTRSM_XG(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
#endif
!
!  -- Reference BLAS level3 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      COMPLEX*16 ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),B(LDB,*)
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: SIDE,UPLO,TRANSA,DIAG,M,N,LDA,LDB
      ATTRIBUTES(DEVICE) ::  A, B, ALPHA
#endif
!     End of ZTRSM_XG
!
      END SUBROUTINE ZTRSM_XG
END INTERFACE ztrsm_interf

