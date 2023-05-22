INTERFACE zlacpy_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZLACPY_XG( UPLO, M, N, A, LDA, B, LDB )
#else
      SUBROUTINE ZLACPY_XG( UPLO, M, N, A, LDA, B, LDB )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * )
!     ..
!
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  UPLO, M, N, LDA, LDB 
      ATTRIBUTES(DEVICE) ::  A, B
#endif
!
!     End of ZLACPY
!
      END SUBROUTINE ZLACPY_XG
END INTERFACE zlacpy_interf
