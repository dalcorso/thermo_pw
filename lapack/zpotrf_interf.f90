INTERFACE zpotrf_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZPOTRF_XG( UPLO, N, A, LDA, INFO )
#else
      SUBROUTINE ZPOTRF_XG( UPLO, N, A, LDA, INFO )
#endif
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  UPLO, N, LDA
      ATTRIBUTES(DEVICE) ::  A, INFO
#endif
!
      END SUBROUTINE ZPOTRF_XG
END INTERFACE zpotrf_interf
