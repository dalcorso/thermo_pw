INTERFACE zlanhe_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) FUNCTION ZLANHE_XG( NORM, UPLO, N, A, LDA, WORK )
#else
      FUNCTION ZLANHE_XG( NORM, UPLO, N, A, LDA, WORK )
#endif

!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   ZLANHE_XG
      CHARACTER          NORM, UPLO
      INTEGER            LDA, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   WORK( * )
      COMPLEX*16         A( LDA, * )
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  NORM, UPLO, N, LDA
      ATTRIBUTES(DEVICE) ::  A, WORK
#endif
!

      END  FUNCTION ZLANHE_XG
END INTERFACE zlanhe_interf
