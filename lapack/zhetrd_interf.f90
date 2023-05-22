INTERFACE zhetrd_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZHETRD_XG( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
#else
      SUBROUTINE ZHETRD_XG( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
#endif
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: UPLO, N, LDA, LWORK
      ATTRIBUTES(DEVICE) ::  D, E, A, TAU, WORK, INFO
#endif

      END  SUBROUTINE ZHETRD_XG
END INTERFACE zhetrd_interf
