INTERFACE zhetd2_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZHETD2_XG( UPLO, N, A, LDA, D, E, TAU, INFO )
#else
      SUBROUTINE ZHETD2_XG( UPLO, N, A, LDA, D, E, TAU, INFO )
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
      DOUBLE PRECISION   D( * ), E( * )
      COMPLEX*16         A( LDA, * ), TAU( * )
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: UPLO, N, LDA
      ATTRIBUTES(DEVICE) :: D, E, A, TAU, INFO
#endif
!
      END SUBROUTINE ZHETD2_XG
END INTERFACE zhetd2_interf
