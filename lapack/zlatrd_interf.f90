INTERFACE zlatrd_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZLATRD_XG( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
#else
      SUBROUTINE ZLATRD_XG( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDW, N, NB
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   E( * )
      COMPLEX*16         A( LDA, * ), TAU( * ), W( LDW, * )
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: UPLO, N, NB, LDA, LDW
      ATTRIBUTES(DEVICE) ::  E, A, TAU, W
#endif
!
      END SUBROUTINE ZLATRD_XG
END INTERFACE zlatrd_interf
