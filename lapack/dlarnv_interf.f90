INTERFACE dlarnv_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE DLARNV_XG( IDIST, ISEED, N, X )
#else
      SUBROUTINE DLARNV_XG( IDIST, ISEED, N, X )
#endif

!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IDIST, N
!     ..
!     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   X( * )
!     ..
!
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  IDIST, N
      ATTRIBUTES(DEVICE) ::  ISEED, X
#endif

!     End of DLARNV
!
      END

END INTERFACE dlarnv_interf
