INTERFACE ieeeck_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) FUNCTION IEEECK_XG( ISPEC, ZERO, ONE )
#else
      FUNCTION IEEECK_XG( ISPEC, ZERO, ONE )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IEEECK_XG
      INTEGER            ISPEC
      REAL               ONE, ZERO
!     ..
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::   ISPEC, ZERO, ONE
#endif
!
      END  FUNCTION IEEECK_XG
END INTERFACE ieeeck_interf

