INTERFACE dlapy3_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) FUNCTION DLAPY3_XG( X, Y, Z )
#else
      FUNCTION DLAPY3_XG( X, Y, Z )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   DLAPY3_XG
      DOUBLE PRECISION   X, Y, Z
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  X, Y, Z
#endif
!
      END  FUNCTION DLAPY3_XG
END INTERFACE dlapy3_interf
