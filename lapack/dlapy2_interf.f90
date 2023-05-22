INTERFACE dlapy2_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) FUNCTION DLAPY2_XG( X, Y )
#else
      FUNCTION DLAPY2_XG( X, Y )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   DLAPY2_XG
      DOUBLE PRECISION   X, Y
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  X, Y
#endif
!
      END FUNCTION DLAPY2_XG
END INTERFACE dlapy2_interf
