INTERFACE dlamch_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) FUNCTION DLAMCH_XG( CMACH )
#else
      FUNCTION DLAMCH_XG( CMACH )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      DOUBLE PRECISION   DLAMCH_XG
!     .. Scalar Arguments ..
      CHARACTER          CMACH
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  CMACH
#endif
!
      END  FUNCTION DLAMCH_XG
END INTERFACE dlamch_interf
