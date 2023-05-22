INTERFACE dlaisnan_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) FUNCTION DLAISNAN_XG( DIN1, DIN2 )
#else
      FUNCTION DLAISNAN_XG( DIN1, DIN2 )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      LOGICAL DLAISNAN_XG
!     .. Scalar Arguments ..
      DOUBLE PRECISION, INTENT(IN) :: DIN1, DIN2
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  DIN1, DIN2
#endif
!  =====================================================================
!
      END FUNCTION DLAISNAN_XG
END INTERFACE dlaisnan_interf
