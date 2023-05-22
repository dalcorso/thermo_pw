INTERFACE disnan_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) FUNCTION DISNAN_XG( DIN )
#else
      FUNCTION DISNAN_XG( DIN )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      LOGICAL :: DISNAN_XG
      DOUBLE PRECISION, INTENT(IN) :: DIN

#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: DIN
#endif

      END FUNCTION DISNAN_XG
END INTERFACE disnan_interf
