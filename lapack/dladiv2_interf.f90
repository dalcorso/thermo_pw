INTERFACE dladiv2_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) FUNCTION DLADIV2_XG( A, B, C, D, R, T )
#else
      FUNCTION DLADIV2_XG( A, B, C, D, R, T )
#endif

!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      DOUBLE PRECISION   DLADIV2_XG
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, R, T
!     ..
!
!  =====================================================================
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: A, B, C, D, R, T
#endif
!

      END FUNCTION DLADIV2_XG
END INTERFACE dladiv2_interf

