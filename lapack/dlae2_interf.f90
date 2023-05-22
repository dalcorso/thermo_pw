INTERFACE dlae2_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE DLAE2_XG( A, B, C, RT1, RT2 )
#else
      SUBROUTINE DLAE2_XG( A, B, C, RT1, RT2 )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, RT1, RT2
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: A, B, C
      ATTRIBUTES(DEVICE) :: RT1, RT2
#endif
!
      END SUBROUTINE DLAE2_XG
END INTERFACE dlae2_interf
