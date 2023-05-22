INTERFACE dladiv_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE DLADIV_XG( A, B, C, D, P, Q )
#else
      SUBROUTINE DLADIV_XG( A, B, C, D, P, Q )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, P, Q
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: A, B, C, D
      ATTRIBUTES(DEVICE) :: P, Q
#endif

!
      END SUBROUTINE DLADIV_XG

END INTERFACE dladiv_interf
