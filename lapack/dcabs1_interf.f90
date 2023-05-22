INTERFACE dcabs1_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) FUNCTION DCABS1_XG(Z)
#else
      FUNCTION DCABS1_XG(Z)
#endif
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      COMPLEX*16 Z
!     ..
!     ..
!  =====================================================================
      DOUBLE PRECISION  DCABS1_XG
!
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) :: Z
#endif

      END FUNCTION DCABS1_XG
END INTERFACE dcabs1_interf
