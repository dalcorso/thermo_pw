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
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
!
!     .. Local Scalars ..
      DOUBLE PRECISION   BR
!     ..
!     .. Executable Statements ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: A, B, C, D, R, T
#endif

      IF( R.NE.ZERO ) THEN
         BR = B * R
         IF( BR.NE.ZERO ) THEN
            DLADIV2_XG = (A + BR) * T
         ELSE
            DLADIV2_XG = A * T + (B * T) * R
         END IF
      ELSE
         DLADIV2_XG = (A + D * (B / C)) * T
      END IF
!
      RETURN
!
!     End of DLADIV2_XG
!
      END


