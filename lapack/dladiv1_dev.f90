#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE DLADIV1_XG( A, B, C, D, P, Q )
#else
      SUBROUTINE DLADIV1_XG( A, B, C, D, P, Q )
#endif

#include<dladiv2_interf.f90>

!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, P, Q
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
!
!     .. Local Scalars ..
      DOUBLE PRECISION   R, T
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLADIV2_XG
!      EXTERNAL           DLADIV2_XG
!     ..
!     .. Executable Statements ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: A, B, C, D
      ATTRIBUTES(DEVICE) :: P, Q
#endif

      R = D / C
      T = ONE / (C + D * R)
      P = DLADIV2_XG(A, B, C, D, R, T)
      A = -A
      Q = DLADIV2_XG(B, A, C, D, R, T)
!
      RETURN
!
!     End of DLADIV1_XG
!
      END


