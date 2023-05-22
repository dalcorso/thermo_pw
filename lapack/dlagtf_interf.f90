INTERFACE dlagtf_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE DLAGTF_XG( N, A, LAMBDA, B, C, TOL, D, IN, INFO )
#else
      SUBROUTINE DLAGTF_XG( N, A, LAMBDA, B, C, TOL, D, IN, INFO )
#endif
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            INFO, N
      DOUBLE PRECISION   LAMBDA, TOL
!     ..
!     .. Array Arguments ..
      INTEGER            IN( * )
      DOUBLE PRECISION   A( * ), B( * ), C( * ), D( * )
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::   N, LAMBDA, TOL 
      ATTRIBUTES(DEVICE) ::  IN, A, B, C, D,INFO
#endif

!     End of DLAGTF
!
      END SUBROUTINE dlagtf_xg
END INTERFACE dlagtf_interf
