INTERFACE dlagts_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE DLAGTS_XG( JOB, N, A, B, C, D, IN, Y, TOL, INFO )
#else
      SUBROUTINE DLAGTS_XG( JOB, N, A, B, C, D, IN, Y, TOL, INFO )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            INFO, JOB, N
      DOUBLE PRECISION   TOL
!     ..
!     .. Array Arguments ..
      INTEGER            IN( * )
      DOUBLE PRECISION   A( * ), B( * ), C( * ), D( * ), Y( * )
!     ..
!

#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::   JOB, N
      ATTRIBUTES(DEVICE) ::  A, B, C, D, IN, Y, TOL, INFO
#endif
!
!     End of DLAGTS
!
      END SUBROUTINE DLAGTS_XG
END INTERFACE dlagts_interf
