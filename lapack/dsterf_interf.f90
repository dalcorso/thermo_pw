INTERFACE dsterf_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE DSTERF_XG( N, D, E, INFO )
#else
      SUBROUTINE DSTERF_XG( N, D, E, INFO )
#endif
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::   N
      ATTRIBUTES(DEVICE) ::  D, E, INFO
#endif
      END SUBROUTINE DSTERF_XG
END INTERFACE dsterf_interf
