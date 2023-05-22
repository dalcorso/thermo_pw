INTERFACE dlanst_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) FUNCTION DLANST_XG( NORM, N, D, E )
#else
      FUNCTION DLANST_XG( NORM, N, D, E )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   DLANST_XG
      CHARACTER          NORM
      INTEGER            N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
!     ..
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  NORM, N
      ATTRIBUTES(DEVICE) ::  D, E
#endif
!
      END FUNCTION DLANST_XG
END INTERFACE dlanst_interf
