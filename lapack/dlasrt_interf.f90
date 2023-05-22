INTERFACE dlasrt_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE DLASRT_XG( ID, N, D, INFO )
#else
      SUBROUTINE DLASRT_XG( ID, N, D, INFO )
#endif
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * )
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: ID, N
      ATTRIBUTES(DEVICE) :: D, INFO
#endif

      END SUBROUTINE DLASRT_XG
END INTERFACE dlasrt_interf
