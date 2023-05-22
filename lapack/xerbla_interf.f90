INTERFACE xerbla_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE XERBLA_XG( SRNAME, INFO )
#else
      SUBROUTINE XERBLA_XG( SRNAME, INFO )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
!      CHARACTER*(*)      SRNAME
      CHARACTER          SRNAME
      INTEGER            INFO
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: SRNAME
      ATTRIBUTES(DEVICE) :: INFO
#endif
      END SUBROUTINE XERBLA_XG
END INTERFACE xerbla_interf
