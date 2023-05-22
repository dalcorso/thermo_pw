INTERFACE zlarfg_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZLARFG_XG( N, ALPHA, X, INCX, TAU )
#else
      SUBROUTINE ZLARFG_XG( N, ALPHA, X, INCX, TAU )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      COMPLEX*16         ALPHA, TAU
!     ..
!     .. Array Arguments ..
      COMPLEX*16         X( * )
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  N, INCX
      ATTRIBUTES(DEVICE) ::  X,ALPHA, TAU
#endif
      END SUBROUTINE ZLARFG_XG
END INTERFACE zlarfg_interf
