INTERFACE zdotc_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) FUNCTION ZDOTC_XG(N,ZX,INCX,ZY,INCY)
#else
      FUNCTION ZDOTC_XG(N,ZX,INCX,ZY,INCY)
#endif
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      COMPLEX*16 ZDOTC_XG
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 ZX(*),ZY(*)
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  N,INCX,INCY
      ATTRIBUTES(DEVICE) ::  ZX, ZY
#endif
!
      END FUNCTION ZDOTC_XG
END INTERFACE zdotc_interf

