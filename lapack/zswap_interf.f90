INTERFACE zswap_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZSWAP_XG(N,ZX,INCX,ZY,INCY)
#else
      SUBROUTINE ZSWAP_XG(N,ZX,INCX,ZY,INCY)
#endif
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
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


      END  SUBROUTINE ZSWAP_XG
END INTERFACE zswap_interf
