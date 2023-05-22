INTERFACE zcopy_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZCOPY_XG(N,ZX,INCX,ZY,INCY)
#else
      SUBROUTINE ZCOPY_XG(N,ZX,INCX,ZY,INCY)
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
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  N,INCX,INCY
      ATTRIBUTES(DEVICE) ::  ZX, ZY
#endif

      END SUBROUTINE ZCOPY_XG
END INTERFACE zcopy_interf
