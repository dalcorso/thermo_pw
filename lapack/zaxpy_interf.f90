INTERFACE zaxpy_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZAXPY_XG(N,ZA,ZX,INCX,ZY,INCY)
#else
      SUBROUTINE ZAXPY_XG(N,ZA,ZX,INCX,ZY,INCY)
#endif
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      COMPLEX*16 ZA
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 ZX(*),ZY(*)
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  N,INCX,INCY
      ATTRIBUTES(DEVICE) ::  ZX, ZY,ZA
#endif


      END SUBROUTINE ZAXPY_XG
END INTERFACE zaxpy_interf

