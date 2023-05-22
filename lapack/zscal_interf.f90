INTERFACE zscal_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZSCAL_XG(N,ZA,ZX,INCX)
#else
      SUBROUTINE ZSCAL_XG(N,ZA,ZX,INCX)
#endif
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      COMPLEX*16 ZA
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 ZX(*)
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  N,INCX
      ATTRIBUTES(DEVICE) ::  ZX,ZA
#endif 
!
      END SUBROUTINE ZSCAL_XG
END INTERFACE zscal_interf

