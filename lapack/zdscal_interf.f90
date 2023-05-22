INTERFACE zdscal_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZDSCAL_XG(N,DA,ZX,INCX)
#else
      SUBROUTINE ZDSCAL_XG(N,DA,ZX,INCX)
#endif
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 ZX(*)
!     ..

#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: N,DA,INCX
      ATTRIBUTES(DEVICE) :: ZX
#endif
      END SUBROUTINE ZDSCAL_XG
END INTERFACE zdscal_interf

