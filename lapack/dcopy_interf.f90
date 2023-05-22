INTERFACE dcopy_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE DCOPY_XG(N,DX,INCX,DY,INCY)
#else
      SUBROUTINE DCOPY_XG(N,DX,INCX,DY,INCY)
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
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!     ..
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: N,INCX,INCY
      ATTRIBUTES(DEVICE) :: DX, DY
#endif

!
!     End of DCOPY
!
      END
END INTERFACE dcopy_interf
