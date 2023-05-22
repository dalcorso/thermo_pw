INTERFACE dscal_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE DSCAL_XG(N,DA,DX,INCX)
#else 
      SUBROUTINE DSCAL_XG(N,DA,DX,INCX)
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
      DOUBLE PRECISION DX(*)
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  N,DA,INCX
      ATTRIBUTES(DEVICE) ::  DX
#endif
!
      END SUBROUTINE DSCAL_XG
END INTERFACE dscal_interf

