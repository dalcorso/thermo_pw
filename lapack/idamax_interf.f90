INTERFACE idamax_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) INTEGER FUNCTION IDAMAX_XG(N,DX,INCX)
#else
      INTEGER FUNCTION IDAMAX_XG(N,DX,INCX)
#endif

!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  N,INCX
      ATTRIBUTES(DEVICE) ::  DX
#endif
!
!     End of IDAMAX
!
      END FUNCTION idamax_xg
!
END INTERFACE idamax_interf
