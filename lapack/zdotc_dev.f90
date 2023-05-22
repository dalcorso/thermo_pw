#if defined(__CUDA)
      ATTRIBUTES(DEVICE) COMPLEX*16 FUNCTION ZDOTC_XG(N,ZX,INCX,ZY,INCY)
#else
      COMPLEX*16 FUNCTION ZDOTC_XG(N,ZX,INCX,ZY,INCY)
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
!  =====================================================================
!
!     .. Local Scalars ..
      COMPLEX*16 ZTEMP
      INTEGER I,IX,IY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCONJG
!     ..

#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  N,INCX,INCY
      ATTRIBUTES(DEVICE) ::  ZX, ZY
#endif

      ZTEMP = (0.0d0,0.0d0)
      ZDOTC_XG = (0.0d0,0.0d0)
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!        code for both increments equal to 1
!
         DO I = 1,N
            ZTEMP = ZTEMP + DCONJG(ZX(I))*ZY(I)
         END DO
      ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            ZTEMP = ZTEMP + DCONJG(ZX(IX))*ZY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      ZDOTC_XG = ZTEMP
      RETURN
!
!     End of ZDOTC_XG
!
      END

