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
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER I,NINCX
!     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D+0)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE, DCMPLX, DIMAG
!     ..

#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: N,DA,INCX
      ATTRIBUTES(DEVICE) :: ZX
#endif

      IF (N.LE.0 .OR. INCX.LE.0 .OR. DA.EQ.ONE) RETURN
      IF (INCX.EQ.1) THEN
!
!        code for increment equal to 1
!
         DO I = 1,N
            ZX(I) = DCMPLX(DA*DBLE(ZX(I)),DA*DIMAG(ZX(I)))
         END DO
      ELSE
!
!        code for increment not equal to 1
!
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            ZX(I) = DCMPLX(DA*DBLE(ZX(I)),DA*DIMAG(ZX(I)))
         END DO
      END IF
      RETURN
!
!     End of ZDSCAL_XG
!
      END

