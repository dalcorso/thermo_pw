INTERFACE zher2_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZHER2_XG(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
#else
      SUBROUTINE ZHER2_XG(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
#endif
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      COMPLEX*16 ALPHA
      INTEGER INCX,INCY,LDA,N
      CHARACTER UPLO
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),X(*),Y(*)
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: UPLO,N,INCX,INCY,LDA
      ATTRIBUTES(DEVICE) ::  A, X, Y,ALPHA
#endif

      END SUBROUTINE ZHER2_XG
END INTERFACE zher2_interf
