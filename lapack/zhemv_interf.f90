INTERFACE zhemv_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZHEMV_XG(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
#else
      SUBROUTINE ZHEMV_XG(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
#endif
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      COMPLEX*16 ALPHA,BETA
      INTEGER INCX,INCY,LDA,N
      CHARACTER UPLO
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),X(*),Y(*)
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: UPLO,N,LDA,INCX,INCY
      ATTRIBUTES(DEVICE) ::  A, X, Y,ALPHA,BETA
#endif


      END SUBROUTINE ZHEMV_XG
END INTERFACE zhemv_interf
