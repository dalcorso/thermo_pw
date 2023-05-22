INTERFACE zgemv_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZGEMV_XG(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
#else
      SUBROUTINE ZGEMV_XG(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
#endif
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      COMPLEX*16 ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),X(*),Y(*)
!     ..
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: TRANS,M,N,LDA,INCX,INCY
      ATTRIBUTES(DEVICE) :: A, X, Y, ALPHA,BETA
#endif
!
!     End of ZGEMV_XG
!
      END SUBROUTINE ZGEMV_XG
END INTERFACE zgemv_interf
