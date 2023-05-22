INTERFACE zgerc_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZGERC_XG(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
#else
      SUBROUTINE ZGERC_XG(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
#endif
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      COMPLEX*16 ALPHA
      INTEGER INCX,INCY,LDA,M,N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),X(*),Y(*)
!     ..
!

#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: M,N,INCX,INCY,LDA
      ATTRIBUTES(DEVICE) :: A, X, Y, ALPHA
#endif


      END SUBROUTINE ZGERC_XG
END INTERFACE zgerc_interf
