INTERFACE zlaset_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZLASET_XG( UPLO, M, N, ALPHA, BETA, A, LDA )
#else
      SUBROUTINE ZLASET_XG( UPLO, M, N, ALPHA, BETA, A, LDA )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      COMPLEX*16         ALPHA, BETA
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  UPLO, M, N, LDA 
      ATTRIBUTES(DEVICE) ::  A,ALPHA, BETA
#endif
!
      END SUBROUTINE ZLASET_XG
END INTERFACE zlaset_interf
