INTERFACE zungtr_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZUNGTR_XG( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
#else
      SUBROUTINE ZUNGTR_XG( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
#endif
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  UPLO, N, LDA, LWORK
      ATTRIBUTES(DEVICE) ::  A, TAU, WORK, INFO
#endif

      END SUBROUTINE ZUNGTR_XG
END INTERFACE zungtr_interf
