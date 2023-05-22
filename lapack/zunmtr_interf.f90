INTERFACE zunmtr_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZUNMTR_XG( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC, &
                         WORK, LWORK, INFO )
#else
      SUBROUTINE ZUNMTR_XG( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC, &
                         WORK, LWORK, INFO )
#endif



!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS, UPLO
      INTEGER            INFO, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  SIDE, UPLO, TRANS, M, N, LDA, LDC, LWORK
      ATTRIBUTES(DEVICE) ::    A, C, TAU, WORK, INFO
#endif

!
!     End of ZUNMTR
!
      END SUBROUTINE ZUNMTR_XG
END INTERFACE zunmtr_interf
