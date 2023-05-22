INTERFACE zunmqr_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZUNMQR_XG( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, LWORK, INFO )
#else
      SUBROUTINE ZUNMQR_XG( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, LWORK, INFO )
#endif
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  SIDE, TRANS, M, N, K, LDA, LDC, LWORK
      ATTRIBUTES(DEVICE) ::  A, C, TAU, WORK, INFO
#endif

!
!     End of ZUNMQR
!
      END SUBROUTINE ZUNMQR_XG
END INTERFACE zunmqr_interf
