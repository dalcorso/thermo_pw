INTERFACE zunm2l_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZUNM2L_XG( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, INFO )
#else
      SUBROUTINE ZUNM2L_XG( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, INFO )
#endif
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: SIDE, TRANS, M, N, K, LDA, LDC
      ATTRIBUTES(DEVICE) ::   A, C, TAU, WORK, INFO
#endif
!
!     End of ZUNM2L
!
      END SUBROUTINE ZUNM2L_XG
END INTERFACE zunm2l_interf
